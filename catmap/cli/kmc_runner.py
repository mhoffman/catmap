#!/usr/bin/env python

import os
import kmos.run
import catmap
import pickle

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import numpy as np


INIT_STEPS = int(1e6)
SAMPLE_STEPS = INIT_STEPS
SEED = 'CO_oxidation'
TEMPERATURE = 500
DIFFUSION_FACTOR = 1e-3

from matplotlib.mlab import griddata

def setup_model(model, data_point=0, data_file='CO_oxidation.pkl'):
    import numpy.random
    import pickle

    with open(data_file) as infile:
        data = pickle.load(infile)
    coverage = data['coverage_map'][data_point][1]

    choices = [
                model.proclist.co,
                model.proclist.o,
                model.proclist.empty,
                ]
    choices_weights = [
        coverage[0],
        coverage[1],
        1 - sum(coverage),
    ]

    config = model._get_configuration()
    X, Y, Z = model.lattice.system_size
    for x in range(X):
        for y in range(Y):
            choice = numpy.random.choice(
                choices,
                p=choices_weights,
            )
            config[x, y, 0, 0] = choice
            #print('{x} {y} {choice}'.format(**locals()))

    model._set_configuration(config)

def set_rate_constants(kmos_model, catmap_data, data_point, diffusion_factor=None):
    """
        A rate constants of a kmos model from a corresponding CatMAP model for
        a given data point (i.e. a tuple of of reactivity descriptors).

        Diffusion factors has a special role of mimicking the effect of diffusion
        in a kMC model as described typically by a mean field model. If the diffusion-factor
        is left at its default value (None) is it simply faithfully set to the corresponding
        CatMAP value which will usually be the normal prefactor (kT/h). If instead it is
        set to a finite float it will be set to diffusion_factor * max_rate_constant
        where max_rate_constant is the fastet rate constant of all non-diffusion
        reaction steps.

    """

    # set rate constant of kMC according to current descriptor tuple
    max_rate_constant = float('-inf')

    for i in range(len(catmap_data['forward_rate_constant_map'][data_point][1])):
        forward_rate_constant = catmap_data['forward_rate_constant_map'][data_point][1][i]
        reverse_rate_constant = catmap_data['reverse_rate_constant_map'][data_point][1][i]

        if hasattr(kmos_model.parameters, 'forward_{i}'.format(**locals())):
            max_rate_constant = max(max_rate_constant, forward_rate_constant)
            setattr(kmos_model.parameters, 'forward_{i}'.format(
                **locals()), forward_rate_constant)
        if hasattr(kmos_model.parameters, 'reverse_{i}'.format(**locals())):
            max_rate_constant = max(max_rate_constant, reverse_rate_constant)
            setattr(kmos_model.parameters, 'reverse_{i}'.format(
                **locals()), reverse_rate_constant)

    for i in range(len(catmap_data['forward_rate_constant_map'][data_point][1])):
        forward_rate_constant = catmap_data['forward_rate_constant_map'][data_point][1][i] if diffusion_factor is None else max_rate_constant * diffusion_factor
        reverse_rate_constant = catmap_data['reverse_rate_constant_map'][data_point][1][i] if diffusion_factor is None else max_rate_constant * diffusion_factor

        if hasattr(kmos_model.parameters, 'diff_forward_{i}'.format(**locals())):
            setattr(kmos_model.parameters, 'diff_forward_{i}'.format(
                **locals()), forward_rate_constant)
        if hasattr(kmos_model.parameters, 'diff_reverse_{i}'.format(**locals())):
            setattr(kmos_model.parameters, 'diff_reverse_{i}'.format(
                **locals()), reverse_rate_constant)

def set_rate_constants_from_descriptors(kmos_model, catmap_model, descriptors, diffusion_factor=None):
    # set rate constant of kMC according to current descriptor tuple
    max_rate_constant = float('-inf')

    rate_constants = catmap_model.get_rate_constants(descriptors)
    n_rate_constants = len(rate_constants)

    forward_rate_constants = rate_constants[: n_rate_constants / 2]
    reverse_rate_constants = rate_constants[n_rate_constants / 2 :]

    for i, (forward_rate_constant, reverse_rate_constant) in \
        enumerate(zip(forward_rate_constants, reverse_rate_constants)):

        if hasattr(kmos_model.parameters, 'forward_{i}'.format(**locals())):
            max_rate_constant = max(max_rate_constant, forward_rate_constant)
            setattr(kmos_model.parameters, 'forward_{i}'.format(
                **locals()), forward_rate_constant)
        if hasattr(kmos_model.parameters, 'reverse_{i}'.format(**locals())):
            max_rate_constant = max(max_rate_constant, reverse_rate_constant)
            setattr(kmos_model.parameters, 'reverse_{i}'.format(
                **locals()), reverse_rate_constant)

    for i, (forward_rate_constant, reverse_rate_constant) in \
        enumerate(zip(forward_rate_constants, reverse_rate_constants)):

        if diffusion_factor is not None:
            forward_rate_constant = max_rate_constant * diffusion_factor
            reverse_rate_constant = max_rate_constant * diffusion_factor

        if hasattr(kmos_model.parameters, 'diff_forward_{i}'.format(**locals())):
            setattr(kmos_model.parameters, 'diff_forward_{i}'.format(
                **locals()), forward_rate_constant)
        if hasattr(kmos_model.parameters, 'diff_reverse_{i}'.format(**locals())):
            setattr(kmos_model.parameters, 'diff_reverse_{i}'.format(
                **locals()), reverse_rate_constant)



def run_model(seed, init_steps, sample_steps):
    data_filename = '{seed}_kMC_output.log'.format(**locals())
    lock_filename = '{seed}.lock'.format(**locals())
    done_filename = '{seed}.done'.format(**locals())

    catmap_model = catmap.ReactionModel(
        setup_file='{seed}.mkm'.format(**locals()))

    catmap_model.output_variables.append('rate_constant')
    catmap_model.output_variables.append('forward_rate_constant')
    catmap_model.output_variables.append('reverse_rate_constant')
    catmap_model.run()

    with open('{seed}.pkl'.format(**locals())) as infile:
        catmap_data = pickle.load(infile)

    print(catmap_data.keys())

    data_point = 1

    if not os.path.exists(lock_filename):
        with open(lock_filename, 'w'):
            pass

    if not os.path.exists(data_filename):
        with open(data_filename, 'w') as outfile:
            with kmos.run.KMC_Model(print_rates=False, banner=True) as kmos_model:
                data_header = kmos_model.get_std_header()[1:]
                outfile.write(
                    '# descriptor0 descriptor1 {data_header}'.format(**locals()))

    for data_point in range(len(catmap_data['forward_rate_constant_map'])):
        descriptors = catmap_data['forward_rate_constant_map'][data_point][0]

        # multi IO mechanism: keep lockfile with one line descriptors of datapoint
        # if datapoint is already in there, skip to next datapoint
        descriptor_string = str(descriptors) + '\n'

        with open(lock_filename, 'r') as lockfile:
            if descriptor_string in lockfile.readlines():
                print('Skipping {descriptor_string}'.format(**locals()))
                continue
        with open(lock_filename, 'a') as lockfile:
            lockfile.write('{descriptor_string}'.format(**locals()))
            lockfile.flush()

        with kmos.run.KMC_Model(print_rates=False, banner=False) as kmos_model:
            print('running DATAPOINT {data_point} DESCRIPTOR {descriptor_string}'.format(**locals()))
            set_rate_constants(kmos_model, catmap_data, data_point, diffusion_factor=DIFFUSION_FACTOR)
            setup_model(kmos_model, data_point)

            print(kmos_model.rate_constants)
            print('INIT STEPS {init_steps} |  SAMPLE STEPS {sample_steps}'.format(**locals()))
            kmos_model.do_steps(init_steps)
            atoms = kmos_model.get_atoms()
            data = kmos_model.get_std_sampled_data(1, sample_steps, tof_method='integ')

            #with open('procstat_{data_point}.log'.format(**locals()), 'w') as outfile:
                #outfile.write(kmos_model.print_procstat(False))

            #with open('rateconstants_{data_point}.log'.format(**locals()), 'w') as outfile:
                #outfile.write(kmos_model.rate_constants)

            with open(data_filename, 'a') as outfile:
                outfile.write(
                    '{descriptors[0]} {descriptors[1]} {data}'.format(**locals()))

        with open(done_filename, 'a') as outfile:
            outfile.write('{descriptor_string}'.format(**locals()))

def contour_plot_data(x, y, z, filename,
                      n_gp=101,
                      m_gp=20,
                      title='',
                      xlabel='',
                      ylabel='',
                      zmin=None,
                      zmax=None,
                      xlabel_unit='',
                      ylabel_unit='',
                      ticks=None,
                      seed=None,
                      catmap_model=None,
                      normalized=False,
                      colorbar_label=None):
    import numpy
    import scipy.interpolate
    fig = plt.figure()

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)


    # with golden ration and the whole shebang ...
    # settings size and font for revtex stylesheet
    fig_width_pt = 2*246.0  # Get this from LaTeX using \showthe\columnwidth
    #fig_width_pt *= 300./72 # convert to 300 dpi
    inches_per_pt = 1.0/72.27               # Convert pt to inches
    #inches_per_pt = 1.0/300               # Convert pt to inches
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean       # height in inches
    fig_size = [fig_width, 1.3 * fig_height]

    font_size = 10
    tick_font_size = 10
    xlabel_pad = 8
    ylabel_pad = 8
    matplotlib.rcParams['ps.usedistiller'] = 'xpdf'
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.serif'] = 'Gill Sans'
    matplotlib.rcParams['font.sans-serif'] = 'Gill Sans'
    matplotlib.rcParams['text.usetex'] = 'true'

    matplotlib.rcParams['lines.linewidth'] = 1.

    fig = plt.figure(figsize=fig_size)

    #x, y = np.linspace(x.min(), x.max(), m_gp), np.linspace(y.min(), y.max(), m_gp)

    print(x, y)
    print(z)
    xi, yi = np.linspace(x.min(), x.max(), n_gp), np.linspace(y.min(), y.max(), n_gp)
    xi, yi = np.meshgrid(xi, yi)

    #print(z)
    rbf = scipy.interpolate.Rbf(x, y, z, function='linear',)
    zi = rbf(xi, yi)
    #zi = griddata(x, y, z, xi, yi, interp='linear')

    if zmin is None :
        zmin = z.min()
    if zmax is None :
        zmax = z.max()

    if normalized:
        levels = np.linspace(0, 1, 11)
        print('NORMALIZED {zmin} {zmax}'.format(**locals()))
        zmax = 1.
        zmin = 0.
    else:
        levels = np.linspace(zmin, zmax, int(zmax - zmin))


    contour_plot = plt.contourf(zi, vmin=zmin, vmax=zmax, origin='lower',
               extent=[x.min(), x.max(), y.min(), y.max()],
               levels=levels,
               extend='both')


    #plt.scatter(x, y, c=z, s=.2)

    cbar = plt.colorbar(contour_plot, ticks=ticks)

    if colorbar_label is None:
        if normalized:
            cbar.set_label(r'${\rm ML}$')
        else:
            cbar.set_label(r'${\rm s}^{-1} {\rm cell}^{-1}$')
    else:
        cbar.set_label(colorbar_label)

    if catmap_model is not None:
        for substrate, (xi, yi) in catmap_model.descriptor_dict.items():
            z = rbf(xi, yi)
            print(substrate, z)
            plt.scatter(xi, yi, c=z, s=35, vmin=z.min(), vmax=z.max(), cmap=cbar.cmap)
            plt.annotate(substrate, xy=(xi +.05, yi + .05), size='small',
             bbox={'facecolor':'white', 'alpha':0.5, 'ec':'white', 'pad':1, 'lw':0 })


    if seed is not None:
        model = catmap.ReactionModel(setup_file='{seed}.mkm'.format(**locals()))
        plt.xlabel(r'${{\rm {} }}$ [{xlabel_unit}]'.format(model.descriptor_names[0], **locals()))
        plt.ylabel(r'${{\rm {} }}$ [{ylabel_unit}]'.format(model.descriptor_names[1], **locals()))
        plt.title(title)

    plt.xlim((x.min(), x.max()))
    plt.ylim((y.min(), y.max()))

    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)

    plt.axis('image')
    plt.xticks(np.arange(x.min(), x.max() + .5, .5))
    plt.yticks(np.arange(y.min(), y.max() + .5, .5))

    plt.savefig(filename, bbox_inches='tight')

def main(options):
    if not options.dontrun:
        run_model(seed=SEED,
             init_steps=INIT_STEPS,
             sample_steps=SAMPLE_STEPS)

    if options.plot:
        data = np.recfromtxt('{SEED}_kMC_output.log'.format(**locals()), names=True)

        for name in data.dtype.names:
            print(name)
            if name in ['descriptor1', 'descriptor0']:
                continue
            if '_2_' in name: # we are plotting a rate
                normalized = False
                plot_data = np.log10(data[name])
                #plot_data = data[name]
                if not np.isfinite(plot_data).any():
                    plot_data[:] = 0.
                else:
                    minimum = np.nanmin(plot_data[np.isfinite(plot_data)])
                    print('MINIMUM {minimum}'.format(**locals()))
                    plot_data[np.logical_or(np.isnan(plot_data),
                                            np.isinf(plot_data))] = minimum
            else: # we are plotting a coverage
                plot_data = data[name]
                normalized = True

            title = r'${{\rm {}}}$'.format(name \
                       .replace('_2_', r' \rightarrow ') \
                       .replace('_n_', ' + ') \
                       .replace('_default', '') \
                       .replace('empty', '*') \
                       .replace('_0', ''))

            catmap_model = catmap.ReactionModel(
                setup_file='{SEED}.mkm'.format(**locals()))

            catmap_model.run()

            zmin = None
            zmax = None
            ticks = None
            xlabel = 'O reactivity [eV]'
            ylabel = 'CO reactivity [eV]'
            if name == 'CO_s_n_O_s_2_empty_s_n_empty_s_0':
                xlabel = 'O reactivity [eV]'
                ylabel = 'CO reactivity [eV]'
                zmin = -48
                zmax = 2
                ticks = range(zmin, zmax+1, 6)

            contour_plot_data(data['descriptor0'],
                              data['descriptor1'],
                              plot_data,
                              'output_{name}.pdf'.format(**locals()),
                              #seed=SEED,
                              catmap_model=catmap_model,
                              normalized=normalized,
                              title=title,
                              zmin=zmin,
                              zmax=zmax,
                              ticks=ticks,
                              xlabel_unit='eV',
                              ylabel_unit='eV',
                              xlabel=xlabel,
                              ylabel=ylabel,
                              )

if __name__ == '__main__':
    main()
