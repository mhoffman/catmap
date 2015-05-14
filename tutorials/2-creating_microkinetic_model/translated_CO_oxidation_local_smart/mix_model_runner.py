#!/usr/bin/env python

import os
import kmos.run
import catmap
import pickle

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import numpy as np


INIT_STEPS = int(1e8)
SAMPLE_STEPS = INIT_STEPS

SEED = 'CO_oxidation'
TEMPERATURE = 500

from matplotlib.mlab import griddata


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
            set_rate_constants(kmos_model, catmap_data, data_point, diffusion_factor=10)
            from setup_model import setup_model
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

def contour_plot_data(x, y, z, filename, n_gp=101, m_gp=20, title='', seed=None, normalized=False):
    import numpy
    import scipy.interpolate
    fig = plt.figure()

    # with golden ration and the whole shebang ...
    # settings size and font for revtex stylesheet
    fig_width_pt = 2*246.0  # Get this from LaTeX using \showthe\columnwidth
    #fig_width_pt *= 300./72 # convert to 300 dpi
    inches_per_pt = 1.0/72.27               # Convert pt to inches
    #inches_per_pt = 1.0/300               # Convert pt to inches
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height =fig_width*golden_mean       # height in inches
    fig_size = [fig_width,fig_height]

    font_size = 10
    tick_font_size = 10
    xlabel_pad = 8
    ylabel_pad = 8
    #matplotlib.rcParams['ps.usedistiller'] = 'xpdf'
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['font.serif'] = 'Computer Modern Roman'
    matplotlib.rcParams['font.sans-serif'] = 'Computer Modern Sans serif'
    matplotlib.rcParams['text.usetex'] = 'true'



    matplotlib.rcParams['lines.linewidth'] = 1.

    fig = plt.figure(figsize=fig_size)



    #x, y = np.linspace(x.min(), x.max(), m_gp), np.linspace(y.min(), y.max(), m_gp)

    xi, yi = np.linspace(x.min(), x.max(), n_gp), np.linspace(y.min(), y.max(), n_gp)
    xi, yi = np.meshgrid(xi, yi)

    #print(z)
    rbf = scipy.interpolate.Rbf(x, y, z, function='linear', )
    zi = rbf(xi, yi)
    #zi = griddata(x, y, z, xi, yi, interp='linear')

    if normalized:
        levels = np.linspace(0, 1, 11)
        zmin = z.min()
        zmax = z.max()
        print('NORMALIZED {zmin} {zmax}'.format(**locals()))
    else:
        levels = np.linspace(zi.min(), zi.max(), 21)

    contour_plot = plt.contourf(zi, vmin=z.min(), vmax=z.max(), origin='lower',
               extent=[x.min(), x.max(), y.min(), y.max()],
               levels=levels,
               extend='both')

    plt.scatter(x, y, c=z, s=1.5)
    cbar = plt.colorbar(contour_plot)
    if normalized:
        cbar.set_label(r'${\rm ML}$')
    else:
        cbar.set_label(r'${\rm s}^{-1} {\rm cell}^{-1}$')

    if seed is not None:
        model = catmap.ReactionModel(setup_file='{seed}.mkm'.format(**locals()))
        plt.xlabel(r'${{\rm {} }}$'.format( model.descriptor_names[0]))
        plt.ylabel(r'${{\rm {} }}$'.format(model.descriptor_names[1]))
        plt.title(title)

    plt.savefig(filename, bbox_inches='tight')

if __name__ == '__main__':
    import optparse


    parser = optparse.OptionParser()

    parser.add_option('-n', '--dont-run', dest='dontrun', action='store_true', default=False)
    parser.add_option('-p', '--plot', dest='plot', action='store_true', default=False)

    options, args = parser.parse_args()

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

            contour_plot_data(data['descriptor0'],
                              data['descriptor1'],
                              plot_data,
                              'output_{name}.pdf'.format(**locals()),
                              seed=SEED,
                              normalized=normalized,
                              title=title)
