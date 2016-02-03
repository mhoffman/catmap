#!/usr/bin/env python

import os
import catmap
import pickle
import copy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.mlab import griddata

import numpy as np

INIT_STEPS = int(1e6)
SAMPLE_STEPS = INIT_STEPS
SEED = None
TEMPERATURE = 500
DIFFUSION_FACTOR = 1e-3


def get_seed_from_path(import_path):
    import sys
    orig_path = copy.copy(sys.path)
    sys.path.insert(0, import_path)

    import kmc_settings
    seed = kmc_settings.model_name

    sys.path = orig_path

    return seed

def setup_model_species(model, species):
    if not species in model.settings.species_tags:
        tags = model.settings.species_tags.keys()
        raise UserWarning("Species '{species}' unknown, choose from {tags}.".format(**locals()))

    fortran_species = species.lower()
    kmos_species = int(eval('model.proclist.{fortran_species}'.format(**locals())))

    # extract the model dimensions
    X, Y, Z = model.lattice.system_size
    N = int(model.lattice.spuck)
    config = model._get_configuration()

    config[:] = kmos_species

    model._set_configuration(config)
    model._adjust_database()

def setup_model_probabilistic(model, data_point=0, majority=False ):
    """Make an educated initial guess for coverages by setting
    adsorbates for each with probabilites according to the
    mean field model result.

    Note that this naive way of implementing this initialization
    might lattices which are locally not realistic.

    If majority==True all sites are two the respective majority-species (a.k.a winner takes it all).
    """

    import numpy.random
    import pickle

    seed = model.settings.model_name

    data = merge_catmap_output(seed=seed)
    coverage = data['coverage_map'][data_point][1]
    catmap_sitenames = list(data['site_names'])
    if 'g' in catmap_sitenames:
        catmap_sitenames.remove('g')

    # extract the model dimensions
    X, Y, Z = model.lattice.system_size
    N = int(model.lattice.spuck)
    S = int(model.proclist.nr_of_species)
    config = model._get_configuration()

    # construct the mapping
    # from kmos species number to catmap coverage of given species
    catmap_coverages = dict(zip(data['adsorbate_names'], map(float, data['coverage_map'][data_point][1])))

    for model_sitename in model.settings.site_names:
        catmap_sitename = model_sitename.split('_')[1]
        model_sitename_index = int(eval('model.lattice.{model_sitename}'.format(**locals())))

        choices = range(S)
        choices_weights = [0.] * S

        for kmos_speciesname in model.settings.species_tags.keys():
            fortran_speciesname = kmos_speciesname.lower()

            n = int(eval('model.proclist.{fortran_speciesname}'.format(**locals())))
            catmap_adsorbatename = '{kmos_speciesname}_{catmap_sitename}'.format(**locals())
            choices_weights[n] = catmap_coverages.get(catmap_adsorbatename, 0.)

        # fill up the 'empty' (= default_species) so that all choices per site sum to 1
        choices_weights[model.proclist.default_species] =  1 - sum(choices_weights)

        choices_weights = map(float, choices_weights)

        if majority:
            argmax = np.argmax(choices_weights)
            choices_weights = [0.] * len(choices_weights)
            choices_weights[argmax] = 1.


        ## DEBUGGING
        import pprint
        print("CATMAP Coverages")
        pprint.pprint(catmap_coverages)
        #print("DEBUGGING")
        print("SITE NAME")
        print(catmap_sitename, model_sitename_index)
        print("CHOICES")
        pprint.pprint(choices)
        print("CHOICES WEIGHTS")
        pprint.pprint(choices_weights)

        # iterate over every lattice site and fill in proportionally weighted species
        for x in range(X):
            for y in range(Y):
                for z in range(Z):
                    choice = numpy.random.choice(choices, p=choices_weights,)
                    config[x, y, z, model_sitename_index - 1] = choice

    model._set_configuration(config)
    model._adjust_database()

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

def run_model(seed, init_steps, sample_steps,
              call_path=None, options=None):
    # a path we need to add to make sure kmc model import works
    if call_path is not None:
        import sys
        orig_path = copy.copy(sys.path)
        sys.path.insert(0, call_path)
    else:
        orig_path = None

    seed = seed or get_seed_from_path(call_path)

    import kmos.run

    data_filename = '{seed}_kMC.log'.format(**locals())
    lock_filename = '{seed}_kMC.lock'.format(**locals())
    done_filename = '{seed}_kMC.done'.format(**locals())

    # Let's first run the CatMAP model again with the
    # forward/back-wards rate constants
    # generating necessary rate constants
    catmap_model = catmap.ReactionModel(
        setup_file='{seed}.mkm'.format(**locals()))

    catmap_model.output_variables.append('rate_constant')
    catmap_model.output_variables.append('forward_rate_constant')
    catmap_model.output_variables.append('reverse_rate_constant')
    catmap_model.run()
    catmap_data = merge_catmap_output(seed=seed)

    # create of lock-file for currently running data-points
    # if it doesn't exist
    if not os.path.exists(lock_filename):
        with open(lock_filename, 'w'):
            pass

    # write out the data file header
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
                #print('Skipping {descriptor_string}'.format(**locals()))
                continue
        with open(lock_filename, 'a') as lockfile:
            lockfile.write('{descriptor_string}'.format(**locals()))
            lockfile.flush()

        # generate the data
        with kmos.run.KMC_Model(print_rates=False, banner=False) as kmos_model:
            print('running DATAPOINT {data_point} DESCRIPTOR {descriptor_string}'.format(**locals()))
            set_rate_constants(kmos_model, catmap_data, data_point, diffusion_factor=DIFFUSION_FACTOR)
            if options.initial_configuration == 'probabilistic':
                setup_model_probabilistic(kmos_model, data_point)
            elif options.initial_configuration.startswith('species:'):
                setup_model_species(kmos_model, options.initial_configuration.split(':')[-1])
            elif options.initial_configuration == 'empty':
                pass
            elif options.initial_configuration == 'majority':
                setup_model_probabilistic(kmos_model, data_point, majority=True)
            else:
                raise UserWarning("Directions for initial configuration '{options.initial_configuration}' can not be processed".format(**locals()))


            # DEBUGGING
            print("kmos coverages")
            kmos_model.print_accum_rate_summation()

            print(kmos_model.rate_constants)
            print('INIT STEPS {init_steps} |  SAMPLE STEPS {sample_steps}'.format(**locals()))

            kmos_model.do_steps(init_steps)
            atoms = kmos_model.get_atoms()
            data = kmos_model.get_std_sampled_data(1, sample_steps, tof_method='integ')

            with open(data_filename, 'a') as outfile:
                outfile.write(
                    '{descriptors[0]} {descriptors[1]} {data}'.format(**locals()))

        with open(done_filename, 'a') as outfile:
            outfile.write('{descriptor_string}'.format(**locals()))

        if options.single_point:
            print("User requested to run only a single-descriptor point, stopping here.")
            break

    else:
        print("Looks like all descriptor points are evaluated. Consider plotting results with 'catmap run_kmc -p'")

    # Restore old path
    if orig_path is not None:
        sys.path = orig_path

def line_plot_data(x, y, filename,
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
    fig = plt.figure()
    x = np.array(x)
    y = np.array(y)


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

    plot = plt.plot(x, y)

    plt.savefig(filename, bbox_inches='tight')
    print("Plotted {filename}".format(**locals()))

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
        zmax = 1.
        zmin = 0.
    else:
        levels = np.linspace(zmin, zmax, max(int(zmax - zmin), 2))


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

def main(options, call_path=None):
    if not options.dontrun:
        init_steps = options.equilibration_steps if options.equilibration_steps else INIT_STEPS
        sample_steps = options.sampling_steps if options.sampling_steps else SAMPLE_steps
        run_model(seed=SEED,
             init_steps=init_steps,
             sample_steps=sample_steps,
             call_path=call_path,
             options=options)

    if options.plot:
        seed = get_seed_from_path(call_path)
        data = np.recfromtxt('{seed}_kMC.log'.format(**locals()), names=True)

        catmap_model = catmap.ReactionModel(
            setup_file='{seed}.mkm'.format(**locals()))

        catmap_model.run()

        for name in data.dtype.names:
            print(name)
            if name in ['descriptor1', 'descriptor0']:
                continue
            if '_2_' in name: # we are plotting a rate
                normalized = False
                # Plot only the log base 10 of rates
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

            diff = lambda x: x[1] - x[0]

            if diff(catmap_model.descriptor_ranges[0]) == 0 and diff(catmap_model.descriptor_ranges[1]) == 0:
                raise UserWarning("Both descriptor ranges are 0, I don't know how to plot that!")
            elif diff(catmap_model.descriptor_ranges[0]) == 0 or diff(catmap_model.descriptor_ranges[1]) == 0:
                if diff(catmap_model.descriptor_ranges[0]) == 0  :
                    iv = independent_variable = 1

                else:
                    iv = independent_variable = 0


                x_data = data['descriptor{iv}'.format(**locals())]

                sort_order = np.argsort(x_data)
                x_data = x_data[sort_order]
                y_data = plot_data[sort_order]

                line_plot_data(x_data,
                               y_data,
                               'plot_kMC_{name}.pdf'.format(**locals()),
                               catmap_model=catmap_model,
                               normalized=normalized,
                               title=title,
                               )

            else:
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
                                  'plot_kMC_{name}.pdf'.format(**locals()),
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

def merge_catmap_output(seed=None, log_filename=None, pickle_filename=None):
    """The entire output from the CatMAP run is distributed
    to a *.log file and a *.pkl file. The log file contains
    all variables that require less than 100 lines, and every
    output variable that is longer is thrown to the Pickle file.

    This function merges both data sources together in one
    dictionary where the log file data takes precedence over
    the pickle data.

    :param seed: The prefix of the *.log and *.pickle files
    :type seed: str
    :param log_filename: Explicit filename of *.log file if it cannot be built from one <seed>.
    :type log_filename: str
    :param pkl_filename: Explicit filename of *.pkl file if it cannot be built from one <seed>.
    :type pkl_filename: str

    """
    import pickle

    if log_filename is None:
        log_filename = '{seed}.log'.format(**locals())
    if pickle_filename is None:
        pickle_filename =  '{seed}.pkl'.format(**locals())

    if os.path.exists(pickle_filename) and os.path.getsize(pickle_filename):
        with open(pickle_filename) as pickle_file:
            pickle_data = pickle.load(pickle_file)
    else:
        pickle_data = {}

    if os.path.exists(log_filename) and os.path.getsize(log_filename):
        log_globals, log_locals = {}, {}
        execfile(log_filename, log_globals, log_locals)
    else:
        log_locals = {}

    overlap = set(pickle_data).intersection(set(log_locals))
    pickle_data.update(log_locals)

    return pickle_data

if __name__ == '__main__':
    main()
