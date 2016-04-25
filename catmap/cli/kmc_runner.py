#!/usr/bin/env python

import os
import sys
import pprint
import catmap
import pickle
import copy
import re
import time
import traceback

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.mlab import griddata
import matplotlib.ticker

import numpy as np

INIT_STEPS = int(1e6)
SAMPLE_STEPS = INIT_STEPS
SEED = None
TEMPERATURE = 500
DIFFUSION_FACTOR = None

def contour_plot_data(x, y, z, filename,
                      n_gp=201,
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
    if os.system('which latex'):
        matplotlib.rcParams['text.usetex'] = 'false'
    else:
        matplotlib.rcParams['text.usetex'] = 'true'

    matplotlib.rcParams['lines.linewidth'] = 1.

    fig = plt.figure(figsize=fig_size)

    #x, y = np.linspace(x.min(), x.max(), m_gp), np.linspace(y.min(), y.max(), m_gp)

    xi, yi = np.linspace(x.min(), x.max(), n_gp), np.linspace(y.min(), y.max(), n_gp)
    xi, yi = np.meshgrid(xi, yi)

    #print(z)
    try:
        rbf = scipy.interpolate.Rbf(x, y, z, function='linear',)
    except Exception as e:
        print("Trouble printing {title}: {e}".format(**locals()))
        return
    zi = rbf(xi, yi)
    #zi = griddata(x, y, z, xi, yi, interp='linear')
    if normalized:
        zi = np.clip(zi, 0, 1)

    if zmin is None :
        zmin = z.min()
    if zmax is None :
        zmax = z.max()

    if normalized:
        zmax = 1
        zmin = 0
        levels = np.linspace(zmin, zmax, 11)
    else:
        #levels = np.linspace(max(-10, zmin), max(2, zmax), min(18, max(int(zmax - zmin), 2)))
        if np.allclose(zmin, zmax, 1e-2):
            levels = np.linspace(zmin * (1-1e-2 * np.sign(zmin)), zmin * (1+1e-2 * np.sign(zmin)), 6)
        else:
            print(zmin, zmax)
            attempt_list = np.abs((np.array(map(lambda x: (zmax - zmin) / x, range(1, int(zmax-zmin) + 1))) - 50))
            if len(attempt_list) > 0 :
                divider = attempt_list.argmin() + 1
            else:
                divider = 1

            levels = np.linspace(round(zmin), round(zmax), round(zmax - zmin + 1))[::divider]
            print(divider, levels)

    print("Levels {levels}".format(**locals()))



    contour_plot = plt.contourf(np.nan_to_num(zi), vmin=zmin, vmax=zmax, origin='lower',
               extent=[x.min(), x.max(), y.min(), y.max()],
               levels=levels,
               extend={False: 'both', True: 'neither'}[normalized],
               )


    ## plot the data point which we actually evaluated
    plt.scatter(x, y, c=z, s=.2)

    def cbar_fmt(x, pos):
        return '{:.3g}'.format(x)


    cbar = plt.colorbar(contour_plot, ticks=ticks, fraction=0.046, pad=0.04, format=matplotlib.ticker.FuncFormatter(cbar_fmt))

    if colorbar_label is None:
        if normalized:
            cbar.set_label(cbar_label)
        else:
            cbar.set_label(cbar_label)
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

        if hasattr(model, 'descriptor_labels'):
            print(model.descriptor_labels)
            plt.xlabel(model.descriptor_labels[0])
            plt.ylabel(model.descriptor_labels[1])
        else:
            plt.xlabel(r'${{\rm {} }}$ [{xlabel_unit}]'.format(model.descriptor_names[0], **locals()))
            plt.ylabel(r'${{\rm {} }}$ [{ylabel_unit}]'.format(model.descriptor_names[1], **locals()))

        print("Setting title {title}".format(**locals()))
        plt.title(title)

    else:
        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)


    #plt.xticks(np.arange(x.min(), x.max(), .5))
    #plt.yticks(np.arange(y.min(), y.max(), .5))
    #plt.autoscale(enable=True, axis='both', tight=True)

    #plt.xticks(np.arange(x.min(), x.max() + .5, .5))
    #plt.yticks(np.arange(y.min(), y.max() + .5, .5))

    #for axis in [fig.gca().xaxis, fig.gca().yaxis]:
        #axis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))

    plt.axis('tight', clip_on=True)
    plt.xlim((x.min(), x.max()))
    plt.ylim((y.min(), y.max()))
    #plt.axis('image')
    #plt.axis('scaled')

    try:
        plt.savefig(filename, bbox_inches='tight')
    except:
        traceback.print_stack()
        print("Had trouble saving {filename}".format(**locals()))


def get_seed_from_path(import_path):
    import sys
    orig_path = copy.copy(sys.path)
    sys.path.insert(0, import_path)

    import kmc_settings
    seed = kmc_settings.model_name

    sys.path = orig_path

    return seed

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
    if os.system('which latex'):
        matplotlib.rcParams['text.usetex'] = 'false'
    else:
        matplotlib.rcParams['text.usetex'] = 'true'
    matplotlib.rcParams['lines.linewidth'] = 1.

    if 'forward' in filename or 'reverse' in filename or 'time' in filename  or '_2_' in filename :
        # we are plotting a rate constant, do it on a log-scale
        #plot = plt.semilogy(x, y)
        plot = plt.plot(x, y)
    else:
        plot = plt.plot(x, y)

    if 'Theta' in title:
        plt.ylabel('coverage [ML]')
    elif 'rightarrow' in title:
        plt.ylabel('rate [$\\rm{cell}^{-1} \\rm{s}^{-1}$]')

    plt.xlabel(xlabel)
    plt.title(title)
    plt.xlim((x.min(), x.max()))


    plt.savefig(filename, bbox_inches='tight')
    print("Plotted {filename}".format(**locals()))

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
        data_filename = 'kMC_run_{seed}.log'.format(**locals())
        data = np.recfromtxt(data_filename, names=True)
        print("Opening {data_filename} for plotting rates and coverages".format(**locals()))

        catmap_model = catmap.ReactionModel(
            setup_file='{seed}.mkm'.format(**locals()))


        catmap_model.output_variables.append('rate_constant')
        catmap_model.output_variables.append('forward_rate_constant')
        catmap_model.output_variables.append('reverse_rate_constant')
        catmap_model.run()

        # plot in reverse to that we start with the coverages
        for name in reversed(data.dtype.names):
            if name.startswith('descriptor') or name.startswith('datapoint') or name == 'T':
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
                colorbar_label = r'$\log({\rm s}^{-1} {\rm cell}^{-1})$'
            elif 'forward' in name or 'reverse' in name: # a rate-constant
                plot_data = np.log10(data[name])
                #plot_data = data[name]
                normalized = False
                colorbar_label = r'$\log({\rm s}^{-1})$'
            elif 'kmc_steps' in name: # kmc_steps or kmc_time
                plot_data = np.log10(data[name])
                #plot_data = data[name]
                normalized = False
                colorbar_label = r'${\rm steps}$'

            else: # we are plotting a coverage
                plot_data = data[name]
                normalized = True
                colorbar_label = r'${\rm ML}$'


            # generate the plot title
            if 'forward' in name or 'reverse' in name:
                import kmc_settings
                for pname, (param, _) in kmc_settings.rate_constants.items():
                    if name == param:
                        pname = '_'.join(pname.split('_')[:-1])
                        break
                else:
                    raise UserWarning("Process corresponding to {name} not found.".format(**locals()))

                pname = r'${{\rm {}}}$'.format(pname \
                           .replace('_2_', r' \rightarrow ') \
                           .replace('_n_', ' + ') \
                           .replace('_default', '') \
                           .replace('empty', '*') \
                           .replace('_0', ''))

                if 'forward' in name:
                    k = r'$k$'
                else:
                    k = r'$k$'

                title = '{k}({pname})'.format(**locals())


            elif name == 'kmc_time':
                title = '$t_{\\rm kMC}$'
            elif name == 'kmc_steps':
                title = '# kMC steps'
            elif name == 'T':
                title = '$T$ [K]'
            elif name == 'simulated_time':
                title = '$t_{\\rm sim.}$'
            elif '_2_' in name :
                title = r'${{\rm {}}}$'.format(name \
                           .replace('_2_', r' \rightarrow ') \
                           .replace('_n_', ' + ') \
                           .replace('_default', '') \
                           .replace('empty', '*') \
                           .replace('_0', ''))
            else:
                title = r'$\Theta({{\rm {}}})$'.format(name \
                           .replace('_2_', r' \rightarrow ') \
                           .replace('_n_', ' + ') \
                           .replace('_default', '') \
                           .replace('empty', '*') \
                           #.replace('_0', ''))
                           )
                title = re.sub('_([0-9]+)', r'\1', title)
                title = re.sub('_([^_() ]+)', r'_{\1}', title)

            print('{name} => {title}'.format(**locals()))
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
                               'kMC_plot_{name}.pdf'.format(**locals()),
                               catmap_model=catmap_model,
                               normalized=normalized,
                               title=title,
                               xlabel=catmap_model.descriptor_names[iv],
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
                                  'kMC_plot_{name}.pdf'.format(**locals()),
                                  seed=seed,
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
                                  colorbar_label = colorbar_label,
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
        with open(pickle_filename, 'rb') as pickle_file:
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


def set_kmc_model_coverage_at_data_point(kmos_model, catmap_data, options, data_point):
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

def run_kmc_model_at_data_point(catmap_data, options, data_point, make_plots=False, log_target=None, bias_threshold=0.10, L=4, alpha=0.01,):
        """
            Evaluate a kMC model at a given data-point using initial bias detection and
            rescaling of fast processes for speeding up the sampling process.

            In pseudo-code the function will do the following

                initialize rate-constants for catmap_model at data_point
                model_has_fast_processes = True
                while model_has_fast_processes:
                    bring model to steady-state
                    sample steady-state
                    model_has_fast_processes <- evaluate equilibration_report

            Parameters:
                :param catmap_data: The data from a finished CatMAP evaluation as provided by <catmap.cli.kmc_runner.merge_catmap_data>
                :type catmap_data: dict
                :param options: Object as generated by optparse.parse_args() that hold further options (here: diffusion_factor, batch_size)
                :type options: object
                :param data_point: The index of the CatMAP data-point (point in descriptor space) at which to evaluate
                :type data_point: int

        """
        import kmos.run
        import kmos.utils
        import kmos.utils.progressbar

        sort_catmap_maps_inplace(catmap_data)

        t_0 = time.time()
        # generate the data
        n_current_point = data_point + 1

        descriptors = catmap_data['forward_rate_constant_map'][data_point][0]
        descriptor_string = str(descriptors)

        with kmos.run.KMC_Model(print_rates=False, banner=False) as kmos_model:
            start_batch = 0
            fast_processes = True
            fast_processes_adaption = 0
            renormalizations = {}
            while fast_processes == True :
                set_rate_constants(kmos_model, catmap_data, data_point, options=options)
                print(kmos_model.rate_constants)
                set_kmc_model_coverage_at_data_point(kmos_model, catmap_data, options, data_point)

                # DEBUGGING
                #f_init_steps = float(init_steps)
                #f_sample_steps = float(sample_steps)
                #print('INIT STEPS {f_init_steps:.2e} |  SAMPLE STEPS {f_sample_steps:.2e}'.format(**locals()))

                t_startup = time.time()
                #progress_bar = kmos.utils.progressbar.ProgressBar()
                #for i in range(100):
                    #kmos_model.do_steps(init_steps/100)
                    #progress_bar.render(i+1, 'Equilibration')
                t_equilibrate = time.time()
                #atoms = kmos_model.get_atoms()
                #print(kmos_model.rate_constants)
                #try:
                    #data = kmos_model.get_std_sampled_data(options.coverage_samples, sample_steps, tof_method='integ', reset_time_overrun=True)
                import kmos.run.steady_state
                data_dict, data = kmos.run.steady_state.sample_steady_state(kmos_model,
                       show_progress=True,
                       make_plots=make_plots,
                       start_batch=start_batch,
                       batch_size=options.batch_size,
                       bias_threshold=bias_threshold,
                       L=L,
                       alpha=alpha,
                       output='both',
                       seed='EWMA_{data_point:04d}'.format(**locals()),
                       )

                if len(data) >= 0:
                    start_batch = int(float(data.split()[-1])) / options.batch_size
                else:
                    start_batch = 0

                #except Exception as e:
                    #print("Warning: Encountered error in sampling: {e}. Make sure this is correct.".format(**locals()))
                    #continue

                t_sample = time.time()


                import kmos.run.steady_state
                with open("procstat_{:04d}.dat".format(data_point), 'a') as procstat_file:
                    if log_target is None:
                        outfile = procstat_file
                    else:
                        outfile = log_target

                    outfile.write("========\n\n")
                    outfile.write("datapoint [descriptor string]\n\n")
                    outfile.write("{data_point} {descriptor_string}\n\n".format(**locals()))
                    outfile.write("\n\nProcstat (number of executed processes)\n")
                    outfile.write(kmos_model.print_procstat(to_stdout=False))
                    outfile.write("\n\nCoverages\n")
                    outfile.write(kmos_model.print_coverages(to_stdout=False))
                    outfile.write('\n\nRate Constants\n')
                    outfile.write(kmos_model.rate_constants())
                    outfile.write("\n\nEquilibration Report\n")
                    equilibration_report, equilibration_data = kmos.run.steady_state.report_equilibration(kmos_model)
                    outfile.write("\nSampled rates and coverages\n")
                    outfile.write(pprint.pformat(data_dict))
                    outfile.write("\n\nRelative rate differences between reversing processes\n")
                    outfile.write(equilibration_report)

                    fast_processes = False
                    outfile.write("\nEvaluating equilibration report\n")
                    for ratio, pn1, pn2 in equilibration_data :
                        outfile.write("{pn1} <=> {pn2} : {ratio}\n".format(**locals()))
                        EQUIB_THRESHOLD = 1e-2
                        if abs(ratio) < EQUIB_THRESHOLD:
                            fast_processes = True
                            for pn in [pn1, pn2]:
                                old_rc = kmos_model.rate_constants.by_name(pn)
                                rc_tuple = kmos_model.settings.rate_constants[pn]
                                #rc_tuple = (rc_tuple[0] + '*.5', rc_tuple[1])
                                rescale_factor = max(EQUIB_THRESHOLD, (abs(ratio) / EQUIB_THRESHOLD))
                                rescale_multiplication = '*%.2e' % rescale_factor
                                rescale_division = '/%.2e' % rescale_factor
                                renormalizations[pn] = renormalizations.get(pn, '1.') + '/{:.2e}'.format(rescale_factor)
                                rc_tuple = (rc_tuple[0] + rescale_multiplication, rc_tuple[1])
                                kmos_model.settings.rate_constants[pn] = rc_tuple
                                new_rc = kmos_model.rate_constants.by_name(pn)
                                #kmos_model.rate_constants.set(pn, new_rc)
                                outfile.write("Found a fast equilibrated process {ratio}: {pn}, reduced rate constant from {old_rc:.2e} to {new_rc:.2e}\n".format(**locals()))
                    print(kmos_model.proclist)
                    for proc in range(kmos_model.proclist.nr_of_proc.max()):
                        kmos_model.base.set_procstat(proc + 1, 0)
                    kmos_model.base.set_kmc_step(0)

                    if fast_processes or fast_processes_adaption == 0:
                        outstring = '{data_point:9d} {descriptors[0]: .5e} {descriptors[1]: .5e} {data}'.format(**locals())

                    # Check if every process has been touch in this round
                    if all([rate < EQUIB_THRESHOLD for (rate, _, _) in equilibration_data]):
                        fast_processes = False
                        outfile.write("Found all processes, to be equilibrated. So further adjustments will not help. Quit.")

                    ## Check if every process has been touched at least once
                    #len_renormalizations = len(renormalizations)
                    #len_equilibration_data = len(equilibration_data)
                    #print("len(renormalizations) = {len_renormalizations};  len(equilibration_data) = {len_equilibration_data}".format(**locals()))
                    #if len_renormalizations >= len_equilibration_data :
                        #fast_processes = False
                        #outfile.write("Every process has been found to be equilibrated at least once, exiting.")
                    # Does not work, let's abandon for now.

                    # Check if ill-defined values in result, if so run again
                    if not(data_dict):
                        outfile.write("Found ill-defined processes, will keep running\n")
                        fast_processes = True

                    # Confirm result
                    if not fast_processes :
                        outfile.write("Found no fast processes after {fast_processes_adaption} adaptions, so this result will count\n".format(**locals()))
                        outfile.write("\nRenormalizations\n")
                        outfile.write(pprint.pformat(renormalizations))

                fast_processes_adaption += 1
            return outstring

def sort_catmap_maps_inplace(data):
    for key in data:
        if key.endswith('_map'):
            data[key] = sorted(data[key], key=lambda x: x[0])


def run_model(seed, init_steps, sample_steps,
              call_path=None, options=None):
    import kmos.run
    # a path we need to add to make sure kmc model import works
    if call_path is not None:
        import sys
        orig_path = copy.copy(sys.path)
        sys.path.insert(0, call_path)
    else:
        orig_path = None

    seed = seed or get_seed_from_path(call_path)

    data_filename = 'kMC_run_{seed}.log'.format(**locals())
    lock_filename = 'kMC_run_{seed}.lock'.format(**locals())
    done_filename = 'kMC_run_{seed}.done'.format(**locals())

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
    sort_catmap_maps_inplace(catmap_data)


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
                    'datapoint descriptor0 descriptor1 {data_header}'.format(**locals()))

    total_points = len(catmap_data['forward_rate_constant_map'])

    for data_point in range(total_points):
        descriptors = catmap_data['forward_rate_constant_map'][data_point][0]

        # multi IO mechanism: keep lockfile with one line descriptors of datapoint
        # if datapoint is already in there, skip to next datapoint
        descriptor_string = str(descriptors) + '\n'

        print('\n\nrunning DATAPOINT {data_point}/{total_points} DESCRIPTOR {catmap_model.descriptor_names} = {descriptor_string}'.format(**locals()))

        with open(lock_filename, 'r') as lockfile:
            if descriptor_string in lockfile.readlines():
                #print('Skipping {descriptor_string}'.format(**locals()))
                continue

        with open(lock_filename, 'a') as lockfile:
            lockfile.write('{descriptor_string}'.format(**locals()))
            lockfile.flush()

        outstring = run_kmc_model_at_data_point(catmap_data, options, data_point)

        if outstring is not None:
            with open(data_filename, 'a') as outfile:
                outfile.write(
                    '{outstring}'.format(**locals()))

            with open(done_filename, 'a') as outfile:
                outfile.write('{descriptor_string}'.format(**locals()))

        if options.single_point:
            print("User requested to run only a single-descriptor point, stopping here.")
            break

    else:
        print("\nLooks like all descriptor points are evaluated. Consider plotting results with 'catmap run_kmc -p'")

    # Restore old path
    if orig_path is not None:
        sys.path = orig_path

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
        print("Initial CatMAP Coverages")
        pprint.pprint(catmap_coverages)

        # iterate over every lattice site and fill in proportionally weighted species
        for x in range(X):
            for y in range(Y):
                for z in range(Z):
                    choice = numpy.random.choice(choices, p=choices_weights,)
                    config[x, y, z, model_sitename_index - 1] = choice

    model._set_configuration(config)
    model._adjust_database()

    print("Initial kmos Coverages")
    model.print_coverages()

def set_rate_constants(kmos_model, catmap_data, data_point, options=None):
    """
        A rate constants of a kmos model from a corresponding CatMAP model for
        a given data point (i.e. a tuple of of reactivity descriptors).

        Diffusion factors has a special role of mimicking the effect of diffusion
        in a kMC model as described typically by a mean field model. If the diffusion-factor
        is left at its default value (None) is it simply faithfully set to the corresponding
        CatMAP value which will usually be the normal prefactor (kT/h). If instead it is
        set to a finite float it will be set to diffusion_factor * max_rate_constant
        where max_rate_constant is the fastest rate constant of all non-diffusion
        reaction steps.

    """

    # set rate constant of kMC according to current descriptor tuple
    max_rate_constant = float('-inf')

    for i in range(len(catmap_data['forward_rate_constant_map'][data_point][1])):
        forward_rate_constant = float(catmap_data['forward_rate_constant_map'][data_point][1][i])
        reverse_rate_constant = float(catmap_data['reverse_rate_constant_map'][data_point][1][i])

        print('{i} Forward {forward_rate_constant:.3e} Reverse {reverse_rate_constant:.3e}'.format(**locals()))

        if hasattr(kmos_model.parameters, 'forward_{i}'.format(**locals())):
            max_rate_constant = max(max_rate_constant, forward_rate_constant)
            setattr(kmos_model.parameters, 'forward_{i}'.format(
                **locals()), forward_rate_constant)
        if hasattr(kmos_model.parameters, 'reverse_{i}'.format(**locals())):
            max_rate_constant = max(max_rate_constant, reverse_rate_constant)
            setattr(kmos_model.parameters, 'reverse_{i}'.format(
                **locals()), reverse_rate_constant)

    # set the rate-constant of diffusion rate-constants
    for i in range(len(catmap_data['forward_rate_constant_map'][data_point][1])):
        forward_rate_constant = catmap_data['forward_rate_constant_map'][data_point][1][i] if getattr(options, 'diffusion_factor', None) is None else max_rate_constant * options.diffusion_factor
        reverse_rate_constant = catmap_data['reverse_rate_constant_map'][data_point][1][i] if getattr(options, 'diffusion_factor', None) is None else max_rate_constant * options.diffusion_factor

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



def find_pairs(project):
    """Find pairs of elementary processes that are reverse processes with respect
    to each others from a kmos.types.Project

    """
    pairs = []
    for p1 in sorted(project.process_list):
        for p2 in sorted(project.process_list):
            if p1.condition_list == p2.action_list and p2.condition_list == p1.action_list:
                if not (p1, p2) in pairs and not (p2, p1) in pairs:
                    pairs.append((p1, p2))
    return pairs

def report_equilibration(model):
    """Iterate over pairs of reverse proceses and print
        rate1 * rho1 / rate2 * rho2

      for each.
    """
    import kmos.types
    import StringIO

    project = kmos.types.Project()
    project.import_ini_file(StringIO.StringIO(model.settings.xml))
    pairs = find_pairs(project)

    atoms = model.get_atoms(geometry=False)

    procstat = dict(zip(sorted(model.settings.rate_constants), atoms.procstat))
    rate_constants = dict(zip(sorted(model.settings.rate_constants), (model.base.get_rate(i+1) for i in range(len(procstat)))))

    integ_rates = dict(zip(sorted(model.settings.rate_constants), atoms.integ_rates))


    report = ''
    for pair in pairs:
        pn1, pn2 = pair[0].name, pair[1].name
        left = integ_rates[pn1] * rate_constants[pn2]
        right = integ_rates[pn2] * rate_constants[pn1]
        ratio = left/right
        eq_score = 4 * left * right / (left + right)**2
        report += ('{pn1} : {pn2} => {left:.2e}/{right:.2e} = {ratio:.2e}, eq_score = {eq_score:.7f}\n'.format(**locals()))
    return report

if __name__ == '__main__':
    main()
