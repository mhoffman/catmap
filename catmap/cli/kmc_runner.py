#!/usr/bin/env python

import os
import pprint
import catmap
import random
import copy
import re
import math
import datetime
import bisect

import catmap.cli.kmc_translation

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.ticker
import matplotlib.colors


INIT_STEPS = int(1e6)
SAMPLE_STEPS = INIT_STEPS
SEED = None
TEMPERATURE = 500
DIFFUSION_FACTOR = None
ROUND_DIGITS = 5

def takeClosest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect.bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return after
    else:
       return before

def process_name_to_latex(pname, arrow=r' \rightarrow '):
    """Translate the name of an elementary processes as used inside the kMC model
    into its LaTeX version.

    """
    pname = r'${{\rm {}}}$'.format(pname
                                   .replace('_2_', arrow)
                                   .replace('_n_', ' + ')
                                   .replace('_default', '')
                                   .replace('empty', '*')
                                   .replace('_0', ''))
    return pname


def get_canonical_process_names(data, empty_species=catmap.cli.kmc_translation.EMPTY_SPECIES):
    """Translate the elementary reactions (elementary_rxns) used in a CatMAP model
    into a form that can be used as a variable name.

    """
    process_names = []
    site_names = data.site_positions.keys()
    for elementary_rxn in data.elementary_rxns:
        step = {}
        surface_intermediates = {}
        # N.B: The general form of an elementary reaction in CatMAP is
        #      A <-> B -> C
        #      where A is the initial state,
        #            C is the final state,
        #            B is the transition state
        #            B may be skipped
        if len(elementary_rxn) == 2:
            step['A'], step['C'] = elementary_rxn
            step['B'] = None
        elif len(elementary_rxn) == 3:
            step['A'], step['B'], step['C'] = elementary_rxn
        for reversible, (X, Y) in [[True, ('A', 'C')], ]:
            if step[X] and step[Y]:
                # add reversible step between A and B
                surface_intermediates[X] = []
                surface_intermediates[Y] = []

                for x in [X, Y]:
                    surface_intermediates[x] = get_canonical_intermediates(step[x], site_names=site_names, empty_species=catmap.cli.kmc_translation.EMPTY_SPECIES)

                forward_name_root, reverse_name_root = catmap.cli.kmc_translation.surface_intermediates_to_process_names(surface_intermediates[X], surface_intermediates[Y])

                process_names.append((forward_name_root, reverse_name_root))
    return process_names


def get_canonical_intermediates(step, site_names=None, empty_species=catmap.cli.kmc_translation.EMPTY_SPECIES):
    """Helper function for  get_canonical_process_names.

    """

    surface_intermediates = []
    # print('step = {step}'.format(**locals()))
    for intermediate in step:
        if '_' in intermediate:
            species, site = intermediate.split('_')
            if any([s.startswith('{site}'.format(**locals())) for s in site_names]):
                surface_intermediates.append([species, site])
        elif any([s.startswith('{intermediate}'.format(**locals())) for s in site_names]):
            surface_intermediates.append(
                [catmap.cli.kmc_translation.EMPTY_SPECIES, intermediate])
        else:
            pass
    return surface_intermediates


class MidpointNormalize(matplotlib.colors.Normalize):
    """Credit goes to

    http://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        if self.midpoint - self.vmin > self.vmax - self.midpoint:
            y = [0, .5, .5 + .5 * float(self.vmax - self.midpoint) / float(self.midpoint - self.vmin)]
        else:
            y = [.5 - .5 * float(self.midpoint - self.vmin) / float(self.vmax - self.midpoint), .5, 1.]
        print("Mid point normalize vmin = {self.vmin}, vmax = {self.vmax}, midpoint = {self.midpoint}, ({x}, {y})".format(**locals()))
        return np.ma.masked_array(np.interp(value, x, y))


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
                      colorbar_label=None,
                      show_evaluated_points=True,
                      cmap=None,):
    """
    Create contour plot from one column (z) over (x, y) as generated by
    kMC runner. Additional arguments can fine-tune appearance (xlabel, ylabel),
    Pass in catmap_model to generate descriptor points specific surfaces.

    """
    import scipy.interpolate
    fig = plt.figure()

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    centered_norm = False
    if cmap is None:
        cmap = plt.cm.jet
    elif type(cmap) is str:
        if cmap in ['bwr', 'seismic']:
            centered_norm = True
        cmap = eval('plt.cm.{cmap}'.format(**locals()))

    # with golden ration and the whole shebang ...
    # settings size and font for revtex stylesheet
    fig_width_pt = 2 * 246.0  # Get this from LaTeX using \showthe\columnwidth
    # fig_width_pt *= 300./72 # convert to 300 dpi
    inches_per_pt = 1.0 / 72.27               # Convert pt to inches
    # inches_per_pt = 1.0/300               # Convert pt to inches
    golden_mean = (np.sqrt(5) - 1.0) / 2.0         # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = fig_width * golden_mean       # height in inches
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

    # x, y = np.linspace(x.min(), x.max(), m_gp), np.linspace(y.min(), y.max(), m_gp)

    xi, yi = np.linspace(x.min(), x.max(), n_gp), np.linspace(y.min(), y.max(), n_gp)
    xi, yi = np.meshgrid(xi, yi)

    # print(z)
    try:
        # rbf = scipy.interpolate.Rbf(x, y, z, function='linear', smooth=0.2)
        rbf = scipy.interpolate.Rbf(x, y, z, function='linear', smooth=0.01)
    except Exception as e:
        import traceback
        print("Trouble printing {title} {colorbar_label}: {e}".format(**locals()))
        traceback.print_stack()
        return
    zi = rbf(xi, yi)
    # zi = griddata(x, y, z, xi, yi, interp='linear')
    if normalized:
        zi = np.clip(zi, 0, 1)

    if zmin is None:
        zmin = z.min()
    if zmax is None:
        zmax = z.max()

    if normalized:
        zmax = 1
        zmin = 0
        levels = np.linspace(zmin, zmax, 11)
    else:
        # levels = np.linspace(max(-10, zmin), max(2, zmax), min(18, max(int(zmax - zmin), 2)))
        if np.allclose(zmin, zmax, 1e-2):
            levels = np.linspace(zmin * (1 - 1e-2 * np.sign(zmin)), zmin * (1 + 1e-2 * np.sign(zmin)), 6)
        else:
            print(zmin, zmax)
            if np.abs(zmax - zmin) > 1e6:
                print("Linear plots goes over range > 1e6, refusing to plot.")
                return

            attempt_list = np.abs((np.array(map(lambda x: (zmax - zmin) / x, range(1, int(zmax - zmin) + 1))) - 50))
            if len(attempt_list) > 0:
                divider = attempt_list.argmin() + 1
            else:
                divider = 1

            levels = np.linspace(round(zmin), round(zmax), max(round(zmax - zmin + 1), 2))[::divider]
            print(divider, levels)

    print("Levels {levels}".format(**locals()))
    if levels[0] == levels[-1]:
        print("Somehow levels, are not increasing, nothing to plot.")
        return

    if centered_norm:
        norm = MidpointNormalize(midpoint=0.)
        print(norm)
        contour_plot = plt.contourf(np.nan_to_num(zi), vmin=zmin, vmax=zmax, origin='lower',
                                    extent=[x.min(), x.max(), y.min(), y.max()],
                                    levels=levels,
                                    extend={False: 'both', True: 'neither'}[normalized],
                                    cmap=cmap,
                                    norm=norm,
                                    )
    else:
        contour_plot = plt.contourf(np.nan_to_num(zi), vmin=zmin, vmax=zmax, origin='lower',
                                    extent=[x.min(), x.max(), y.min(), y.max()],
                                    levels=levels,
                                    extend={False: 'both', True: 'neither'}[normalized],
                                    cmap=cmap,
                                    )

    if show_evaluated_points:
        # plot the data point which we actually evaluated
        plt.scatter(x, y, c=z, s=.2)

    def cbar_fmt(x, pos):
        return '{:.3g}'.format(x)

    cbar = plt.colorbar(contour_plot, ticks=ticks, fraction=0.046, pad=0.04, format=matplotlib.ticker.FuncFormatter(cbar_fmt))

    if colorbar_label is not None:
        cbar.set_label(colorbar_label)

    if catmap_model is not None:
        for substrate, (xi, yi) in catmap_model.descriptor_dict.items():
            _z = rbf(xi, yi)
            print(substrate, _z)
            plt.scatter(xi, yi, c=_z, s=35, vmin=_z.min(), vmax=_z.max(), cmap=cbar.cmap)
            plt.annotate(substrate, xy=(xi + .05, yi + .05), size='small',
                         bbox={'facecolor': 'white',
                               'alpha': 0.5,
                               'ec': 'white',
                               'pad': 1,
                               'lw': 0})
        if seed is None:
            seed = catmap_model.model_name

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

    plt.axis('tight', clip_on=True)
    plt.xlim((x.min(), x.max()))
    plt.ylim((y.min(), y.max()))

    if title == 'data point':
        for _x, _y, _label in zip(x, y, z):
            plt.annotate(_label, xy=(_x, _y), size='small',)

    try:
        plt.savefig(filename, bbox_inches='tight')
    except:
        import traceback
        traceback.print_stack()
        print("Had trouble saving {filename}".format(**locals()))


def get_seed_from_path(import_path):
    """Return the seed (=internal name) of kMC model that is stored in `import_path` by inspecting the `kmc_settings` module.

    """
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
    fig_width_pt = 2 * 246.0  # Get this from LaTeX using \showthe\columnwidth
    # fig_width_pt *= 300./72 # convert to 300 dpi
    inches_per_pt = 1.0 / 72.27               # Convert pt to inches
    # inches_per_pt = 1.0/300               # Convert pt to inches
    golden_mean = (np.sqrt(5) - 1.0) / 2.0         # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches
    fig_height = fig_width * golden_mean       # height in inches
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

    if 'forward' in filename or 'reverse' in filename or 'time' in filename or '_2_' in filename:
        # we are plotting a rate constant, do it on a log-scale
        # plot = plt.semilogy(x, y)
        plot = plt.plot(x, y)
        plot = plt.plot(x, y, 'o')
    else:
        plot = plt.plot(x, y)
        plot = plt.plot(x, y, 'o')

    if 'Theta' in title:
        plt.ylabel('coverage [ML]')
    elif 'rightarrow' in title:
        plt.ylabel('rate [$\\rm{cell}^{-1} \\rm{s}^{-1}$]')

    plt.xlabel(xlabel)
    plt.title(title)
    plt.xlim((x.min(), x.max()))

    plt.savefig(filename, bbox_inches='tight')
    print("Plotted {filename}".format(**locals()))


def plot_mft_coverages(catmap_model, kmos_data, seed=None):
    for i, adsorbate_name in enumerate(catmap_model.adsorbate_names):
        xs, ys, zs = [], [], []
        zs_log = []
        zraw = []
        zMFT_dict = {}

        mft_signal = False

        # plot MFT rates
        for (x, y), rates in catmap_model.coverage_map:
            xs.append(x)
            ys.append(y)
            zraw.append(float(rates[i]))
            zs.append((float(rates[i])))
            zs_log.append(np.log10(float(rates[i])))
            zMFT_dict.setdefault(round(x, ROUND_DIGITS), {})[round(y, ROUND_DIGITS)] = float(rates[i])
            if float(rates[i]) != 0.:
                mft_signal = True
        zs_MFT_max = max(zraw)
        pname = process_name_to_latex(adsorbate_name, arrow=r' \rightleftharpoons ')
        title = '{pname}'.format(**locals())
        colorbar_label = '{{\\rm ML}}'
        contour_plot_data(xs, ys, zs,
                          'kMC_plot_MFT_coverage_{i}.pdf'.format(**locals()),
                          colorbar_label=colorbar_label,
                          title=title,
                          normalized=True,
                          catmap_model=catmap_model,
                          )

        contour_plot_data(xs, ys, zs_log,
                          'kMC_plot_MFT_coverage_log_{i}.pdf'.format(**locals()),
                          colorbar_label=colorbar_label,
                          title=title,
                          normalized=False,
                          catmap_model=catmap_model,
                          )


def plot_mft_kmc_differences(catmap_model, kmos_data, seed=None):
    process_names = (get_canonical_process_names(catmap_model))
    for i, elementary_rxn in enumerate(catmap_model.elementary_rxns):
        print(elementary_rxn)
        xs, ys, zs = [], [], []
        zraw = []
        zMFT_dict = {}

        mft_signal = False

        # plot MFT rates
        for (x, y), rates in catmap_model.rate_map:
            xs.append(x)
            ys.append(y)
            zraw.append(float(rates[i]))
            zs.append(np.log10(float(rates[i])))
            zMFT_dict.setdefault(round(x, ROUND_DIGITS), {})[round(y, ROUND_DIGITS)] = float(rates[i])
            if float(rates[i]) != 0.:
                mft_signal = True
        zs_MFT_max = max(zraw)
        pname = process_name_to_latex(process_names[i][0], arrow=r' \rightleftharpoons ')
        title = '(R_{{\\rm MFT}})$ ({pname})'.format(**locals())
        colorbar_label = '$\\log_{{10}}({{\\rm TOF}})$ [s$^{{-1}}$ cell$^{{-1}}$]'.format(**locals())
        contour_plot_data(xs, ys, zs,
                          'kMC_plot_MFT_rate_{i}.pdf'.format(**locals()),
                          colorbar_label=colorbar_label,
                          title=title,
                          catmap_model=catmap_model,
                          )

        print("PLOTTING DELTA")
        xs, ys, zs = [], [], []
        zs_lin = []
        zs_forward = []
        zs_reverse = []
        zs_MFT = []
        x_kMC, y_kMC, z_kMC = [], [], []
        zs_kMC = np.array(kmos_data[process_names[i][0]]) - np.array(kmos_data[process_names[i][1]])

        zs_kMC_max = max(zs_kMC)
        z_ratio = zs_MFT_max / zs_kMC_max

        #print('zMFT dict')
        #print(zMFT_dict)
        #print(zs_MFT_max, zs_kMC_max, z_ratio)

        print(process_names[i])
        for x, y, z0, z1 in zip(kmos_data['descriptor0'], kmos_data['descriptor1'], kmos_data[process_names[i][0]], kmos_data[process_names[i][1]]):

            x, y = round(x, ROUND_DIGITS), round(y, ROUND_DIGITS)
            if np.isfinite(np.log10(np.abs(z0 - z1))):
                x_kMC.append(x)
                y_kMC.append(y)
                z_kMC.append(np.log10(np.abs(z0 - z1)))
            elif np.isfinite(np.log10(- (z0 - z1))):
                x_kMC.append(x)
                y_kMC.append(y)
                z_kMC.append(np.log10(-(z0 - z1)))
            else:
                process_names_i_ = process_names[i]
                print("Check {process_names_i_} {i} {z0} {z1}".format(**locals()))

            x_closest = takeClosest(sorted(zMFT_dict.keys()), x)
            y_closest = takeClosest(sorted(zMFT_dict[x_closest].keys()), y)
            #print('    {x} => {x_closest}'.format(**locals()))
            #print('    {y} => {y_closest}'.format(**locals()))
            zMFT_xy = zMFT_dict[x_closest][y_closest]
            _correction = 0.
            _correction = 1.
            if mft_signal:
                ztest = - np.log10((z0 - z1) / zMFT_xy)
                if (z0 - z1) > zMFT_xy :
                    ztest_lin = - (z0 - z1) / zMFT_xy + _correction
                else:
                    ztest_lin =  zMFT_xy / (z0 - z1) - _correction
            else:
                ztest = - np.log10((z0 - z1))
                ztest_lin = - (z0 - z1)

            if not math.isnan(ztest):
                xs.append(x)
                ys.append(y)
                zs.append(ztest)
                zs_lin.append(ztest_lin)

            zs_forward.append(z0)
            zs_reverse.append(z1)
            zs_MFT.append(zMFT_xy)
        zs = np.array(zs)
        zs[np.logical_not(np.isfinite(zs))] = 0.

        zs_lin = np.array(zs_lin)
        zs_lin[np.logical_not(np.isfinite(zs_lin))] = 0.

        print(zs)
        print(zs_lin)

        pname = process_name_to_latex(process_names[i][0], arrow=r' \rightleftharpoons ')
        colorbar_label = '$\\log_{{10}}(R_{{\\rm MFT}}/R_{{\\rm kMC}})$ ({pname})'.format(**locals())
        colorbar_label_lin = '$(R_{{\\rm MFT}}/R_{{\\rm kMC}})$ ({pname})'.format(**locals())
        print("Colorbar label {colorbar_label}".format(**locals()))

        def diff(x):
            return x[1] - x[0]

        if diff(catmap_model.descriptor_ranges[0]) == 0 or diff(catmap_model.descriptor_ranges[1]) == 0:
            if diff(catmap_model.descriptor_ranges[0]) == 0:
                iv = independent_variable = 1
                x_data = ys
            else:
                iv = independent_variable = 0
                x_data = xs
            y_data = zs
            y_lin_data = zs_lin

            print('DELTA {i}'.format(**locals()))
            print(x_data)
            print(y_data)
            print(y_lin_data)
            print('Zs forward')
            print(zs_forward)
            print('Zs reverse')
            print(zs_reverse)
            print('Zs MFT')
            print(zs_MFT)

            line_plot_data(x_data,
                           y_data,
                           'kMC_plot_delta_{i}.pdf'.format(**locals()),
                           catmap_model=catmap_model,
                           normalized=False,
                           title='log net rate delta {i}'.format(**locals()),
                           xlabel=catmap_model.descriptor_names[iv],
                           )

            line_plot_data(x_data,
                           y_lin_data,
                           'kMC_plot_delta_lin_{i}.pdf'.format(**locals()),
                           catmap_model=catmap_model,
                           normalized=False,
                           title='lin net rate delta {i}'.format(**locals()),
                           xlabel=catmap_model.descriptor_names[iv],
                           )
        else:
            try:
                contour_plot_data(xs, ys, zs,
                                  'kMC_plot_delta_kMC_MFT_{i}.pdf'.format(**locals()),
                                  colorbar_label=colorbar_label, catmap_model=catmap_model, seed=seed, cmap='seismic')
                contour_plot_data(xs, ys, zs_lin,
                                  'kMC_plot_delta_kMC_MFT_lin_{i}.pdf'.format(**locals()),
                                  colorbar_label=colorbar_label_lin, catmap_model=catmap_model, seed=seed, cmap='seismic')
                contour_plot_data(xs, ys, np.clip(zs_lin, -10, 10),
                                  'kMC_plot_delta_kMC_MFT_lin_clipped_{i}.pdf'.format(**locals()),
                                  colorbar_label=colorbar_label_lin, catmap_model=catmap_model, seed=seed, cmap='seismic')
            except:
                process_name = process_names[i]
                print("Trouble plotting delta {process_name}".format(**locals()))
                print("PLOTTED DELTA KMC_MFT {process_names}".format(**locals()))

        pname = process_name_to_latex(process_names[i][0], arrow=r' \rightleftharpoons ')
        colorbar_label = '$\\log_{{10}}(R_{{\\rm kMC}})$ ({pname})'.format(**locals())
        print(z_kMC)
        try:
            contour_plot_data(x_kMC, y_kMC, z_kMC,
                              'kMC_plot_kMC_{i}.pdf'.format(**locals()),
                              colorbar_label=colorbar_label,
                              catmap_model=catmap_model,
                              seed=seed,
                              )
        except:
            process_name = process_names[i]
            print("Trouble plotting kMC {process_name}".format(**locals()))


def main(options, call_path=None):
    """Main function of kMC runner, that can either run the kMC model or plots its results

    """
    if not options.dontrun:
        # init_steps = options.batch_size if options.batch_size else INIT_STEPS
        # sample_steps = options.batch_size if options.batch_size else SAMPLE_steps
        run_model(seed=SEED,
                  # init_steps=init_steps,
                  # sample_steps=sample_steps,
                  call_path=call_path,
                  options=options)

    if options.plot:
        seed = get_seed_from_path(call_path)
        data_filename = 'kMC_run_{seed}.log'.format(**locals())
        data = np.recfromtxt(data_filename, names=True)
        data, inverse_indices = np.unique(data, return_inverse=True)  # remove duplicate rows
        if len(inverse_indices) > len(data):
            print("Warning: Data file {data_filename} contained duplicate lines, I ignore them for now but you want to check it out.".format(**locals()))
        print("Opening {data_filename} for plotting rates and coverages".format(**locals()))

        catmap_model = catmap.ReactionModel(
            setup_file='{seed}.mkm'.format(**locals()))

        catmap_model.output_variables.append('rate_constant')
        catmap_model.output_variables.append('forward_rate_constant')
        catmap_model.output_variables.append('reverse_rate_constant')
        catmap_model.run()

        # possible move further down
        plot_mft_kmc_differences(catmap_model, data, seed=seed)
        plot_mft_coverages(catmap_model, data, seed=seed)

        # plot in reverse to that we start with the coverages
        for name in reversed(data.dtype.names):
            print("COLUMN {name}".format(**locals()))
            if name.startswith('descriptor') or name == 'T':
                continue
            elif name.startswith('datapoint'):
                plot_data = data[name]
            elif '_2_' in name:  # we are plotting a rate
                normalized = False
                # Plot only the log base 10 of rates
                plot_data = np.log10(data[name])
                if not np.isfinite(plot_data).any():
                    plot_data[:] = 0.
                else:
                    minimum = np.nanmin(plot_data[np.isfinite(plot_data)])
                    print('MINIMUM {minimum}'.format(**locals()))
                    plot_data[np.logical_or(np.isnan(plot_data),
                                            np.isinf(plot_data))] = minimum
                colorbar_label = r'$\log({\rm s}^{-1} {\rm cell}^{-1})$'
            elif 'forward' in name or 'reverse' in name:  # a rate-constant
                plot_data = np.log10(data[name])
                # plot_data = data[name]
                normalized = False
                colorbar_label = r'$\log({\rm s}^{-1})$'
            elif 'kmc_steps' in name:  # kmc_steps or kmc_time
                plot_data = np.log10(data[name])
                normalized = False
                colorbar_label = r'${\rm steps}$'

            else:  # we are plotting a coverage
                plot_data = data[name]
                normalized = True
                colorbar_label = r'${\rm ML}$'

            # generate the plot title
            if 'forward' in name or 'reverse' in name:
                import kmc_settings
                for pname, (param, _) in kmc_settings.rate_constants.items():
                    if name == param.split('/')[0]:
                        pname = '_'.join(pname.split('_')[:-1])
                        break
                else:
                    raise UserWarning("Process corresponding to {name} not found.".format(**locals()))

                pname = process_name_to_latex(pname)

                if 'forward' in name:
                    k = r'$k$'
                else:
                    k = r'$k$'

                title = '{k}({pname})'.format(**locals())

            elif name == 'datapoint':
                title = 'data point'

            elif name == 'kmc_time':
                title = '$t_{\\rm kMC}$'
            elif name == 'kmc_steps':
                title = '# kMC steps'
            elif name == 'T':
                title = '$T$ [K]'
            elif name == 'simulated_time':
                title = '$t_{\\rm sim.}$'
            elif '_2_' in name:
                title = r'${{\rm {}}}$'.format(name
                                               .replace('_2_', r' \rightarrow ')
                                               .replace('_n_', ' + ')
                                               .replace('_default', '')
                                               .replace('empty', '*')
                                               .replace('_0', ''))
            else:
                title = r'$\Theta({{\rm {}}})$' \
                        .format(name
                                .replace('_2_', r' \rightarrow ')
                                .replace('_n_', ' + ')
                                .replace('_default', '')
                                .replace('empty', '*')
                                # .replace('_0', ''))
                                )
                title = re.sub('_([0-9]+)', r'\1', title)
                title = re.sub('_([^_() ]+)', r'_{\1}', title)

            print('{name} => {title}'.format(**locals()))

            def diff(x):
                return x[1] - x[0]

            if diff(catmap_model.descriptor_ranges[0]) == 0 and diff(catmap_model.descriptor_ranges[1]) == 0:
                raise UserWarning("Both descriptor ranges are 0, I don't know how to plot that!")
            elif diff(catmap_model.descriptor_ranges[0]) == 0 or diff(catmap_model.descriptor_ranges[1]) == 0:
                if diff(catmap_model.descriptor_ranges[0]) == 0:
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
                    ticks = range(zmin, zmax + 1, 6)

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
                                  colorbar_label=colorbar_label,
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
        pickle_filename = '{seed}.pkl'.format(**locals())

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


def setup_edged_model_at_datapoint(model, data_point, reset_configuration=False):
    """Convenience function for carrying out all steps necessary to setup a
    compiled kMC model corresponding to a CatMAP model.

    """
    set_mft_parameters_at_data_point(model, data_point)
    set_rate_constants(model, data_point)
    if reset_configuration:
        setup_model_probabilistic(model, data_point)
    if hasattr(model.proclist, 'mft_'):
        setup_mft_edges_2d(model)
    if hasattr(model.parameters, 'N_sites'):
        model.parameters.N_sites = model.lattice.system_size.prod()

def setup_model_species(model, species):
    if species not in model.settings.species_tags:
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
    """Dispatch function for initializing kMC model in one of several ways
    more or less corresponding to the (mean-field) coverages of a CatMAP model.

    """
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

    if hasattr(kmos_model.proclist, 'mft_'):
        setup_mft_edges_2d(kmos_model)
        print("Added MFT boundary conditions")
        kmos_model.print_coverages()


def run_kmc_model_at_data_point(catmap_data, options, data_point,
                                make_plots=None,
                                log_target=None,
                                L=4,
                                alpha=0.05,
                                coverage_tolerance=0.1,
                                tof_method='procrates',
                                ):
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

        # generate the data
        n_current_point = data_point + 1

        if make_plots is None:
            if hasattr(options, 'make_plots'):
                make_plots = options.make_plots
            else:
                make_plots = False

        descriptors = catmap_data['forward_rate_constant_map'][data_point][0]
        # descriptor_string = str(descriptors)
        descriptor_string = '{data_point}\t{descriptors}'.format(**locals())

        log_filename = "procstat_{:04d}.dat".format(data_point)

        # important to change the lattice size before the Fortran object is instantiated
        import kmc_settings
        print(options)
        print(dir(options))
        kmc_settings.simulation_size = options.simulation_size

        with kmos.run.KMC_Model(print_rates=False, banner=False) as kmos_model:
            # -1 is needed to go form 1 based Fortran to 0 based C.
            elementary_process_index = dict([(i, eval('kmos_model.proclist.' + i.lower()) - 1) for i in sorted(kmos_model.settings.rate_constants)])
            start_batch = 0
            fast_processes = True
            fast_processes_adaption = 0
            renormalizations = {}
            coverages_history = {}
            rates_history = {}
            sampled_rates = {}
            increased_batch = False

            numeric_renormalizations = np.ones([kmos_model.proclist.nr_of_proc])

            old_nondiff = float('-inf')
            setup_edged_model_at_datapoint(kmos_model, data_point, reset_configuration=True)
            start_time = datetime.datetime.now()
            with open(log_filename, 'a') as procstat_file:
                if log_target is None:
                    outfile = procstat_file
                else:
                    outfile = log_target

                mft_rates = ('\n'.join(map(lambda x: str(float(x)), sorted(catmap_data['rate_map'])[data_point][1])))
                outfile.write('\nMFT Rates Benchmark\n{mft_rates}\n\n'.format(**locals()))
                outfile.write("\n\nInitial coverages, lattice-size {kmos_model.lattice.system_size}\n".format(**locals()))
                outfile.write("Sampling started at {}.\n\n".format(start_time.isoformat()))
                outfile.write(kmos_model.print_coverages(to_stdout=False))

            ####################################################################
            # Main sampling/convergence loop
            ####################################################################
            while fast_processes:
                update_outstring = True

                # set_rate_constants(kmos_model, data_point, catmap_data=catmap_data, options=options)
                # print(kmos_model.rate_constants)
                # set_kmc_model_coverage_at_data_point(kmos_model, catmap_data, options, data_point)

                setup_edged_model_at_datapoint(kmos_model, data_point)

                # reset kmc_time and kmc_steps before going any further
                kmos_model.base.set_kmc_time(0.)
                kmos_model.base.set_kmc_step(0)

                import kmos.run.steady_state

                ####################################################################
                # 1. Generate steady-state data
                ####################################################################

                full_data = kmos.run.steady_state.sample_steady_state(kmos_model,
                                                                      show_progress=True,
                                                                      make_plots=make_plots,
                                                                      start_batch=start_batch,
                                                                      batch_size=options.batch_size,
                                                                      bias_threshold=options.bias_threshold,
                                                                      L=L,
                                                                      alpha=alpha,
                                                                      output='both',
                                                                      tof_method='procrates',
                                                                      seed='EWMA_{data_point:04d}'.format(**locals()),
                                                                      renormalizations=numeric_renormalizations,
                                                                      log_filename=log_filename,
                                                                      sub_batches=1,
                                                                      )

                #data_dict, data = full_data['integ']
                data_dict, data = full_data

                normalize_coverage_data(kmos_model, data_dict)
                if len(data.strip()) > 0:
                    start_batch = int(float(data.split()[-1])) / options.batch_size
                else:
                    start_batch = 0

                import socket
                import os
                import kmos.run.steady_state

                batch_size0 = options.batch_size

                hostname = socket.getfqdn()
                working_directory = os.getcwd()
                script_path = os.path.dirname(os.path.abspath(__file__))
                with open(log_filename, 'a') as procstat_file:
                    if log_target is None:
                        outfile = procstat_file
                    else:
                        outfile = log_target

                    ####################################################################
                    # 2. Retrieve and log other data
                    ####################################################################
                    outfile.write("========\n\n")
                    outfile.write("Running at {hostname}/{working_directory}\n".format(**locals()))
                    outfile.write("\tscript located at {script_path}\n".format(**locals()))
                    outfile.write("datapoint [descriptor string]\n\n")
                    outfile.write("{data_point} {descriptor_string}\n\n".format(**locals()))
                    outfile.write("fast-process adaption step {fast_processes_adaption}\n".format(**locals()))
                    outfile.write("\n\nProcstat (number of executed processes)\n")
                    outfile.write(kmos_model.print_procstat(to_stdout=False))
                    outfile.write("\n\nCoverages, lattice-size {kmos_model.lattice.system_size}\n".format(**locals()))
                    outfile.write(kmos_model.print_coverages(to_stdout=False))
                    outfile.write("\n\nkMC state\n")
                    outfile.write(kmos_model.print_kmc_state(to_stdout=False))
                    outfile.write("\n\nProcstat (number of sampled processes)\n")
                    outfile.write(kmos_model.print_procstat(to_stdout=False, integ=True))
                    outfile.write('\n\nRate Constants\n')
                    outfile.write(kmos_model.rate_constants())
                    outfile.write('\n\nParameters\n')
                    outfile.write(kmos_model.parameters())
                    outfile.write("\n\nEquilibration Report\n")
                    equilibration_report, equilibration_data = kmos.run.steady_state.report_equilibration(kmos_model, tof_method='procrates')
                    integ_equilibration_report, integ_equilibration_data = kmos.run.steady_state.report_equilibration(kmos_model, tof_method='integ')
                    outfile.write("\n\nNumpy Renormalizations\n{numeric_renormalizations}\n\n".format(**locals()))
                    outfile.write("\nSampled rates and coverages\n")
                    outfile.write(pprint.pformat(data_dict))
                    #outfile.write("\nProcrates Sampled rates and coverages\n")
                    #outfile.write(pprint.pformat(full_data['procrates'][0]))
                    outfile.write("\n\nStatistically best sampled rates\n")
                    outfile.write(pprint.pformat(sampled_rates))
                    outfile.write("\n\nRelative rate differences between reversing processes\n")
                    outfile.write(equilibration_report)
                    #outfile.write("\n\nDEBUG\n\n")
                    #outfile.write(pprint.pformat(equilibration_data))

                    # update coverages history if we have obtained meaningful sampling
                    _current_coverages = []
                    for key, value in data_dict.items():
                        if 'time' not in key and 'forward' not in key and 'reverse' not in key and key != 'T' and 'steps' not in key and '_2_' not in key:
                            _current_coverages.append(value)

                    if not sum(_current_coverages) == 0.:


                        for key, value in data_dict.items():
                            if 'time' not in key and 'forward' not in key and 'reverse' not in key and key != 'T' and 'steps' not in key and '_2_' not in key:
                                coverages_history.setdefault(key, []).append(value)

                    EQUIB_THRESHOLD = options.equilibration_threshold
                    SAMPLE_MIN = options.sampling_min
                    STAT_MIN = SAMPLE_MIN * 10

                    sums0 = [s for r, _, s, _, _, _ in equilibration_data if s >= SAMPLE_MIN and abs(r) < EQUIB_THRESHOLD]
                    ratios = [r for r, _, s, _, _, _ in equilibration_data if s >= SAMPLE_MIN and abs(r) < EQUIB_THRESHOLD]
                    sums = [EQUIB_THRESHOLD * STAT_MIN / s for r, _, s, _, _, _ in equilibration_data if s >= SAMPLE_MIN and abs(r) < EQUIB_THRESHOLD]

                    rescale_factor = 1. / float(options.lowering_factor)

                    outfile.write("\nrescale factor calculation threshold {EQUIB_THRESHOLD}, SAMPLE_MIN. min. {SAMPLE_MIN} STAT_MIN {STAT_MIN}\n".format(**locals()))
                    outfile.write("ratios {ratios}\n".format(**locals()))
                    outfile.write("original sums {sums0}\n".format(**locals()))
                    outfile.write("rescaled sums {sums}\n".format(**locals()))

                    outfile.write("\nGlobal rescale factor {rescale_factor}".format(**locals()))

                    outfile.write("\nEvaluating equilibration report\n")

                    fastest_nondiff_unsampled_rconstant = float('-inf')
                    fastest_nondiff_unsampled_pname = ''

                    fastest_nondiff_rconstant = float('-inf')
                    fastest_nondiff_pname = ''

                    ####################################################################
                    # 3. Determine if we have sufficient data from every reaction step
                    ####################################################################
                    # Extra step: if least_sampled_pair > 0, increase batch size but don't adjust rate-constants
                    least_sampled_pair = float('inf')
                    least_sampled_pair_name = ''
                    if not all([s >= SAMPLE_MIN and data_dict[pn1] > 0. for r, pn1, _, pair, s, _ in equilibration_data if
                        'mft' not in pn1 and
                        #'mft' not in pn2 and
                        '_1p_' not in pn1 and
                        #'_1p_' not in pn2 and
                        not pair[0].rate_constant.startswith('diff')]):
                        for r, pn1, _, pair, s, _ in equilibration_data:
                            if 'mft' in pn1 or \
                               '_1p_' in pn1 or \
                               pair[0].rate_constant.startswith('diff'):
                                continue
                            if s < least_sampled_pair:
                                least_sampled_pair = s
                                least_sampled_pair_name = pn1

                        #outfile.write("\nLeast sampled pair is {least_sampled_pair_name} with {least_sampled_pair} events.\n".format(**locals()))

                    # Extra step: if least_sampled_pair > 0, increase batch size but don't adjust rate-constants
                    integ_least_sampled_pair = float('inf')
                    integ_least_sampled_pair_name = ''
                    if not all([s >= SAMPLE_MIN and data_dict[pn1] > 0. for r, pn1, _, pair, s, _ in integ_equilibration_data if
                        'mft' not in pn1 and
                        #'mft' not in pn2 and
                        '_1p_' not in pn1 and
                        #'_1p_' not in pn2 and
                        not pair[0].rate_constant.startswith('diff')]):
                        for r, pn1, _, pair, s, _ in integ_equilibration_data:
                            if 'mft' in pn1 or \
                               '_1p_' in pn1 or \
                               pair[0].rate_constant.startswith('diff'):
                                continue
                            if s < integ_least_sampled_pair:
                                integ_least_sampled_pair = s
                                integ_least_sampled_pair_name = pn1

                        outfile.write("\nLeast sampled integ pair is {integ_least_sampled_pair_name} with {integ_least_sampled_pair} events.\n".format(**locals()))
                        #if 0 < integ_least_sampled_pair and least_sampled_pair < SAMPLE_MIN :
                            ##options.lowering_factor = math.sqrt(float(options.lowering_factor))
                            #options.lowering_factor = 1.
                            #outfile.write(("Ok, we have an integ signal from every proces.\n"
                                           #"Will reduce lowering factor by applying np.sqrt to "
                                           ##"{options.lowering_factor}.\n\n"
                                    ##).format(**locals()))

                    # determine most sampled non-diff elementary process
                    # to determine if it is justified to adjust diffusion processes
                    most_sampled_pair = float('-inf')
                    most_sampled_pair_name = ''
                    for ratio, pn1, left_right_sum, pair, s, _ in equilibration_data:
                        pn = pair[0].name
                        if (not kmos_model.settings.rate_constants[pn][0].startswith('diff')
                            or 'mft' in pn
                               ):
                            if left_right_sum > most_sampled_pair:
                                most_sampled_pair = left_right_sum
                                most_sampled_pair_name = pn
                    outfile.write('\nMost sampled pair is {most_sampled_pair_name} with {most_sampled_pair}.\n\n'.format(**locals()))


                    ####################################################################
                    # 4. If necessary adjust rate constants
                    ####################################################################
                    if least_sampled_pair < options.freezing_fraction * SAMPLE_MIN:

                        ####################################################################
                        # 4.b  ... of diff events
                        ####################################################################
                        touched_diff_events = False
                        if True or fastest_nondiff_pname == '':
                            outfile.write("\n\nChecking for alternative fastest non-diff rate-constants\n")
                            for ratio, pn1, left_right_sum, pair, s, _ in integ_equilibration_data:
                                if s > SAMPLE_MIN * options.batch_size :
                                    # outfile.write("Could be {pn1} or {pn2}".format(**locals()))
                                    # outfile.write("{pn1} <=> {pn2} : {ratio}\n".format(**locals()))

                                    #pn = pn1 if (kmos_model.rate_constants.by_name(pn1) < kmos_model.rate_constants.by_name(pn2)) else pn2
                                    #pn = pn1 if (kmos_model.rate_constants.by_name(pn1) > kmos_model.rate_constants.by_name(pn2)) else pn2
                                    #pn = pn1
                                    pn = pair[0].name

                                    new_rc = kmos_model.rate_constants.by_name(pn)
                                    if new_rc > fastest_nondiff_rconstant:
                                        if 'mft' not in pn and \
                                                not kmos_model.settings.rate_constants[pn][0].startswith('diff'):
                                            outfile.write(' - found k({pn}) = {new_rc}\n'.format(**locals()))
                                            fastest_nondiff_rconstant = new_rc
                                            fastest_nondiff_pname = pn

                        outfile.write("\n\nFastest non-diff process k({fastest_nondiff_pname}) = {fastest_nondiff_rconstant:.2e}\n".format(**locals()))

                        # if non-diff processes where not sampled at all, reduce all diffusion rate-constants consistently
                        if True or fastest_nondiff_pname == '':
                            outfile.write("\n\nchecking alternatively for fastest diff rate-constants\n")
                            for ratio, pn1, left_right_sum, pair, _, _ in integ_equilibration_data:
                                if left_right_sum > 0:
                                    #outfile.write("could be {pn1} or {pn2}".format(**locals()))
                                    #outfile.write("{pn1} <=> {pn2} : {ratio}\n".format(**locals()))
                                    #for pn in [pn1, pn2]:
                                    #for pn in [pn1]:
                                    pn = pair[0].name
                                    new_rc = kmos_model.rate_constants.by_name(pn)
                                    if new_rc > fastest_nondiff_rconstant:
                                        if 'mft' not in pn and  \
                                           kmos_model.settings.rate_constants[pn][0].startswith('diff'):
                                            outfile.write(' - found k({pn}) = {new_rc}\n'.format(**locals()))
                                            fastest_nondiff_rconstant = new_rc
                                            fastest_nondiff_pname = pn

                        outfile.write("\n\nfastest diff process k({fastest_nondiff_pname}) = {fastest_nondiff_rconstant:.2e}\n".format(**locals()))

                        # bracket the maximum change of the fastest non-diff rate-constant
                        if np.isfinite(old_nondiff):
                            if fastest_nondiff_rconstant < old_nondiff / 10.:
                                fastest_nondiff_rconstant = old_nondiff / 10.
                                outfile.write('\n\nthrottling non-diff reduction, to a maximum factor of 10. to {fastest_nondiff_rconstant:.3e}\n'.format(**locals()))
                            if fastest_nondiff_rconstant >=  old_nondiff:
                                fastest_nondiff_rconstant = old_nondiff / float(options.lowering_factor)
                                outfile.write("\n\nreducing maximum non-diff factor by minimum of factor {options.lowering_factor} to {fastest_nondiff_rconstant:.3e}\n".format(**locals()))

                        # loop again for preserving the fast diffusion rate constants
                        if fastest_nondiff_rconstant > float('-inf') and options.diffusion_factor is not None:
                            # options.batch_size = batch_size0
                            for pn in kmos_model.settings.rate_constants:
                                if kmos_model.settings.rate_constants[pn][0].startswith('diff'):
                                    #for ratio, pn1, left_right_sum, pair, s, _ in integ_equilibration_data:
                                    for ratio, pn1, left_right_sum, pair, s, _ in equilibration_data:
                                        #if pn in [pair[0].name, pair[1].name] and (left_right_sum >= 1000 * sample_min or fast_processes_adaption < 1) :
                                        if (pn in [pair[0].name, pair[1].name]
                                            #and left_right_sum >= options.batch_size * SAMPLE_MIN
                                            #and left_right_sum >= most_sampled_pair
                                            # EXPERIMENTAL
                                            #and ((left_right_sum >= most_sampled_pair and left_right_sum > 0)  or fast_processes_adaption < 1)
                                            and ((left_right_sum >= 1. * most_sampled_pair and left_right_sum > 0)  or fast_processes_adaption < 1)
                                            and not '_1p_' in pn
                                            and not 'mft' in pn
                                            ):
                                            outfile.write('\t- adapting {pn}, bc. ((lr_sum = {left_right_sum})>(ms_pair={most_sampled_pair}))\n'.format(**locals()))
                                            # note: the fast_process_adaption < 1 condition
                                            # makes sure that the diffusion factor gets applied at least
                                            # once in the begining, even if the corresponding process is
                                            # not sufficiently sampled
                                            rc_tuple = kmos_model.settings.rate_constants[pn]
                                            diff_const = rc_tuple[0].split('*')[0]
                                            #theta_terms = '*'.join([term for term in rc_tuple[0].split('*') if term.startswith('theta_')])
                                            # we need the `re` version of split to make sure we only split on '*' but not on '**'
                                            theta_terms = '*'.join([term for term in re.split("(?<![*])[*](?![*])", rc_tuple[0])
                                                if 'theta_' in term
                                                or 'n_sites' in term
                                                ])
                                            safe_rc = (options.diffusion_factor * fastest_nondiff_rconstant)
                                            diff_rconst = '{diff_const}*({diff_const})**(-1.)*{options.diffusion_factor}*{fastest_nondiff_rconstant}'.format(**locals())
                                            if theta_terms:
                                                diff_rconst += '*' + theta_terms
                                            # we filter the diffusion processes by checking which rate constants start with 'diff', thus
                                            # we need this weird way of writing "1*" to keep the 'diff' prefix throughout adaptations
                                            ndiff_rconst = options.diffusion_factor * fastest_nondiff_rconstant
                                            rc_tuple = (diff_rconst, rc_tuple[1])
                                            # experimental
                                            #lowering_factor = math.pow(float(options.lowering_factor), 0.1)
                                            lowering_factor = float(options.lowering_factor)
                                            rc_tuple = ('{rc_tuple[0]}/{lowering_factor:5.3e}'.format(**locals()), rc_tuple[1])
                                            outfile.write("\t- diff reset k({pn}) = {diff_rconst} = {rc_tuple}\n".format(**locals()))
                                            touched_diff_events = True
                                            kmos_model.settings.rate_constants[pn] = rc_tuple
                                            break
                                        else:
                                            pass
                                            #outfile.write('\t\t- not touching {pn}, bc. not ((lr_sum = {left_right_sum})>(ms_pair={most_sampled_pair}))\n'.format(**locals()))

                        else:
                            options.batch_size *= 2
                            outfile.write("\n\ndidn't sample non-diff processes, doubling batch-size to {options.batch_size:.3e}\n".format(**locals()))

                        # outfile.write("\n\nrate constants after all adjustments\n")
                        # outfile.write(kmos_model.rate_constants())

                        if old_nondiff == fastest_nondiff_rconstant:
                            options.batch_size *= 2
                            outfile.write("fastest non-diff rate-constant unchanged, doubling batch size to {options.batch_size:.3e}\n".format(**locals()))

                        old_nondiff = fastest_nondiff_rconstant


                        ####################################################################
                        # 4.a  ... of non-diff events
                        ####################################################################
                        # First loop: test for equilibrated pairs of elementary processes
                        # that have been sampled many times and adjust those rate constants
                        if touched_diff_events == True:
                            outfile.write('\nJust adapted diffusion rate constant. Will skip non-diff adaption for the rest of this round.\n')
                        else:
                            for ratio, pn1, left_right_sum, pair, _, _ in equilibration_data:
                                # Minimum number of events, to produce statistically meaningful results
                                if abs(ratio) < EQUIB_THRESHOLD and left_right_sum >= 10 * SAMPLE_MIN:
                                #if abs(ratio) < EQUIB_THRESHOLD and left_right_sum >= 100 * SAMPLE_MIN:
                                #if abs(ratio) < EQUIB_THRESHOLD and left_right_sum >= 1000 * SAMPLE_MIN:
                                #if abs(ratio) < EQUIB_THRESHOLD and left_right_sum >= options.batch_size * SAMPLE_MIN:
                                #if abs(ratio) < EQUIB_THRESHOLD and left_right_sum >= 2 * options.batch_size * SAMPLE_MIN:
                                    fast_processes = True
                                    pn = pair[0].name
                                    if (
                                            not kmos_model.settings.rate_constants[pn][0].startswith('diff')
                                            or 'mft' in pn
                                            or '1p' in pn

                                            #or ('mft' in pn and '1p' not in pn)
                                            #(not kmos_model.settings.rate_constants[pn][0].startswith('diff')
                                            #or 'mft' in pn)
                                            #and '1p' not in pn
                                            ):
                                        old_rc = kmos_model.rate_constants.by_name(pn)
                                        rc_tuple = kmos_model.settings.rate_constants[pn]
                                        # rc_tuple = (rc_tuple[0] + '*.5', rc_tuple[1])
                                        # rescale_factor = max(EQUIB_THRESHOLD, (abs(ratio) / EQUIB_THRESHOLD))
                                        rescale_multiplication = '*%.2e' % rescale_factor
                                        rescale_division = '/%.2e' % rescale_factor
                                        renormalizations[pn] = renormalizations.get(pn, '1.') + '/{:.2e}'.format(rescale_factor)
                                        # Hackish way of producing the same normalizations for 2 significant digits
                                        numeric_renormalizations[elementary_process_index[pn]] /= eval('{:.2e}'.format(rescale_factor))
                                        rc_tuple = (rc_tuple[0] + rescale_multiplication, rc_tuple[1])
                                        kmos_model.settings.rate_constants[pn] = rc_tuple
                                        new_rc = kmos_model.rate_constants.by_name(pn)
                                        kmos_model.rate_constants.set(pn, new_rc)
                                        outfile.write(("Found a fast equilibrated process {ratio}:"
                                                       "\t{pn}, reduced rate constant from {old_rc:.2e} to {new_rc:.2e}\n"
                                                       "\tconst. sample-min {SAMPLE_MIN}, left-right sum {left_right_sum}\n"
                                                       "\tconst. equib-threshol {EQUIB_THRESHOLD}, ratio {ratio}\n\n"

                                            ).format(**locals()))

                                    #pn = pn1 if (kmos_model.rate_constants.by_name(pn1) > kmos_model.rate_constants.by_name(pn2)) else pn2
                                    #pn = pn1
                                    if not kmos_model.settings.rate_constants[pn][0].startswith('diff'):
                                        new_rc = kmos_model.rate_constants.by_name(pn)
                                        if new_rc > fastest_nondiff_rconstant:
                                            if 'mft' not in pn:
                                                fastest_nondiff_rconstant = new_rc
                                                fastest_nondiff_pname = pn
                                else:
                                    pn = pair[0].name
                                    outfile.write("Not touching {pn}, because ratio = {ratio} < {EQUIB_THRESHOLD} or left_right_sum = {left_right_sum} < {SAMPLE_MIN}\n".format(**locals()))


                    ####################################################################
                    # 5. or else just make sampling longer
                    ####################################################################
                    else:
                        options.batch_size *= 2
                        #increased_batch = True
                        outfile.write("Least sampled pair signal from {least_sampled_pair_name} is {least_sampled_pair}. (SAMPLEMIN = {SAMPLE_MIN})\n".format(**locals()))
                        outfile.write('Ok, we have a signal from every processes pair, will just increase batch-size to {options.batch_size:.3f}\n'.format(**locals()))

                    ####################################################################
                    # 6. Final checks if we are done sampling, and prepare next loop
                    ####################################################################
                    if False:
                        # Randomize coverage a bit to avoid dead-locks
                        random_bits = 0
                        special_species = [getattr(kmos_model.proclist, species_name.lower())
                                           for species_name in kmos_model.settings.representations
                                           if species_name.endswith('_')]
                        for x in range(kmos_model.lattice.system_size[0]):
                            for y in range(kmos_model.lattice.system_size[1]):
                                for z in range(kmos_model.lattice.system_size[2]):
                                    for n in range(kmos_model.lattice.spuck):
                                        site = (x, y, z, n + 1)
                                        species = kmos_model.lattice.get_species(site)
                                        if species not in special_species:
                                            if random.random() < 1./10:
                                                kmos_model._put(site, kmos_model.proclist.default_species)
                                                random_bits += 1
                        outfile.write('\nRandomized configuration on {random_bits} sites.\n\n'.format(**locals()))
                        kmos_model._adjust_database()
                        outfile.write('\n\nPeturbed coverages\n')
                        outfile.write(kmos_model.print_coverages(to_stdout=False))



                    # Reset procstat and kmc steps
                    for proc in range(kmos_model.proclist.nr_of_proc.max()):
                        kmos_model.base.set_procstat(proc + 1, 0)
                        kmos_model.base.set_integ_rate(proc + 1, 0)
                        kmos_model.base.set_integ_event(proc + 1, 0)
                    kmos_model.base.set_kmc_step(0)

                    #outfile.write("\n\nReset Procstat (number of executed processes)\n")
                    #outfile.write(kmos_model.print_procstat(to_stdout=False))
                    #outfile.write(kmos_model.print_procstat(integ=True, to_stdout=False))

                    outfile.write("\nRenormalizations\n")
                    outfile.write(pprint.pformat(renormalizations))

                    # Check if every process has been touched in this round
                    if all([rate < EQUIB_THRESHOLD for (rate, _, _, _, _, _) in equilibration_data]):
                        fast_processes = False
                        update_outstring = True
                        outfile.write("\nFound all processes, to be equilibrated. So further adjustments will not help. Quit.\n")

                    # Check if we have sufficient sampling for every process pair
                    if all([s >= SAMPLE_MIN for r, pn1, _, pair, s, _ in integ_equilibration_data if
                        'mft' not in pn1 and
                        #'mft' not in pn2 and
                        '_1p_' not in pn1 and
                        #'_1p_' not in pn2 and
                        not pair[0].rate_constant.startswith('diff')]):
                        fast_processes = False
                        update_outstring = True
                        outfile.write('\n\nFinal pair-sampling statistic\n\n')
                        for r, pn1, _, _, s, _ in integ_equilibration_data:
                            outfile.write('\n\t- {s} events for ({pn1})'.format(**locals()))
                        outfile.write("\n\nObtained well-sampled statistics for every process-pair, no further sampling needed. Done.\n")
                    else:
                        outfile.write('\n\nProcess pairs that are not sufficiently sampled :\n')
                        for r, pn1, _, pair, s, _ in integ_equilibration_data:
                            if 'mft' in pn1 or \
                               '_1p_' in pn1 or \
                               pair[0].rate_constant.startswith('diff'):
                                continue
                            if s < SAMPLE_MIN:
                                outfile.write('\n\t- only {s} events for ({pn1})'.format(**locals()))
                        outfile.write('\n')

                    # Check if we have obtained meaningful data at all
                    if sum(_current_coverages) == 0.:
                        outfile.write("\nCould not obtain coverage data at all, will not update result.\n")
                        update_outstring = False

                    # Check if one or more coverages has become sensitive to adaptations
                    for key, values in coverages_history.items():
                        if len(values) >= 2:
                            _vm1 = values[-1]
                            _vm2 = values[-2]
                            if abs(_vm2 - _vm1) > coverage_tolerance:
                                outfile.write(kmos_model.print_coverages(to_stdout=False))
                                fast_processes = False
                                update_outstring = False
                                outfile.write("\nCoverage {key} changed from {_vm2} to {_vm1}, critical. Exiting!\n".format(**locals()))

                    # Check if time is overrun
                    if (math.isinf(data_dict['kmc_time'])
                            or math.isnan(data_dict['kmc_time'])
                            or math.isinf(data_dict['simulated_time'])
                            or math.isnan(data_dict['simulated_time'])
                        ):
                        options.batch_size *= .3
                        outfile.write("kMC time overrun. Will give it another round and reduce batch-size to {options.batch_size:.3e}.\n".format(**locals()))
                        fast_processes = True
                        update_outstring = False

                    #if fast_processes_adaption <= 0:
                        #outfile.write("Not so fast. Let's do at least two rounds.\n".format(**locals()))
                        #fast_processes = True
                        #update_outstring = False

                    if update_outstring:
                        outfile.write('\n\nUpdating outstring\n')
                        best_sampled_rates = {}
                        for key, value in sampled_rates.items():
                            best_sampled_rates[key] = value['rate']
                        data = ' '.join(format(data_dict[key.replace('#', '')], '.5e') for key in kmos_model.get_std_header().split()) + '\n'
                        outstring = _get_outstring(data_point, descriptors, data)
                        #outstring = _get_outstring(data_point, descriptors, full_data['procrates'][1])
                        outfile.write('\n\n-outstring\n')
                        outfile.write(outstring)

                    isoformat = datetime.datetime.now().isoformat()
                    outfile.write('\n\nLocal time: {isoformat}\n'.format(**locals()))
                    outfile.write("\nWalltime spent until now : {}\n\n".format(str(datetime.datetime.now() - start_time)))

                fast_processes_adaption += 1
            return outstring


def normalize_coverage_data(model, data):
    """When the kMC model uses a background MFT species to avoid
    extremely large kMC jumps this function normalizes the sites
    filled with "MFT" out.
    """
    for site in model.settings.site_names:
        site_sum = sum([value
                        for key, value in data.items()
                        if key.endswith(site) and
                        not key.startswith('MFT_')
                        ])
        if site_sum > 0:
            for key, value in data.items():
                if key.endswith(site):
                    if key.startswith('MFT_'):
                        data[key] = 0.
                    else:
                        data[key] /= site_sum


def update_mft_parameters(model, data):
    for site in model.settings.site_names:
        site_name = '_'.join(site.split('_')[1:])
        for species in model.settings.species_tags:
            if species == 'MFT':
                continue

            setattr(model.parameters,
                    'Theta_{site_name}_{species}'.format(**locals()),
                    data['{species}_{site}'.format(**locals())]
                    )


def _get_outstring(data_point, descriptors, data):
    if data.strip():
        return '{data_point:9d} {descriptors[0]: .5e} {descriptors[1]: .5e} {data}'.format(**locals())
    else:
        return '# {data_point:9d} {descriptors[0]: .5e} {descriptors[1]: .5e} EMPTY DATA RETURN CHECK procstat_{data_point:04d}.dat TO SEE WHAT WENT WRONG\n'.format(**locals())


def sort_catmap_maps_inplace(data):
    for key in data:
        if key.endswith('_map'):
            data[key] = sorted(data[key], key=lambda x: x[0])


def run_model(seed, call_path=None, options=None):
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

    for data_point in reversed(range(total_points)):
        descriptors = catmap_data['forward_rate_constant_map'][data_point][0]

        # multi IO mechanism: keep lockfile with one line descriptors of datapoint
        # if datapoint is already in there, skip to next datapoint
        # descriptor_string = str(descriptors) + '\n'
        descriptor_string = '{data_point}\t{descriptors}\n'.format(**locals())

        print('\n\nrunning DATAPOINT {data_point}/{total_points} DESCRIPTOR {catmap_model.descriptor_names} = {descriptor_string}'.format(**locals()))

        with open(lock_filename, 'r') as lockfile:
            if descriptor_string in lockfile.readlines():
                # print('Skipping {descriptor_string}'.format(**locals()))
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


def set_mft_parameters_at_data_point(model, data_point):
    seed = model.settings.model_name
    data = merge_catmap_output(seed=seed)
    sort_catmap_maps_inplace(data)
    coverage = data['coverage_map'][data_point][1]
    print("CatMAP coverage {coverage}".format(**locals()))
    site_names = [site for site in data['site_names'] if not site == 'g']
    adsorbate_names = data['adsorbate_names']
    coverage = dict(zip(adsorbate_names, coverage))
    print("CatMAP coverage {coverage}".format(**locals()))
    for site_name in site_names:
        empty_species = catmap.cli.kmc_translation.EMPTY_SPECIES
        coverage['{empty_species}_{site_name}'.format(**locals())] =  \
            1 - sum([_theta
                     for (_species, _theta)
                     in coverage.items()
                     if _species.endswith('_{site_name}'.format(**locals()))])
    print("CatMAP coverage {coverage}".format(**locals()))

    for species_site, value in coverage.items():
        for n in range(10):  # there shouldn't more equivalent sites per unit cell than this ...
            species, site = species_site.split('_')
            if hasattr(model.parameters, 'Theta_{site}_{n}_{species}'.format(**locals())):
                setattr(model.parameters,
                        'Theta_{site}_{n}_{species}'.format(**locals()),
                        value
                        )


def setup_model_probabilistic(model, data_point=0, majority=False):
    """Make an educated initial guess for coverages by setting
    adsorbates for each with probabilites according to the
    mean field model result.

    Note that this naive way of implementing this initialization
    might lattices which are locally not realistic.

    If majority==True all sites are two the respective majority-species (a.k.a winner takes it all).
    """

    import numpy.random

    seed = model.settings.model_name

    data = merge_catmap_output(seed=seed)
    sort_catmap_maps_inplace(data)
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
        choices_weights[model.proclist.default_species] = 1 - sum(choices_weights)

        choices_weights = map(float, choices_weights)

        if majority:
            argmax = np.argmax(choices_weights)
            choices_weights = [0.] * len(choices_weights)
            choices_weights[argmax] = 1.

        # DEBUGGING
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

    print("kMC lattice size {model.lattice.system_size}".format(**locals()))
    print("Initial kmos Coverages")
    model.print_coverages()

def setup_mft_edges_2d(model, grid_size=None):
    """To avoid extremely slow transitions (i.e. 1e-30 seconds)
    one can encase the active lattice with a layer of MFT species.
    MFT species aims to circumvent these extremely large kMC jumps
    by offering a reservoir of background species. This function
    sets a up a boundary of MFT along the edges (for 2d).
    """
    X, Y, Z = model.lattice.system_size
    if grid_size is None:
        grid_size = X
    #old version
    #for x in range(X):
        #model._put([x, Y - 1, 0, 1], model.proclist.mft_)
        #model._put([x, 0, 0, 1], model.proclist.mft_)
#
    #for y in range(Y):
        #model._put([0, y, 0, 1], model.proclist.mft_)
        #model._put([X - 1, y, 0, 1], model.proclist.mft_)

    # newer version: create finer grid of MFT processes
    for x in  range(X):
        for y in range(Y):
            if x % grid_size == 0 or y % grid_size == 0:
                model._put([x, y, 0, 1], model.proclist.mft_)
                pass

    model._adjust_database()

def set_rate_constants(kmos_model, data_point, catmap_data=None, options=None):
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
    if catmap_data is None:
        seed = kmos_model.settings.model_name
        catmap_data = merge_catmap_output(seed=seed)
        sort_catmap_maps_inplace(catmap_data)

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
    reverse_rate_constants = rate_constants[n_rate_constants / 2:]

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


def set_rate_constants_from_procstat_logfile(model, log_filename, step=0):
    """Utility function to support debugging: set the (possibly rescaled)
    rate constants from a generated procstat_XXXX.dat file.

    """

    with open(log_filename) as log_file:
        step0 = step
        while step >= 0:
            line = log_file.readline()
            if line.strip() == 'Rate Constants':
                step += -1
            if line == '':
                print("Warning: Logfile {log_filename} ended before step {step0} was reached, pick an earlier step".format(**locals()))
                return

        for i in range(model.proclist.nr_of_proc):
            line = log_file.readline().replace(':', '')
            elements = line.split()
            process_name = elements[1]
            rate_constant = float(elements[4])
            model.rate_constants.set(process_name, rate_constant)


def get_rates_history_from_procstat_logfile(log_filename):
    with open(log_filename) as infile:
        rates_history = {}
        while True:
            line = infile.readline()
            if line == '':
                break

            if line.startswith('datapoint'):
                infile.readline()
                datapoint_line = infile.readline()
                data_point, descriptor0, descriptor1 = map(float,
                                                           datapoint_line.strip()
                                                           .replace('[', '')
                                                           .replace(']', '')
                                                           .replace(',', '')
                                                           .split())
                data_point = int(data_point)
            if line.startswith('rates history: updated'):
                raw_rates_history = ''
                while not line.strip() == '':
                    line = infile.readline()
                    raw_rates_history += line
                rates_history = eval(raw_rates_history)

        return rates_history


def get_rates_from_procstat_logfile(model, log_filename, output='str'):
    rates_list = []
    with open(log_filename) as infile:
        while True:
            line = infile.readline()

            if line.startswith('datapoint'):
                infile.readline()
                datapoint_line = infile.readline()
                data_point, descriptor0, descriptor1 = map(float,
                                                           datapoint_line.strip()
                                                           .replace('[', '')
                                                           .replace(']', '')
                                                           .replace(',', '')
                                                           .split())
                data_point = int(data_point)

            if line == '':
                break
            if 'Sampled rates and coverages' in line:
                rates = ''
                while True:
                    line = infile.readline()
                    if line.strip() == '':
                        break
                    rates += line
                rates = eval(rates)
                rates['data_point'] = data_point
                rates['descriptor0'] = descriptor0
                rates['descriptor1'] = descriptor1
                if output == 'dict':
                    rates_list.append(rates)
                elif output == 'str':
                    format_str = ' '.join(map(lambda x: '{{{x:s}:.5e}}'.format(x=x), model.get_std_header()[1:].split()))
                    data = format_str.format(**rates)
                    descriptors = [descriptor0, descriptor1]

                    rates_list.append(_get_outstring(data_point, descriptors, data) + '\n')
                else:
                    raise UserWarning("Don't know this format, should be either str or dict")

    return rates_list


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
    rate_constants = dict(zip(sorted(model.settings.rate_constants), (model.base.get_rate(i + 1) for i in range(len(procstat)))))

    integ_rates = dict(zip(sorted(model.settings.rate_constants), atoms.integ_rates))

    report = ''
    for pair in pairs:
        pn1, pn2 = pair[0].name, pair[1].name
        left = integ_rates[pn1] * rate_constants[pn2]
        right = integ_rates[pn2] * rate_constants[pn1]
        ratio = left / right
        eq_score = 4 * left * right / (left + right)**2
        report += ('{pn1} : {pn2} => {left:.2e}/{right:.2e} = {ratio:.2e}, eq_score = {eq_score:.7f}\n'.format(**locals()))
    return report

if __name__ == '__main__':
    main()
