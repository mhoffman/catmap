#!/usr/bin/env python

import os
import kmos.run
import catmap
import pickle

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import numpy as np


INIT_STEPS = int(1e5)
SAMPLE_STEPS = int(1e5)
SEED = 'CO_oxidation'


def set_rate_constants(kmos_model, catmap_data, data_point):
    # set rate constant of kMC according to current descriptor tuple
    for i in range(len(catmap_data['forward_rate_constant_map'][data_point][1])):
        setattr(kmos_model.parameters, 'forward_{i}'.format(
            **locals()), catmap_data['forward_rate_constant_map'][data_point][1][i])
        setattr(kmos_model.parameters, 'reverse_{i}'.format(
            **locals()), catmap_data['reverse_rate_constant_map'][data_point][1][i])


def run_model(seed, init_steps, sample_steps):
    data_filename = '{seed}_DATA.log'.format(**locals())
    lock_filename = '{seed}.lock'.format(**locals())

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
            with kmos.run.KMC_Model(print_rates=False, banner=False) as kmos_model:
                data_header = kmos_model.get_std_header()[1:]
                outfile.write(
                    '# descriptor0 descriptor1 {data_header}'.format(**locals()))

    for data_point in range(len(catmap_data['forward_rate_constant_map'])):
        descriptors = catmap_data['forward_rate_constant_map'][data_point][0]

        # multi IO mechanism: keep lockfile with one line descriptors of datapoint
        # if datapoint is already in there, skip to next datapoint
        descriptor_string = str(descriptors)

        with open(lock_filename, 'r') as lockfile:
            if descriptor_string in lockfile.readlines():
                continue
        with open(lock_filename, 'a') as lockfile:
            lockfile.write('{descriptor_string}\n'.format(**locals()))

        with kmos.run.KMC_Model(print_rates=False) as kmos_model:
            set_rate_constants(kmos_model, catmap_data, data_point)

            print('DESCRIPTORS {descriptor_string} DATAPOINT {data_point}'
                   .format(**locals()))
            print(kmos_model.rate_constants)

            # run model (hopefully) to steady state and evaluate
            kmos_model.do_steps(init_steps)
            data = kmos_model.get_std_sampled_data(1, sample_steps, verbose=True)

            with open(data_filename, 'a') as outfile:
                outfile.write(
                    '{descriptors[0]} {descriptors[1]} {data}'.format(**locals()))


def contour_plot_data(x, y, z, filename, n_gp=15):
    import numpy
    import scipy.interpolate
    print(z)

    fig = plt.figure()

    xi, yi = np.linspace(x.min(), x.max(), n_gp), np.linspace(y.min(), y.max(), n_gp)
    xi, yi = np.meshgrid(xi, yi)

    rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    zi = rbf(xi, yi)

    plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
               extent=[x.min(), x.max(), y.min(), y.max()])

    plt.scatter(x, y, c=z)
    plt.colorbar()
    plt.savefig(filename)



if __name__ == '__main__':
    run_model(seed=SEED,
         init_steps=INIT_STEPS,
         sample_steps=SAMPLE_STEPS)
