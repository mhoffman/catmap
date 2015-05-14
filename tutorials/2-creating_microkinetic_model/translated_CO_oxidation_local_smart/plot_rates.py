#!/usr/bin/env python

import mix_model_runner
import numpy as np
import catmap


import optparse

parser = optparse.OptionParser()

options, args = parser.parse_args()

if len(args) >= 1:
    seed = args[1]
else:
    seed = 'CO_oxidation'

model = catmap.ReactionModel(setup_file='{seed}.mkm'.format(**locals()))

data = np.recfromtxt('{seed}_kMC_output.log'.format(**locals()), names=True)

for name in data.dtype.names:
    print(name)
    if name in ['descriptor1', 'descriptor0']:
        continue
    if '_2_' in name:
        plot_data = np.log(data[name])
        normalized = False
    else:
        plot_data = data[name]
        normalized = True

    plot_data[np.logical_not(np.isfinite(plot_data))] = np.nanmin(plot_data)
    mix_model_runner.contour_plot_data(data['descriptor0'],
                                       data['descriptor1'],
                                       plot_data,
                                       filename='output_{name}.pdf'.format(**locals()),
                                       title=name,
                                       normalized=normalized,
                                      )
