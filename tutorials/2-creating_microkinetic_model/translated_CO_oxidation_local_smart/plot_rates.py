#!/usr/bin/env python

import mix_model_runner
import numpy as np

data = np.recfromtxt('CO_oxidation_DATA.log', names=True)

for name in data.dtype.names:
    print(name)
    if name in ['descriptor1', 'descriptor0']:
        continue
    mix_model_runner.contour_plot_data(data['descriptor0'],
                                       data['descriptor1'],
                                       data[name],
                                       'output_{name}.pdf'.format(**locals()))
