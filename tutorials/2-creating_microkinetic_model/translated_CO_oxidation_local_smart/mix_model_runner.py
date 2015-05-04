#!/usr/bin/env python

import kmos.run
import catmap
import pickle


seed = 'CO_oxidation'
kmos_model = kmos.run.KMC_Model(print_rates=False)
catmap_model = catmap.ReactionModel(setup_file='{seed}.mkm'.format(**locals()))

catmap_model.output_variables.append('rate_constant')
catmap_model.output_variables.append('forward_rate_constant')
catmap_model.output_variables.append('reverse_rate_constant')
catmap_model.run()

with open('{seed}.pkl'.format(**locals())) as infile:
    data = pickle.load(infile)

print(data.keys())

print(data['forward_rate_constant_map'][0][1])

data_point = 1


print('descriptors = ', data['forward_rate_constant_map'][data_point][0])
for i in range(len(data['forward_rate_constant_map'][data_point][1])):
    setattr(kmos_model.parameters, 'forward_{i}'.format(**locals()), data['forward_rate_constant_map'][data_point][1][i])
    setattr(kmos_model.parameters, 'reverse_{i}'.format(**locals()), data['reverse_rate_constant_map'][data_point][1][i])

print(kmos_model.rate_constants)
