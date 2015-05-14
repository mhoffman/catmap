import numpy.random
import pickle

def setup_model(model, data_point=0, data_file='CO_oxidation.pkl'):
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
