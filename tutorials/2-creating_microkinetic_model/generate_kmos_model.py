#!/usr/bin/env python

from kmos.types import Project, Site, Condition, Action
from catmap import ReactionModel
import numpy as np


# TODO : write function which generates n-nearest sites of certain type
# TODO : write function which finds nearest neighbors shell for adsorbate interaction

def find_nth_neighbor_sites(central_site, site, neighbor_site, n):
    distance_dict = {}
    neighbor_sites = []

    for site in
        if site.name == neighbor_site:
            distance_list.append(np.linalg.norm(
            CONTINUE HERE


def catmap2kmos(cm_model,
    unit_cell=None,
    site_positions=None,
    diffusion_barriers=None,
    species_representations=None,
    surface_representation=None,
    ):

    EMPTY_SPECIES = 'empty'

    # initialize model
    pt = Project()

    # add meta information
    # TODO : should we add corresponding information
    #        in CatMAP ?
    pt.set_meta(author='CatMAP USER',
                email='mkm-developers-request@lists.stanford.edu',
                model_name='CatMAP_translated_model',
                model_dimension='2',   # default ...
                debug=0,  # let's be pessimistic
                )
    
    # add unit cell
    layer = pt.add_layer(name='default')
    pt.layer_list.cell = unit_cell

    # add site positions
    for site_name in cm_model.site_names:
        if not site_name == 'g':
            layer.sites.append(Site(name=site_name, pos=site_positions[site_name]))

    Coord = pt.layer_list.generate_coord



    # add species
    pt.add_species(name=EMPTY_SPECIES)
    for species_definition in cm_model.species_definitions.keys():
        if '_' in species_definition:
            species_name, site = species_definition.split('_')
            if not species_name == '*' \
                and not species_name == 's' \
                and not species_name in [x.name for x in pt.species_list]:
                species_name = species_name.replace('-', '_')
                pt.add_species(name=species_name)
            

    # add parameters

    # add processes
    site_names = [x.name for x in pt.layer_list[0].sites]
    for ri, elementary_rxn in enumerate(cm_model.elementary_rxns):
        step = {}
        surface_intermediates = {}
        if len(elementary_rxn) == 2:
            step['A'] = None
            step['B'], step['C'] = elementary_rxn
        elif len(elementary_rxn) == 3:
            step['A'], step['B'], step['C'] = elementary_rxn
        print('Steps {step}'.format(**locals()))


        print(surface_intermediates)
        for reversible, (X, Y) in [[True, ('A', 'B')],
                                   [False, ['B' , 'C']]]:
            if step[X] and step[Y]:
                # add reversible step between A and B
                surface_intermediates[X] = []
                surface_intermediates[Y] = []

                for x in [X, Y]:
                    print(x)
                    for intermediate in step[x]:
                        if '_' in intermediate:
                            species, site = intermediate.split('_')
                            if site in site_names :
                                surface_intermediates[x].append([species, site])
                        elif intermediate in site_names:
                            surface_intermediates[x].append([EMPTY_SPECIES, intermediate])
                print('Elementary Rxn: {elementary_rxn}, Surface intermediates {surface_intermediates}'.format(**locals()))

                conditions = []
                actions = []
                for i, (species, site) in enumerate(surface_intermediates[X]):
                    conditions.append(Condition(species=species.replace('-', '_'), coord=Coord('{site}.(0, 0, {i})'.format(**locals()))))
                for i, (species, site) in enumerate(surface_intermediates[Y]):
                    actions.append(Condition(species=species.replace('-', '_'), coord=Coord('{site}.(0, 0, {i})'.format(**locals()))))

                if reversible :
                    prefix = 'reversible'
                else:
                    prefix = 'irreversible'
                pt.add_process(name='{prefix}_process_{ri}'.format(**locals()),
                               conditions=conditions,
                               actions=actions,
                               rate_constant='1.')

                if reversible:
                    pt.add_process(name='{prefix}_process_{ri}_back'.format(**locals()),
                                   conditions=actions,
                                   actions=conditions,
                                   rate_constant='1.')



        if step['B'] and step['C']:
            pass
            # add irreversible from B to C
                



    return pt

if __name__ == '__main__':
    catmap_model = ReactionModel(setup_file='CO_oxidation.mkm')
    catmap_model.run()

    # additional information needed for CatMAP/kMOS translation
    unit_cell = np.diag([3.1, 3.1, 18.])

    site_positions = {
        's': [0, 0, 0],
    }

    kmos_model = catmap2kmos(catmap_model,
                             unit_cell=unit_cell,
                             site_positions=site_positions,
                             )

    #print(kmos_model)
    kmos_model.print_statistics()
    for suffix in ['ini', 'xml']:
        kmos_model.save('translated_CatMAP_model.{suffix}'.format(**locals()))
