#!/usr/bin/env python

import os.path
import copy
import pprint
import itertools

from kmos.types import Project, Site, Condition, Action
import kmos.utils
import ase.atoms
from ase.atoms import Atoms

from catmap import ReactionModel
import numpy as np


# TODO : write function which finds nearest neighbors shell for adsorbate interaction
# test for more than one site per unit cell
# test for more more than one identical site per unit cell (e.g. bridge in
# Pd(1000) cell)


# def find_nth_neighbor_sites(central_site, site, neighbor_site, n):
#distance_dict = {}
#neighbor_sites = []
#
# for site in:
# if site.name == neighbor_site:
# print(site)
# distance_list.append(np.linalg.norm(
# CONTINUE HERE


def catmap2kmos(cm_model,
                unit_cell=None,
                site_positions=None,
                diffusion_barriers=None,
                species_representations=None,
                surface_representation='None',
                adsorbate_interaction=2,
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

    if unit_cell is None:
        unit_cell = np.diag([1., 1., 1.]).tolist()

    # Need dummy atom (here 'H') so that ase.atoms.Atoms doesn't puke further
    # down
    if hasattr(cm_model, 'background_representation') and hasattr(cm_model, 'background_representation'):
        print("Warning: unit cell in background_representation will override the given background representation")
    background_representation = getattr(cm_model,
                                        'background_representation',
                                        'Atoms("H", [[0., 0., -0.1]], cell={unit_cell})'.format(**locals())
    )
    pt.layer_list.representation = '[{background_representation}]'.format(**locals())


    # add site positions
    for site_name in cm_model.site_names:
        if not site_name == 'g':
            print(cm_model.site_positions)
            for i, site_position in enumerate(cm_model.site_positions[site_name]):
                layer.sites.append(Site(name='{site_name}_{i}'.format(**locals()),
                                        pos=site_position))

    Coord = pt.layer_list.generate_coord

    # add species

    # for generating random colors
    import random

    pt.add_species(name=EMPTY_SPECIES, color='#ffffff')
    species_names = set()
    for species_definition in cm_model.species_definitions.keys():
        color = '#%02X%02X%02X' % (random.randint(0, 255),
                                   random.randint(0, 255),
                                   random.randint(0, 255))

        print('SPECIES DEFINITION {species_definition}'.format(**locals()))
        if '_' in species_definition:
            species_name, site = species_definition.split('_')
            species_name = species_name.replace('-', '_')

            species_representation = getattr(cm_model, 'species_representation', {}).get(species_name, "Atoms()")
            if not species_name == '*' \
                    and not species_name == 's' \
                    and not '_' in species_name \
                    and not species_name in [x.name for x in pt.species_list]:
                pt.add_species(name=species_name,
                               representation=species_representation,
                               color=color,
                              )

    # figure out which adsorbates can be one which site
    # extract all existing reaction terms from catmap model
    # this will be needed in case we later want to auto-generate adsorbate-adsorbate interaction
    reaction_terms = list(cm_model.elementary_rxns)
    # flatten the reaction terms list (twice)
    reaction_terms = [item for sublist in reaction_terms for item in sublist]
    reaction_terms = [item for sublist in reaction_terms for item in sublist]
    # build up site -> species mapping
    site_species = {}
    for term in reaction_terms :
        if not '-' in term: # skip terms with a transition state
            if '_' in term:
                species, site = term.split('_')
            else:
                species, site = 'empty', term
            if not species in site_species.get(site, []):
                site_species.setdefault(site, []).append(species)

    # add parameters
    # TODO : let's avoid this for now and just accept rate constants from
    # catmap

    # add processes
    site_names = [x.name for x in pt.layer_list[0].sites]
    print('SITE NAMES {site_names}'.format(**locals()))
    for ri, elementary_rxn in enumerate(cm_model.elementary_rxns):
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
        print('\n\n\nSteps {step}'.format(**locals()))

        print(surface_intermediates)
        # for reversible, (X, Y) in [[True, ('A', 'B')],
        #[False, ['B' , 'C']]]:
        # DEBUGGING: make everything reversible for now
        #            since we want a non-crashing model first
        for reversible, (X, Y) in [[True, ('A', 'C')], ]:
            if step[X] and step[Y]:
                # add reversible step between A and B
                surface_intermediates[X] = []
                surface_intermediates[Y] = []

                for x in [X, Y]:
                    print('x = {x}'.format(**locals()))
                    for intermediate in step[x]:
                        print(
                            'intermediate = {intermediate}'.format(**locals()))
                        if '_' in intermediate:
                            species, site = intermediate.split('_')
                            print(
                                'SPECIES {species}, SITE {site}, SITE_NAMES {site_names}'.format(**locals()))
                            print(
                                [s.startswith('{site}_'.format(**locals())) for s in site_names])
                            if any([s.startswith('{site}_'.format(**locals())) for s in site_names]):
                                surface_intermediates[
                                    x].append([species, site])
                        elif any([s.startswith('{intermediate}_'.format(**locals())) for s in site_names]):
                            surface_intermediates[x].append(
                                [EMPTY_SPECIES, intermediate])
                        else:
                            print('NOTHING MATCHED!!!')

                print('Elementary Rxn: {elementary_rxn}, Surface intermediates {surface_intermediates}'.format(
                    **locals()))

                # some validation checks to raise better error messages
                for letter_step, surface_intermediate in surface_intermediates.items():
                    if len(surface_intermediate) > 2:
                        raise UserWarning(
                            "Elementary reaction steps with more than two intermediates cannot be automatically translated into a kMC model: '{surface_intermediate}".format(**locals()))
                if len(surface_intermediates[X]) != len(surface_intermediates[Y]):
                    raise UserWarning(
                        "Number of surface different for elementary reaction: {surface_intermediates}".format(**locals()))

                for _, siteX in surface_intermediates[X]:
                    sitesY = [s for _, s in surface_intermediates[Y]]
                    if not siteX in sitesY:
                        raise UserWarning(
                            'Site {siteX} is mentioned on one side of the equation but not on the other: {surface_intermediates}'.format(**locals()))

                # first populate conditions and actions with one site
                print(surface_intermediates, X)
                condition_species, condition_site_name = surface_intermediates[
                    X][0]
                for cs_i, condition_site in enumerate([cs for cs in site_names if cs.startswith('{condition_site_name}_'.format(**locals()))]):
                    condition_coord = pt.layer_list.generate_coord(
                        '{condition_site}.(0, 0, 0)'.format(**locals()))
                    conditions = [Condition(species=condition_species.replace(
                        '-', '_').replace('-', '_'), coord=condition_coord)]

                    action_species, action_site = surface_intermediates[Y][0]
                    action_coord = pt.layer_list.generate_coord(
                        '{condition_site}.(0, 0, 0)'.format(**locals()))

                    actions = [
                        Action(species=action_species.replace('-', '_'), coord=action_coord)]

                    # if not action_site == condition_site:
                    # raise UserWarning(
                    #'Positions of sites seems to have changed: {surface_intermediates}.'.format(**locals()))

                    condition_string = '_n_'.join([species.replace(
                        '-', '_') + '_' + site for (species, site) in surface_intermediates[X]])
                    action_string = '_n_'.join([species.replace(
                        '-', '_') + '_' + site for (species, site) in surface_intermediates[Y]])

                    forward_name_root = '{condition_string}_2_{action_string}_{cs_i}'.format(
                        **locals())
                    reverse_name_root = '{action_string}_2_{condition_string}_{cs_i}'.format(
                        **locals())

                    # then create all second (auxiliary) sites which have
                    # the same nearest distance

                    # this is the cut-off with which positional equality is
                    # tested
                    # if a process is geometrically degenerate (more than one direction)
                    # it is important that distance is identical within the
                    # cut-off
                    dist_tol = 1.e-3

                    # geometrically complex case: two-sites

                    two_site_process = False
                    if len(surface_intermediates[X]) == 2:
                        two_site_process = True
                        other_condition_species, other_condition_site = surface_intermediates[
                            X][1]
                        other_action_species, other_action_site = surface_intermediates[
                            Y][1]

                        if not other_condition_site == other_action_site:
                            raise UserWarning(
                                'Positions of sites seem to have changed {surface_intermediates}'.format(**locals()))

                        auxiliary_coords = set(
                            pt.layer_list.generate_coord_set([2, 2, 2]))
                        # first drop the initial site and sites with other
                        # labels
                        for auxiliary_site in copy.copy(auxiliary_coords):
                            if np.linalg.norm(action_coord.pos - auxiliary_site.pos) < dist_tol:
                                auxiliary_coords.discard(auxiliary_site)

                            if not auxiliary_site.name.startswith('{other_condition_site}_'.format(**locals())):
                                auxiliary_coords.discard(auxiliary_site)

                        # find the shortest distance to one of the remaining
                        # sites
                        min_dist = min([
                            np.linalg.norm(aux_site.pos - condition_coord.pos)
                            for aux_site in auxiliary_coords
                        ])

                        # produce reaction steps with all auxiliary sites with-in dist_tol of that
                        # minimum distance
                        aux_counter = 0

                        # if the process is a diffusion process
                        # mark it as such in the rate constant
                        # so that we can later fine tune its rate-constant
                        # easier
                        if other_condition_species == EMPTY_SPECIES \
                           and action_species == EMPTY_SPECIES \
                           and condition_species == other_action_species:
                            diff_prefix = 'diff_'
                        else:
                            diff_prefix = ''


                        for auxiliary_coord in auxiliary_coords:
                            aux_dist = np.linalg.norm(
                                auxiliary_coord.pos - condition_coord.pos)
                            close = abs(aux_dist - min_dist) < dist_tol
                            if close:
                                aux_conditions = conditions + \
                                    [Condition(
                                        species=other_condition_species.replace('-', '_'), coord=auxiliary_coord)]
                                aux_actions = actions + \
                                    [Action(
                                        species=other_action_species.replace('-', '_'), coord=auxiliary_coord)]

                                process_name = '{forward_name_root}_{aux_counter}'.format(
                                    **locals())

                                process = pt.add_process(name=process_name,
                                               conditions=aux_conditions,
                                               actions=aux_actions,
                                               rate_constant='{diff_prefix}forward_{ri}'.format(
                                                   **locals()),
                                               tof_count={forward_name_root: 1})

                                if reversible:  # swap conditions and actions
                                    process_name = '{reverse_name_root}_{aux_counter}'.format(
                                        **locals())
                                    process_r = pt.add_process(name=process_name,
                                                   conditions=aux_actions,
                                                   actions=aux_conditions,
                                                   rate_constant='{diff_prefix}reverse_{ri}'.format(
                                                       **locals()),
                                                   tof_count={forward_name_root: -1})

                                aux_counter += 1

                    else:  # trivial case: single-site processes
                        diff_prefix = ''
                        process_name = forward_name_root
                        process = pt.add_process(name=process_name,
                                       conditions=conditions,
                                       actions=actions,
                                       tof_count={forward_name_root: 1},
                                       rate_constant='forward_{ri}'.format(**locals()))

                        if reversible:
                            # swap conditions and actions
                            process_name = reverse_name_root
                            process_r = pt.add_process(name=process_name,
                                           conditions=actions,
                                           actions=conditions,
                                           tof_count={forward_name_root: -1},
                                           rate_constant='reverse_{ri}'.format(**locals()))
                    if adsorbate_interaction > 0 :
                        # Collect the nearest-neighbor sites up to a certain cut-off
                        coordinate_set = pt.layer_list.generate_coord_set([5, 5, 1])
                        # regenerate all action coordinates via generation string to set the
                        # absolute position of the coord correctly

                        action_coords = [pt.layer_list.generate_coord(action.coord._get_genstring()) for action in process.action_list]
                        min_dists = [min(map(lambda x: np.linalg.norm(x.pos-c.pos), action_coords)) for c in coordinate_set]

                        coordinate_set = sorted(zip(min_dists, coordinate_set))
                        print(min_dists)
                        print(coordinate_set)
                        curr_dist = 0.
                        neighbor_shell = 0
                        n_neighbors = {}
                        for dist, coord in coordinate_set:
                            if abs(curr_dist - dist) > dist_tol:
                                curr_dist = dist
                                neighbor_shell += 1
                            n_neighbors.setdefault(neighbor_shell, []).append(coord)
                        pprint.pprint(n_neighbors)

                        interacting_coords = []
                        for i in range(adsorbate_interaction):
                            interacting_coords.extend(n_neighbors[i+1])

                        pprint.pprint(interacting_coords)

                        species_options = []
                        for interacting_coord in interacting_coords:
                            site_name = coord.name.split('_')[0]
                            species_options.append(site_species[site_name])
                        print(species_options)
                        species_sets = list(itertools.product(*species_options))

                        pprint.pprint(species_sets)

                        for i, species_set in enumerate(species_sets):
                            condition_list = process.condition_list
                            action_list = process.action_list
                            rate_constant = process.rate_constant
                            name = process.name
                            tof_count = process.tof_count

                            aux_conditions = [
                                kmos.types.ConditionAction(
                                    coord=coord,
                                    species=species,
                                ) for (species, coord) in zip(species_set, interacting_coords)
                            ]
                            pt.add_process(name='{name}_{i}'.format(**locals()),
                                           action_list=action_list,
                                           condition_list=condition_list + aux_conditions,
                                           rate_constant='{rate_constant}_{i}'.format(**locals()),
                                           tof_count=tof_count)

                        pt.process_list.remove(process)
                        # TODO: Add similar extra process for reverse reaction






                        # Remove the original process


                        # Add adsorbate interaction processes with all necessary combinations
                        # or nearest neighbor sites

                        # Add sensible markers to adsorption rates to calculate
                        # adsorbate interactions from DBMI model


                pt.add_parameter(name='{diff_prefix}forward_{ri}'.format(**locals()), value=1.)
                pt.add_parameter(name='{diff_prefix}reverse_{ri}'.format(**locals()), value=1.)




    return pt

if __name__ == '__main__':
    import optparse

    parser = optparse.OptionParser('Usage: %script *.mkm')

    options, args = parser.parse_args()

    if not len(args) >= 1:
        raise UserWarning('*.mkm filename expected!')

    model_filename = args[0]

    seed, _ = os.path.splitext(model_filename)

    catmap_model = ReactionModel(setup_file=model_filename)
    catmap_model.run()

    kmos_model = catmap2kmos(catmap_model,)

    # print(kmos_model)
    kmos_model.print_statistics()
    for suffix in ['xml']:
        kmos_model.save('translated_{seed}.{suffix}'.format(**locals()))
