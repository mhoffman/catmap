#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import itertools

import itertools

def itertools_product_no_repetition(*vectors):
    n = len(vectors)
    for out in itertools.product(*vectors):
        if len(set(out)) == n:
            yield out

def sum_edge_length_metric(sites):
    import itertools
    import numpy.linalg
    return sum([
        numpy.linalg.norm(site_A.pos - site_B.pos)**2
        for (site_A, site_B) in itertools.combinations(sites, 2)
    ])

def sum_edge_length_squared_metric(sites):
    import itertools
    import numpy.linalg
    return sum([
        numpy.linalg.norm(site_A.pos - site_B.pos)**2
        for (site_A, site_B) in itertools.combinations(sites, 2)
    ])

def get_color(string):
 """Generate a color from any string using the hexdigest of
 the md5 hash
 """
 import md5
 return '#{}'.format(md5.md5(string).hexdigest()[:6])

def translate_model_file(mkm_filename, options):
    import catmap
    import os.path
    seed, _ = os.path.splitext(mkm_filename)
    catmap_model = catmap.ReactionModel(setup_file=mkm_filename)
    catmap_model.run() # Running the model once is needed to initialize all model values
    kmos_model = catmap2kmos(catmap_model,
                 options=options,
                 model_name=seed,
                 )
    kmos_model.print_statistics()
    kmos_model.save('{seed}_kmc.ini'.format(**locals()))



def catmap2kmos(cm_model,
                unit_cell=None,
                site_positions=None,
                diffusion_barriers=None,
                species_representations=None,
                surface_representation='None',
                model_name='CatMAP_translated_model',
                options=None,
                ):
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

    import copy
    import pprint

    import numpy as np

    from kmos.types import Project, Site, Condition, Action
    import kmos.utils
    import ase.atoms
    from ase.atoms import Atoms

    EMPTY_SPECIES = 'empty'

    # initialize model
    pt = Project()

    # add meta information
    # TODO : should we add corresponding information
    #        in CatMAP ?
    pt.set_meta(author=options.author_name,
                email=options.author_email,
                model_name=model_name,
                model_dimension='2',   # default ...
                debug=0,  # let's be pessimistic
                )

    # add unit cell
    layer = pt.add_layer(name='default')

    if unit_cell is None:
        unit_cell = np.diag([1., 1., 1.]).tolist()

    # Need dummy atom (here 'H') so that ase.atoms.Atoms doesn't puke further
    # down
    if hasattr(cm_model, 'background_representation') and hasattr(cm_model, 'unit_cell'):
        print("Warning: unit cell in background_representation will override the given unit cell {unit_cell}".format(**locals()))
    background_representation = getattr(cm_model,
                                        'background_representation',
                                        'Atoms("H", [[0., 0., -0.1]], cell={unit_cell})'.format(**locals())
    )
    pt.layer_list.representation = '[{background_representation}]'.format(**locals())


    # add site positions
    for site_name in sorted(cm_model.site_names):
        if not site_name == 'g':
            print(cm_model.site_positions)
            for i, site_position in enumerate(sorted(cm_model.site_positions[site_name])):
                layer.sites.append(Site(name='{site_name}_{i}'.format(**locals()),
                                        pos=site_position))

    Coord = pt.layer_list.generate_coord

    # add species


    pt.add_species(name=EMPTY_SPECIES, color='#ffffff')
    species_names = set()
    for species_definition in sorted(cm_model.species_definitions.keys()):

        print('SPECIES DEFINITION {species_definition}'.format(**locals()))
        if '_' in species_definition:
            species_name, site = species_definition.split('_')
            species_name = species_name.replace('-', '_')


            species_representation = getattr(cm_model, 'species_representation', {}).get(species_name, None)
            if species_representation is None:
                species_representation = "Atoms()"
                color = get_color(species_name)
            elif type(species_representation) is str:
                species_representation = species_representation
                color = get_color(species_name)
            elif type(species_representation) is dict:
                species_representation, color = species_representation['geometry'], species_representation['color']
            else:
                raise UserWarning(("Specification for species representation {species_representation} not understood!\n\n\n"
                                   "The species representation should a string in the form of an ASE contructor (\"Atoms(...)\")\n"
                                   "or a dictionary with keys 'geometry' and 'color' where the geometry hold the ASE constructor\n"
                                   "and the color is a HTML hex string representing the color for the 'kmos edit' GUI\n"
                                   "like {{'color': '#ff0000', 'geometry': 'Atoms(\"CO\", ... )}}"
                ).format(**locals()))

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
    for ri, elementary_rxn in enumerate(sorted(cm_model.elementary_rxns)):
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


        # Now let's try to bring terms into the same order
        # Note: if two sites have the same name, they
        # are assumed to appear in the same order
        # This latter case is important to describe
        # diffusion processes between equal sites.
        term_order_dict = {}
        for i, term in enumerate(step['A']):
            site = term.split('_')[-1]
            term_order_dict.setdefault(site, []).append(i)
        step_C_tmp = {}
        for i, term in enumerate(step['C']):
            site = term.split('_')[-1]
            if site in term_order_dict:
                step_C_tmp[term_order_dict[site].pop(0)] = term
            else:
                # special case e.g.  A_a + B_b -> *_a + *_b + AB_g
                # special treatment for gas phase 'site' if site does show up on left side of reaction
                step_C_tmp[len(step['A']) + i] = term
        step['C'] = [_x[1] for _x in sorted(step_C_tmp.items())]

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
                if len(surface_intermediates[X]) != len(surface_intermediates[Y]):
                    raise UserWarning(
                        "Number of surface sites different for elementary reaction: {elementary_rxn} equiv. to {surface_intermediates}".format(**locals()))

                for _, siteX in surface_intermediates[X]:
                    sitesY = [s for _, s in surface_intermediates[Y]]
                    if not siteX in sitesY:
                        raise UserWarning(
                            'Site {siteX} is mentioned on one side of the equation but not on the other: {surface_intermediates}'.format(**locals()))

                for letter_step, surface_intermediate in surface_intermediates.items():
                    n = len(surface_intermediate)
                    print('State {letter_step} Number of surface intermediates {n}: {surface_intermediate}'.format(**locals()))
                    if not letter_step == 'A': # only generat sites set from initial step
                        continue
                    sites_vectors = []
                    for i, (_ads, site_name) in enumerate(surface_intermediate):
                        if i == 0:
                            sites_vectors.append(pt.layer_list.generate_coord_set(size=[1, 1, 1], site_name=site_name))
                        else:
                            sites_vectors.append(pt.layer_list.generate_coord_set(size=[2, 2, 2], site_name=site_name))

                sites_list = itertools_product_no_repetition(*sites_vectors)

                dist_tol =  1e-3
                metric_site_list = sorted(map(lambda x: (sum_edge_length_metric(x), x), sites_list))
                min_dist = metric_site_list[0][0]
                #sites_list = [sites for (dist, sites) in metric_site_list ]
                sites_list = [sites for (dist, sites) in metric_site_list if np.abs(dist - min_dist) < dist_tol]
                n_sets = len(sites_list)
                print("N {n_sets} Sites list {sites_list} Min dist {min_dist}".format(**locals()))

                for s_i, sites in enumerate(sorted(sites_list)):
                    ads_initial = [ads for (ads, site) in surface_intermediates['A']]
                    ads_final = [ads for (ads, site) in surface_intermediates['C']]

                    condition_string = '_n_'.join([species.replace( '-', '_') + '_' + site for (species, site) in surface_intermediates[X]])
                    action_string = '_n_'.join([species.replace( '-', '_') + '_' + site for (species, site) in surface_intermediates[Y]])

                    forward_name_root = '{condition_string}_2_{action_string}'.format( **locals())
                    reverse_name_root = '{action_string}_2_{condition_string}'.format( **locals())

                    actions = [Action(species=species, coord=coord) for (species, coord) in zip(ads_final, sites)]
                    conditions = [Condition(species=species, coord=coord) for (species, coord) in zip(ads_initial, sites)]
                    diff_prefix = ''

                    # if the process is a diffusion process
                    # mark it as such in the rate constant
                    # so that we can later fine tune its rate-constant
                    # easier
                    if len(sites_vectors) == 2 :
                        si = surface_intermediates
                        if ((si['A'][0][0] == EMPTY_SPECIES) and (si['C'][1][0] == EMPTY_SPECIES) and (si['A'][1][0] == si['C'][0][0])) \
                        or ((si['A'][1][0] == EMPTY_SPECIES) and (si['C'][0][0] == EMPTY_SPECIES) and (si['A'][0][0] == si['C'][1][0])) :
                            diff_prefix = 'diff_'
                        else:
                            diff_prefix = ''

                    process = pt.add_process(name='{forward_name_root}_{s_i}'.format(**locals()),
                                   conditions=conditions,
                                   actions=actions,
                                   rate_constant='{diff_prefix}forward_{ri}'.format(
                                       **locals()),
                                   tof_count={forward_name_root: 1})

                    reverse_process = pt.add_process(name='{reverse_name_root}_{s_i}'.format(**locals()),
                                   conditions=actions,
                                   actions=conditions,
                                   rate_constant='{diff_prefix}reverse_{ri}'.format(
                                       **locals()),
                                   tof_count={forward_name_root: -1})

                    if options.interaction > 0 :
                        r = 1 * options.interaction + 1
                        # Collect the nearest-neighbor sites up to a certain cut-off
                        coordinate_set = pt.layer_list.generate_coord_set([r, r, 1])
                        # regenerate all action coordinates via generation string to set the
                        # absolute position of the coord correctly

                        action_coords = [pt.layer_list.generate_coord(action.coord._get_genstring()) for action in process.action_list]
                        min_dists = [min(map(lambda x: np.linalg.norm(x.pos-c.pos), action_coords)) for c in coordinate_set]

                        coordinate_set = sorted(zip(min_dists, coordinate_set))

                        curr_dist = 0.
                        neighbor_shell = 0
                        n_neighbors = {}
                        for dist, coord in coordinate_set:
                            if abs(curr_dist - dist) > dist_tol:
                                curr_dist = dist
                                neighbor_shell += 1
                            n_neighbors.setdefault(neighbor_shell, []).append(coord)
                        print(process.name)
                        print("==========> Neighbor Shells <=============")
                        pprint.pprint(n_neighbors)

                        interacting_coords = []
                        for i in range(options.interaction):
                            interacting_coords.extend(n_neighbors[i+1])

                        pprint.pprint(interacting_coords)

                        species_options = []
                        for interacting_coord in interacting_coords:
                            site_name = coord.name.split('_')[0]
                            species_options.append(site_species[site_name])
                        print(species_options)

                        bystander_list = []
                        for allowed_species, interacting_coord in zip(species_options, interacting_coords):
                            _X, _Y, _ = interacting_coord.offset
                            flag = '{interacting_coord.name}_{_X}_{_Y}'.format(**locals())
                            flag = flag.replace('-', 'm')
                            try:
                                bystander_list.append(kmos.types.Bystander(
                                    coord=interacting_coord,
                                    allowed_species=allowed_species,
                                    flag=flag
                                ))
                            except AttributeError as e:
                                raise type(e)(('{e.message}: adsorbate-adsorbate interaction using bystanders not support in this kmos version.'
                                               'Please try a again from a branch that supported the otf backend ("kmos export -botf ....").'
                                ).format(**locals()))

                        print('Process {process} Bystanders {bystander_list}'.format(**locals()))
                        process.bystander_list = bystander_list
                    reverse_process.bystander_list = bystander_list

            pt.add_parameter(name='{diff_prefix}forward_{ri}'.format(**locals()), value=1.)
            pt.add_parameter(name='{diff_prefix}reverse_{ri}'.format(**locals()), value=1.)



    return pt
