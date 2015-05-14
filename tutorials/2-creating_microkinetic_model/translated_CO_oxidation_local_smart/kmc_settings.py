model_name = 'CatMAP_translated_model'
simulation_size = 20
random_seed = 1

def setup_model(model):
    """Write initialization steps here.
       e.g. ::
    model.put([0,0,0,model.lattice.default_a], model.proclist.species_a)
    """
    from setup_model import setup_model
    setup_model(model)

# Default history length in graph
hist_length = 30

parameters = {
    "diff_forward_3":{"value":"1.0", "adjustable":False, "min":"0.0", "max":"0.0","scale":"linear"},
    "diff_forward_4":{"value":"1.0", "adjustable":False, "min":"0.0", "max":"0.0","scale":"linear"},
    "diff_reverse_3":{"value":"1.0", "adjustable":False, "min":"0.0", "max":"0.0","scale":"linear"},
    "diff_reverse_4":{"value":"1.0", "adjustable":False, "min":"0.0", "max":"0.0","scale":"linear"},
    "forward_0":{"value":"1.0", "adjustable":False, "min":"0.0", "max":"0.0","scale":"linear"},
    "forward_1":{"value":"1.0", "adjustable":False, "min":"0.0", "max":"0.0","scale":"linear"},
    "forward_2":{"value":"1.0", "adjustable":False, "min":"0.0", "max":"0.0","scale":"linear"},
    "reverse_0":{"value":"1.0", "adjustable":False, "min":"0.0", "max":"0.0","scale":"linear"},
    "reverse_1":{"value":"1.0", "adjustable":False, "min":"0.0", "max":"0.0","scale":"linear"},
    "reverse_2":{"value":"1.0", "adjustable":False, "min":"0.0", "max":"0.0","scale":"linear"},
    }

rate_constants = {
    "CO_s_2_empty_s_0":("reverse_0", True),
    "CO_s_n_O_s_2_empty_s_n_empty_s_0_0":("forward_2", True),
    "CO_s_n_O_s_2_empty_s_n_empty_s_0_1":("forward_2", True),
    "CO_s_n_O_s_2_empty_s_n_empty_s_0_2":("forward_2", True),
    "CO_s_n_O_s_2_empty_s_n_empty_s_0_3":("forward_2", True),
    "CO_s_n_O_s_2_empty_s_n_empty_s_0_4":("forward_2", True),
    "CO_s_n_O_s_2_empty_s_n_empty_s_0_5":("forward_2", True),
    "CO_s_n_empty_s_2_empty_s_n_CO_s_0_0":("diff_forward_3", True),
    "CO_s_n_empty_s_2_empty_s_n_CO_s_0_1":("diff_forward_3", True),
    "CO_s_n_empty_s_2_empty_s_n_CO_s_0_2":("diff_forward_3", True),
    "CO_s_n_empty_s_2_empty_s_n_CO_s_0_3":("diff_forward_3", True),
    "CO_s_n_empty_s_2_empty_s_n_CO_s_0_4":("diff_forward_3", True),
    "CO_s_n_empty_s_2_empty_s_n_CO_s_0_5":("diff_forward_3", True),
    "O_s_n_O_s_2_empty_s_n_empty_s_0_0":("reverse_1", True),
    "O_s_n_O_s_2_empty_s_n_empty_s_0_1":("reverse_1", True),
    "O_s_n_O_s_2_empty_s_n_empty_s_0_2":("reverse_1", True),
    "O_s_n_O_s_2_empty_s_n_empty_s_0_3":("reverse_1", True),
    "O_s_n_O_s_2_empty_s_n_empty_s_0_4":("reverse_1", True),
    "O_s_n_O_s_2_empty_s_n_empty_s_0_5":("reverse_1", True),
    "O_s_n_empty_s_2_empty_s_n_O_s_0_0":("diff_forward_4", True),
    "O_s_n_empty_s_2_empty_s_n_O_s_0_1":("diff_forward_4", True),
    "O_s_n_empty_s_2_empty_s_n_O_s_0_2":("diff_forward_4", True),
    "O_s_n_empty_s_2_empty_s_n_O_s_0_3":("diff_forward_4", True),
    "O_s_n_empty_s_2_empty_s_n_O_s_0_4":("diff_forward_4", True),
    "O_s_n_empty_s_2_empty_s_n_O_s_0_5":("diff_forward_4", True),
    "empty_s_2_CO_s_0":("forward_0", True),
    "empty_s_n_CO_s_2_CO_s_n_empty_s_0_0":("diff_reverse_3", True),
    "empty_s_n_CO_s_2_CO_s_n_empty_s_0_1":("diff_reverse_3", True),
    "empty_s_n_CO_s_2_CO_s_n_empty_s_0_2":("diff_reverse_3", True),
    "empty_s_n_CO_s_2_CO_s_n_empty_s_0_3":("diff_reverse_3", True),
    "empty_s_n_CO_s_2_CO_s_n_empty_s_0_4":("diff_reverse_3", True),
    "empty_s_n_CO_s_2_CO_s_n_empty_s_0_5":("diff_reverse_3", True),
    "empty_s_n_O_s_2_O_s_n_empty_s_0_0":("diff_reverse_4", True),
    "empty_s_n_O_s_2_O_s_n_empty_s_0_1":("diff_reverse_4", True),
    "empty_s_n_O_s_2_O_s_n_empty_s_0_2":("diff_reverse_4", True),
    "empty_s_n_O_s_2_O_s_n_empty_s_0_3":("diff_reverse_4", True),
    "empty_s_n_O_s_2_O_s_n_empty_s_0_4":("diff_reverse_4", True),
    "empty_s_n_O_s_2_O_s_n_empty_s_0_5":("diff_reverse_4", True),
    "empty_s_n_empty_s_2_CO_s_n_O_s_0_0":("reverse_2", True),
    "empty_s_n_empty_s_2_CO_s_n_O_s_0_1":("reverse_2", True),
    "empty_s_n_empty_s_2_CO_s_n_O_s_0_2":("reverse_2", True),
    "empty_s_n_empty_s_2_CO_s_n_O_s_0_3":("reverse_2", True),
    "empty_s_n_empty_s_2_CO_s_n_O_s_0_4":("reverse_2", True),
    "empty_s_n_empty_s_2_CO_s_n_O_s_0_5":("reverse_2", True),
    "empty_s_n_empty_s_2_O_s_n_O_s_0_0":("forward_1", True),
    "empty_s_n_empty_s_2_O_s_n_O_s_0_1":("forward_1", True),
    "empty_s_n_empty_s_2_O_s_n_O_s_0_2":("forward_1", True),
    "empty_s_n_empty_s_2_O_s_n_O_s_0_3":("forward_1", True),
    "empty_s_n_empty_s_2_O_s_n_O_s_0_4":("forward_1", True),
    "empty_s_n_empty_s_2_O_s_n_O_s_0_5":("forward_1", True),
    }

site_names = ['default_s_0']
representations = {
    "CO":"""Atoms('CO', [[0., 0., 0.], [0., 0., 1.7]])""",
    "CO2":"""Atoms('OCO', [[1., 1., 0.], [0., 0., 0.,], [-1., 1., 0]])""",
    "O":"""Atoms('O')""",
    "O2":"""Atoms()""",
    "empty":"""""",
    }

lattice_representation = """[Atoms(symbols='H',
          pbc=np.array([False, False, False], dtype=bool),
          cell=np.array(
      [[  2.68700577,   0.        ,   0.        ],
       [  1.34350288,   2.32701526,   0.        ],
       [  0.        ,   0.        ,  26.58179307]]),
          scaled_positions=np.array(
      [[0.0, 0.0, -0.0037619734581734905]]),
),]"""

species_tags = {
    "CO":"""""",
    "CO2":"""""",
    "O":"""""",
    "O2":"""""",
    "empty":"""""",
    }

tof_count = {
    "CO_s_2_empty_s_0":{'empty_s_2_CO_s_0': -1},
    "CO_s_n_O_s_2_empty_s_n_empty_s_0_0":{'CO_s_n_O_s_2_empty_s_n_empty_s_0': 1},
    "CO_s_n_O_s_2_empty_s_n_empty_s_0_1":{'CO_s_n_O_s_2_empty_s_n_empty_s_0': 1},
    "CO_s_n_O_s_2_empty_s_n_empty_s_0_2":{'CO_s_n_O_s_2_empty_s_n_empty_s_0': 1},
    "CO_s_n_O_s_2_empty_s_n_empty_s_0_3":{'CO_s_n_O_s_2_empty_s_n_empty_s_0': 1},
    "CO_s_n_O_s_2_empty_s_n_empty_s_0_4":{'CO_s_n_O_s_2_empty_s_n_empty_s_0': 1},
    "CO_s_n_O_s_2_empty_s_n_empty_s_0_5":{'CO_s_n_O_s_2_empty_s_n_empty_s_0': 1},
    "CO_s_n_empty_s_2_empty_s_n_CO_s_0_0":{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': 1},
    "CO_s_n_empty_s_2_empty_s_n_CO_s_0_1":{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': 1},
    "CO_s_n_empty_s_2_empty_s_n_CO_s_0_2":{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': 1},
    "CO_s_n_empty_s_2_empty_s_n_CO_s_0_3":{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': 1},
    "CO_s_n_empty_s_2_empty_s_n_CO_s_0_4":{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': 1},
    "CO_s_n_empty_s_2_empty_s_n_CO_s_0_5":{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': 1},
    "O_s_n_O_s_2_empty_s_n_empty_s_0_0":{'empty_s_n_empty_s_2_O_s_n_O_s_0': -1},
    "O_s_n_O_s_2_empty_s_n_empty_s_0_1":{'empty_s_n_empty_s_2_O_s_n_O_s_0': -1},
    "O_s_n_O_s_2_empty_s_n_empty_s_0_2":{'empty_s_n_empty_s_2_O_s_n_O_s_0': -1},
    "O_s_n_O_s_2_empty_s_n_empty_s_0_3":{'empty_s_n_empty_s_2_O_s_n_O_s_0': -1},
    "O_s_n_O_s_2_empty_s_n_empty_s_0_4":{'empty_s_n_empty_s_2_O_s_n_O_s_0': -1},
    "O_s_n_O_s_2_empty_s_n_empty_s_0_5":{'empty_s_n_empty_s_2_O_s_n_O_s_0': -1},
    "O_s_n_empty_s_2_empty_s_n_O_s_0_0":{'O_s_n_empty_s_2_empty_s_n_O_s_0': 1},
    "O_s_n_empty_s_2_empty_s_n_O_s_0_1":{'O_s_n_empty_s_2_empty_s_n_O_s_0': 1},
    "O_s_n_empty_s_2_empty_s_n_O_s_0_2":{'O_s_n_empty_s_2_empty_s_n_O_s_0': 1},
    "O_s_n_empty_s_2_empty_s_n_O_s_0_3":{'O_s_n_empty_s_2_empty_s_n_O_s_0': 1},
    "O_s_n_empty_s_2_empty_s_n_O_s_0_4":{'O_s_n_empty_s_2_empty_s_n_O_s_0': 1},
    "O_s_n_empty_s_2_empty_s_n_O_s_0_5":{'O_s_n_empty_s_2_empty_s_n_O_s_0': 1},
    "empty_s_2_CO_s_0":{'empty_s_2_CO_s_0': 1},
    "empty_s_n_CO_s_2_CO_s_n_empty_s_0_0":{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': -1},
    "empty_s_n_CO_s_2_CO_s_n_empty_s_0_1":{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': -1},
    "empty_s_n_CO_s_2_CO_s_n_empty_s_0_2":{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': -1},
    "empty_s_n_CO_s_2_CO_s_n_empty_s_0_3":{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': -1},
    "empty_s_n_CO_s_2_CO_s_n_empty_s_0_4":{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': -1},
    "empty_s_n_CO_s_2_CO_s_n_empty_s_0_5":{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': -1},
    "empty_s_n_O_s_2_O_s_n_empty_s_0_0":{'O_s_n_empty_s_2_empty_s_n_O_s_0': -1},
    "empty_s_n_O_s_2_O_s_n_empty_s_0_1":{'O_s_n_empty_s_2_empty_s_n_O_s_0': -1},
    "empty_s_n_O_s_2_O_s_n_empty_s_0_2":{'O_s_n_empty_s_2_empty_s_n_O_s_0': -1},
    "empty_s_n_O_s_2_O_s_n_empty_s_0_3":{'O_s_n_empty_s_2_empty_s_n_O_s_0': -1},
    "empty_s_n_O_s_2_O_s_n_empty_s_0_4":{'O_s_n_empty_s_2_empty_s_n_O_s_0': -1},
    "empty_s_n_O_s_2_O_s_n_empty_s_0_5":{'O_s_n_empty_s_2_empty_s_n_O_s_0': -1},
    "empty_s_n_empty_s_2_CO_s_n_O_s_0_0":{'CO_s_n_O_s_2_empty_s_n_empty_s_0': -1},
    "empty_s_n_empty_s_2_CO_s_n_O_s_0_1":{'CO_s_n_O_s_2_empty_s_n_empty_s_0': -1},
    "empty_s_n_empty_s_2_CO_s_n_O_s_0_2":{'CO_s_n_O_s_2_empty_s_n_empty_s_0': -1},
    "empty_s_n_empty_s_2_CO_s_n_O_s_0_3":{'CO_s_n_O_s_2_empty_s_n_empty_s_0': -1},
    "empty_s_n_empty_s_2_CO_s_n_O_s_0_4":{'CO_s_n_O_s_2_empty_s_n_empty_s_0': -1},
    "empty_s_n_empty_s_2_CO_s_n_O_s_0_5":{'CO_s_n_O_s_2_empty_s_n_empty_s_0': -1},
    "empty_s_n_empty_s_2_O_s_n_O_s_0_0":{'empty_s_n_empty_s_2_O_s_n_O_s_0': 1},
    "empty_s_n_empty_s_2_O_s_n_O_s_0_1":{'empty_s_n_empty_s_2_O_s_n_O_s_0': 1},
    "empty_s_n_empty_s_2_O_s_n_O_s_0_2":{'empty_s_n_empty_s_2_O_s_n_O_s_0': 1},
    "empty_s_n_empty_s_2_O_s_n_O_s_0_3":{'empty_s_n_empty_s_2_O_s_n_O_s_0': 1},
    "empty_s_n_empty_s_2_O_s_n_O_s_0_4":{'empty_s_n_empty_s_2_O_s_n_O_s_0': 1},
    "empty_s_n_empty_s_2_O_s_n_O_s_0_5":{'empty_s_n_empty_s_2_O_s_n_O_s_0': 1},
    }

xml = """<?xml version="1.0" ?>
<kmc version="(0, 2)">
    <meta author="CatMAP USER" debug="0" email="mkm-developers-request@lists.stanford.edu" model_dimension="2" model_name="CatMAP_translated_model"/>
    <species_list default_species="empty">
        <species color="#6029FD" name="CO" representation="Atoms('CO', [[0., 0., 0.], [0., 0., 1.7]])" tags=""/>
        <species color="#5BB41B" name="CO2" representation="Atoms('OCO', [[1., 1., 0.], [0., 0., 0.,], [-1., 1., 0]])" tags=""/>
        <species color="#B59FC2" name="O" representation="Atoms('O')" tags=""/>
        <species color="#FBCAA9" name="O2" representation="Atoms()" tags=""/>
        <species color="#ffffff" name="empty" representation="" tags=""/>
    </species_list>
    <parameter_list>
        <parameter adjustable="False" max="0.0" min="0.0" name="diff_forward_3" scale="linear" value="1.0"/>
        <parameter adjustable="False" max="0.0" min="0.0" name="diff_forward_4" scale="linear" value="1.0"/>
        <parameter adjustable="False" max="0.0" min="0.0" name="diff_reverse_3" scale="linear" value="1.0"/>
        <parameter adjustable="False" max="0.0" min="0.0" name="diff_reverse_4" scale="linear" value="1.0"/>
        <parameter adjustable="False" max="0.0" min="0.0" name="forward_0" scale="linear" value="1.0"/>
        <parameter adjustable="False" max="0.0" min="0.0" name="forward_1" scale="linear" value="1.0"/>
        <parameter adjustable="False" max="0.0" min="0.0" name="forward_2" scale="linear" value="1.0"/>
        <parameter adjustable="False" max="0.0" min="0.0" name="reverse_0" scale="linear" value="1.0"/>
        <parameter adjustable="False" max="0.0" min="0.0" name="reverse_1" scale="linear" value="1.0"/>
        <parameter adjustable="False" max="0.0" min="0.0" name="reverse_2" scale="linear" value="1.0"/>
    </parameter_list>
    <lattice cell_size="2.68700577 0.0 0.0 1.34350288 2.32701526 0.0 0.0 0.0 26.58179307" default_layer="default" representation="[Atoms(symbols='H',
          pbc=np.array([False, False, False], dtype=bool),
          cell=np.array(
      [[  2.68700577,   0.        ,   0.        ],
       [  1.34350288,   2.32701526,   0.        ],
       [  0.        ,   0.        ,  26.58179307]]),
          scaled_positions=np.array(
      [[0.0, 0.0, -0.0037619734581734905]]),
),]" substrate_layer="default">
        <layer color="#ffffff" name="default">
            <site default_species="default_species" pos="0.0 0.0 0.0" tags="" type="s_0"/>
        </layer>
    </lattice>
    <process_list>
        <process enabled="True" name="CO_s_2_empty_s_0" rate_constant="reverse_0" tof_count="{'empty_s_2_CO_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
        </process>
        <process enabled="True" name="CO_s_n_O_s_2_empty_s_n_empty_s_0_0" rate_constant="forward_2" tof_count="{'CO_s_n_O_s_2_empty_s_n_empty_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="empty"/>
        </process>
        <process enabled="True" name="CO_s_n_O_s_2_empty_s_n_empty_s_0_1" rate_constant="forward_2" tof_count="{'CO_s_n_O_s_2_empty_s_n_empty_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="empty"/>
        </process>
        <process enabled="True" name="CO_s_n_O_s_2_empty_s_n_empty_s_0_2" rate_constant="forward_2" tof_count="{'CO_s_n_O_s_2_empty_s_n_empty_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="empty"/>
        </process>
        <process enabled="True" name="CO_s_n_O_s_2_empty_s_n_empty_s_0_3" rate_constant="forward_2" tof_count="{'CO_s_n_O_s_2_empty_s_n_empty_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="empty"/>
        </process>
        <process enabled="True" name="CO_s_n_O_s_2_empty_s_n_empty_s_0_4" rate_constant="forward_2" tof_count="{'CO_s_n_O_s_2_empty_s_n_empty_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="empty"/>
        </process>
        <process enabled="True" name="CO_s_n_O_s_2_empty_s_n_empty_s_0_5" rate_constant="forward_2" tof_count="{'CO_s_n_O_s_2_empty_s_n_empty_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="empty"/>
        </process>
        <process enabled="True" name="CO_s_n_empty_s_2_empty_s_n_CO_s_0_0" rate_constant="diff_forward_3" tof_count="{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="CO"/>
        </process>
        <process enabled="True" name="CO_s_n_empty_s_2_empty_s_n_CO_s_0_1" rate_constant="diff_forward_3" tof_count="{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="CO"/>
        </process>
        <process enabled="True" name="CO_s_n_empty_s_2_empty_s_n_CO_s_0_2" rate_constant="diff_forward_3" tof_count="{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="CO"/>
        </process>
        <process enabled="True" name="CO_s_n_empty_s_2_empty_s_n_CO_s_0_3" rate_constant="diff_forward_3" tof_count="{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="CO"/>
        </process>
        <process enabled="True" name="CO_s_n_empty_s_2_empty_s_n_CO_s_0_4" rate_constant="diff_forward_3" tof_count="{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="CO"/>
        </process>
        <process enabled="True" name="CO_s_n_empty_s_2_empty_s_n_CO_s_0_5" rate_constant="diff_forward_3" tof_count="{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="CO"/>
        </process>
        <process enabled="True" name="O_s_n_O_s_2_empty_s_n_empty_s_0_0" rate_constant="reverse_1" tof_count="{'empty_s_n_empty_s_2_O_s_n_O_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="empty"/>
        </process>
        <process enabled="True" name="O_s_n_O_s_2_empty_s_n_empty_s_0_1" rate_constant="reverse_1" tof_count="{'empty_s_n_empty_s_2_O_s_n_O_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="empty"/>
        </process>
        <process enabled="True" name="O_s_n_O_s_2_empty_s_n_empty_s_0_2" rate_constant="reverse_1" tof_count="{'empty_s_n_empty_s_2_O_s_n_O_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="empty"/>
        </process>
        <process enabled="True" name="O_s_n_O_s_2_empty_s_n_empty_s_0_3" rate_constant="reverse_1" tof_count="{'empty_s_n_empty_s_2_O_s_n_O_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="empty"/>
        </process>
        <process enabled="True" name="O_s_n_O_s_2_empty_s_n_empty_s_0_4" rate_constant="reverse_1" tof_count="{'empty_s_n_empty_s_2_O_s_n_O_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="empty"/>
        </process>
        <process enabled="True" name="O_s_n_O_s_2_empty_s_n_empty_s_0_5" rate_constant="reverse_1" tof_count="{'empty_s_n_empty_s_2_O_s_n_O_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="empty"/>
        </process>
        <process enabled="True" name="O_s_n_empty_s_2_empty_s_n_O_s_0_0" rate_constant="diff_forward_4" tof_count="{'O_s_n_empty_s_2_empty_s_n_O_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="O"/>
        </process>
        <process enabled="True" name="O_s_n_empty_s_2_empty_s_n_O_s_0_1" rate_constant="diff_forward_4" tof_count="{'O_s_n_empty_s_2_empty_s_n_O_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="O"/>
        </process>
        <process enabled="True" name="O_s_n_empty_s_2_empty_s_n_O_s_0_2" rate_constant="diff_forward_4" tof_count="{'O_s_n_empty_s_2_empty_s_n_O_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="O"/>
        </process>
        <process enabled="True" name="O_s_n_empty_s_2_empty_s_n_O_s_0_3" rate_constant="diff_forward_4" tof_count="{'O_s_n_empty_s_2_empty_s_n_O_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="O"/>
        </process>
        <process enabled="True" name="O_s_n_empty_s_2_empty_s_n_O_s_0_4" rate_constant="diff_forward_4" tof_count="{'O_s_n_empty_s_2_empty_s_n_O_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="O"/>
        </process>
        <process enabled="True" name="O_s_n_empty_s_2_empty_s_n_O_s_0_5" rate_constant="diff_forward_4" tof_count="{'O_s_n_empty_s_2_empty_s_n_O_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="O"/>
        </process>
        <process enabled="True" name="empty_s_2_CO_s_0" rate_constant="forward_0" tof_count="{'empty_s_2_CO_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
        </process>
        <process enabled="True" name="empty_s_n_CO_s_2_CO_s_n_empty_s_0_0" rate_constant="diff_reverse_3" tof_count="{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="empty"/>
        </process>
        <process enabled="True" name="empty_s_n_CO_s_2_CO_s_n_empty_s_0_1" rate_constant="diff_reverse_3" tof_count="{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="empty"/>
        </process>
        <process enabled="True" name="empty_s_n_CO_s_2_CO_s_n_empty_s_0_2" rate_constant="diff_reverse_3" tof_count="{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="empty"/>
        </process>
        <process enabled="True" name="empty_s_n_CO_s_2_CO_s_n_empty_s_0_3" rate_constant="diff_reverse_3" tof_count="{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="empty"/>
        </process>
        <process enabled="True" name="empty_s_n_CO_s_2_CO_s_n_empty_s_0_4" rate_constant="diff_reverse_3" tof_count="{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="empty"/>
        </process>
        <process enabled="True" name="empty_s_n_CO_s_2_CO_s_n_empty_s_0_5" rate_constant="diff_reverse_3" tof_count="{'CO_s_n_empty_s_2_empty_s_n_CO_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="empty"/>
        </process>
        <process enabled="True" name="empty_s_n_O_s_2_O_s_n_empty_s_0_0" rate_constant="diff_reverse_4" tof_count="{'O_s_n_empty_s_2_empty_s_n_O_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="empty"/>
        </process>
        <process enabled="True" name="empty_s_n_O_s_2_O_s_n_empty_s_0_1" rate_constant="diff_reverse_4" tof_count="{'O_s_n_empty_s_2_empty_s_n_O_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="empty"/>
        </process>
        <process enabled="True" name="empty_s_n_O_s_2_O_s_n_empty_s_0_2" rate_constant="diff_reverse_4" tof_count="{'O_s_n_empty_s_2_empty_s_n_O_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="empty"/>
        </process>
        <process enabled="True" name="empty_s_n_O_s_2_O_s_n_empty_s_0_3" rate_constant="diff_reverse_4" tof_count="{'O_s_n_empty_s_2_empty_s_n_O_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="empty"/>
        </process>
        <process enabled="True" name="empty_s_n_O_s_2_O_s_n_empty_s_0_4" rate_constant="diff_reverse_4" tof_count="{'O_s_n_empty_s_2_empty_s_n_O_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="empty"/>
        </process>
        <process enabled="True" name="empty_s_n_O_s_2_O_s_n_empty_s_0_5" rate_constant="diff_reverse_4" tof_count="{'O_s_n_empty_s_2_empty_s_n_O_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="empty"/>
        </process>
        <process enabled="True" name="empty_s_n_empty_s_2_CO_s_n_O_s_0_0" rate_constant="reverse_2" tof_count="{'CO_s_n_O_s_2_empty_s_n_empty_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="O"/>
        </process>
        <process enabled="True" name="empty_s_n_empty_s_2_CO_s_n_O_s_0_1" rate_constant="reverse_2" tof_count="{'CO_s_n_O_s_2_empty_s_n_empty_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="O"/>
        </process>
        <process enabled="True" name="empty_s_n_empty_s_2_CO_s_n_O_s_0_2" rate_constant="reverse_2" tof_count="{'CO_s_n_O_s_2_empty_s_n_empty_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="O"/>
        </process>
        <process enabled="True" name="empty_s_n_empty_s_2_CO_s_n_O_s_0_3" rate_constant="reverse_2" tof_count="{'CO_s_n_O_s_2_empty_s_n_empty_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="O"/>
        </process>
        <process enabled="True" name="empty_s_n_empty_s_2_CO_s_n_O_s_0_4" rate_constant="reverse_2" tof_count="{'CO_s_n_O_s_2_empty_s_n_empty_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="O"/>
        </process>
        <process enabled="True" name="empty_s_n_empty_s_2_CO_s_n_O_s_0_5" rate_constant="reverse_2" tof_count="{'CO_s_n_O_s_2_empty_s_n_empty_s_0': -1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="CO"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="O"/>
        </process>
        <process enabled="True" name="empty_s_n_empty_s_2_O_s_n_O_s_0_0" rate_constant="forward_1" tof_count="{'empty_s_n_empty_s_2_O_s_n_O_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 0 0" species="O"/>
        </process>
        <process enabled="True" name="empty_s_n_empty_s_2_O_s_n_O_s_0_1" rate_constant="forward_1" tof_count="{'empty_s_n_empty_s_2_O_s_n_O_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="1 -1 0" species="O"/>
        </process>
        <process enabled="True" name="empty_s_n_empty_s_2_O_s_n_O_s_0_2" rate_constant="forward_1" tof_count="{'empty_s_n_empty_s_2_O_s_n_O_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 -1 0" species="O"/>
        </process>
        <process enabled="True" name="empty_s_n_empty_s_2_O_s_n_O_s_0_3" rate_constant="forward_1" tof_count="{'empty_s_n_empty_s_2_O_s_n_O_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 0 0" species="O"/>
        </process>
        <process enabled="True" name="empty_s_n_empty_s_2_O_s_n_O_s_0_4" rate_constant="forward_1" tof_count="{'empty_s_n_empty_s_2_O_s_n_O_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="-1 1 0" species="O"/>
        </process>
        <process enabled="True" name="empty_s_n_empty_s_2_O_s_n_O_s_0_5" rate_constant="forward_1" tof_count="{'empty_s_n_empty_s_2_O_s_n_O_s_0': 1}">
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="empty"/>
            <condition coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="empty"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 0 0" species="O"/>
            <action coord_layer="default" coord_name="s_0" coord_offset="0 1 0" species="O"/>
        </process>
    </process_list>
    <output_list/>
</kmc>
"""
