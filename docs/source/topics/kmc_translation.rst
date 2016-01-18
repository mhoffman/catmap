

.. warning::
    This feature is highly experimental and this documentation intendend
    for early-adsopters as the feature is evolving.

    Main limitations include but are not limited to:
    
        * limited support of adsorbate-adsorbate interactions
        * significantly manual testing for convergence and steady-state required


.. note::
    For users on the Stanford Sherlock Cluster a light-weight version can
    be used after adding the following two lines to your `~/.bashrc` ::

        export PYTHONPATH=/home/maxjh/env/kmos/lib/python2.7/site-packages:${PYTHONPATH}
        export PATH=/home/vossj/suncat/bin:/home/maxjh/env/kmos/bin/${PATH}

Evaluation of a CatMAP model using kinetic Monte Carlo
==========================================================



Short Overview
-------------------------------------------

In order to execute a given CatMAP model as a kMC model you have to do the following steps

* add aditional geometrical information to your *\*.mkm* file
* translate the *\*.mkm* file to a kMC model using `catmap to_kmc *.mkm`
* compile the resulting *translated_\*.mkm* file using `kmos translate_*.mkm`
* evaluate the compiled kMC model using `catmap run_kmc` to generate a kMC evaluation similar to the catmap job
  (beware of various options in this step you will find be consulting `catmap help run_kmc`)


Long Version
-------------------------------------------

In order to evaluate your CatMAP model as a kMC model we use kmos as dependency. Check out it's
`documentation <http://kmos.readthedocs.org/en/latest/tutorials/index.html#installation`_ for how
to install it.


Adapt the \*.mkm file
-------------------------------------------

Second we will have to add minimal geometrical information to the *mean-field* CatMAP model in order
to enable the kMC translation feature figure out what nearest-neighbor site and such means.
The information we need to included is the *unit cell* of our surface model and the positions
of the *adsorption sites*.

It is recommended that you add this information at the end of the \*.mkm file. Normally specifing the
CatMAP requires "no programming", however the \*.mkm file is interpreted by the Python interpreter,
so nothing keeps us from adding a little of Python code here. No worries, this shouldn't hurt.

If you have a ASE \*.traj file of your empty surface somewhere now would be a good time to copy
that into the same directory as your \*.mkm file. Let's pretend this file is called `slab.py`
and it is stored somewhere on a computer `cluster` under `/path/to`. Then you would `cd`
to the \*.mkm file and run ::

    scp <username>@cluster:/path/to/slab.traj .

Once this done we'll use a bit of Python to include it in the \*.mkm file. So open your \*.mkm file
and add the end  ::

    .
    .
    .
    ### kMC specific settings
    import ase.io
    import kmos.utils

    background_representaton = kmos.utils.get_ase_constructor(ase.io.read("slab.traj"))


This will achieve to things at the same time:
(1) the translator will extract the unit cell from the ASE Atoms object which is important to extract nearest neighbor
information. Consider for example the hollow sites on a (111) slab, a (100) slab and a (110) slab: we
have 3 nearest neighbor sites, 4 nearest neighbor sites, and (2) neighbor sites respectively.

The second piece of crucial information are the actual positions of the surface sites as those cannot
be reliably extracted from the surface slab geometry alone. To this end you will have to define
a dictionary `site_positions = {...}` which contains the names of the sites as keys and the
positions of the site as a list of lists in fractional coordinates. Note: this is a list
of lists since a surface could have several identical sites. E.g. the `fcc(100)` has typically
two bridge sites.

The names of the sites should correspond to those used in the reaction definition. If you haven't
used any explicit site names, catmap defaults to `s`. Thus a minimal setting could be ::

    ### kMC specific settings
    .
    .
    .
    site_positions = {
        's': [[0., 0., 0.]
              ],
    }

It should be straightforward to expand this example to more complex surfaces. For instance
for a `fcc(100)` surface the site settings could look like ::

    ### kMC specific settings
    .
    .
    .
    site_positions = {
        'bridge': [
          [.0, .5., .5], [.5, .0, .5]
          ],
        'hollow': [
          [.5, .5, .5]
    ]}


There are is another key that allows to set the geometries of individual molecules. This is mostly
if you want to make really spiffy visualizations. Just to give an example, this could look like ::

    ### kMC specific settings
    .
    .
    .
    species_representation = {
        'O': "Atoms('O')",
        'CO': "Atoms('CO', [[0., 0., 0.], [0., 0., 1.7]])",
        'CO2': "Atoms('OCO', [[1., 1., 0.], [0., 0., 0.,], [-1., 1., 0]])",
    }

And of course you might have to adapt this to the species in your system. But we will not worry about this for now.


Translate to model to kmc
-------------------------------------------

Having such an adapted \*.mkm file we can go on an attempt to translate it into a kMC model. Let's pretend the
\*.mkm file is called `model.mkm`. You should then run::

    catmap to_kmc model.mkm

Depending on the model this can take a while. The translation step will also evaluat the model as a CatMAP model
to ensure that all variables are set and the rate constants are evaluated as stored in the \*.pkl file.
The kMC model will use the values of those exact rate constants in its run later on. If this step
crashes or does not seem to finish try with a smaller first. Please feel free to file a bug
report if your model does not translate.

Once this step is done you should see a file named `model_kmc.ini` in the same directory. You can have either
open this file with a text editor or inspect it with kmos' graphical editor by running ::

    kmos edit model_kmc.ini

If you are satisfied with the result you should go ahead and compile the model by issuing ::


    kmos export model_kmc.ini

If this takes a very long time you could try some other backend, like e.g.::

    kmos export -blat_int  model_kmc.ini


.. note::
    If you get an f2py error complaining that it cannot compile a certain submodule, it
    might help to explicit set an environment variable for the Fortran compiler you are using, like
    e.g ::
    
        exporrt F2PY_FCOMPILER=intelem

    Though the actual value may depend on your environment.

For more details please consult the `kmos documentation <http://kmos.rtfd.org>`_ . 

Afterwards you should see a new directory named `model_kmc_local_smart` or
something along the line. The kmc model (settings and source code) are all contained
within that folder. 

You should now copy the original `model.mkm` into the new folder `model_kmc_local_smart`
as well as other files required by the `*.mkm` such as `energies.txt` ::

    cp model.mkm energies.txt model_kmc_local_smart

Evaluating the kMC model
-------------------------------------------

From now one there are two routes forward: (1) you can either evaluate the model
for comparison with the mean-field model. To this end it is a good idea
to copy the files resulting from the successful CatMAP run into the kmc folder::

    cp model.mkm model.log model.pkl model_kmc_local_smart

For a simple straightforward comparison you should `cd` into the directory
and run ::

    catmap run_kmc


which should generate more files and finally plots. Or you could directly run ::


    kmos run


The default model runner can evaluate every descriptor tuple in a different process. A file
name `model.lock` keeps track of which data points are currently evaluated. In a cluster environemnt
you can therefore start many instances of the same evaluation process. Each of them will
keep running until all descriptor points have been evaluated. A simple submission script could
look as follows ::


    #!/bin/bash -eux
    #SBATCH -p slac
    #exclusive means the nodes shouldn't be shared with
    #other users
    #don't specify when running jobs that don't use all
    #cores on the nodes
    #SBATCH --exclusive
    #################
    #set a job name
    #SBATCH --job-name=myjob
    #################
    #a file for job output, you can check job progress
    #SBATCH --output=myjob.out
    #################
    # a file for errors from the job
    #SBATCH --error=myjob.err
    #################
    #time you think you need; default is one hour
    #in minutes in this case
    #SBATCH --time=08:00:00
    #################
    #number of nodes you are requesting
    #SBATCH --nodes=1
    #################
    #SBATCH --mem-per-cpu=4000
    #################
    #get emailed about job BEGIN, END, and FAIL
    #SBATCH --mail-type=ALL
    #################
    #who to send email to; please change to your email
    #SBATCH  --mail-user=maxjh@stanford.edu
    #################
    #task to run per node; each node has 16 cores
    #SBATCH --ntasks-per-node=16
    #################

    for i in {1..16}
    do
        catmap run_kmc -E 100000000 -S 100000000 &
        sleep 2
    done

Please change the email address and make the number of equilibration steps and sampling steps are sufficient.
