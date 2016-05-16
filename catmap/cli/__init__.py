import os
import pprint
import shutil

usage = {}

SCRIPT = 'catmap'

usage['all'] = """{SCRIPT} help all
    Display documentation for all commands.
                """.format(**locals())

usage['help'] = """{SCRIPT} help <command>
    Print usage information for the given command.
                """.format(**locals())
usage['import'] = """{SCRIPT} import <mkm-file>
    Open a *.mkm project file and work with it interactively.
""".format(**locals())

usage['to_kmc'] = """{SCRIPT} to_kmc <mkm-file> -i <N>
    Translate a given *.mkm to a kmos kMC model.
    Options :
      -i, --interaction N (default: 0)
            How many shells of nearest neighbor interactions are considered.
            Default is 0, beware of computational effort of increasing
            interaction range. Usually 1 or 2 should be sufficient
     -l, --validate (default: True)
            Validate the kmos kMC model before writing it to INI
""".format(**locals())

usage['run_kmc'] = """{SCRIPT} run_kmc -E <N> -S <N>
    Options :
        -E, --equilibration-steps (default: 1e8)
            number of equilibration steps before sampling
            CHECK CAREFULLY IF SUFFICIENT

        -S, --sampling-steps (default 1e8)
           number of sampling steps to calculate rates
           and coverages.
           CHECK CAREFULLY IF SUFFICIENT

        -s, --single-point
          Force 'catmap run_kmc'  to only evalute one descriptor point
          and finish evaluation.
          Otherwise the model runner will keep evaluating descriptor
          points until they are all exhausted.

        -n, --dont-run
         If 'catmap run_kmc' should only plot results and
         not evaluate any descriptor points

        -p, --print
         If 'catmap run_kmc' should plot results

        -c, --initial-configuration
         Pick how the kmc model runner initializes the lattice.
         Currently implemented are:
         probabilistic [default], empty, species:<species_name>, and majority.


""".format(**locals())


def get_options(args=None, get_parser=False):
    import optparse
    import os
    from glob import glob
    import catmap

    parser = optparse.OptionParser(
        'Usage: %prog [help] ('
        + '|'.join(sorted(usage.keys()))
        + ') [options]',
        version=catmap.__version__)

    parser.add_option('-E', '--equilibration-steps', type='float', dest='equilibration_threshold', default=1e-2, help="The number of kmc steps before it starts sampling averages. If < 100, will be interpreted as 10**x.")
    parser.add_option('-S', '--sampling-min', type='int', dest='sampling_min', default=int(1e3), help="The number of kmc steps used to calculate averages. If < 100, will be interpreted as 10**x.")
    parser.add_option('-a', '--author-name', dest='author_name', default='CatMAP User', help="Specify your name for the translated kMC model (catmap to_kmc -a 'Joe Blow' ..")
    parser.add_option('-b', '--batch-size', dest='batch_size', default=int(1e6), type='float', help="The default number of steps per batch for determinining the steady-state and sampling averages of rates and coverages")
    parser.add_option('-C', '--coverage-samples', dest='coverage_samples', default=100, help="Set how many sample are taken of the coverage during the sampling run [100]", type='int')
    parser.add_option('-c', '--initial-configuration', dest='initial_configuration', default='probabilistic', help="Pick how the kmc model runner initializes the lattice. Currently implemented are: probabilistic [default], empty, species:<species_name>, and majority")
    parser.add_option('-d', '--diffusion-factor', dest='diffusion_factor', default=None, help="Parameter to tune the rate-constant of diffusion processes. The rationale is that for transition metal surfaces diffusion barriers are typically quite low, leading to fast surface diffusion processes. If no factor is specified the diffusion barrier is derived from whatever scaling relation is specified. If a diffusion factor X is specified all diffusion rate constants are set to k_max*X where k_max is the faster rate-constant among all non-diffusion elementary processes. If X < 10, the factor is interpreted as 10**X")
    parser.add_option('-e', '--author-email', dest='author_email', default='mkm-developers-request@lists.stanford.edu', help="Specify your email address for the kmc translated model (catmap to_kmc -e ...)")
    parser.add_option('-i', '--interaction', dest='interaction', type='int', default=0)
    parser.add_option('-l', '--validate', dest='validate', action='store_false', default=True, help="Validate the kmos kMC model before writing it to INI")
    parser.add_option('-m', '--simulation-size', dest='simulation_size', default=20, type='int', help="The lattice size of the kMC lattice")
    parser.add_option('-n', '--dont-run', dest='dontrun', action='store_true', default=False, help="If 'catmap run_kmc' should only plot results")
    parser.add_option('-p', '--plot', dest='plot', action='store_true', default=False, help="If 'catmap run_kmc' should plot results")
    parser.add_option('-s', '--single-point', dest='single_point', action='store_true', default=False, help="Force 'catmap run_kmc'  to only evalute one descriptor point and finish evaluation. Otherwise the model runner will keep evaluating descriptor points until they are all exhausted.")
    parser.add_option('-t', '--bias-threshold', dest='bias_threshold', default=0.1, type='float', help="Sample trajectory until it has been determined that this amount or less is in transient part")
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', default=False, help="Make catmap CLI verbose about what it is doing")
    parser.add_option('-w', '--ewma-plots', dest='make_plots', action='store_true', default=False, help="Have the steady-state runner create EWMA plots at every evaluation step.")


    if args is not None:
        options, args = parser.parse_args(args.split())
    else:
        options, args = parser.parse_args()

    # Assume that if we plot, we don't want to
    # run at the same time
    if options.plot:
        options.dontrun = True

    ## Interpret sampling steps and equilibration steps < 100 as 10^x

    if options.batch_size < 100 :
        options.batch_size = int(10.**options.batch_size)
    #else:
        #options.sampling_steps = int(options.sampling_steps)
    #if options.sampling_steps < 100 :
        #options.sampling_steps = int(10.**options.sampling_steps)
    #else:
        #options.sampling_steps = int(options.sampling_steps)

    #if options.equilibration_steps < 100 :
        #options.equilibration_steps = int(10.**options.equilibration_steps)
    #else:
        #options.equilibration_steps = int(options.equilibration_steps)

    if options.diffusion_factor is not None:
        options.diffusion_factor = float(options.diffusion_factor)
        # to make batch calculations simpler assume that any negative
        # diffusion factor means that the original diffusion rate
        # constant is not affected
        if options.diffusion_factor < -100:
            options.diffusion_factor = None
        elif abs(options.diffusion_factor) < 100:
            options.diffusion_factor = 10.**float(options.diffusion_factor)

    if len(args) < 1:
        parser.error('Command expected')
    if get_parser:
        return options, args, parser
    else:
        return options, args


def match_keys(arg, usage, parser):
    """Try to match part of a command against
       the set of commands from usage. Throws
       an error if not successful.

    """
    possible_args = [key for key in usage if key.startswith(arg)]
    if len(possible_args) == 0:
        parser.error('Command "%s" not understood.' % arg)
    elif len(possible_args) > 1:
        parser.error(('Command "%s" ambiguous.\n'
                      'Could be one of %s\n\n') % (arg, possible_args))
    else:
        return possible_args[0]


def main(args=None):
    """The CLI main entry point function.

    The optional argument args, can be used to
    directly supply command line argument like

    $ catmap <args>

    otherwise args will be taken from STDIN.

    """
    from glob import glob

    options, args, parser = get_options(args, get_parser=True)

    if options.verbose:
        print("CatMAP CLI options")
        pprint.pprint(options)

    if not args[0] in usage.keys():
        args[0] = match_keys(args[0], usage, parser)
        print(args)

    elif args[0] == 'help':
        if len(args) < 2:
            parser.error('Which help do you  want?')
        if args[1] == 'all':
            for command in sorted(usage):
                print(usage[command])
        elif args[1] in usage:
            print('Usage: %s\n' % usage[args[1]])
        else:
            arg = match_keys(args[1], usage, parser)
            print('Usage: %s\n' % usage[arg])

    elif args[0] == 'import':
        if len(args) < 2:
            parser.error('mkm filename expected.')

        from catmap import ReactionModel
        mkm_file = args[1]
        global model
        model = ReactionModel(setup_file=mkm_file)
        sh(banner='Note: model = catmap.ReactionModel(setup_file=\'%s\')\n# do model.run()\nfor a fully initialized model.' %
           args[1])

    elif args[0] == 'to_kmc':
        import catmap.cli.kmc_translation
        mkm_file = args[1]
        catmap.cli.kmc_translation.translate_model_file(mkm_file, options)

    elif args[0] == 'run_kmc':
        call_path = os.path.abspath(os.getcwd())
        import catmap.cli.kmc_runner
        catmap.cli.kmc_runner.main(options, call_path=call_path)

    elif args[0] == 'version':
        import catmap
        print(catmap.__version__)
    else:
        parser.error('Command "%s" not understood.' % args[0])


def sh(banner):
    """Wrapper around interactive ipython shell
    that factors out ipython version depencies.

    """

    from distutils.version import LooseVersion
    import IPython
    if hasattr(IPython, 'release'):
        try:
            from IPython.terminal.embed import InteractiveShellEmbed
            InteractiveShellEmbed(banner1=banner)()

        except ImportError:
            try:
                from IPython.frontend.terminal.embed \
                    import InteractiveShellEmbed
                InteractiveShellEmbed(banner1=banner)()

            except ImportError:
                from IPython.Shell import IPShellEmbed
                IPShellEmbed(banner=banner)()
    else:
        from IPython.Shell import IPShellEmbed
        IPShellEmbed(banner=banner)()
