#!/usr/bin/env python
#
# Copyright (c) 2010-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
#  nuskell/framework.py
#  NuskellCompilerProject
#
import logging
log = logging.getLogger(__name__)

import os
import sys
import argparse
import pkg_resources

from . import __version__
from .dsdcompiler import translate
from .dsdenumerator import enumerate_solution, enumerate_modules, interpret_species
from .crnverifier import verify, verify_modules
from .ioutils import (natural_sort, 
                      write_pil,
                      load_pil,
                      get_strands,
                      write_vdsd)
from .crnutils import (parse_crn_string, Reaction, 
                       split_reversible_reactions, 
                       assign_crn_species,
                       removeSpecies,
                       removeTrivial,
                       genCRN)


class InvalidSchemeError(Exception):
    """Raise Error: Cannot find translation scheme."""
    def __init__(self, ts_file, builtin=None):
        self.message = "Cannot find translation scheme: {}\n".format(ts_file)

        if builtin:
            self.message += "You may want to use one of the built-in schemes instead:\n"
            self.message += "Schemes in {}:\n".format(builtin)
            for s in sorted(os.listdir(builtin)):
                self.message += " * {}\n".format(s)

        super(InvalidSchemeError, self).__init__(self.message)

class colors:
    RED = '\033[91m'
    YELLOW = '\033[93m'
    GREEN = '\033[92m'
    BLUE = '\033[94m'
    PINK = '\033[95m'
    CYAN = '\033[96m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    colors = [RED, YELLOW, GREEN, CYAN, BLUE, PINK]

    @staticmethod
    def color(string):
        pass

    @staticmethod
    def legend(keys=None):
        if keys is None:
            l = enumerate(colors.colors)
        else:
            l = zip(keys, colors.colors)
        return "\n".join([(c + str(i) + colors.ENDC) for i, c in l])

class ColorFormatter(logging.Formatter):
    def __init__(self, msg, use_color=True):
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color
        self.COLORS = {
            'DEBUG': colors.CYAN,
            'INFO': colors.BLUE,
            'WARNING': colors.YELLOW,
            'ERROR': colors.RED,
            'Exception': colors.PINK,
        }
        self.RESET = colors.ENDC

    def format(self, record):
        levelname = record.levelname
        if self.use_color:
            record.levelname = self.COLORS[levelname] + \
                levelname + self.RESET

        return super(ColorFormatter, self).format(record)

def get_peppercorn_args(parser):
    """ Selected arguments for the peppercorn interface. """
    peppercorn = parser.add_argument_group('Peppercorn Reaction Enumerator Arguments')
    peppercorn.add_argument('--max-complex-size', default=50, type=int, metavar='<int>',
            help="""Maximum number of strands allowed in a complex (used to prevent polymerization)""")
    peppercorn.add_argument('--max-complex-count', default=1000, type=int, metavar='<int>',
            help="""Maximum number of complexes that may be enumerated before the enumerator halts.""")
    peppercorn.add_argument('--max-reaction-count', default=10000, type=int, metavar='<int>',
            help="""Maximum number of reactions that may be enumerated before the enumerator halts.""")

    peppercorn.add_argument('--reject-remote', action='store_true',
            help="Discard remote toehold mediated 3-way and 4-way branch migration reactions.")
    peppercorn.add_argument('--ignore-branch-3way', action='store_true',
            help="Ignore 3-way branch migration events during enumeration.")
    peppercorn.add_argument('--ignore-branch-4way', action='store_true',
            help="Ignore 4-way branch migration events during enumeration.")

    # TODO: explain these options in more detail!
    peppercorn.add_argument('--release-cutoff-1-1', type=int, default=6, metavar='<int>',
            help="""Maximum number of bases that will be released spontaneously in a 1-1 `open` reaction""")
    peppercorn.add_argument('--release-cutoff-1-n', type=int, default=6, metavar='<int>',
            help="""Maximum number of bases that will be released spontaneously in a 1-n `open` reaction.""")
    peppercorn.add_argument('--release-cutoff', type=int, default=None, metavar='<int>',
            help="""Maximum number of bases that will be released spontaneously
            in an `open` reaction, for either 1-1 or 1-n reactions (equivalent
            to setting --release-cutoff-1-1 and --release-cutoff-1-n to the
            same value)""")

    peppercorn.add_argument('--no-max-helix', action='store_true',
            help="""Do not apply 'max helix at a time' semantics to 3-way branch migration reactions.""")

    # NOTE: Output formatting: this option is not directly passed on to peppercorn
    peppercorn.add_argument('--enum-detailed', action='store_true',
                            help="Do not condense reactions into only resting complexes")

    # NOTE: The option --no-rates was removed, because peppercorn always computes
    # rates, but you may choose to ignore them for condensed reaction graphs.
    # k-fast enables to prune the condensed network, leaving the default of 0 M/s, has
    # the same effect.

    peppercorn.add_argument('--k-slow', default=0.0, type=float, metavar='<flt>',
                            help="Unimolecular reactions slower than this rate will be discarded")

    peppercorn.add_argument('--k-fast', default=0.0, type=float, metavar='<flt>',
                            help="Unimolecular reactions slower than this rate will be marked as slow")

    return parser

def get_nuskell_args(parser):
    """ A collection of arguments for Nuskell """
    #input = parser.add_mutually_exclusive_group(required=True)
    input = parser.add_argument_group('Nuskell Input Arguments (required - choose one)')
    default = parser.add_argument_group('Nuskell Output Arguments')
    verify = parser.add_argument_group('Nuskell Verification Arguments')
    simulate = parser.add_argument_group('Nuskell Simulation Arguments')

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument("--schemes", action='store_true', 
            help="Print a list of available schemes and exit.")
    parser.add_argument("-v", "--verbose", action='count', default=0,
            help="Print logging output. (-vv increases verbosity.)")
    parser.add_argument('--logfile', default='', action='store', metavar='<str>',
        help="""Redirect logging information to a file.""")

    # Options for the translation mode of Nuskell
    input.add_argument("--ts", action='store', metavar='</path/to/file>',
            help="""Specify a translation scheme. Either choose from the
            default schemes which are installed with Nuskell, or specify the
            full path to a not-installed scheme. Omitting both --ts and
            --readpil will print a list of available schemes and exit.""")

    input.add_argument("--readpil", action='store', metavar='</path/to/file>',
            help="""Read a domain-level strand displacement system from a file
            in PIL format. (Modular bisimulation verification is not available
            through this mode.)""")

    input.add_argument("--enumerated-pil", action='store', metavar='</path/to/file>',
            help="""Read a domain-level strand displacement system from a file
            in PIL format. Use this option to provide a custom set of reactions
            and supress built-in reaction enumeration. (Modular bisimulation
            verification is not available through this mode.)""")

    default.add_argument("-o", "--output", default='domainlevel', action='store', metavar='<str>', 
            help="Specify basename of output files.")

    default.add_argument("--pilfile", action='store_true',
            help="""Print (intermediate) results in the Pepper Internal
            Language (*.pil) file format.""")

    default.add_argument("--dnafile", action='store_true',
            help="Print (intermediate) results in the VisualDSD (*.dna) file format.")

    default.add_argument("-u", "--concentration-units", default='nM', action='store',
            choices=('M', 'mM', 'uM', 'nM', 'pM'),
            help="""Specify default concentration units when writing results to 
            ouptut files (reaction rates, initial concentrations). """)

    default.add_argument("--enumerate", action='store_true',
            help="""Enumerate the DSD system. This is turned on automatically
            when using the argument --verify in combination with --ts or --readpil.""")

    # Choose a verification method.
    verify.add_argument("--verify", nargs = '+', default = [], action = 'store',
            choices = ('crn-bisimulation', 
                       'crn-bisimulation-ls', 
                       'crn-bisimulation-rs', 
                       'modular-crn-bisimulation', 
                       'modular-crn-bisimulation-ls', 
                       'modular-crn-bisimulation-rs', 
                       'pathway-decomposition', 
                       'compositional-hybrid',
                       'integrated-hybrid'), metavar = '<str>', 
            help="""Specify verification methods. Choose one or more from:
            crn-bisimulation, crn-bisimulation-ls, crn-bisimulation-gs, 
            modular-crn-bisimulation, 
            modular-crn-bisimulation-ls, modular-crn-bisimulation-gs, 
            pathway-decomposition, integrated-hybrid, compositional-hybrid.""")

    verify.add_argument("--modular", action = 'store_true',
            #help="""After enumeration of the full system, enumerate individual
            #CRN modules separately, to identify crosstalk between reactions.
            #This is turned on automatically when using bisimulation
            #verification.""")
            help=argparse.SUPPRESS)

    verify.add_argument("--verify-timeout", type = int, default = 30, metavar = '<int>',
            help="Specify time in seconds to wait for verification to complete.")

    return parser

def header(msg):
    return ' {} '.format(msg).center(80, '*')

def set_handle_verbosity(h, v):
    if v == 0:
        h.setLevel(logging.WARNING)
    elif v == 1:
        h.setLevel(logging.INFO)
    elif v == 2:
        h.setLevel(logging.DEBUG)
    elif v >= 3:
        h.setLevel(logging.NOTSET)

def simulate_me(icrn):
    # A wrapper for enumeration? -> dsdenumerate.py
    raise NotImplementedError

def main():
    """ The Nuskell compiler.

      - translate formal CRNs into domain-level strand displacement systems.
      - verify the equivalence between formal CRN and implementation CRN.

    Output:
      - Domain-level DSD circuits printed into .pil and/or .dna files
      - Verbose information to STDOUT
      - verification results
      - simulatior scripts
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = """Nuskell: Compile a formal CRN to a DSD system specification.""")
    parser = get_nuskell_args(parser)
    parser = get_peppercorn_args(parser)
    args = parser.parse_args()

    # ~~~~~~~~~~~~~
    # Logging Setup 
    # ~~~~~~~~~~~~~
    title = "Nuskell Domain-level System Compiler {}".format(__version__)
    if args.verbose >= 3: # Get root logger.
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        lformat = '[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s'
    else:
        logger = logging.getLogger('nuskell')
        logger.setLevel(logging.DEBUG)
        lformat = '%(levelname)s %(message)s'

    if args.logfile:
        fh = logging.FileHandler(args.logfile)
        formatter = logging.Formatter(lformat)
        set_handle_verbosity(fh, args.verbose)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    else:
        ch = logging.StreamHandler()
        formatter = ColorFormatter(lformat, use_color = True)
        set_handle_verbosity(ch, args.verbose)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
    logger.info(header(title))
    logger.info("")

    if args.schemes:
        schemedir = pkg_resources.resource_filename('nuskell', 'schemes')
        print("Listing schemes in: {}".format(schemedir))
        for s in sorted(os.listdir(schemedir)):
            print("   --ts {}".format(s))
        raise SystemExit

    # ~~~~~~~~~~~~~~~~~~~
    # Argument processing
    # ~~~~~~~~~~~~~~~~~~~
    comppil = args.output + '_sys.pil' if args.pilfile else None
    enumpil = args.output + '_enum.pil' if args.pilfile else None
    dnafile = args.output + '.dna' if args.dnafile else None

    if not args.modular:
        args.modular = any(map(lambda x: 'modular' in x, args.verify))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Parse and process input CRN
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    input_crn = sys.stdin.readlines()
    input_crn = "".join(input_crn)
    fcrn, fsc = parse_crn_string(input_crn)

    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Initialize the TestTube
    # ~~~~~~~~~~~~~~~~~~~~~~~
    if args.ts:  # Translate CRN using a translation scheme
        logger.info(header("Translating"))

        if not os.path.isfile(args.ts):
            builtin = 'schemes/' + args.ts
            try:
                args.ts = pkg_resources.resource_filename('nuskell', builtin)
            except KeyError:
                schemedir = pkg_resources.resource_filename('nuskell', 'schemes')
                raise InvalidSchemeError(args.ts, schemedir)

        solution, modules = translate(input_crn, args.ts,
                                      modular = args.modular)
    #elif args.readpil:  # Parse information from a PIL file (cannot give modules!)
    #    logger.info(header("Parsing file {}".format(args.readpil)))
    #    solution = TestTube()
    #    TestTubeIO(solution).load_pil(args.readpil, is_file = True)
    #    # Assign fuel complexes and signal complexes.

    #    # signals = species that correspond to formal species in the formal CRN
    #    # fuels   = all non-signal species that have concentration > 0 nM
    #    for cplx in solution.complexes:
    #        if cplx.name in fsc:
    #            solution.nodes[cplx]['ctype'] = 'signal'
    #        elif cplx.concentration is not None:
    #            solution.nodes[cplx]['ctype'] = 'fuel'
    else:
        # At some point Nuskell should choose translation schemes automatically,
        # but we are not there yet ... use nuskellCMP for such things.
        logger.error("Please specify a translation scheme, see option --ts.")
        schemedir = pkg_resources.resource_filename('nuskell', 'schemes')
        logger.error("For exaple the schemes in: {}".format(schemedir))
        for s in sorted(os.listdir(schemedir)):
            print("   --ts {}".format(s))
        raise SystemExit

    fuels = [x for x in solution.values() if x.name[0] == 'f']
    wastes = [x for x in solution.values() if x.name[0] == 'w']
    intermediates = [x for x in solution.values() if x.name[0] == 'i']
    signals = [x for x in solution.values() if x.name[0] not in ('f', 'i', 'w')]

    if args.ts and intermediates:
        raise SystemExit('EXIT: solution contains intermediate species.')

    if signals == []:
        raise SystemExit('EXIT: solution does not contain signals.')

    logger.info(f"Formal species: {', '.join(natural_sort(fsc))}")
    logger.info("Signal Complexes:")
    for cplx in natural_sort(signals):
        logger.info('   {} = {}'.format(cplx.name, cplx.kernel_string))
    logger.info("Fuel Complexes:")
    for cplx in natural_sort(fuels):
        logger.info('   {} = {}'.format(cplx.name, cplx.kernel_string))

    if dnafile:
        with open(dnafile, 'w') as dna:
            write_vdsd(solution, 
                       fh = dna, 
                       molarity = args.concentration_units, 
                       crn = fcrn,
                       fsc = fsc,
                       ts = args.ts if args.ts else None)
        print("Wrote file: {}".format(dnafile))

    if args.pilfile:
        with open(comppil, 'w') as pil:
            write_pil(solution, None,
                      fh = pil,
                      molarity = args.concentration_units, 
                      crn = fcrn, 
                      fsc = fsc,
                      ts = args.ts if args.ts else None)
        print("Compilation successfull. Wrote file: {}".format(comppil))
    else:
        print("Compilation successfull. Use --pilfile to inspect or modify.")

    print(" - signal species: {}\n - fuel species: {}\n".format(
            ' '.join(natural_sort(map(str, signals))), 
            ' '.join(natural_sort(map(str, fuels)))))

    if args.verify or args.enumerate:
        logger.info(header("Enumerating reaction pathways"))
        complexes, reactions = enumerate_solution(solution, args)

        if not reactions:
            raise SystemExit('No DSD reactions have been enumerated.')

        logger.info(f"After enumeration: {len(complexes)} species, {len(reactions)} reactions")

        # Only for debugging.
        logger.debug('\n' + write_pil(complexes, reactions, fh = None,
                               molarity = args.concentration_units))

        logger.info("Removing unnecessary complexes with history domains.")
        # History domains within constant and intermediate species will not get replaced.
        interpretation, complexes, reactions = interpret_species(complexes, 
                                                                reactions,
                                                                fsc.keys(),
                                                                prune = True)
        
        # Update species assignments
        fuels = [x for x in complexes.values() if x.name[0] == 'f']
        wastes = [x for x in complexes.values() if x.name[0] == 'w']
        intermediates = [x for x in complexes.values() if x.name[0] == 'i']
        signals = [x for x in complexes.values() if x.name[0] not in ('f', 'i', 'w')]

        logger.info(f"Enumerated CRN: \n  " + \
                    '\n  '.join([rxn.full_string() for rxn in reactions]))

        if args.pilfile:
            with open(enumpil, 'w') as pil:
                write_pil(complexes, reactions,
                          fh = pil,
                          molarity = args.concentration_units, 
                          crn = fcrn, 
                          fsc = fsc,
                          ts = args.ts if args.ts else None)
            print("Enumeration successfull. Wrote file: {}".format(enumpil))
        else:
            print("Enumeration successfull. Use --pilfile to inspect or modify.")

        print(" - {:3d} species\n - {:3d} reactions".format(len(complexes), len(reactions)))
        print(" - {:3d} signal species".format(len(signals)))
        print(" - {:3d} fuel species".format(len(fuels)))
        print(" - {:3d} intermediate species".format(len(intermediates + wastes)))
        print(' - {:3d} distinct strands in the system'.format(len(get_strands(complexes))))
        print(' - {:3d} nucleotides to be designed\n'.format(sum(
            sum(map(lambda d: d.length, s)) for s in get_strands(complexes))))

        if args.modular:
            logger.info("")
            logger.info("Modular network enumeration ...")
            mcomplexes, mreactions = enumerate_modules(modules, 
                                                       interpretation,
                                                       complexes,
                                                       reactions,
                                                       args)
            print(f"Split reaction enumeration into {len(mcomplexes)} modules:")
            for e in range(len(mcomplexes)):
                print(f' - module {e+1}: {len(mcomplexes[e])} complexes and {len(mreactions[e])} reactions.')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Verify equivalence of CRNs
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    if args.verify: # Needs genCRN, removeSpecies, assign_crn_species, split_rev_reactions, ...
        # TODO: clean_crn!
        logger.info(header("Verification using: {}".format(args.verify)))
        formals = set(fsc.keys())
        fuels = set([x.name for x in fuels])

        logger.info(f"Formal CRN with {len(formals)} species:\n  " + \
                    '\n  '.join(genCRN(fcrn, reversible = True)))

        icrn = [Reaction([str(x) for x in rxn.reactants], 
                         [str(x) for x in rxn.products], 
                         rxn.rate.constant, 0) for rxn in reactions]
        logger.debug(f"Implementation CRN:\n  " + \
                    '\n  '.join(natural_sort(genCRN(icrn, reversible = True))))

        icrn = removeSpecies(icrn, fuels)
        intermediates, wastes, _ = assign_crn_species(icrn, set(map(str, signals)))
        icrn = removeTrivial(removeSpecies(icrn, wastes))
        logger.info(f"Implementation CRN (no fuels, no wastes, no rates): \n  " + \
                    '\n  '.join(natural_sort(genCRN(icrn, reversible = True, rates = False))))
        logger.info("Partial interpretation:\n  " + \
                    '\n  '.join(f"{k} => {', '.join(v)}" \
                        for k, v in natural_sort(interpretation.items())))

        if args.modular:
            fcrns = [[m] for m in fcrn]
            icrns = []
            for e, module in enumerate(mreactions, 1):
                mcrn = [Reaction([str(x) for x in rxn.reactants], 
                                 [str(x) for x in rxn.products], 
                                 rxn.rate.constant, 0) for rxn in module]
                mcrn = removeTrivial(removeSpecies(mcrn, fuels | wastes))
                icrns.append(mcrn)
                logger.info(f"Implementation Module {e}:\n  " + \
                            '\n  '.join(natural_sort(
                                        genCRN(mcrn, reversible = True, rates = False))))

        for meth in args.verify:
            logger.info(header("Verification method: {}".format(meth)))
            if 'modular-' in meth and len(fcrns) > 1:
                v, i = verify_modules(fcrns, icrns, formals, meth[8:], 
                                      interpretation = interpretation, 
                                      timeout = args.verify_timeout)
            else:
                if 'modular' in meth: meth = meth[8:]
                v, i = verify(fcrn, icrn, formals, meth, 
                              interpretation = interpretation, 
                              timeout = args.verify_timeout)

            if i:
                logger.info(f"Returned interpretation for {meth}:\n  " + \
                            '\n  '.join(f"{k} => {', '.join(v)}" \
                                for k, v in natural_sort(i.items())))
                logger.info(f"Interpreted CRN: \n  " + \
                    '\n  '.join(natural_sort(genCRN(icrn, 
                                                    reversible = True, 
                                                    rates = False,
                                                    interpretation = i))))
            if v is True:
                print(f"Verification result: {v}.",
                      f"The implementation CRN is correct according to {meth}.")
            elif v is False:
                print(f"Verification result: {v}.",
                      f"The implementation CRN is *not* correct according to {meth}.")
            elif v is None:
                print(f"No verification result for {meth}.", 
                      f"Verification did not terminate within {args.verify_timeout} seconds.")

if __name__ == '__main__':
   main()

