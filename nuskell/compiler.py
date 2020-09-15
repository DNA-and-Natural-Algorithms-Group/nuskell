#!/usr/bin/env python
#
# Copyright (c) 2010-2020 Caltech. All rights reserved.
# Written by Seung Woo Shin (seungwoo.theory@gmail.com).
#            Stefan Badelt (stefan.badelt@gmail.com)
#
#  nuskell/compiler.py
#  NuskellCompilerProject
#
import logging

import os
import sys
import argparse
import pkg_resources
from dsdobjects.utils import natural_sort

from nuskell import __version__

from nuskell.dsdcompiler import parse_ts_file, interpret
from nuskell.objects import TestTube, TestTubeIO
from nuskell.crnutils import (Reaction, 
                              parse_crn_string, 
                              assign_crn_species,
                              split_reversible_reactions, 
                              removeSpecies,
                              genCRN, genCON)
from nuskell.verifier import verify, modular_bisimulation

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

def translate(input_crn, ts_file, modular = False):
    """CRN-to-DSD translation wrapper function.

    A formal chemical reaction network (CRN) is translated into a domain-level
    strand displacement (DSD) system. The translation-scheme and the CRN are
    parsed into low-level instructions using the **nuskell.parser** module,
    passed on to the **nuskell.interpreter** and returned in form of a
    **nuskell.objects.TestTube()** object.

    Args:
      input_crn (str): An input string representation of the formal CRN.
      ts_file (str): The input file name of a translation scheme.
      modular (bool, optional): Split CRN into modules.

    Returns:
      [:obj:`TestTube()`,...]: A list of TestTube objects.
      The first object contains signal and fuel species of the full DSD
      system, followed by the modular system specifications.
    """
    if not os.path.isfile(ts_file):
        builtin = 'schemes/' + ts_file
        try:
            ts_file = pkg_resources.resource_filename('nuskell', builtin)
        except KeyError:
            schemedir = pkg_resources.resource_filename('nuskell', 'schemes')
            raise InvalidSchemeError(ts_file, schemedir)

    ts = parse_ts_file(ts_file)
    crn, fs = parse_crn_string(input_crn)

    solution, modules = interpret(ts, crn, fs, modular = modular)
    return solution, modules

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

    default.add_argument("-u", "--concentration-units", default='M', action='store',
            choices=('M', 'mM', 'uM', 'nM', 'pM'),
            help="""Specify default concentration units when writing results to 
            ouptut files (reaction rates, initial concentrations). """)

    default.add_argument("--enumerate", action='store_true',
            help="""Enumerate the DSD system. This is turned on automatically
            when using the argument --verify in combination with --ts or --readpil.""")

    # Choose a verification method.
    verify.add_argument("--verify", nargs = '+', default = [], action = 'store',
            choices=('bisimulation', 'pathway', 'integrated', 'modular-bisimulation',
                'bisim-loop-search', 'bisim-depth-first', 'bisim-whole-graph',
                'modular-bisim-loop-search', 'modular-bisim-depth-first',
                'modular-bisim-whole-graph'), metavar = '<str>', 
            help="""Specify verification methods. Choose one or more from:
            bisimulation, pathway, integrated, modular-bisimulation,
            bisim-loop-search, bisim-depth-first, bisim-whole-graph,
            modular-bisim-loop-search, modular-bisim-depth-first,
            modular-bisim-whole-graph.""")

    verify.add_argument("--modular", action = 'store_true',
            #help="""After enumeration of the full system, enumerate individual
            #CRN modules separately, to identify crosstalk between reactions.
            #This is turned on automatically when using modular-bisimulation
            #verification.""")
            help = argparse.SUPPRESS)

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
    fcrn, fs = parse_crn_string(input_crn)

    # ~~~~~~~~~~~~~~~~~~~~~~~
    # Initialize the TestTube
    # ~~~~~~~~~~~~~~~~~~~~~~~
    if args.ts:  # Translate CRN using a translation scheme
        logger.info(header("Translating"))
        solution, modules = translate(input_crn, args.ts,
                                      modular = args.modular)
    elif args.readpil:  # Parse information from a PIL file (cannot give modules!)
        logger.info(header("Parsing file {}".format(args.readpil)))
        solution = TestTube()
        TestTubeIO(solution).load_pil(args.readpil, is_file = True)
        # Assign fuel complexes and signal complexes.

        # signals = species that correspond to formal species in the formal CRN
        # fuels   = all non-signal species that have concentration > 0 nM
        for cplx in solution.complexes:
            if cplx.name in fs:
                solution.nodes[cplx]['ctype'] = 'signal'
            elif cplx.concentration is not None:
                solution.nodes[cplx]['ctype'] = 'fuel'
    else:
        # At some point Nuskell should choose translation schemes automatically,
        # but we are not there yet ... use nuskellCMP for such things.
        logger.error("Please specify a translation scheme, see option --ts.")
        schemedir = pkg_resources.resource_filename('nuskell', 'schemes')
        logger.error("For exaple the schemes in: {}".format(schemedir))
        for s in sorted(os.listdir(schemedir)):
            print("   --ts {}".format(s))
        raise SystemExit

    fuels = solution.fuel_complexes
    signals = solution.signal_complexes

    if args.ts and solution.unspecified_complexes != []:
        raise SystemExit('EXIT: solution contains intermediate species.')

    if signals == []:
        raise SystemExit('EXIT: solution does not contain signals.')

    logger.info(f"Formal species: {', '.join(natural_sort(fs))}")
    logger.info("Signal Complexes:")
    for cplx in natural_sort(signals):
        logger.info('   {} = {}'.format(cplx.name, cplx.kernel_string))
    logger.info("Fuel Complexes:")
    for cplx in natural_sort(fuels):
        logger.info('   {} = {}'.format(cplx.name, cplx.kernel_string))

    if dnafile:
        with open(dnafile, 'w') as dna:
            TestTubeIO(solution).write_dnafile(dna, signals = signals,
                                               crn = fcrn, ts = args.ts if args.ts else None)
        print("Wrote file: {}".format(dnafile))

    if args.pilfile:
        with open(comppil, 'w') as pil:
            TestTubeIO(solution).write_pil(fh = pil,
                                           unit = args.concentration_units, 
                                           crn = fcrn, fs = fs,
                                           ts = args.ts if args.ts else None)
        print("Compilation successfull. Wrote file: {}".format(comppil))
    else:
        print("Compilation successfull. Use --pilfile to inspect or modify.")
    print(" - signal species: {}\n - fuel species: {}\n".format(
            ' '.join(natural_sort(map(str, signals))), 
            ' '.join(natural_sort(map(str, fuels)))))

    if args.verify or args.enumerate:
        logger.info(header("Enumerating reaction pathways"))
        solution.enumerate_reactions(args, prefix = 'i', condensed = not args.enum_detailed)

        if not solution.reactions:
            raise SystemExit('No DSD reactions have been enumerated.')
        logger.info("After enumeration: {:3d} species, {:3d} reactions".format(
            len(solution.complexes), len(solution.reactions)))

        # Only for debugging.
        #pilstring = TestTubeIO(solution).write_pil(fh = None,
        #                                           unit = args.concentration_units, 
        #                                           crn = fcrn, fs = fs,
        #                                           ts = args.ts if args.ts else None)

        logger.info("Replacing regular-expression complexes ...")
        # History domains within constant and intermediate species will not get replaced.
        interpret = solution.interpret_species(fs.keys(), prune = True)
        signals = solution.signal_complexes

        # Return the implementation CRN. (= enumerated CRN)
        icrn = []
        for r in natural_sort(solution.reactions):
            logger.info("reaction {}".format(r.full_string()))
            rxn = Reaction(map(str, r.reactants), map(str, r.products), r.rate.constant, 0)
            icrn.append(rxn)

        if args.pilfile:
            with open(enumpil, 'w') as pil:
                TestTubeIO(solution).write_pil(fh = pil,
                                               unit = args.concentration_units, 
                                               crn = fcrn, fs = fs,
                                               ts = args.ts if args.ts else None)
            print("Enumeration successfull. Wrote file: {}".format(enumpil))
        else:
            print("Enumeration successfull. Use --pilfile to inspect or modify.")
        print(" - {:3d} species\n - {:3d} reactions".format(
            len(solution.complexes), len(solution.reactions)))
        print(" - {:3d} signal species".format(
            len([c for c in solution.complexes if c.name in interpret.keys()])))
        print(" - {:3d} fuel species".format(
            len([c for c in solution.complexes if c.name in map(str, fuels)])))
        print(" - {:3d} intermediate species".format(
            len([c for c in solution.complexes if c.name not in (
                list(interpret.keys()) + list(map(str, fuels)))])))
        print(' - {:3d} distinct strands in the system'.format(len(solution.strands)))
        print(' - {:3d} nucleotides to be designed\n'.format(sum(sum(
            map(lambda d: d.length, s)) for s in solution.strands.values())))

        if args.modular:
            # TODO(SB): cross-check verification:
            # - enumeration detects crosstalk between modules and append it as
            #   the last iCRN which has no corresponding fCRN.
            logger.info("")
            logger.info("Modular network enumeration ...")

            seen_reactions = set()
            all_reactions = set(solution.reactions)

            module_crns = []
            for e, module in enumerate(modules):
                # first, replace history complexes with their interpretation!
                for cplx in module.complexes:
                    # TODO quite inefficient loops
                    for k, v in interpret.items():
                        if (cplx.name in v) and k != cplx.name:
                            [newc] = solution.selected_complexes([k])
                            module.add_complex(newc, solution.get_complex_concentration(newc))
                            if module.has_complex(cplx):
                                module.rm_complex(cplx)

                module.enumerate_reactions(args, prefix = 'tmp', 
                        condensed = not args.enum_detailed)

                # after enumeration, make sure there are no new 'tmp' species present.
                for cplx in module.complexes:
                    assert cplx.name[:3] != 'tmp'

                # append the CRN
                logger.info("Module {}:".format(e))
                mcrn = []
                for rxn in module.reactions:
                    assert rxn in all_reactions
                    seen_reactions.add(rxn)
                    logger.info("reaction {}".format(rxn.full_string()))
                    rxn = Reaction([str(r) for r in rxn.reactants],
                                   [str(r) for r in rxn.products],
                                   rxn.rate.constant, 0)
                    mcrn.append(rxn)
                module_crns.append(mcrn)
                logger.info("")

            # last, identify crosstalk as the last implementation CRN 
            # (without formal correspondence)
            logger.info("Crosstalk:")
            mcrn = [] # crosstalk
            for r in all_reactions :
                if r not in seen_reactions :
                    logger.info("reaction {}".format(r.full_string()))
                    rxn = Reaction(list(map(str, r.reactants)), list(map(str, r.products)), r.rate.constant, 0)
                    mcrn.append(rxn)

            # append the CRN
            if mcrn:
                module_crns.append(mcrn)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Verify equivalence of CRNs
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    fcrm = list(map(lambda x: split_reversible_reactions([x]), fcrn))  # formal chemical reaction module
    fcrn = split_reversible_reactions(fcrn)

    if args.verify:
        logger.info(header("Verification using: {}".format(args.verify)))
        formals = set(fs.keys())

        logger.info(f"Formal CRN with {len(formals)} species:")
        [logger.info('    {}'.format(r)) for r in genCRN(fcrn, reversible = True)]
        logger.info("")

        vcrn = removeSpecies(icrn, list(map(str, fuels)))
        intermediates, wastes, _ = assign_crn_species(vcrn, set(map(str, signals)))
        vcrn = removeSpecies(vcrn, wastes)

        logger.info("Verification CRN (no fuels, no wastes):")
        # TODO: remove trival reactions and rates! Use clean_crn?
        [logger.info('    {}'.format(r)) for r in genCRN(vcrn, reversible = True)]
        logger.info("")

        logger.info("Partial interpretation:")
        for impl, formal in sorted(interpret.items()):
            logger.info("    {} => {}".format(impl, ', '.join(formal.elements())))
        logger.info("")

        if args.modular:
            icrns = list(map(lambda x: removeSpecies(x, set(map(str, fuels)) | wastes), module_crns))
            for e, m in enumerate(icrns):
                if e < len(fcrm):
                    logger.info("CRN Module {}:".format(e + 1))
                    logger.info("--")
                    [logger.info('    {}'.format(r)) for r in genCRN(fcrm[e], reversible = True)]
                else:
                    logger.info("CROSSTALK:")
                logger.info("--")
                [logger.info('    {}'.format(r)) for r in genCRN(m, reversible = True)]
                logger.info("")

        for meth in args.verify:
            logger.info(header("Verification method: {}".format(meth)))
            if 'modular-' in meth:
                if len(fcrm) > 1:
                    import copy # Temporary to fix a bug in testModules
                    backup = copy.deepcopy(interpret) if interpret else None
                    v, i = modular_bisimulation(fcrm, icrns, formals,
                                                interpret = backup, 
                                                method = meth[8:], 
                                                timeout = args.verify_timeout)
                else :
                    assert fcrn == fcrm[0]
                    v, i = verify(fcrn, vcrn, formals,
                                  interpret = interpret, 
                                  method = meth[8:], 
                                  timeout = args.verify_timeout)
            else:
                v, i = verify(fcrn, vcrn, formals,
                              interpret = interpret, 
                              method = meth, 
                              timeout = args.verify_timeout)

            if i and args.verbose:
                if not v:
                    i = i[2 if 'modular-' in meth and len(fcrm) > 1 else 0]
                logger.info("Returned CRN:")
                [logger.info('    {}'.format(r)) for r in genCRN(vcrn, 
                    reversible = True, rates = False, interpretation = i)]
                logger.info("Returned CRN {}:".format(meth))
                for impl, formal in sorted(i.items()):
                    logger.info("    {} => {}".format(impl, ', '.join(formal.elements())))

            logger.info("")
            if v is True:
                print("Verification result: {} - CRNs are {} equivalent.".format(v, meth))
            elif v is False:
                print("Verification result: {} - CRNs are not {} equivalent.".format(v, meth))
            elif v is None:
                print("No verification result: {} verification did not terminate within {} seconds.".format(meth, args.verify_timeout))

if __name__ == '__main__':
   main()

