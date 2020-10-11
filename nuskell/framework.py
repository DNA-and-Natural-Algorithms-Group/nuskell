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

from . import __version__
from .dsdcompiler import translate, get_builtin_schemes
from .dsdenumerator import enumerate_solution, enumerate_modules, interpret_species
from .crnverifier import verify, verify_modules
from .ioutils import (natural_sort, 
                      write_pil,
                      load_pil,
                      get_strands,
                      write_vdsd)
from .crnutils import (parse_crn_string, Reaction, 
                       assign_crn_species,
                       remove_species,
                       cleanup_rxns,
                       genCRN)

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

def header(msg):
    return ' {} '.format(msg).center(80, '*')

def get_peppercorn_args(parser):
    """ Selected arguments for the peppercorn interface. """
    peppercorn = parser.add_argument_group('Peppercorn Reaction Enumerator Arguments')
    peppercorn.add_argument('--max-complex-size', default=10, type=int, metavar='<int>',
            help="""Maximum number of strands allowed in a complex (to prevent polymerization).""")
    peppercorn.add_argument('--max-complex-count', default=5_000, type=int, metavar='<int>',
            help="""Maximum number of complexes that may be enumerated before the enumerator halts.""")
    peppercorn.add_argument('--max-reaction-count', default=10_000, type=int, metavar='<int>',
            help="""Maximum number of reactions that may be enumerated before the enumerator halts.""")

    peppercorn.add_argument('--reject-remote', action='store_true',
            help="Discard remote toehold mediated 3-way and 4-way branch migration reactions.")
    peppercorn.add_argument('--ignore-branch-3way', action='store_true',
            help="Ignore 3-way branch migration events during enumeration.")
    peppercorn.add_argument('--ignore-branch-4way', action='store_true',
            help="Ignore 4-way branch migration events during enumeration.")

    peppercorn.add_argument('--release-cutoff-1-1', type=int, default=7, metavar='<int>',
            help="""Maximum number of bases that will be released spontaneously in a 1-1 `open` reaction.""")
    peppercorn.add_argument('--release-cutoff-1-2', type=int, default=7, metavar='<int>',
            help="""Maximum number of bases that will be released spontaneously in a 1-2 `open` reaction.""")
    peppercorn.add_argument('--release-cutoff', type=int, default=None, metavar='<int>',
            help="""Maximum number of bases that will be released spontaneously
            in an `open` reaction, for either 1-1 or 1-2 reactions (equivalent
            to setting --release-cutoff-1-1 and --release-cutoff-1-2 to the
            same value).""")

    peppercorn.add_argument('--no-max-helix', action='store_true',
            help="""Do not apply 'max helix at a time' semantics to 3-way branch migration reactions.""")

    peppercorn.add_argument('--enum-detailed', action='store_true',
            help="Do not condense the reaction network.")

    peppercorn.add_argument('--k-slow', default=0.0, type=float, metavar='<flt>',
            help="Unimolecular reactions with a rate constant lower than k-slow will be discarded.")

    peppercorn.add_argument('--k-fast', default=0.0, type=float, metavar='<flt>',
            help="Unimolecular reactions with a rate constant lower than k-fast will be slow.")

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
                       'crn-bisimulation-bf', 
                       'modular-crn-bisimulation', 
                       'modular-crn-bisimulation-ls', 
                       'modular-crn-bisimulation-bf', 
                       'pathway-decomposition', 
                       'compositional-hybrid',
                       'integrated-hybrid'), metavar = '<str>', 
            help="""Specify verification methods. Choose one or more from:
            crn-bisimulation, crn-bisimulation-ls, crn-bisimulation-bf, 
            modular-crn-bisimulation, 
            modular-crn-bisimulation-ls, modular-crn-bisimulation-bf, 
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

def set_handle_verbosity(h, v):
    if v == 0:
        h.setLevel(logging.WARNING)
    elif v == 1:
        h.setLevel(logging.INFO)
    elif v == 2:
        h.setLevel(logging.DEBUG)
    elif v >= 3:
        h.setLevel(logging.NOTSET)

def assign_species(complexes):
    """ The nuskell naming standard for all types of species. """
    fuels = [x for x in complexes.values() if x.name[0] == 'f']
    wastes = [x for x in complexes.values() if x.name[0] == 'w']
    intermediates = [x for x in complexes.values() if x.name[0] == 'i']
    signals = [x for x in complexes.values() if x.name[0] not in ('f', 'i', 'w')]
    return fuels, wastes, intermediates, signals

def get_verification_crn(reactions, fuels, signals):
    # Prepare the verification CRN - Step 1: 
    fuels = set([x.name for x in fuels]) 
    signals = set([x.name for x in signals])
    icrn = [Reaction([str(x) for x in rxn.reactants], 
                     [str(x) for x in rxn.products], 
                     rxn.rate.constant, 0) for rxn in reactions]
    log.debug(f"Implementation CRN:\n  " + \
                '\n  '.join(natural_sort(genCRN(icrn, reversible = True))))

    # Prepare the verification CRN - Step 2: 
    icrn = remove_species(icrn, fuels)
    intermediates, wastes, _ = assign_crn_species(icrn, signals)
    icrn = cleanup_rxns(remove_species(icrn, wastes))
    return icrn, fuels, wastes

def get_verification_modules(fcrn, mreactions, fuels, wastes):
    fcrns = [[m] for m in fcrn]
    icrns = []
    for e, module in enumerate(mreactions, 1):
        mcrn = [Reaction([str(x) for x in rxn.reactants], 
                         [str(x) for x in rxn.products], 
                         rxn.rate.constant, 0) for rxn in module]
        mcrn = cleanup_rxns(remove_species(mcrn, fuels | wastes))
        icrns.append(mcrn)
    return fcrns, icrns

def main():
    """ The nuskell compiler framework commandline interface.

    Translates formal chemical reaction networks (CRNs) into domain-level
    strand displacement (DSD) systems. Enumerates DSD systems using the
    domain-level reaction enumerator "peppercorn". Verifies the correctness of
    the traslation using the notions of "pathway decomposition", "crn
    bisimulation", or hybrid notions thereof using the Python package
    "crnverifier". Use nuskell --help for options. 

    Nuskell can report every stage of the compilation by printing the
    corresponding *.pil file. After initial translation, this file only
    contains signal and fuel species. After enumeration it contains also
    the (potentially processed) enumerated reaction network.
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = """Nuskell: Compile a formal CRN to a DSD system specification.""")
    parser = get_nuskell_args(parser)
    parser = get_peppercorn_args(parser)
    args = parser.parse_args()

    # ~~~~~~~~~~~~~ #
    # Logging Setup #
    # ~~~~~~~~~~~~~ #
    title = f"Nuskell Domain-level System Compiler {__version__}"
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
        schemes = get_builtin_schemes()
        for k, v in schemes.items():
            print(f'Listing installed schemes: "{k}"')
            for s in v:
                print(f"   --ts {s}")
        raise SystemExit

    # ~~~~~~~~~~~~~~~~~~~ #
    # Argument processing #
    # ~~~~~~~~~~~~~~~~~~~ #
    comppil = args.output + '_sys.pil' if args.pilfile else None
    enumpil = args.output + '_enum.pil' if args.pilfile else None
    dnafile = args.output + '.dna' if args.dnafile else None

    if not args.modular:
        args.modular = any(map(lambda x: 'modular' in x, args.verify))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Parse and process input CRN #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    input_crn = "".join(sys.stdin.readlines())
    fcrn, fsc = parse_crn_string(input_crn)

    # ~~~~~~~~~~~~~~~~~~ #
    # Do the translation #
    # ~~~~~~~~~~~~~~~~~~ #
    if args.ts:  # Translate CRN using a translation scheme
        log.info(header(f"Translating using scheme {args.ts}"))
        solution, modules = translate(input_crn, args.ts, modular = args.modular)
    elif args.readpil:  # Parse information from a PIL file 
        if args.modular:
            raise NotImplementedError('Modular verification cannot be used in combiation with --readpil input.')
        log.info(header(f"Parsing file {args.readpil}"))
        dom, solution, rms, det, con = load_pil(args.readpil, is_file = True)
        if det:
            log.warning(f'Ignoring {len(det)} detailed reactions in {args.readpil}.')
        if con:
            log.warning(f'Ignoring {len(con)} condensed reactions in {args.readpil}.')
        if rms:
            log.warning(f'Ignoring {len(rms)} resting macrostates in {args.readpil}.')
    else:
        log.error("Please specify a translation scheme, see option --ts.")
        schemes = get_builtin_schemes()
        for k, v in schemes.items():
            log.error(f'Listing installed schemes: "{k}"')
            for s in v:
                log.error(f"   --ts {s}")
        raise SystemExit

    fuels, wastes, intermediates, signals = assign_species(solution)

    if args.ts and intermediates:
        raise SystemExit('EXIT: solution contains intermediate species.')
    if signals == []:
        raise SystemExit('EXIT: solution does not contain signals.')

    log.info(f"Formal species: {', '.join(natural_sort(fsc))}")
    log.info("Signal Complexes:\n" + '\n'.join(
        [f'   {cplx.name} = {cplx.kernel_string}' for cplx in natural_sort(signals)]))
    log.info("Fuel Complexes:\n" + '\n'.join(
        [f'   {cplx.name} = {cplx.kernel_string}' for cplx in natural_sort(fuels)]))

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
        log.info(header("Enumerating reaction pathways."))
        complexes, reactions = enumerate_solution(solution, args, molarity = args.concentration_units)

        if not reactions:
            raise SystemExit('No DSD reactions have been enumerated.')

        log.info("After enumeration: " + \
                f"{len(complexes)} species, {len(reactions)} reactions")

        # Only for debugging.
        log.debug('\n' + write_pil(complexes, reactions, fh = None,
                               molarity = args.concentration_units))

        log.info("Removing unnecessary complexes with history domains.")
        # History domains within constant and intermediate species will not get replaced.
        interpretation, complexes, reactions = interpret_species(complexes, 
                                                                reactions,
                                                                fsc.keys(),
                                                                prune = True)
        # Update species assignments
        fuels, wastes, intermediates, signals = assign_species(complexes)
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

        print(" - {:3d} species\n - {:3d} reactions".format(len(complexes), 
                                                            len(reactions)))
        print(" - {:3d} signal species".format(len(signals)))
        print(" - {:3d} fuel species".format(len(fuels)))
        print(" - {:3d} intermediate species".format(len(intermediates + wastes)))
        print(' - {:3d} distinct strands in the system'.format(
                                                        len(get_strands(complexes))))
        print(' - {:3d} nucleotides to be designed\n'.format(sum(
            sum(map(lambda d: d.length, s)) for s in get_strands(complexes))))

        if args.modular:
            log.info("")
            log.info("Modular network enumeration ...")
            mcomplexes, mreactions = enumerate_modules(modules, 
                                                       interpretation,
                                                       complexes,
                                                       reactions,
                                                       args)
            print(f"Split reaction enumeration into {len(mcomplexes)} modules:")
            for e in range(len(mcomplexes)):
                print(f' - module {e+1}: {len(mcomplexes[e])} complexes',
                                   f'and {len(mreactions[e])} reactions.')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Verify correctness of implementation CRN #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    if args.verify:
        logger.info(header("Verification using: {}".format(args.verify)))
        formals = set(fsc.keys())
        log.info(f"Formal CRN with {len(formals)} species:\n  " + \
                    '\n  '.join(genCRN(fcrn, reversible = True)))
        icrn, fuels, wastes = get_verification_crn(reactions, fuels, signals)
        log.info(f"Implementation CRN (no fuels, no wastes, no rates): \n  " + \
            '\n  '.join(natural_sort(genCRN(icrn, reversible = True, rates = False))))
        log.info("Partial interpretation:\n  " + \
                '\n  '.join(f"{k} => {', '.join(v)}" \
                    for k, v in natural_sort(interpretation.items())))

        if args.modular:
            fcrns, icrns = get_verification_modules(fcrn, mreactions, fuels, wastes)
            for e, mcrn in enumerate(icrns, 1):
                log.info(f"Implementation Module {e}:\n  " + \
                        '\n  '.join(natural_sort(genCRN(mcrn, 
                            reversible = True, rates = False))))
        for meth in args.verify:
            log.info(header("Verification method: {}".format(meth)))
            if 'modular-' in meth and len(fcrns) > 1:
                v, i = verify_modules(fcrns, icrns, formals, meth[8:], 
                                      interpretation = interpretation, 
                                      timeout = args.verify_timeout)
            else:
                if 'modular' in meth: meth = meth[8:]
                v, i = verify(fcrn, icrn, formals, meth, 
                              interpretation = interpretation, 
                              timeout = args.verify_timeout)

            if v:
                log.info(f"Returned interpretation for {meth}:\n  " + \
                            '\n  '.join(f"{k} => {', '.join(v)}" \
                                for k, v in natural_sort(i.items())))
                log.info(f"Interpreted CRN: \n  " + \
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

