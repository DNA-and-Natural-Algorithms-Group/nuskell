#!/usr/bin/env python
#
#  compare_schemes.py
#  EnumeratorProject
#
import logging
log = logging.getLogger(__name__)

import os
import sys
import argparse
import pkg_resources
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import nan

from dsdobjects import clear_memory, DSDDuplicationError
from peppercornenumerator import (CondensationError, 
                                  PolymerizationError)

from . import __version__
from .framework import (get_peppercorn_args, 
                        set_handle_verbosity,
                        get_verification_crn,
                        get_verification_modules,
                        find_scheme_file,
                        enumerate_solution,
                        enumerate_modules,
                        verify,
                        verify_modules,
                        interpret_species,
                        assign_species)
from .dsdcompiler import translate, NuskellExit
from .ioutils import (natural_sort, 
                      write_pil,
                      load_pil,
                      get_strands,
                      write_vdsd)
from .crnutils import (parse_crn_string, Reaction, 
                       parse_crn_file, 
                       split_reversible_reactions, 
                       assign_crn_species,
                       removeSpecies,
                       removeTrivial,
                       genCRN)

def process_directories(crndir, schemedir):
    crns = []
    if crndir is None:
        log.info("Compiling single CRN:")
        input_crn = sys.stdin.readlines()
        input_crn = "".join(input_crn)
        crns.append((input_crn, input_crn))
        log.info("")
    else:
        if crndir[-1] != '/': 
            crndir += '/'
        log.info("Compiling CRNs in: {}".format(crndir))
        for crnfile in [x for x in os.listdir(crndir) if x[-4:] == '.crn']:
            log.info(' *{}'.format(crnfile))
            with open(crndir + crnfile) as cf:
                input_crn = ''.join(cf.readlines())
            crns.append((crnfile, input_crn))
        log.info("")

    if schemedir is None:
        log.info("Comparing default schemes:")
        schemedir = pkg_resources.resource_filename('nuskell', 'schemes') + '/'
    else:
        if schemedir[-1] != '/': schemedir += '/'
        log.info("Comparing schemes in: {}".format(schemedir))

    schemes = []
    for ts in [s for s in sorted(os.listdir(schemedir)) if s[-3:] == '.ts']:
        log.info('  {}'.format(ts))
        schemes.append(schemedir + ts)
    log.info("")
    return crns, schemes

def compare_schemes(crns, schemes, args = None):
    """ Compare different schemes for different CRNs.

    Args:
        crns (list[str]): A list of CRN stings.
        schemes (list[str]): A path to a directory with translation schemes.
        args (argparse(), optional): An object that contains arguments for peppercorn.
        verify (list, optional): An list of correctness notions.
        verify_timeout (int, optional): Timeout of verification [seconds].

    Returns:
        A dataframe.
    """
    plotdata = []  # Scheme, CRN, Cost, Speed
    for tsn in schemes:
        ts = find_scheme_file(tsn)
        for (name, input_crn) in crns:
            clear_memory()
            current = [name, tsn] # Make named tuple
            log.info(f"Compiling CRN {name} using translation scheme {ts}.")

            fcrn, fsc = parse_crn_string(input_crn)

            # TRANSLATE
            try:
                solution, modules = translate(input_crn, ts, modular = args.modular)
                fuels, wastes, intermediates, signals = assign_species(solution)
            except NuskellExit as e:
                log.error(e)
                log.error(f"Exiting translation for {name} and {ts}.")
                current.extend([None] * len(vnotions) + [None, None, None])
                plotdata.append(current)
                continue

            # ENUMERATE
            try:
                semantics = ['d'] if args.enum_detailed else ['c']
                complexes, reactions = enumerate_solution(solution, args)
                interpretation, complexes, reactions = interpret_species(complexes, 
                                                                reactions,
                                                                fsc.keys(),
                                                                prune = True)
                fuels, wastes, intermediates, signals = assign_species(complexes)
                if args.modular:
                    mcomplexes, mreactions = enumerate_modules(modules, 
                                                               interpretation,
                                                               complexes,
                                                               reactions,
                                                               args)
            except PolymerizationError as e:
                # NOTE: Many PolymerizationErrors are easiest to fix with
                # --reject-remote enumeration semantics, but that is not a
                # guarantee it that it will work.
                log.warning(f"Changing enumeration parameters to reject-remote ({name=}, {ts=}).")
                try:
                    args.reject_remote = True
                    complexes, reactions = enumerate_solution(solution, args)
                    interpretation, complexes, reactions = interpret_species(complexes, 
                                                                    reactions,
                                                                    fsc.keys(),
                                                                    prune = True)
                    fuels, wastes, intermediates, signals = assign_species(complexes)
                    if args.modular:
                        mcomplexes, mreactions = enumerate_modules(modules, 
                                                                   interpretation,
                                                                   complexes,
                                                                   reactions,
                                                                   args)
                    semantics.append('rr')
                    args.reject_remote = False
                except PolymerizationError as e:
                    log.error(e)
                    current.extend([None] * len(vnotions) + [None, None, None])
                    plotdata.append(current)
                    continue

            if not reactions:
                log.error("No DSD reactions have been enumerated ({name=}, {ts=}).")
                current.extend([None] * len(vnotions) + [None, None, None])
                plotdata.append(current)
                continue
            current.append(semantics)

            # VERIFY
            formals = set(fsc.keys())
            icrn, fuels, wastes = get_verification_crn(reactions, fuels, signals)
            if args.modular:
                fcrns, icrns = get_verification_modules(fcrn, mreactions, fuels, wastes)

            for meth in args.verify:
                if 'modular-' in meth and len(fcrns) > 1:
                    v, i = verify_modules(fcrns, icrns, formals, meth[8:], 
                                          interpretation = interpretation, 
                                          timeout = args.verify_timeout)
                else:
                    if 'modular' in meth: meth = meth[8:]
                    v, i = verify(fcrn, icrn, formals, meth, 
                                  interpretation = interpretation, 
                                  timeout = args.verify_timeout)
                current.append(v)

            cost = sum(sum(map(lambda d: d.length, s)) for s in get_strands(complexes))
            current.append(cost)
            speed = len(icrn)
            current.append(speed)
            plotdata.append(current)
    return plotdata

def get_nuskellCMP_args(parser):
    """ A collection of arguments for Nuskell """
    parser.add_argument('--version', action = 'version', version = '%(prog)s ' + __version__)
    parser.add_argument("-v", "--verbose", action = 'count', default = 0,
            help = "Print logging output. -vv increases verbosity level.")

    input = parser.add_argument_group('NuskellCMP Input Arguments')

    input.add_argument("--ts-dir", action='store', metavar='<path/to/dir>',
            help="""Specify path to the translation scheme directory. Only
            files that have a *.ts ending will be compared.""")

    input.add_argument("--crn-dir", action='store', metavar='<path/to/dir>',
            help="""Specify path to a CRN directory. Only files that have a
            *.crn ending will be compared.""")

    input.add_argument("--reference", action='store', metavar='<path/to/file>',
            help="Specify a translation scheme that serves as a reference.")

    input.add_argument("--from-csv", action='store', metavar='<path/to/file>',
            help="Read results from a *.CSV file. All other inputs are ignored.")

    default = parser.add_argument_group('NuskellCMP Output Arguments')

    default.add_argument("--pyplot", default='', action='store', metavar='<path/to/file>',
            help="Specify name of plot file. Choose from fileformats *.pdf or *.png")

    default.add_argument("--to-csv", action='store', default='', metavar='<path/to/file>',
            help="Print results to a *.CSV file.")

    # Choose a verification method.
    default.add_argument("--verify", nargs = '+', default = [], action = 'store',
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

    default.add_argument("-u", "--concentration-units", default='nM', action='store',
            choices=('M', 'mM', 'uM', 'nM', 'pM'),
            help="""Specify default concentration units when writing results to 
            ouptut files (reaction rates, initial concentrations). """)

    default.add_argument("--modular", action = 'store_true',
            #help="""After enumeration of the full system, enumerate individual
            #CRN modules separately, to identify crosstalk between reactions.
            #This is turned on automatically when using bisimulation
            #verification.""")
            help=argparse.SUPPRESS)

    default.add_argument("--verify-timeout", type = int, default = 30, metavar = '<int>',
            help="Specify time in seconds to wait for verification to complete.")

    return parser

def main():
    """Compare multiple tranlation schemes for a given CRN. """
    parser = argparse.ArgumentParser()
    get_nuskellCMP_args(parser)
    get_peppercorn_args(parser)
    args = parser.parse_args()

    # ~~~~~~~~~~~~~
    # Logging Setup 
    # ~~~~~~~~~~~~~
    title = "NuskellCMP - Comparison of translation schemes {}".format(__version__)
    if args.verbose >= 3: # Get root logger.
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
                '[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s')
    else:
        logger = logging.getLogger('nuskell')
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(levelname)s %(message)s')

    ch = logging.StreamHandler()
    set_handle_verbosity(ch, args.verbose)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info(title)
    logger.info("")

    # *************** #
    # Check arguments #
    # _______________ #
    if not args.modular:
        args.modular = any(map(lambda x: 'modular' in x, args.verify))

    # ***************** #
    # Process CSV input #
    # _________________ #
    if args.from_csv:
        logger.info('# Parsing data from file ... ')
        df = pd.DataFrame().from_csv(args.from_csv)
    else:
        crns, schemes = process_directories(args.crn_dir, args.ts_dir)

        # ********* #
        # MAIN LOOP #
        # _________ #

        plotdata = compare_schemes(crns, schemes, args = args)

        # Results:
        for row in plotdata:
            print(row)

        idx = list(zip(list(zip(*plotdata))[0], list(zip(*plotdata))[1]))
        df = pd.DataFrame(plotdata, index = idx,
                          columns = ['Translation scheme', 
                                     'CRN', 
                                     'enumerated'] + args.verify + [
                                     'number of nucleotides', 
                                     'reactions in condensed network'])
    # Save to portable format:
    if args.to_csv:
        df.to_csv(path_or_buf=args.to_csv)

    def equiv(x):
        if True in x:
            return True
        elif (nan in x) or (None in x):
            return 'timeout'
        else:
            return False

    # Add column that combines equivalence notions.
    e = list(map(equiv, zip(*map(lambda x: df[x], args.verify))))
    df['equivalent'] = e
    df = df.sort_values(['Translation scheme', 'CRN'], ascending = [True, True])
    print(df.to_string(index=False, justify='left'))

    if args.pyplot:
        single_plot(df, pfile=args.pyplot)

if __name__ == '__main__':
   main()

