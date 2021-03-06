#!/usr/bin/env python
#
#  compare_schemes.py
#  EnumeratorProject
#
import logging
log = logging.getLogger(__name__)

import gc
import sys
import argparse
from natsort import natsorted
from itertools import chain
from peppercornenumerator import PolymerizationError
from peppercornenumerator.objects import show_memory as pepper_memory

from . import __version__
from .objects import NuskellDomain, NuskellComplex, show_memory
from .ioutils import get_strands
from .crnutils import parse_crn_string 
from .dsdcompiler import translate, get_builtin_schemes, NuskellExit
from .framework import (get_peppercorn_args, 
                        set_handle_verbosity,
                        enumerate_solution,
                        enumerate_modules,
                        interpret_species,
                        get_verification_crn,
                        get_verification_modules,
                        verify,
                        verify_modules,
                        assign_species)

def process_input(crns, schemes):
    mycrns = []
    if crns is None:
        log.info("Reading single CRN from STDIN.")
        input_crn = sys.stdin.readlines()
        input_crn = "".join(input_crn)
        mycrns.append((input_crn, input_crn))
    else:
        for crnfile in crns:
            log.info(f'Processing {crnfile=}.')
            with open(crnfile) as cf:
                input_crn = ''.join(cf.readlines())
            mycrns.append((crnfile, input_crn))
    log.info("")

    if schemes is None:
        schemes = get_builtin_schemes()
        schemes = list(chain(*(schemes.values())))
    return mycrns, schemes

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
    plotdata = [['scheme', 'CRN', 'enumerated', '# nuc', '# rxns'
                ] + args.verify + ['equivalent']]

    for ts in schemes:
        for (name, input_crn) in crns:
            assert list(show_memory()) == []
            assert list(pepper_memory()) == []
            NuskellDomain.ID = 1
            NuskellComplex.ID = 1
            log.info(f"Compiling CRN {name=} using translation scheme {ts=}.")
            current = [ts, name] # TODO: make named tuple?
            fcrn, fsc = parse_crn_string(input_crn)

            try: # TRANSLATE
                solution, modules = translate(input_crn, ts, modular = args.modular)
            except NuskellExit as e:
                log.error(e)
                log.error(f"Exiting translation for {name} and {ts}.")
                # Garbage collection to clear singleton objects.
                gc.collect()
                # Return empty results
                current.extend([None, None, None] + [None] * (len(args.verify)+1))
                plotdata.append(current)
                continue
            
            fuels, wastes, intermediates, signals = assign_species(solution)

            log.debug(f"Formal species: {', '.join(natsorted(fsc))}")
            log.debug("Signal Complexes:\n" + '\n'.join(
                [f'   {cplx.name} = {cplx.kernel_string}' for cplx in natsorted(signals, lambda x: x.name)]))
            log.debug("Fuel Complexes:\n" + '\n'.join(
                [f'   {cplx.name} = {cplx.kernel_string}' for cplx in natsorted(fuels, lambda x: x.name)]))

            log.info(f"Enumerating CRN {name=} using translation scheme {ts=}.")
            assert args.reject_remote is False
            semantics = ['d'] if args.enum_detailed else ['c']
            complexes = dict()
            reactions = set()
            try: # ENUMERATE
                complexes, reactions = enumerate_solution(solution, args)
            except PolymerizationError as e:
                complexes.clear()
                reactions.clear()
                # NOTE: Many PolymerizationErrors are easiest to fix with
                # --reject-remote enumeration semantics, but that is not a
                # guarantee it that it will work.
                log.warning(f"Changing enumeration parameters to reject-remote ({name=}, {ts=}).")
                args.reject_remote = True
                semantics.append('rr')
                try:
                    complexes, reactions = enumerate_solution(solution, args)
                except PolymerizationError as e:
                    log.error(e)
                    log.error(f"Exiting enumeration for {name} and {ts}.")
                    # Garbage collection to clear singleton objects.
                    complexes.clear()
                    reactions.clear()
                    solution.clear()
                    [m.clear() for m in modules]
                    del fuels, wastes, intermediates, signals
                    gc.collect()
                    assert list(show_memory()) == []
                    # Append results
                    args.reject_remote = False
                    current.extend([None, None, None] + [None] * (len(args.verify)+1))
                    plotdata.append(current)
                    continue

            if not reactions:
                log.error("No DSD reactions have been enumerated ({name=}, {ts=}).")
                # Garbage collection to clear singleton objects.
                solution.clear()
                [m.clear() for m in modules]
                del fuels, wastes, intermediates, signals
                complexes.clear()
                gc.collect()
                assert list(show_memory()) == []
                current.extend([None, None, None] + [None] * (len(args.verify)+1))
                plotdata.append(current)
                continue

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
            # Restore to default, just in case.
            args.reject_remote = False

            current.append(':'.join(semantics))
            cost = sum(sum(map(lambda d: d.length, s)) for s in get_strands(complexes))
            current.append(cost)
            speed = len(reactions)
            current.append(speed)
 
            # VERIFY
            formals = set(fsc.keys())
            icrn, fuels, wastes = get_verification_crn(reactions, fuels, signals)
            if args.modular:
                fcrns, icrns = get_verification_modules(fcrn, mreactions, fuels, wastes)

            equiv = False
            for meth in args.verify:
                if 'modular-' in meth and len(fcrns) > 1:
                    try:
                        v, i = verify_modules(fcrns, icrns, formals, meth[8:], 
                                          interpretation = interpretation, 
                                          timeout = args.verify_timeout)
                    except NotImplementedError:
                        v = i = '-'
                else:
                    if 'modular' in meth: meth = meth[8:]
                    v, i = verify(fcrn, icrn, formals, meth, 
                                  interpretation = interpretation, 
                                  timeout = args.verify_timeout)
                current.append(v)
                if equiv is True:
                    pass
                elif v is True:
                    equiv = True
                elif v is None:
                    equiv = 'timeout'
            current.append(equiv)
            plotdata.append(current)

            solution.clear()
            [m.clear() for m in modules]
            del fuels, wastes, intermediates, signals
            complexes.clear()
            reactions.clear()
            interpretation.clear()
            if args.modular:
                [m.clear() for m in mreactions]
                [m.clear() for m in mcomplexes]
            gc.collect()
            assert list(show_memory()) == []
            log.info(f"Done with CRN {name=} using translation scheme {ts=}. Moving on ...")
    return plotdata

def get_nuskellCMP_args(parser):
    """ A collection of arguments for NuskellCMP """
    parser.add_argument('--version', action = 'version', 
                        version = '%(prog)s ' + __version__)
    parser.add_argument("-v", "--verbose", action = 'count', default = 0,
            help = "Print logging output. -vv increases verbosity level.")

    inp = parser.add_argument_group('NuskellCMP Input Arguments')
    inp.add_argument("--schemes", action='store', nargs = '+', 
            help="""A list of one or more translation schemes.""")
    inp.add_argument("--crns", action='store', nargs = '+', 
            help="""A list of one or more CRNs.""")

    # Choose a verification method.
    out = parser.add_argument_group('NuskellCMP Output Arguments')
    out.add_argument("--verify", nargs = '+', action = 'store',
            default = ['crn-bisimulation', 
                       'modular-crn-bisimulation', 
                       'pathway-decomposition', 
                       'compositional-hybrid',
                       'integrated-hybrid'], 
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
    out.add_argument("--modular", action = 'store_true',
            help=argparse.SUPPRESS)
    out.add_argument("--verify-timeout", type = int, default = 30, metavar = '<int>',
            help="Specify time in seconds to wait for verification to complete.")
    return parser

def parse_args(args):
    parser = argparse.ArgumentParser()
    get_nuskellCMP_args(parser)
    get_peppercorn_args(parser)
    args = parser.parse_args(args)
    if not args.modular:
        args.modular = any(map(lambda x: 'modular' in x, args.verify))
    return args

def main():
    """Compare multiple tranlation schemes for a given CRN. """
    args = parse_args(sys.argv[1:])
    # ~~~~~~~~~~~~~
    # Logging Setup 
    # ~~~~~~~~~~~~~
    title = f"NuskellCMP - Comparison of translation schemes {__version__}"
    if args.verbose >= 3: # Get root logger.
        args.verbose -= 2
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

    # ***************** #
    # Process CSV input #
    # _________________ #
    crns, schemes = process_input(args.crns, args.schemes)

    # ********* #
    # MAIN LOOP #
    # _________ #

    plotdata = compare_schemes(crns, schemes, args = args)

    # Results:
    for row in plotdata:
        print(row)

if __name__ == '__main__':
   main()

