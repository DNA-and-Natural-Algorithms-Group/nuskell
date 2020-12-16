#!/usr/bin/env python
#
#  tests/test_schemes.py
#  NuskellCompilerProject
#
"""
nuskellCMP is a script used to compare schemes that are distributed with
nuskell. These tests generate snapshots, comparing enumeration and
verification of *all* builtin schemes for some selected CRNs in
"tests/crns/". Outputs are written and compared to "tests/snapshots/".

This is not a regular unittest, run overnight to check consistency of results
before every release.
"""
import logging
# logger = logging.getLogger('nuskell')
# logger.setLevel(logging.INFO)
# fh = logging.FileHandler('tests/test_schemes.log')
# fh.setLevel(logging.INFO)
# formatter = logging.Formatter('[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s')
# fh.setFormatter(formatter)
# logger.addHandler(fh)

import os
import unittest
from itertools import chain

from peppercornenumerator.objects import clear_memory as clear_pepper_memory
import nuskell.dsdcompiler.compiler as comp
from nuskell.objects import clear_memory
from nuskell.dsdcompiler import get_builtin_schemes
from nuskell.compare_schemes import (parse_args, 
                                     process_input,
                                     compare_schemes)

SKIP = False # These are only skipped for debugging, or if pandas is missing.
SKIP_SLOW = True # These here are slow tests, which should be skipped by default!

try: 
    import pandas as pd
except ImportError as err:
    log.warning(f'Unittest: {__name__} needs pandas installed.')
    SKIP = True

def compare_snapshots(cmp_file, new_file, crns, schemes, args = None):
    """ A helper function to compare compilations.

    This function does not assert False if the results differ. Instead, it
    shows the results whenever there is a problem.
    """
    calculate = True
    display = True

    if calculate:
        if args is None:
            args = parse_args(['--max-complex-size', '50',
                               '--max-complex-count', '10000',
                               '--max-reaction-count', '50000'])
        crns, schemes = process_input(crns, schemes)
        output = compare_schemes(crns, schemes, args)
        dfhead = output[0]
        dfdata = output[1:]
        df = pd.DataFrame(dfdata, columns = dfhead)
        df.to_csv(path_or_buf = new_file)

    if display:
        df_cmp = pd.read_csv(cmp_file, index_col = 0, na_values = None)
        df_new = pd.read_csv(new_file, index_col = 0, na_values = None)
        same = df_cmp.equals(df_new)
        if not same:
            print()
            print('Previous:')
            print(df_cmp.where(df_cmp.notnull(), None).to_string(index=False, justify='left'))
            print('New:')
            print(df_new.where(df_new.notnull(), None).to_string(index=False, justify='left'))

@unittest.skipIf(SKIP, "missing pandas requirement.")
class QuickSnapshotCMP(unittest.TestCase):
    def tearDown(self):
        clear_memory()
        clear_pepper_memory()

    def test_small_nuskellCMP(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/small_test.csv'
        new_file = 'tests/snapshots/small_test.new'
        # Input 
        schemes = ['soloveichik2010.ts', 'qian2011_3D_var1.ts']
        crns = ['tests/crns/binary/irr/irr_bin_14.crn', 
                'tests/crns/binary/rev/rev_bin_14.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    @unittest.skipIf(SKIP_SLOW, "slow tests are disabled by default")
    def test_bigger_nuskellCMP(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/bigger_test.csv'
        new_file = 'tests/snapshots/bigger_test.new'
        # Input 
        schemes = ['lakin2016_2D_3I.ts', 'cardelli2013_2D_3I.ts']
        crns = ['tests/crns/oscillator_01.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

@unittest.skipIf(SKIP or SKIP_SLOW, "slow tests are disabled by default")
class SystemSnapshotTests(unittest.TestCase):
    def setUp(self):
        comp.SCHEME_DIRS = ['schemes/literature/'] 
        self.lit = list(chain(*get_builtin_schemes().values()))
        comp.SCHEME_DIRS = ['schemes/variants/'] 
        self.var = list(chain(*get_builtin_schemes().values()))
        comp.SCHEME_DIRS = ['schemes/literature/', 'schemes/variants/'] 

    def tearDown(self):
        comp.SCHEME_DIRS = ['schemes/literature/', 'schemes/variants/'] 
        clear_memory()
        clear_pepper_memory()

    def test_oscillator_01_lit(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/oscillator_01_lit.csv'
        new_file = 'tests/snapshots/oscillator_01_lit.new'
        # Input 
        schemes = self.lit
        crns = ['tests/crns/oscillator_01.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_oscillator_01_var(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/oscillator_01_var.csv'
        new_file = 'tests/snapshots/oscillator_01_var.new'
        # Input 
        schemes = self.var
        crns = ['tests/crns/oscillator_01.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_oscillator_02_lit(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/oscillator_02_lit.csv'
        new_file = 'tests/snapshots/oscillator_02_lit.new'
        # Input 
        schemes = self.lit
        crns = ['tests/crns/oscillator_02.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_oscillator_02_var(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/oscillator_02_var.csv'
        new_file = 'tests/snapshots/oscillator_02_var.new'
        # Input 
        schemes = self.var
        crns = ['tests/crns/oscillator_02.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_roessler_01_lit(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/roessler_01_lit.csv'
        new_file = 'tests/snapshots/roessler_01_lit.new'
        # Input 
        schemes = self.lit
        crns = ['tests/crns/roessler_01.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_roessler_01_var(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/roessler_01_var.csv'
        new_file = 'tests/snapshots/roessler_01_var.new'
        # Input 
        schemes = self.var
        crns = ['tests/crns/roessler_01.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_bin_counter_01_lit(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/bin_counter_01_lit.csv'
        new_file = 'tests/snapshots/bin_counter_01_lit.new'
        # Input
        schemes = self.lit
        crns = ['tests/crns/bin_counter_01.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_bin_counter_01_var(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/bin_counter_01_var.csv'
        new_file = 'tests/snapshots/bin_counter_01_var.new'
        # Input
        schemes = self.var
        crns = ['tests/crns/bin_counter_01.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

@unittest.skipIf(SKIP or SKIP_SLOW, "slow tests are disabled by default")
class IrreversibleRxnSnapshotTest(unittest.TestCase):
    def setUp(self):
        comp.SCHEME_DIRS = ['schemes/literature/'] 
        self.lit = list(chain(*get_builtin_schemes().values()))
        comp.SCHEME_DIRS = ['schemes/variants/'] 
        self.var = list(chain(*get_builtin_schemes().values()))
        comp.SCHEME_DIRS = ['schemes/literature/', 'schemes/variants/'] 

    def tearDown(self):
        comp.SCHEME_DIRS = ['schemes/literature/', 'schemes/variants/'] 
        clear_memory()
        clear_pepper_memory()

    def test_basic_lit(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/basic_irr_lit.csv'
        new_file = 'tests/snapshots/basic_irr_lit.new'
        # Input 
        crndir = 'tests/crns/basic/irr/'
        schemes = self.lit
        crns = [crndir + x for x in sorted(os.listdir(crndir)) if x[-4:] == '.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_basic_var(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/basic_irr_var.csv'
        new_file = 'tests/snapshots/basic_irr_var.new'
        # Input 
        crndir = 'tests/crns/basic/irr/'
        schemes = self.var
        crns = [crndir + x for x in sorted(os.listdir(crndir)) if x[-4:] == '.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_binary_lit(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/binary_irr_lit.csv'
        new_file = 'tests/snapshots/binary_irr_lit.new'
        # Input 
        crndir = 'tests/crns/binary/irr/'
        schemes = self.lit
        crns = [crndir + x for x in sorted(os.listdir(crndir)) if x[-4:] == '.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_binary_var(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/binary_irr_var.csv'
        new_file = 'tests/snapshots/binary_irr_var.new'
        # Input 
        crndir = 'tests/crns/binary/irr/'
        schemes = self.var
        crns = [crndir + x for x in sorted(os.listdir(crndir)) if x[-4:] == '.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

@unittest.skipIf(SKIP or SKIP_SLOW, "slow tests are disabled by default")
class ReversibleRxnSnapshotTest(unittest.TestCase):
    def setUp(self):
        comp.SCHEME_DIRS = ['schemes/literature/'] 
        self.lit = list(chain(*get_builtin_schemes().values()))
        comp.SCHEME_DIRS = ['schemes/variants/'] 
        self.var = list(chain(*get_builtin_schemes().values()))
        comp.SCHEME_DIRS = ['schemes/literature/', 'schemes/variants/'] 

    def tearDown(self):
        comp.SCHEME_DIRS = ['schemes/literature/', 'schemes/variants/'] 
        clear_memory()
        clear_pepper_memory()

    def test_basic_lit(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/basic_rev_lit.csv'
        new_file = 'tests/snapshots/basic_rev_lit.new'
        # Input 
        crndir = 'tests/crns/basic/rev/'
        schemes = self.lit
        crns = [crndir + x for x in sorted(os.listdir(crndir)) if x[-4:] == '.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_basic_var(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/basic_rev_var.csv'
        new_file = 'tests/snapshots/basic_rev_var.new'
        # Input 
        crndir = 'tests/crns/basic/rev/'
        schemes = self.var
        crns = [crndir + x for x in sorted(os.listdir(crndir)) if x[-4:] == '.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_binary_lit(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/binary_rev_lit.csv'
        new_file = 'tests/snapshots/binary_rev_lit.new'
        # Input 
        crndir = 'tests/crns/binary/rev/'
        schemes = self.lit
        crns = [crndir + x for x in sorted(os.listdir(crndir)) if x[-4:] == '.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

    def test_binary_var(self):
        # Exisiting & generated data files.
        cmp_file = 'tests/snapshots/binary_rev_var.csv'
        new_file = 'tests/snapshots/binary_rev_var.new'
        # Input 
        crndir = 'tests/crns/binary/rev/'
        schemes = self.var
        crns = [crndir + x for x in sorted(os.listdir(crndir)) if x[-4:] == '.crn']
        # Output
        compare_snapshots(cmp_file, new_file, crns, schemes)

if __name__ == '__main__':
    unittest.main()
