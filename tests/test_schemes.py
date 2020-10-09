#!/usr/bin/env python
#
#  test_schemes.py
#  NuskellCompilerProject
#
import warnings
import os
import filecmp
import unittest
import subprocess as sub

#
# nuskellCMP is a script used to compare schemes that are distributed with
# nuskell. These tests generate snapshots, comparing enumeration and
# verification of *all* schemes in "schemes/literature/" and
# "schemes/variants/" for some selected CRNs in "tests/crns/". Outputs
# are written and compared to "tests/snapshots/".
#
# This is not a regular unittest, run overnight to check consistency of results
# before every release.
#

try: 
    import pandas as pd
except ImportError as err:
    log.warning('')
SKIP = False

    #if args.from_csv:
    #    logger.info('# Parsing data from file ... ')
    #    df = pd.read_csv(args.from_csv, index_col = 0, na_values = None)
    #else:
    #df = pd.DataFrame(plotdata, columns = dfheader)

    ## Save to portable format:
    #if args.to_csv:
    #    df.to_csv(path_or_buf=args.to_csv)

    #print(df.to_string(index=False, justify='left'))


@unittest.skipIf(SKIP, "slow tests are disabled by default")
class SinlgeSnapshotCMP(unittest.TestCase):
    def setUp(self):
        self.exe = 'nuskellCMP'
        self.lit = 'schemes/literature'
        self.var = 'schemes/variants'

        self.call = [self.exe]
        self.call.extend(['--verbose'])
        self.call.extend(['--verify', 'pathway', 'integrated', 'bisimulation', 'modular-bisimulation'])
        self.call.extend(['--verify-timeout',  str(30)])

    def test_oscillator_01_lit(self):
        # Exisiting Data
        cmp_file = 'tests/snapshots/oscillator_01_lit.cmp'

        # Generated Data
        new_file = 'tests/snapshots/oscillator_01_lit.new'
        err_file = 'tests/snapshots/oscillator_01_lit.err'

        # Input Schemes
        ts_dir = self.lit

        # Input CRNs
        crn_file = 'tests/crns/oscillator_01.crn'

        # Options
        call = self.call
        call.extend(['--ts-dir', ts_dir])
        call.extend(['--max-complex-size',      str(10)])
        call.extend(['--max-complex-count',   str(5000)])
        call.extend(['--max-reaction-count', str(10000)])

        print(crn_file)
        print(call)

        with open(crn_file) as crn, \
                open(new_file, 'w') as out, \
                open(err_file, 'w') as err :

            proc = sub.Popen(call, stdin=crn, stdout=out, stderr=err)
            proc.communicate(None)
            self.assertEqual(proc.returncode, 0)

        self.assertTrue(filecmp.cmp(cmp_file, new_file, shallow=False))
        os.remove(new_file)
        os.remove(err_file)

    def test_oscillator_01_var(self):
        # Exisiting Data
        cmp_file = 'tests/snapshots/oscillator_01_var.cmp'

        # Generated Data
        new_file = 'tests/snapshots/oscillator_01_var.new'
        err_file = 'tests/snapshots/oscillator_01_var.err'

        # Input Schemes
        ts_dir = self.var

        # Input CRNs
        crn_file = 'tests/crns/oscillator_01.crn'

        # Options
        call = self.call
        call.extend(['--ts-dir', ts_dir])
        call.extend(['--max-complex-size',      str(20)])
        call.extend(['--max-complex-count',   str(5000)])
        call.extend(['--max-reaction-count', str(10000)])

        print(crn_file)
        print(call)

        with open(crn_file) as crn, \
                open(new_file, 'w') as out, \
                open(err_file, 'w') as err :

            proc = sub.Popen(call, stdin=crn, stdout=out, stderr=err)
            proc.communicate(None)
            self.assertEqual(proc.returncode, 0)

        self.assertTrue(filecmp.cmp(cmp_file, new_file, shallow=False))
        os.remove(new_file)
        os.remove(err_file)

    def test_oscillator_02_lit(self):
        # Exisiting Data
        cmp_file = 'tests/snapshots/oscillator_02_lit.cmp'

        # Generated Data
        new_file = 'tests/snapshots/oscillator_02_lit.new'
        err_file = 'tests/snapshots/oscillator_02_lit.err'

        # Input Schemes
        ts_dir = self.lit

        # Input CRNs
        crn_file = 'tests/crns/oscillator_02.crn'

        # Options
        call = self.call
        call.extend(['--ts-dir', ts_dir])
        call.extend(['--max-complex-size',      str(20)])
        call.extend(['--max-complex-count',   str(5000)])
        call.extend(['--max-reaction-count', str(10000)])

        print(crn_file)
        print(call)

        with open(crn_file) as crn, \
                open(new_file, 'w') as out, \
                open(err_file, 'w') as err :

            proc = sub.Popen(call, stdin=crn, stdout=out, stderr=err)
            proc.communicate(None)
            self.assertEqual(proc.returncode, 0)

        self.assertTrue(filecmp.cmp(cmp_file, new_file, shallow=False))
        os.remove(new_file)
        os.remove(err_file)

    def test_oscillator_02_var(self):
        # Exisiting Data
        cmp_file = 'tests/snapshots/oscillator_02_var.cmp'

        # Generated Data
        new_file = 'tests/snapshots/oscillator_02_var.new'
        err_file = 'tests/snapshots/oscillator_02_var.err'

        # Input Schemes
        ts_dir = self.var

        # Input CRNs
        crn_file = 'tests/crns/oscillator_02.crn'

        # Options
        call = self.call
        call.extend(['--ts-dir', ts_dir])
        call.extend(['--max-complex-size',      str(20)])
        call.extend(['--max-complex-count',   str(5000)])
        call.extend(['--max-reaction-count', str(10000)])

        print(crn_file)
        print(call)

        with open(crn_file) as crn, \
                open(new_file, 'w') as out, \
                open(err_file, 'w') as err :

            proc = sub.Popen(call, stdin=crn, stdout=out, stderr=err)
            proc.communicate(None)
            self.assertEqual(proc.returncode, 0)

        self.assertTrue(filecmp.cmp(cmp_file, new_file, shallow=False))
        os.remove(new_file)
        os.remove(err_file)

    def test_roessler_01_lit(self):
        # Exisiting Data
        cmp_file = 'tests/snapshots/roessler_01_lit.cmp'

        # Generated Data
        new_file = 'tests/snapshots/roessler_01_lit.new'
        err_file = 'tests/snapshots/roessler_01_lit.err'

        # Input Schemes
        ts_dir = self.lit

        # Input CRNs
        crn_file = 'tests/crns/roessler_01.crn'

        # Options
        call = self.call
        call.extend(['--ts-dir', ts_dir])
        call.extend(['--max-complex-size',      str(20)])
        call.extend(['--max-complex-count',   str(5000)])
        call.extend(['--max-reaction-count', str(10000)])

        print(crn_file)
        print(call)

        with open(crn_file) as crn, \
                open(new_file, 'w') as out, \
                open(err_file, 'w') as err :

            proc = sub.Popen(call, stdin=crn, stdout=out, stderr=err)
            proc.communicate(None)
            self.assertEqual(proc.returncode, 0)

        self.assertTrue(filecmp.cmp(cmp_file, new_file, shallow=False))
        os.remove(new_file)
        os.remove(err_file)

    def test_roessler_01_var(self):
        # Exisiting Data
        cmp_file = 'tests/snapshots/roessler_01_var.cmp'

        # Generated Data
        new_file = 'tests/snapshots/roessler_01_var.new'
        err_file = 'tests/snapshots/roessler_01_var.err'

        # Input Schemes
        ts_dir = self.var

        # Input CRNs
        crn_file = 'tests/crns/roessler_01.crn'

        # Options
        call = self.call
        call.extend(['--ts-dir', ts_dir])
        call.extend(['--max-complex-size',      str(20)])
        call.extend(['--max-complex-count',   str(5000)])
        call.extend(['--max-reaction-count', str(10000)])

        print(crn_file)
        print(call)

        with open(crn_file) as crn, \
                open(new_file, 'w') as out, \
                open(err_file, 'w') as err :

            proc = sub.Popen(call, stdin=crn, stdout=out, stderr=err)
            proc.communicate(None)
            self.assertEqual(proc.returncode, 0)

        self.assertTrue(filecmp.cmp(cmp_file, new_file, shallow=False))
        os.remove(new_file)
        os.remove(err_file)

    def test_bin_counter_01_lit(self):
        # Exisiting Data
        cmp_file = 'tests/snapshots/bin_counter_01_lit.cmp'

        # Generated Data
        new_file = 'tests/snapshots/bin_counter_01_lit.new'
        err_file = 'tests/snapshots/bin_counter_01_lit.err'

        # Input Schemes
        ts_dir = self.lit

        # Input CRNs
        crn_file = 'tests/crns/bin_counter_01.crn'

        # Options
        call = self.call
        call.extend(['--ts-dir', ts_dir])
        call.extend(['--max-complex-size',      str(20)])
        call.extend(['--max-complex-count',   str(5000)])
        call.extend(['--max-reaction-count', str(10000)])

        print(crn_file)
        print(call)

        with open(crn_file) as crn, \
                open(new_file, 'w') as out, \
                open(err_file, 'w') as err :

            proc = sub.Popen(call, stdin=crn, stdout=out, stderr=err)
            proc.communicate(None)
            self.assertEqual(proc.returncode, 0)

        self.assertTrue(filecmp.cmp(cmp_file, new_file, shallow=False))
        os.remove(new_file)
        os.remove(err_file)

    def test_bin_counter_01_var(self):
        # Exisiting Data
        cmp_file = 'tests/snapshots/bin_counter_01_var.cmp'

        # Generated Data
        new_file = 'tests/snapshots/bin_counter_01_var.new'
        err_file = 'tests/snapshots/bin_counter_01_var.err'

        # Input Schemes
        ts_dir = self.var

        # Input CRNs
        crn_file = 'tests/crns/bin_counter_01.crn'

        # Options
        call = self.call
        call.extend(['--ts-dir', ts_dir])
        call.extend(['--max-complex-size',      str(20)])
        call.extend(['--max-complex-count',   str(5000)])
        call.extend(['--max-reaction-count', str(10000)])

        print(crn_file)
        print(call)

        with open(crn_file) as crn, \
                open(new_file, 'w') as out, \
                open(err_file, 'w') as err :

            proc = sub.Popen(call, stdin=crn, stdout=out, stderr=err)
            proc.communicate(None)
            self.assertEqual(proc.returncode, 0)

        self.assertTrue(filecmp.cmp(cmp_file, new_file, shallow=False))
        os.remove(new_file)
        os.remove(err_file)

@unittest.skipIf(SKIP, "slow tests are disabled by default")
class MultiSnapshotCMP(unittest.TestCase):
    def setUp(self):
        self.exe = 'scripts/nuskellCMP'
        self.lit = 'schemes/literature'
        self.var = 'schemes/variants'

        self.call = [self.exe]
        self.call.extend(['--verbose'])
        self.call.extend(['--verify', 'pathway', 'integrated', 'bisimulation', 'modular-bisimulation'])
        self.call.extend(['--verify-timeout',  str(30)])

    def test_directory_bimolecular_lit(self):
        # Exisiting Data
        cmp_file = 'tests/snapshots/bimolecular_lit.cmp'

        # Generated Data
        new_file = 'tests/snapshots/bimolecular_lit.new'
        err_file = 'tests/snapshots/bimolecular_lit.err'

        # Input Schemes
        ts_dir = self.lit

        # Input CRNs
        crn_dir = 'tests/crns/bimol/'

        # Options
        call = self.call
        call.extend(['--ts-dir', ts_dir])
        call.extend(['--crn-dir', crn_dir])
        call.extend(['--max-complex-size',      str(50)])
        call.extend(['--max-complex-count',   str(5000)])
        call.extend(['--max-reaction-count', str(10000)])

        print(call)

        with open(new_file, 'w') as out, \
                open(err_file, 'w') as err :

            proc = sub.Popen(call, stdin=None, stdout=out, stderr=err)
            proc.communicate(None)
            self.assertEqual(proc.returncode, 0)

        self.assertTrue(filecmp.cmp(cmp_file, new_file, shallow=False))
        os.remove(new_file)
        os.remove(err_file)

    def test_directory_bimolecular_var(self):
        # Exisiting Data
        cmp_file = 'tests/snapshots/bimolecular_var.cmp'

        # Generated Data
        new_file = 'tests/snapshots/bimolecular_var.new'
        err_file = 'tests/snapshots/bimolecular_var.err'

        # Input Schemes
        ts_dir = self.var

        # Input CRNs
        crn_dir = 'tests/crns/bimol/'

        # Options
        call = self.call
        call.extend(['--ts-dir', ts_dir])
        call.extend(['--crn-dir', crn_dir])
        call.extend(['--max-complex-size',      str(50)])
        call.extend(['--max-complex-count',   str(5000)])
        call.extend(['--max-reaction-count', str(10000)])

        print(call)

        with open(new_file, 'w') as out, \
                open(err_file, 'w') as err :

            proc = sub.Popen(call, stdin=None, stdout=out, stderr=err)
            proc.communicate(None)
            self.assertEqual(proc.returncode, 0)

        self.assertTrue(filecmp.cmp(cmp_file, new_file, shallow=False))
        os.remove(new_file)
        os.remove(err_file)

    def test_directory_basic_lit(self):
        # Exisiting Data
        cmp_file = 'tests/snapshots/basic_lit.cmp'

        # Generated Data
        new_file = 'tests/snapshots/basic_lit.new'
        err_file = 'tests/snapshots/basic_lit.err'

        # Input Schemes
        ts_dir = self.lit

        # Input CRNs
        crn_dir = 'tests/crns/basic/'

        # Options
        call = self.call
        call.extend(['--ts-dir', ts_dir])
        call.extend(['--crn-dir', crn_dir])
        call.extend(['--max-complex-size',      str(50)])
        call.extend(['--max-complex-count',   str(5000)])
        call.extend(['--max-reaction-count', str(10000)])

        print(call)

        with open(new_file, 'w') as out, \
                open(err_file, 'w') as err :

            proc = sub.Popen(call, stdin=None, stdout=out, stderr=err)
            proc.communicate(None)
            self.assertEqual(proc.returncode, 0)

        self.assertTrue(filecmp.cmp(cmp_file, new_file, shallow=False))
        os.remove(new_file)
        os.remove(err_file)

    def test_directory_basic_var(self):
        # Exisiting Data
        cmp_file = 'tests/snapshots/basic_var.cmp'

        # Generated Data
        new_file = 'tests/snapshots/basic_var.new'
        err_file = 'tests/snapshots/basic_var.err'

        # Input Schemes
        ts_dir = self.var

        # Input CRNs
        crn_dir = 'tests/crns/basic/'

        # Options
        call = self.call
        call.extend(['--ts-dir', ts_dir])
        call.extend(['--crn-dir', crn_dir])
        call.extend(['--max-complex-size',      str(50)])
        call.extend(['--max-complex-count',   str(5000)])
        call.extend(['--max-reaction-count', str(10000)])

        print(call)

        with open(new_file, 'w') as out, \
                open(err_file, 'w') as err :

            proc = sub.Popen(call, stdin=None, stdout=out, stderr=err)
            proc.communicate(None)
            self.assertEqual(proc.returncode, 0)

        self.assertTrue(filecmp.cmp(cmp_file, new_file, shallow=False))
        os.remove(new_file)
        os.remove(err_file)

if __name__ == '__main__':
    unittest.main()
