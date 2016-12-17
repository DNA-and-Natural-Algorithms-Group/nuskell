#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
  readme = f.read()

with open('LICENSE') as f:
  license = f.read()

# Dynamically figure out the version
version = __import__('nuskell').__version__

setup(
    name='nuskell',
    version=version,
    description='Nucleic acid strand displacement compiler',
    long_description=readme,
    author='Seung Woo Shin, Qing Dong, Robert Johnson, Stefan Badelt, Erik Winfree',
    author_email='winfree@caltech.edu',
    #url='http://www.dna.caltech.edu/nuskell/',
    data_files=[('nuskell/schemes', 
        ['schemes/generalized/soloveichik2010_v1.ts',
        #'schemes/generalized/soloveichik2010_v2.ts', 
         'schemes/original/qian2011.ts', 
         'schemes/generalized/qian2011_gen.ts',
         'schemes/original/cardelli2011_NM.ts',
         'schemes/original/cardelli2011_NM_noGC.ts',
         'schemes/original/cardelli2011_FJ.ts',
         'schemes/original/cardelli2011_FJ_noGC.ts',
         'schemes/original/cardelli2013_2D.ts'])],
    license=license,
    install_requires=['pyparsing>=1.5.5', 'argparse>=1.2.1'],
    test_suite='tests',
    packages=['nuskell', 'nuskell.parser', 'nuskell.interpreter', 'nuskell.verifier', 
        'nuskell.include', 'nuskell.include.peppercorn'],
    scripts=['scripts/nuskell']
)

