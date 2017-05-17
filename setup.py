#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
  readme = f.read()

with open('LICENSE') as f:
  license = f.read()

# Dynamically figure out the version
version = __import__('nuskell').__version__

# Canonical: These schemes should be installed by default ...
canonical = [ 'schemes/canonical/soloveichik2010_v1.ts', 
              #'schemes/canonical/soloveichik2010_v1.ts', # optimized version?
              'schemes/canonical/qian2011_3D_v1.ts', # SWS
              #'schemes/canonical/qian2011_v2.ts', # as shown in the revised paper
              #'schemes/canonical/cardelli2011_NM_noGC.ts',
              #'schemes/canonical/cardelli2013_2D_3I_noGC.ts',
              #'schemes/literature/chen2013_2D_JF_3I_noGC.ts',
              'schemes/canonical/srinivas2015.ts']
              #'schemes/literature/lakin2016_2D_3I_noGC.ts',

# Literature: These schemes should be optional (educational) ...
literature = ['schemes/literature/soloveichik2010.ts',
              'schemes/literature/qian2011_3D.ts',
              'schemes/literature/cardelli2011_FJ.ts',
              'schemes/literature/cardelli2011_FJ_noGC.ts',
              'schemes/literature/cardelli2011_NM.ts',
              'schemes/literature/cardelli2011_NM_noGC.ts',
              'schemes/literature/lakin2012_3D.ts',
              'schemes/literature/cardelli2013_2D.ts',
              'schemes/literature/cardelli2013_2D_noGC.ts',
              'schemes/literature/cardelli2013_2D_2TGC.ts',
              'schemes/literature/cardelli2013_2D_3I.ts',
              'schemes/literature/cardelli2013_2D_3I_noGC.ts',
              'schemes/literature/chen2013_2D_JF.ts',
              'schemes/literature/srinivas2015.ts',
              'schemes/literature/lakin2016_2D_3I.ts']

install_schemes = canonical
if True:
  install_schemes += literature 

setup(
    name='nuskell',
    version=version,
    description='Nucleic acid strand displacement compiler',
    long_description=readme,
    author='Seung Woo Shin, Qing Dong, Robert Johnson, Stefan Badelt, Erik Winfree',
    author_email='winfree@caltech.edu',
    url='http://www.github.com/DNA-and-Natural-Algorithms-Group/nuskell/',
    data_files=[('nuskell/schemes', install_schemes)],
    license=license,
    install_requires=['pyparsing>=1.5.5', 'argparse>=1.2.1'],
    test_suite='tests',
    packages=['nuskell', 'nuskell.parser', 'nuskell.interpreter', 'nuskell.verifier'],
    scripts=['scripts/nuskell', 'scripts/nuskellCMP']
)

