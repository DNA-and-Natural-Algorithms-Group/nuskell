#!/usr/bin/env python

# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
  readme = f.read()

with open('LICENSE') as f:
  license = f.read()

# Literature: These schemes are implemented as described in a publication
literature = ['schemes/literature/soloveichik2010.ts',
              'schemes/literature/qian2011_3D.ts',
              'schemes/literature/cardelli2011_FJ.ts',
              'schemes/literature/cardelli2011_NM.ts',
              'schemes/literature/lakin2012_3D.ts',
              'schemes/literature/cardelli2013_2D.ts',
              'schemes/literature/cardelli2013_2D_2TGC.ts',
              'schemes/literature/cardelli2013_2D_3I.ts',
              'schemes/literature/chen2013_2D_JF.ts',
              'schemes/literature/srinivas2015.ts',
              'schemes/literature/lakin2016_2D_3I.ts']
              #'schemes/literature/srinivas2017.ts',

# Variants: These schemes are variants of schemes in literature.
variants = [  'schemes/variants/qian2011_3D_var1.ts', # SWS version
              'schemes/variants/lakin2012_3D_var1.ts',
              'schemes/variants/cardelli2011_FJ_noGC.ts',
              'schemes/variants/cardelli2011_NM_noGC.ts',
              'schemes/variants/cardelli2013_2D_noGC.ts',
              'schemes/variants/cardelli2013_2D_3I_noGC.ts',
              'schemes/variants/chen2013_2D_JF_var1.ts',
              'schemes/variants/chen2013_2D_JF_var2.ts']
              #'schemes/variants/soloveichik2010_var2.ts', # optimized version
              #'schemes/variants/qian2011_3D_var2.ts', # as shown in the revised paper
              #'schemes/variants/lakin2016_2D_3I_var1.ts', # ML email version

# Optionally, select a list of canonical schemes that are installed with
# nuskell. These schemes will be used by nuskellCMP to compare translations.
canonical = [ 
      'schemes/literature/soloveichik2010.ts',
      'schemes/variants/qian2011_3D_var1.ts', 
      'schemes/literature/cardelli2011_NM.ts',
      'schemes/variants/cardelli2011_NM_noGC.ts',
      'schemes/variants/lakin2012_3D_var1.ts',
      'schemes/variants/cardelli2013_2D_3I_noGC.ts',
      'schemes/variants/chen2013_2D_JF_var1.ts',
      'schemes/literature/srinivas2015.ts']

if False:
  install_schemes = canonical
else :
  install_schemes = literature + variants

setup(
    name='nuskell',
    version='0.5',
    description='Nucleic acid strand displacement compiler',
    long_description=readme,
    author='Seung Woo Shin, Qing Dong, Robert Johnson, Stefan Badelt, Erik Winfree',
    author_email='winfree@caltech.edu',
    url='http://www.github.com/DNA-and-Natural-Algorithms-Group/nuskell/',
    data_files=[('nuskell/schemes', install_schemes)],
    license=license,
    install_requires=[
        'pyparsing>=1.5.5', 
        'networkx>=1.10',
        'seaborn>=0.7.1', # nuskellCMP
        'pandas>=0.19.1', # nuskellCMP
        'numpy>=1.11.0',  # nuskellCMP
        'sympy>=0.7.6.1', 
        'peppercornenumerator==0.5',
        'dsdobjects==0.5',
        'crnsimulator==0.2'],
    dependency_links=[
        'https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator/tarball/development#egg=peppercornenumerator-0.5',
        'https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects/tarball/master#egg=dsdobjects-0.5',
        'https://github.com/bad-ants-fleet/crnsimulator/tarball/master#egg=crnsimulator-0.2'],
    test_suite='tests',
    packages=['nuskell', 'nuskell.parser', 'nuskell.interpreter', 'nuskell.verifier'],
    scripts=['scripts/nuskell', 'scripts/nuskellCMP'],
)

