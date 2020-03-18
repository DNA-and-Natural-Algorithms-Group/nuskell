#!/usr/bin/env python

# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

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


LONG_DESCRIPTION="""
``Nuskell`` compiles formal chemical reaction networks (CRNs) into domain-level
strand displacement (DSD) systems. It provides a library of ``translation
schemes`` (i.e. variations of CRN-to-DSD translations) to explore the diversity
of DSD systems implementing the same formal CRN.

In order to proof/disproof the correctness of a particular translation,
``Nuskell`` includes the domain-level reaction enumerator ``Peppercorn`` [Badelt
et al. (2020)] to find intended and unintended reaction pathways and then
provides two notions of stochastic trajectory-type CRN equivalence:
bisimulation [Johnson et al. (2019)] and pathway decomposition [Shin et al. (2019)].
"""

setup(
    name='nuskell',
    version='0.6',
    description='Nucleic acid strand displacement compiler',
    long_description=LONG_DESCRIPTION,
    author='Stefan Badelt, Seung Woo Shin, Robert Johnson, Qing Dong, Erik Winfree',
    author_email='winfree@caltech.edu',
    url='http://www.github.com/DNA-and-Natural-Algorithms-Group/nuskell/',
    data_files=[('nuskell/schemes', install_schemes)],
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        ],
    install_requires=[
        'future', 
        'pyparsing', 
        'networkx>=2.2',
        'seaborn', # nuskellCMP
        'pandas', # nuskellCMP
        'numpy',  # nuskellCMP
        'peppercornenumerator>=0.8',
        'dsdobjects>=0.7.1',
        'crnsimulator>=0.6'],
    dependency_links=[
        'https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator/archive/v0.7.1.tar.gz#egg=peppercornenumerator-0.7.1'],
    test_suite='tests',
    packages=['nuskell', 'nuskell.parser', 'nuskell.interpreter', 'nuskell.verifier'],
    scripts=['scripts/nuskell', 'scripts/nuskellCMP'],
)

