#!/usr/bin/env python
from setuptools import setup

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

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name = 'nuskell',
    version = '0.8',
    description = 'Domain-level strand displacement compiler',
    long_description = LONG_DESCRIPTION,
    long_description_content_type = 'text/markdown',
    author = 'Stefan Badelt, Seung Woo Shin, Robert Johnson, Qing Dong, Erik Winfree',
    author_email = 'winfree@caltech.edu',
    maintainer = 'Stefan Badelt',
    maintainer_email = 'bad-ants-fleet@posteo.eu',
    url = 'http://www.github.com/DNA-and-Natural-Algorithms-Group/nuskell/',
    data_files = [('nuskell/schemes', install_schemes)],
    license = 'MIT',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        ],
    install_requires = [
        'pyparsing', 
        #'networkx>=2.4', 
        #'seaborn', # nuskellCMP
        #'pandas', # nuskellCMP
        #'numpy',  # nuskellCMP
        'peppercornenumerator>=1.0',
        'dsdobjects>=0.8',
        'crnsimulator>=0.6'],
    test_suite = 'tests',
    packages = ['nuskell', 'nuskell.dsdcompiler'],
    entry_points = {
        'console_scripts': [
            'nuskell=nuskell.framework:main',
            'nuskellCMP=nuskell.compare_schemes:main'
            ],
        }
)

