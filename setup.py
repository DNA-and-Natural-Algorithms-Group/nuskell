#!/usr/bin/env python
from setuptools import setup

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
    license = 'MIT',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        ],
    python_requires = '>=3.8',
    install_requires = [
        'natsort', 
        'pyparsing', 
        'dsdobjects>=0.8',
        'peppercornenumerator>=1.1',
        'crnverifier>=0.2'],
    packages = ['nuskell', 'nuskell.dsdcompiler'],
    include_package_data = True,
    package_data = {'nuskell.dsdcompiler': ['schemes/literature/*.ts', 
                                            'schemes/variants/*.ts']},
    test_suite = 'tests',
    entry_points = {
        'console_scripts': [
            'nuskell=nuskell.framework:main',
            'nuskellCMP=nuskell.compare_schemes:main'
            ],
        }
)

