# Nuskell: nucleic acid strand displacement compiler

[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/DNA-and-Natural-Algorithms-Group/nuskell)](https://github.com/DNA-and-Natural-Algorithms-Group/nuskell/tags)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/DNA-and-Natural-Algorithms-Group/nuskell?include_prereleases)](https://github.com/DNA-and-Natural-Algorithms-Group/nuskell/releases)
[![PyPI version](https://badge.fury.io/py/nuskell.svg)](https://badge.fury.io/py/nuskell)
[![PyPI - License](https://img.shields.io/pypi/l/nuskell)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.com/DNA-and-Natural-Algorithms-Group/nuskell.svg?branch=master)](https://travis-ci.com/github/DNA-and-Natural-Algorithms-Group/nuskell)
[![Codecov branch](https://img.shields.io/codecov/c/github/DNA-and-Natural-Algorithms-Group/nuskell/master)](https://codecov.io/gh/DNA-and-Natural-Algorithms-Group/nuskell)

``Nuskell`` is software framework to compile formal chemical reaction networks
(CRNs) into domain-level strand displacement (DSD) systems.  As there are many
ways to do such a translation, we provide a library of [translation schemes]
(i.e. variations of CRN-to-DSD translations) that the user can select from. 

In order to proof/disproof the correctness of a particular translation,
``Nuskell`` includes the domain-level reaction enumeration package
[Peppercorn][] [[Badelt et al. (2020)]] and the CRN verification package
[crnverifier][].  
Peppercorn finds intended and potentially unintended reaction pathways, the
crnverifier then checks if the implementation CRN is a correct implementation
of the formal CRN using the stochastic trajectory-type CRN correctness notions
of bisimulation [[Johnson et al. (2019)]] and pathway decomposition 
[[Shin et al.  (2019)]].

### Examples

Implement a formal CRN using a particular translation-scheme:

```
  $ echo "A + B <=> X + Y; X -> A" | nuskell --ts srinivas2017.ts --verify crn-bisimulation
```
for options see:
```
  $ nuskell --help
```
## Translation Schemes
Detailed information about translation schemes can be found in the [translation
schemes] directory.
 
## Installation
```
  $ python setup.py install
```

## Documentation
A preview of the documentation for release v1.0 can be found at: [documentation].

## Version
0.8 -- basically a complete rewrite, python>=3.8 only.
  * nuskell.dsdcompiler is now a subpackage to compile from CRN to DSD.
  * crnverifier is now an independent package and therefore a dependency.
  * enumeration interface updated to peppercornenumerator-v1.1.
  * nuskell now uses the prototype objects provided by the dsdobjects library.

### Authors
  - Stefan Badelt
  - Seung Woo Shin
  - Robert Johnson
  - Qing Dong
  - Erik Winfree

### License
MIT

## Cite
Stefan Badelt, Seung Woo Shin, Robert F. Johnson, Qing Dong, Chris Thachuk, and Erik Winfree (2017)
"A General-Purpose CRN-to-DSD Compiler with Formal Verification, Optimization, and Simulation Capabilities"
[[Badelt et al. (2017)]].


[//]: References
[Peppercorn]: <https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator>
[crnverifier]: <https://github.com/DNA-and-Natural-Algorithms-Group/crnverifier>
[translation schemes]: <https://github.com/DNA-and-Natural-Algorithms-Group/nuskell/tree/master/nuskell/dsdcompiler/schemes>
[Badelt et al. (2017)]: <https://doi.org/10.1007/978-3-319-66799-7_15>
[Badelt et al. (2020)]: <https://doi.org/10.1098/rsif.2019.0866>
[Shin et al. (2019)]: <https://doi.org/10.1016/j.tcs.2017.10.011> 
[Johnson et al. (2019)]: <https://doi.org/10.1016/j.tcs.2018.01.002>
[documentation]: <http://dna.caltech.edu/~badelt/nuskell/index.html>

