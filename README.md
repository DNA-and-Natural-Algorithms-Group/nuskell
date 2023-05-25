# Nuskell: a nucleic acid strand displacement compiler

[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/DNA-and-Natural-Algorithms-Group/nuskell)](https://github.com/DNA-and-Natural-Algorithms-Group/nuskell/tags)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/DNA-and-Natural-Algorithms-Group/nuskell?include_prereleases)](https://github.com/DNA-and-Natural-Algorithms-Group/nuskell/releases)
![build](https://github.com/DNA-and-Natural-Algorithms-Group/nuskell/actions/workflows/python-package.yml/badge.svg)
[![Codecov](https://img.shields.io/codecov/c/github/dna-and-natural-algorithms-group/nuskell)](https://codecov.io/gh/dna-and-natural-algorithms-group/nuskell)

**Nuskell** is compiler framework to translate formal chemical reaction
networks (CRNs) into domain-level strand displacement (DSD) systems. 
To support the diversity of proposed CRN-to-DSD translation methods 
from literature, as well as potential new CRN-to-DSD translation
methods that are yet to be developed, Nuskell provides a *domain-specific
programming language*. In this language, CRN-to-DSD translations can be 
formulated as algorithms: so-called *translation schemes*.  We provide a
library of selected existing [translation schemes] that users can use 
without understanding the details of the Nuskell programming language. 

A notion of correctness for a particular translation is established on a
case-by-case basis using the rate-independent stochastic-level theories of
pathway decomposition equivalence [[Shin et al.  (2019)]] and/or CRN
bisimulation [[Johnson et al. (2019)]].
In order to verify a translation, Nuskell includes the domain-level reaction
enumeration package [Peppercorn][] [[Badelt et al. (2020)]] and the CRN
verification package [crnverifier][].  

Peppercorn finds intended and potentially unintended reaction pathways, the
crnverifier then checks if the implementation CRN is a correct implementation
of the formal CRN using the stochastic trajectory-type CRN correctness notions
of CRN bisimulation [[Johnson et al. (2019)]] and pathway decomposition 
[[Shin et al.  (2019)]].

Nuskell is a first step to integrate biophysical modeling of nucleic acids with
rigorous abstraction hierarchies of modern compilers to design and characterize
DSD systems. For more details, see [[Badelt et al. (2017)]].

## Installation
Nuskell must be installed directly from this repository, we recommend clining
the repository and then using:
```bash
$ pip install .
```

For debugging, or if you are planning a contribution to the repository, please
install the development version and make sure all tests pass:
``` 
$ pip install .[dev]
$ pytest 
```

## Quickstart: the nuskell executable 
Upon installation, the package provides an executable `nuskell`, which is the
main interface combining CRN-to-DSD translation, DSD enumeration and CRN
verification.  For example, to implement the formal CRN
```
A + B <=> X + Y
X -> A
```
using the translation-scheme from [Srinivas (2015)], and to verify that
the translation is correct with CRN bisimulation, use the command line call:
```
  $ echo "A + B <=> X + Y; X -> A" | nuskell --ts srinivas2015.ts --verify crn-bisimulation
```
New users may also appreciate the `-v` flag to get more detailed information on
the individual compilation steps, as well as the option `--pilfile` to print the 
DNA complexes generated in different stages of the compilation.

For more options see:
```
  $ nuskell --help
```

### Translation Schemes
Detailed information about existing translation schemes can be found in the
[translation schemes] directory.
 
### Documentation
A (preliminary) documentation can be found online: [documentation].  Apart from 
auto-generated API documentation, the documentation mainly provides details on
how to write new translation schemes.  Suggestions and contributions that 
improve the documentation are very welcome.

## Version
0.8 -- basically a complete rewrite, python>=3.8 only.
  * nuskell.dsdcompiler is now a subpackage to compile from CRN to DSD.
  * crnverifier is now an independent package and therefore a dependency.
  * enumeration interface updated to peppercornenumerator-v1.1.
  * nuskell now uses the prototype objects provided by the dsdobjects library.

## Authors
  - Stefan Badelt
  - Seung Woo Shin
  - Hope Amber Johnson
  - Qing Dong
  - Erik Winfree

## Cite
Stefan Badelt, Seung Woo Shin, Robert F. Johnson, Qing Dong, Chris Thachuk, and Erik Winfree (2017)
"A General-Purpose CRN-to-DSD Compiler with Formal Verification, Optimization, and Simulation Capabilities"
[[Badelt et al. (2017)]].

[//]: References
[Peppercorn]: <https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator>
[crnverifier]: <https://github.com/DNA-and-Natural-Algorithms-Group/crnverifier>
[translation schemes]: <https://github.com/DNA-and-Natural-Algorithms-Group/nuskell/tree/master/nuskell/dsdcompiler/schemes>
[Srinivas (2015)]: <https://github.com/DNA-and-Natural-Algorithms-Group/nuskell/tree/master/nuskell/dsdcompiler/schemes/literature/srinivas2015.ts>
[Badelt et al. (2017)]: <https://doi.org/10.1007/978-3-319-66799-7_15>
[Badelt et al. (2020)]: <https://doi.org/10.1098/rsif.2019.0866>
[Shin et al. (2019)]: <https://doi.org/10.1016/j.tcs.2017.10.011> 
[Johnson et al. (2019)]: <https://doi.org/10.1016/j.tcs.2018.01.002>
[documentation]: <http://dna.caltech.edu/~badelt/nuskell/index.html>

