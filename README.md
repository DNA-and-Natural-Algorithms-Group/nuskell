# Nuskell: nucleic acid strand displacement compiler

[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/DNA-and-Natural-Algorithms-Group/nuskell)](https://github.com/DNA-and-Natural-Algorithms-Group/nuskell/tags)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/DNA-and-Natural-Algorithms-Group/nuskell?include_prereleases)](https://github.com/DNA-and-Natural-Algorithms-Group/nuskell/releases)
![build](https://github.com/DNA-and-Natural-Algorithms-Group/nuskell/actions/workflows/python-package.yml/badge.svg)
[![Codecov](https://img.shields.io/codecov/c/github/dna-and-natural-algorithms-group/nuskell)](https://codecov.io/gh/dna-and-natural-algorithms-group/nuskell)

``Nuskell`` is software framework to compile formal chemical reaction networks
(CRNs) into domain-level strand displacement (DSD) systems.  As there are many
ways to do such a translation, we provide a library of [translation schemes]
(i.e. variations of CRN-to-DSD translations) that the user can select from. 

In order to proof the correctness of a particular translation,
``Nuskell`` includes the domain-level reaction enumeration package
[Peppercorn][] [[Badelt et al. (2020)]] and the CRN verification package
[crnverifier][].  
Peppercorn finds intended and potentially unintended reaction pathways, the
crnverifier then checks if the implementation CRN is a correct implementation
of the formal CRN using the stochastic trajectory-type CRN correctness notions
of bisimulation [[Johnson et al. (2019)]] and pathway decomposition 
[[Shin et al.  (2019)]].

## Installation
Nuskell must be installed from the git repository, we recommend using:
```bash
$ pip install .
```

For debugging, or if you plan to contribute to the repository, please install
the development version and make sure all tests pass:
``` 
$ pip install .[dev]
$ pytest 
```

### Quick start
Upon installation, the package provides an executable `nuskell`: an interface
combining CRN-to-DSD compilation, DSD enumeration and CRN verification.  For
example, to implement the formal CRN
```
A + B <=> X + Y
X -> A
```
using the translation-scheme from [Srinivas (2015)], and then verify that
it is correct with CRN bisimulation, use the command line call:
```
  $ echo "A + B <=> X + Y; X -> A" | nuskell --ts srinivas2015.ts --verify crn-bisimulation
```

For more options see:
```
  $ nuskell --help
```
New users may appreciate the `-v` flag to get more detailed information on the individual 
compilation steps.

## Translation Schemes
Detailed information about translation schemes can be found in the [translation
schemes] directory.
 
## Documentation
A preview of the documentation for release v1.0 can be found at:
[documentation].  Note that this documentation is not fully up to date with the
latest code release. Questions and/or suggestions for improving the
documentation are very welcome.

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

