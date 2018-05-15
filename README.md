# Nuskell: Nucleic acid strand displacement compiler

``Nuskell`` compiles formal chemical reaction networks (CRNs) into domain-level
strand displacement (DSD) systems. It provides a library of ``translation
schemes`` (i.e. variations of CRN-to-DSD translations) to exploit the diversity
of DSD circuits implementing the same CRN.

In order to proof/disproof the correctness of a particular translation,
``Nuskell`` includes the domain-level reaction enumerator ``Peppercorn`` [Grun
et al. (2014)] to find intended and unintended reaction pathways and then
provides two notions of stochastic trajectory-type CRN equivalence:
bisimulation [Johnson et al. (2016)] and pathway decomposition [Shin et al. (2014)].

The domain-level reactions and their approximate rates can be exported in form
of an ODE system to simulate the dynamics of the compiled DSD network.

### Examples

Implement a formal CRN using a particular translation-scheme:

```
  $ echo "A + B <=> X + Y; X -> A" | nuskell --ts scheme.ts --verify modular-bisimulation
```
for options see:
```
  $ nuskell --help
```
## Translation Schemes
Detailed information about translation schemes can be found in the ``/schemes`` directory.
 
## Installation
```
  $ python setup.py install
```
or
```
  $ python setup.py install --user
```

## Documentation
A preview of the documentation for release v1.0 can be found at: [documentation].

## Version
0.5.1

### Authors
  - Seung Woo Shin
  - Qing Dong
  - Robert Johnson
  - Stefan Badelt
  - Erik Winfree

### License
MIT

## Cite
Stefan Badelt, Seung Woo Shin, Robert F. Johnson, Qing Dong, Chris Thachuk, and Erik Winfree (2017)
"A General-Purpose CRN-to-DSD Compiler with Formal Verification, Optimization, and Simulation Capabilities"
[[Badelt et al. (2017)]].


[//]: References
[Badelt et al. (2017)]: <https://doi.org/10.1007/978-3-319-66799-7_15>
[Grun et al. (2014)]: <https://arxiv.org/abs/1505.03738>
[Shin et al. (2014)]: <http://dna.caltech.edu/DNAresearch_publications.html#PathwayDecomposition>
[Johnson et al. (2016)]: <http://dna.caltech.edu/DNAresearch_publications.html#CRN-Bisimulation>
[documentation]: <http://dna.caltech.edu/~badelt/nuskell/index.html>

