# Nuskell - Checklist for Release

# Issues:
  - write generalized version for lakin scheme (roessler)
  - Check whether there is a reason to keep qian2011_gen.ts separate from qian2011_fixed.ts
  - structure the scheme directory into {original, generalized, optimized, implementation}
  - original only makes sense for schemes that do not verify, otherwise they are the same
    as generalized, i.e. genralized by using A+B->X+Y using fuel-strands. ... maybe we can
    make a directory 'verfied' for schemes that are known to be correct on every input?

  - is there a problem when the enumerator gets history-reduced systems in the first place?
  - identify rates for leaky remote-toehold binding in the soloveichik scheme (A+A->A)
  - reproduce results from QD bisimulation on soloveichik and roessler.

## Release 0.2:

### Coding:
  - adapt pathway-decomposition interface
  - adapt integrated-hybrid interface
  - add basic sphinx documentation 
  - add PIL I/O
  - add kernel I/O
  - add domain I/O
  - add VisualDSD I/O

### Testing :
  - add schemes/ directory with well documented examples
  - unittesting for translation
  - unittesting for pathway-decomposition

## Publication/Release 1.0:

### Coding:
  - crn-simulation [TestTube(object)]
  - optimization module
  - case-by-case analysis for implementation CRNs

### Documentation:
  - history domains
  - commandline calls
  - naming conventions: 
    * Nuskell CRN-to-DSD verifying compiler
      - i.e. it compiles CRN language into the DSD language
    * programming languages:
      - CRN is a programming language
      - Domain-level strand displacement (DSD) programming language
      - Nuskell translation scheme language

### Publication content:
  - Describe the .ts language
    * philosopy 
    * difference to VisualDSD
    * functional design choice (reading, debugging)
    * case-by-case examples from literature (Zhang, Song, ...)
       - analog formula
       - turing machine 
       - finite state machine
       - ...

  - performance of `Nuskell`:
    - translation
    - verification
    - enumeration

  - performance of schemes:
    - collect interesting formal CRNs and compare different implementation CRNs
      * size of enumerated crn and parameters for enumeration
      * verify pathway & bisimulation
      * simulation time of enumerated crn
      * nucleotides * concentration (=cost)
      * \# of strands
      * \# of complexes
      * max(strand/complex)
      * \# of domain-IDs 
      * free energy gain, e.g. in terms of toehold bindings
      * leak reduction (maybe)
      * toehold context, length of domains ...
    - improve published schemes
    - discuss garbage collection problems
    - why do verifications disagree
    - preformance on random CRNs (spurious pathways)

  - discuss optimizations on mutiple levels:
    - discuss general criteria for optimization
    - opti at the transation scheme
    - opti at the compiler level
      - remove redundant species (CRN level, including proof)
      - try different schemes
    - write an objective function for DSD design and evaluate multiple schemes

### Figures for publication:
  - schematic workflow using `Nuskell`
  - comparison of translation schemes
  - verification using `Nuskell`
  - optimization using `Nuskell`

