# Nuskell - Checklist for Release

## Release 0.2.2 (for DNA22):
  - include Robert's bisimulation '
  - add basic sphinx documentation 
  - unittesting for bisimulation

### Clarify:
  - history domains
  - commandline calls
  - naming conventions: 
    * Nuskell CRN-to-DSD verifying compiler
      - i.e. it compiles CRNs into the Peppercorn DSD language
    * programming languages:
      - CRN is a programming language
      - Domain-level strand displacement (DSD) programming language
      - Nuskell translation scheme language

## Publication/Release 1.0:
  - add crn-simulation
  - add nuskell.optimization module
  - add schemes/ directory with well documented examples
  - case-by-case analysis for non formal-CRN-type schemes
  - unittesting for translation (as much as possible, but no more)
  - unittesting for pathway-decomposition

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
      * verify pathway & bisimulation 
      * nucleotides * concentration (=cost)
      * \# of strands
      * \# of domain-IDs 
      * maximum count
      * free energy gain, e.g. in terms of toehold bindings
      * leak reduction (maybe)
      * size of enumerated crn
      * simulate and see how fast the computation is
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

### Figures for publication:
  - schematic workflow using `Nuskell`
  - comparison of translation schemes
  - verification using `Nuskell`
  - optimization using `Nuskell`

