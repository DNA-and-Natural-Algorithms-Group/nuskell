# Nuskell -- Checklist for Release

## DNA22:
  - include Robert's bisimulation '
  - add basic shinx documentation 

## Publication/Release 1.0:
  - add crn-simulation
  - add nuskell.optimization module
  - add schemes/ directory with well documented examples
  - case-by-case analysis for non formal-CRN-type schemes

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

### Clarify:
  - history domains
  - naming conventions: 
    * Nuskell compiler, nuskell executable/library
    * CRN is programming language
    * translation schemes are instructrions written in ...?
    * CRN-to-DSD, DOM, DNA, TSD?

  - commandline calls
  - Authors in setup and licence, choice of licence, ...

