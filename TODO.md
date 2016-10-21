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
      - i.e. it compiles CRN language into the Peppercorn DSD language
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

## Okt, 3rd -- Meeting Erik & Stefan
  All schemes used for automatic nuskell translation should be *generalized*,
  i.e. applicable to every possible input CRN. That means we have often two 
  schemes:
    - Name2016.ts (original version from the paper)
    - Name2016_gen.ts (generalized version)
    - if these are the *same* then we have a softlink from *.ts -> *_gen.ts

  For example, the Cardelli-Schemes with carbage collection are in
  Cardelli2011.ts, but the _gen schemes additionally implement " -> A" and do
  not use garbage collection. Claim: garbage collection is not as helpful as
  one might think, because toehold acclusion exists anyways from the fuel
  strands.

  Also, the _gen.ts schemes should implement general rules to extrapolate from
  bimolecular to trimolecular reactions.

    - 0->X via f->X where f=fuel (can we call them "instant" reactions?)
        
        .. p = if len(r.products) == 0 then [formal(0)] else r.products;

    - A+B+C->X via A+B<=>i; i+C->X
    - A-> W+X+Y+Z via A->W+i1; i1->X+i2; i2->Y+Z

  Some schemes want reactions to be irreversible, but others save species with
  reversible reactions. That means it is possible to optimize the formal input
  CRN. The translation scheme should specify how the CRN should be compressed.

  Check Schemes:
    - Seelig lab: Programmable chemical controlers made from DNA 
        -- variant of published 2-domain scheme, see what they did.

    - qian_revfixed should be the *correct* implementation of the paper.
      According to the enumerator, the irreversible reaction is always correctly
      implemented, but that is confusing.. the system for A->B and A<=>B has the
      same size. Need to check what is really going on!

    - send case examples for assert errors and for soloveichik enumeration

  We also have third category of implementation schemes for particular
  reactions. They are hacks that enable comparison to VisualDSD style.
    - Name2016_catalyst.ts

