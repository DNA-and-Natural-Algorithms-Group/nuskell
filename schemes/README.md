# Nuskell CRN-to-DSD translation schemes

## Explore
A collection of translation schemes. We distinguish three types of schemes:

  - **literature**: Schemes implemented *exactly* as described in literature.
    If schemes are only described for a subset of reactions (e.g. only
    bimolecular reactions), then all other reactions are implemented as
    combinations of these reactions, with fuel species (f) or intermediate
    species (i) that look like signal species, for example:
    ```
    - { -> X} is implemented as {f -> X} 
    - {X -> } is implemented as {X -> f} 
    - {X -> Y} is implemented as {X + f -> f + Y}

    - {W + X + Y -> Z} is implemented as {W + X <=> i; i + Y -> Z} 
    - {Z -> W + X + Y} is implemented as  {Z -> i + W; i -> X + Y}
    ```
    Note: there are exceptions where schemes are not generalized on the CRN
    level, but they are generalized on the DNA level. Either, because it was
    pointed out in the publication, or, because the generalization was obvious. 

  - **canonical**: recommended schemes for CRN-to-DSD compilation. These can either
    be identical copies of schemes found in ''literature'', (i.e. schemes with
    the same name), or they are *corrected* or *optimized variants* of those
    schemes. These schemes are globally installed with the Nuskell library.

  - **experimental**: new schemes, or variations of schemes that haven't been
    thoroughly tested.

  - **implementations**: Alternative input format (*.PIL) for case-by-case testing
    of pre-compiled systems as described in literature. 

## Contribute
Most importantly, reference the corresponding literature and yourself as author 
of the translation scheme. Make sure to emphasize potential optimizations on
the formal CRN input. For example, answer the following questions: 

  - *Does the order of reactants/products matter?* A scheme might implement the
    catalytic reactions more efficiently if they are specified as ```x + c <=>
    c + y``` than if they are specified as ```c + x <=> c + y```. 

  - *Is the scheme optimized?* A translation scheme might implement the
    *delayed-choice* optimization: ```{A + B -> C; A + B -> D + E}``` where
    ```A + B``` consumption is implemented only once:
    ```{A + B -> i; i -> C; i -> D + E}``` 

  - *Should reversible reactions be specified as reversible reactions?* Some
    schemes require reactions to be irreversible, others are more efficient for
    reversible reactions. A scheme can make a choice to combine irreversible
    reactions into reversible ones, but it doesn't have to!

## Overview

### ./literature
  * `soloveichik2010.ts`, *DNA as a universal substrate for chemical kinetics*. [Soloveichik et al. (2010)]
    ```
    #
    # Soloveichik, Seelig, Winfree: "DNA as a universal substrate for chemical
    # kinetics", Proceedings of the National Academy of Sciences, 107: 5393-5398,
    # 2010.
    #
    # Note:   * implements Figure 2 (X1 -> X2 + X3) 
    #         * implements Figure 3 (X1 + X2 -> X3) 
    #         * implements (X1 + X2 -> X3 + X4) as combination of the above.
    #         * DNA level generalization for higher order reactions.
    #         * CRN level generalization for {->X; X->}
    #
    # Coded by Seung Woo Shin (seungwoo.theory@gmail.com)
    #          Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `cardelli2011_FJ.ts` *Strand Algebras for DNA Computing*. [Cardelli (2011)]

    ```
    #
    # Luca Cardelli "Strand Algebras for DNA Computing", Natural Computing, 10:
    # 407-428, 2011.
    #
    #   Note: * This implements the `fork' and `join' gates from the paper,
    #           - Figure 4 (Annihilator): {X -> }
    #           - Figure 5 (Transducer):  {X -> Y}
    #           - Figure 6 (2-way Fork):  {X -> Y + Z}
    #           - Figure 7 (2-way Join):  {X + Y -> Z}
    #           - Figure 8 (2-way Join cleanup): garbage collection gates.
    #
    #         * Generalized on the CRN level for { -> X} and trimolecular or higher
    #           order reactions.
    #
    # Coded by Seung Woo Shin (seungwoo.theory@gmail.com) 
    #          Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `cardelli2011_FJ_noGC.ts` *Strand Algebras for DNA Computing*. [Cardelli (2011)]

    ```
    #
    # Luca Cardelli "Strand Algebras for DNA Computing", Natural Computing, 10:
    # 407-428, 2011.
    #
    #   Note: * This implements the `fork' and `join' gates from the paper,
    #           - Figure 4 (Annihilator): {X -> }
    #           - Figure 5 (Transducer):  {X -> Y}
    #           - Figure 6 (2-way Fork):  {X -> Y + Z}
    #           - Figure 7 (2-way Join):  {X + Y -> Z}
    #
    #         * Does not implement Figure 8: garbage collection.
    #
    #         * Generalized on the CRN level for { -> X} and trimolecular or higher
    #           order reactions.
    #
    # Coded by Seung Woo Shin (seungwoo.theory@gmail.com) 
    #          Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `cardelli2011_NM.ts` *Strand Algebras for DNA Computing*. [Cardelli (2011)]

    ```
    #
    # Luca Cardelli "Strand Algebras for DNA Computing", Natural Computing, 10:
    # 407-428, 2011.
    #
    #   Note: * Implements the scheme shown in Figures 9 and 10. Note that the 
    #           gate species at the bottom seems like one molecule, but it is 
    #           in fact two molecules because of an arrowhead separating them.
    #           - Figure 9 (n x m gate): {X1 + ... + Xn -> Y1 + ... + Ym}
    #           - Figure 10 (1 x 1 gate): {X1 > Y1}
    #
    #         * Generalized on the CRN level for { -> X}
    #
    # Coded by Seung Woo Shin (seungwoo.theory@gmail.com) 
    #          Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `cardelli2011_NM_noGC.ts` *Strand Algebras for DNA Computing*. [Cardelli (2011)]
    ```
    #
    # Luca Cardelli "Strand Algebras for DNA Computing", Natural Computing, 10:
    # 407-428, 2011.
    #
    #   Note: * Implements the scheme shown in Figures 9 and 10. Note that the 
    #           gate species at the bottom seems like one molecule, but it is 
    #           in fact two molecules because of an arrowhead separating them.
    #           - Figure 9 (n x m gate): {X1 + ... + Xn -> Y1 + ... + Ym}
    #             *modified* to exclude garbage collection.
    #           - Figure 10 (1 x 1 gate): {X1 > Y1}
    #
    #         * Generalized on the CRN level for { -> X}
    #
    # Coded by Seung Woo Shin (seungwoo.theory@gmail.com) 
    #          Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `qian2011_3D.ts`, *Efficient Turing-universal computation with DNA polymers*. [Qian et al. (2011)]
    ```
    #
    # Qian, Soloveichik, Winfree: "Efficient Turing-universal computation with DNA
    # polymers", DNA Computing and Molecular Programming 16, 2011.
    #
    # Note:   * implements Figure 1 (X + Y -> A + B) 
    #         * implements Figure 2 (X + Y <=> A + B)
    #
    #         * generalized on the DNA level
    #         * automatically combines corresponding irreversible reactions into 
    #           one reversible reaction
    #
    #         * pathway and bisimulation incorrect for irreversible reactions.
    #
    # Coded by Seung Woo Shin (seungwoo.theory@gmail.com)
    #          Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `lakin2012_3D.ts`, *Abstractions for DNA circuit design*. [Lakin et al. (2012)]

    ```
    #
    # Lakin, Youssef, Cardelli, Phillips "Abstractions for DNA circuit design."
    # J. R. Soc. Interface 9 (68) (2012) 470-486 (2012)
    #
    # Note: * Implements Figure 5: {A + B -> C + D}
    #
    #       * no fuel species to reverse the release of second reactant
    #
    #       * generalized on the CRN level
    #
    # Coded by Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `cardelli2013_2D.ts` *Two-Domain DNA Strand Displacement*. [Cardelli (2013)]
    ```
    #
    # Luca Cardelli "Two-Domain DNA Strand Displacement", Mathematical Structures
    # in Computer Science. (2013)
    #
    #   Note: * The scheme implements Figures 4-10 (transducer, fork, 
    #           catalyst, join, 3-way join and the implicit n-way join)
    #           - Figure 3-5 (Transducer): {X -> Y}
    #           - Figure 6 (2-way Fork):  {X -> Y + Z}
    #           - Figure 7 (Catalyst):  {X + Y -> Y + Z}
    #           - Figure 8,9 (2-way Join): {X + Y -> Z}
    #           - Figure 10 (3-way Join): {W + X + Y -> Z}
    #
    #         * Generalized on the CRN level for { -> X}
    #
    #         * Generalized on DNA level for higher order reactions, but supports
    #           only one catalyst optimization. e.g. X + Y + Z -> Z + Y + W will 
    #           only optimize for Z, but not for Y.
    #         * Note, the Catalyst case does not apply to {C -> C + A}, as it is
    #           not obvious from the paper.
    #
    # Coded by Seung Woo Shin (seungwoo.theory@gmail.com).
    #          Stefan Badelt (badelt@caltech.edu)
    #

    ```

  * `cardelli2013_2D_noGC.ts` *Two-Domain DNA Strand Displacement*. [Cardelli (2013)]
    ```
    #
    # Luca Cardelli "Two-Domain DNA Strand Displacement", Mathematical Structures
    # in Computer Science. (2013)
    #
    #   Note: * The scheme implements Figures 4-10 (transducer, fork, 
    #           catalyst, join, 3-way join and the implicit n-way join)
    #           - Figure 3-5 (Transducer): {X -> Y}
    #           - Figure 6 (2-way Fork):  {X -> Y + Z}
    #           - Figure 7 (Catalyst):  {X + Y -> Y + Z}
    #           - Figure 8,9 (2-way Join): {X + Y -> Z}
    #           - Figure 10 (3-way Join): {W + X + Y -> Z}
    #
    #         * All reactios are implemented *without* garbage collection. This 
    #           includes the cooperative binding complexes and additional domains
    #           on the produce complexes.
    #
    #         * Intermediate ('garbage') species on react complexes are fuels 
    #           in this implementation.
    #
    #         * Generalized on the CRN level for { -> X}
    #
    #         * Generalized on DNA level for higher order reactions.
    #
    # Coded by Seung Woo Shin (seungwoo.theory@gmail.com).
    #          Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `cardelli2013_2D_2TGC.ts` *Two-Domain DNA Strand Displacement*. [Cardelli (2013)]

    ```
    #
    # Luca Cardelli "Two-Domain DNA Strand Displacement", Mathematical Structures
    # in Computer Science. (2013)
    #
    #   Note: * The scheme implements Figures 4-6,8-10 (transducer, fork, 
    #           join, 3-way join and the implicit n-way join) but it uses 
    #           garbage collection as described in Figures 11, 12 and 14, 
    #           i.e. avoids cooperative binding by using a second toehold.
    #           - Figure 3-5 (Transducer): {X -> Y}
    #           - Figure 6 (2-way Fork):  {X -> Y + Z}
    #           - Figure 11,12 (2-toehold-2-way Join): {X + Y -> Z}
    #           - Figure 14 (3-input join collector)
    #           - Generalized n-way join collector as described in the text
    #
    #         * Does not implement Figure 7 (Catalyst):  {X + Y -> Y + Z}
    #
    #         * Generalized on the CRN level for { -> X}
    #
    #         * Generalized on the DNA level for higher order reactions.
    #
    # Coded by Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `cardelli2013_2D_3I.ts` *Two-Domain DNA Strand Displacement*. [Cardelli (2013)]
    ```
    #
    # Luca Cardelli "Two-Domain DNA Strand Displacement", Mathematical Structures
    # in Computer Science. (2013)
    #
    #   Note: * The scheme implements Figures 4-10 (transducer, fork, 
    #           catalyst, join, 3-way join and the implicit n-way join)
    #           but it uses a 3-domain strand to introduce an irreversilbe
    #           step after all reactants have been consumed (Figure 15).
    #           - Figure 3-5 (Transducer): {X -> Y}
    #           - Figure 6 (2-way Fork):  {X -> Y + Z}
    #           - Figure 7 (Catalyst):  {X + Y -> Y + Z}
    #           - Figure 8,9 (2-way Join): {X + Y -> Z}
    #           - Figure 10 (3-way Join): {W + X + Y -> Z}
    #           - Figure 15 (3-domain transducer modification)
    #             + same 3-domain mechanism for every other reation.
    #
    #         * Generalized on the CRN level for { -> X}
    #
    #         * Uses garbage collection with cooperative binding, which cannot be 
    #           enumerated by peppercorn -> pathway equivalence fails.
    #
    #         * Generalized on DNA level for higher order reactions, but supports
    #           only one catalyst optimization. e.g. X + Y + Z -> Z + Y + W will 
    #           only optimize for Z, but not for Y.
    #
    #         * Note, the Catalyst case does not apply to {C -> C + A}, as it is
    #           not obvious from the paper.
    #
    #         * Needs to be enumerated using --reject-remote
    #
    # Coded by Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `cardelli2013_2D_3I_noGC.ts` *Two-Domain DNA Strand Displacement*. [Cardelli (2013)]

    ```
    #
    # Luca Cardelli "Two-Domain DNA Strand Displacement", Mathematical Structures
    # in Computer Science. (2013)
    #
    #   Note: * The scheme implements Figures 4-10 (transducer, fork, 
    #           catalyst, join, 3-way join and the implicit n-way join)
    #           but it uses a 3-domain strand to introduce an irreversilbe
    #           step after all reactants have been consumed (Figure 15).
    #           - Figure 3-5 (Transducer): {X -> Y}
    #           - Figure 6 (2-way Fork):  {X -> Y + Z}
    #           - Figure 7 (Catalyst):  {X + Y -> Y + Z}
    #           - Figure 8,9 (2-way Join): {X + Y -> Z}
    #           - Figure 10 (3-way Join): {W + X + Y -> Z}
    #           - Figure 15 (3-domain transducer modification)
    #             + same 3-domain mechanism for every other reation.
    #
    #         * All reactions are implemented *without* garbage collection. This 
    #           includes the cooperative binding complexes and additional domains
    #           on the produce complexes.
    #
    #         * Generalized on the CRN level for { -> X}
    #
    #         * Generalized on DNA level for higher order reactions, but supports
    #           only one catalyst optimization. e.g. X + Y + Z -> Z + Y + W will 
    #           only optimize for Z, but not for Y.
    #         * Note, the Catalyst case does not apply to {C -> C + A}, as it is
    #           not obvious from the paper.
    #
    #         * Needs to be enumerated using --reject-remote
    #
    # Coded by Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `chen2013_2D_FJ.ts` *Programmable chemical controllers made from DNA*. [Chen et al. (2013)]
    ```
    #
    # Chen, Dalchau, Srinivas, Phillips, Cardelli, Soloveichik, Seelig
    # "Programmable chemical controllers made from DNA", Nature Nanotech., 2013.  
    #
    #   Note: * Implements supplementary Figure S7.6 fork and join modules.
    #           {A->R; R->B; A+B->R; R->B+C; A+B+C->R; R->B+C+D}
    #
    #         * generalized on the DNA level higher order reactions and for {X->}
    #
    #         * generalized on the CRN level for {->X}
    #
    #         * includes the delayed choice optimization suggested by the authors
    #
    # Coded by Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `srinivas2015.ts` *Programming chemical kinetics: engineering dynamic
    reaction networks with DNA strand displacement*. [Srinivas (2015)]

    ```
    #
    # Niranjan Srinivas PhD Thesis: "Programming chemical kinetics: engineering
    # dynamic reaction networks with DNA strand displacement.", Caltech (2015)
    #
    # Note: * Implements Figures from Chapter 3: Adventures in programming and
    #         debugging molecular reaction networks:
    #         - Figure 3.4: {B + A -> X + Y}
    #         - Figure 3.5: {B -> X}
    #         - Figure 3.6: {B -> }
    #         - Figure 3.7: { -> X}
    #       * Generalized on the DNA level for higher order reactions.
    #
    # Coded by Stefan Badelt (badelt@caltech.edu)
    #
    ```

  * `lakin2016_2D_3I.ts`, *Modular verification of chemical reaction network
    encodings via serializability analysis*. [Lakin et al. (2016)]
    ```
    #
    # Lakin, Stefanovic, Phillips "Modular verification of chemical reaction
    # network encodings via serializability analysis" (2016), TCS
    #
    #   Note: * The scheme implements Figure 8 from the above publication. 
    #           {x + y -> y + b}
    #
    #         * generalized on the CRN level for all but bimolecular reactions.
    #
    #         * requires autocatalytic format, otherwise incorrect, BUT:
    #           as e.g. {A->B} is implemented as {A+f->f+B} it is correct for 
    #           a wider class of CRNs!
    #
    #         * enumerate using --reject-remote !!!
    #
    # Coded by Stefan Badelt (badelt@caltech.edu)
    #
    ```
  
### ./canonical
  * `soloveichik2010_v1.ts`, *DNA as a universal substrate for chemical kinetics*. 
    [Soloveichik et al. (2010)]

  * `qian2011_3D_v1.ts`, *Efficient Turing-universal computation with DNA polymers*. 
    [Qian et al. (2011)]

  * `srinivas2015.ts` *Programming chemical kinetics: engineering dynamic
    reaction networks with DNA strand displacement*. [Srinivas (2015)]

### ./implementations
  * `song2016-add.ts` Figure 3-5 of [Song et al. (2016)]

  * `song2016-subtract.ts` Figure 6-8 of [Song et al. (2016)]

### Last Update
April, 27th, 2017

[//]: References
[Turberfield et al. (2003)]: <http://dx.doi.org/10.1103/PhysRevLett.90.118102>
[Zhang et al. (2007)]: <http://dx.doi.org/10.1126/science.1148532>
[Soloveichik et al. (2010)]: <http://dx.doi.org/10.1073/pnas.0909380107>
[Qian et al. (2011)]: <http://dx.doi.org/10.1007/978-3-642-18305-8_12>
[Qian & Winfree (2011)]: <http://dx.doi.org/10.1126/science.1200520>
[Lakin et al. (2012)]: <http://dx.doi.org/10.1098/rsif.2011.0343>
[Lakin et al. (2016)]: <http://dx.doi.org/10.1016/j.tcs.2015.06.033>
[Cardelli (2011)]: <http://dx.doi.org/10.1007/s11047-010-9236-7>
[Cardelli (2013)]: <http://dx.doi.org/10.1017/S0960129512000102>
[Chen et al. (2013)]: <http://dx.doi.org/10.1038/NNANO.2013.189>
[Thachuk et al. (2015)]: <http://dx.doi.org/10.1007/978-3-319-21999-8_9>
[Srinivas (2015)]: <http://www.dna.caltech.edu/Papers/NiranjanSrinivas_2015_thesis.pdf>
[Song et al. (2016)]: <http://dx.doi.org/10.1021/acssynbio.6b00144>

