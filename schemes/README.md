# Nuskell CRN-to-DSD translation schemes

A collection of translation schemes taken from literature.  The schemes are
named after the first author on the publication, some of them have been
extended to capture a bigger variety of CRNs

## Example:

Assume you want to compile the CRN {A+B <=> C+D, D+X -> Y} using the
translation scheme from [Soloveichik et al. (2010)], then the respective
commandline call would be:

```
nuskell --ts soloveichik2010.ts "A + B <=> C + D; D + X -> Y"
```

## Schemes
### published schemes:
  * `soloveichik2010.ts`, *DNA as a universal substrate for chemical kinetics*. [Soloveichik et al. (2010)]

    A generalized scheme for transating arbitrary CRNs into DNA strand displacement systems. 

    **TODO**: Check if this needs to be downgraded.

    Coded by Seung Woo Shin (seungwoo.theory@gmail.com).

  * `cardelli2011_FJ.ts` *Strand Algebras for DNA Computing*. [Cardelli (2011)]

    The 2-way fork and join gates described in Figures 4-7. 

    **TODO**: Garbage collection (Figure 8) is currently not implemented.
    
    Coded by Seung Woo Shin (seungwoo.theory@gmail.com), 
      modified by Stefan Badelt (badelt@caltech.edu)

  * `cardelli2011_NM.ts` *Strand Algebras for DNA Computing*. [Cardelli (2011)]

    A generalization from 2-way join and fork gates to a gate with n inputs and
    m outputs. The implementation includes garbage collection as described in the paper.

    Coded by Seung Woo Shin (seungwoo.theory@gmail.com).

  * `qian2011.ts`, *Efficient Turing-universal computation with DNA polymers*. [Qian et al. (2011)]
  
    A generalized scheme for transating arbitrary CRNs into DNA strand displacement systems. 

    Coded by Seung Woo Shin (seungwoo.theory@gmail.com).

  * `srinivas2015_gen.ts` *Programming chemical kinetics: engineering dynamic reaction networks with DNA strand displacement*. [Srinivas (2015)]

    Niranjan Srinivas' adaptation of the `soloveichik2010.ts` scheme to design the "displacillator".

    Coded by Stefan Badelt (badelt@caltech.edu).

### generalized schemes:
  * `soloveichik2010_gen.ts`, *DNA as a universal substrate for chemical kinetics*. [Soloveichik et al. (2010)]

    A generalized scheme for transating arbitrary CRNs into DNA strand displacement systems. 
    `soloveichik2010.ts` and `soloveichik2010_gen.ts` are identical.

    Coded by Seung Woo Shin (seungwoo.theory@gmail.com).

  * `cardelli2011_FJ_gen.ts` *Strand Algebras for DNA Computing*. [Cardelli (2011)]

    The generalized version of `cardelli2011_FJ.ts` implements instant reactions, e.g. {-> A + B} that have not been discussed in the paper.
    
    Coded by Seung Woo Shin (seungwoo.theory@gmail.com), 
      modified by Stefan Badelt (badelt@caltech.edu)

  * `cardelli2011_NM_gen.ts` *Strand Algebras for DNA Computing*. [Cardelli (2011)]

    A generalization from 2-way join and fork gates to a gate with n inputs and
    m outputs. The scheme implements instant reactions, e.g. {-> A + B} that
    have not been discussed in the paper. Garbage collection is not implemented
    for this scheme.

    Coded by Seung Woo Shin (seungwoo.theory@gmail.com).


  * `qian2011_gen.ts`, *Efficient Turing-universal computation with DNA polymers*. [Qian et al. (2011)]
  
    A generalized scheme for transating arbitrary CRNs into DNA strand displacement systems. 
    `qian2011.ts` and `qian2011_gen.ts` are identical.

    **NOTE:** This scheme implements the reactions from the above paper, and it
    does not verify correct for irreversible reactions. There are also *fixed*
    versions of this scheme, which are unpublished. 

    Coded by Seung Woo Shin (seungwoo.theory@gmail.com).

  * `srinivas2015_gen.ts` *Programming chemical kinetics: engineering dynamic reaction networks with DNA strand displacement*. [Srinivas (2015)]

    Niranjan Srinivas' adaptation of the `soloveichik2010.ts` scheme to design the "displacillator".
    `srinivas2015.ts` and `srinivas2015_gen.ts` are identical.

    Coded by Stefan Badelt (badelt@caltech.edu).

### work in progress:

  * `soloveichik2010_opt.ts`
    
    A variant of soloveichik2010.ts that reduces the number of strands for certain CRNs.

    Coded by Seung Woo Shin (seungwoo.theory@gmail.com).

  * `lakin2011.ts`, *Abstractions for DNA circuit design*. [Lakin et al. (2011)]

    The translation scheme for unbuffered gates, e.g. Figure 5.

    Coded by: Stefan Badelt (badelt@caltech.edu).

  * `cardelli2013_2domain.ts` *Two-Domain DNA Strand Displacement*. [Cardelli (2013)]

    **NOTE:** There is a version from SWS with the suffix '_fixed' as well, but
    what is the difference?

    Coded by Seung Woo Shin (seungwoo.theory@gmail.com).

### incomplete schemes:
  * `thachuk2015.ts` *Leakless DNA strand displacement systems* [Thachuk et al. (2015)]

    The section 6, figure 8 "cooperative hybridization" CRN implementation.
    Clamps are not implemented.
    
    In order to enumerate cooperative hybridization without considering all toehold release steps to be slow,
    we use two toehold sizes, and just the "right-side" toeholds are "slow".  Our solution involves enumeration
    using the options --release-cutoff-1-1 4 --release-cutoff-1-n 8 --k-fast 0.1
    and thus the verifying compiler invocation will be

    Coded by Erik Winfree (winfree@caltech.edu).

  * `metha2015_3domain.ts` *unpublished*

    An unpublished "three domain" scheme by Jenish Mehta, 2015.  Each signal
    species is a single long domain flanked by a toehold on each side.

    Coded by Erik Winfree (winfree@caltech.edu).

    **Note:** There appears to be an undesirably 4-way branch migration pathway
    (not yet investigated).  Also, the waste from the irreversible pathway is
    not inert, which causes problems for pathway decomposition, and therefore
    only bisimulation deems it to be correct.

### missing schemes:
  * `qian2011_seesaw.ts` *Scaling Up Digital Circuit Computation with DNA Strand Displacement Cascades*. [Qian & Winfree (2011)]
  * `lakin2011_buffered.ts` *Abstractions for DNA circuit design*. [Lakin et al. (2011)]
  * `chen2013.ts` *Programmable chemical controllers made from DNA*. [Chen et al. (2013)]

### single-reaction-schemes:
  * `turberfield2003-motor.ts` Figure 4 of [Turberfield et al. (2003)] **missing**

  * `zhang2007-catalyst.ts` Figure 1 of [Zhang et al. (2007)]

  Coded by: Stefan Badelt (badelt@caltech.edu).

  * `zhang2007-autocatalyst.ts` Figure 4 of [Zhang et al. (2007)]

  Coded by: Stefan Badelt (badelt@caltech.edu).

  * `song2016-add.ts` Figure 3-5 of [Song et al. (2016)] **missing**

  * `song2016-subtract.ts` Figure 6-8 of [Song et al. (2016)]

  Coded by: Stefan Badelt (badelt@caltech.edu).

  * `song2016-multiply.ts` Figure 9-14 of [Song et al. (2016)] **missing**

  * `song2016-2amplify.ts` Figure 15-16 of [Song et al. (2016)] **missing**
  
### Last Update
July, 27th, 2016

[//]: References
[Turberfield et al. (2003)]: <http://dx.doi.org/10.1103/PhysRevLett.90.118102>
[Zhang et al. (2007)]: <http://dx.doi.org/10.1126/science.1148532>
[Soloveichik et al. (2010)]: <http://dx.doi.org/10.1073/pnas.0909380107>
[Qian et al. (2011)]: <http://dx.doi.org/10.1007/978-3-642-18305-8_12>
[Qian & Winfree (2011)]: <http://dx.doi.org/10.1126/science.1200520>
[Lakin et al. (2011)]: <http://dx.doi.org/10.1098/rsif.2011.0343>
[Cardelli (2011)]: <http://dx.doi.org/10.1007/s11047-010-9236-7>
[Cardelli (2013)]: <http://dx.doi.org/10.1017/S0960129512000102>
[Chen et al. (2013)]: <http://dx.doi.org/10.1038/NNANO.2013.189>
[Thachuk et al. (2015)]: <http://dx.doi.org/10.1007/978-3-319-21999-8_9>
[Srinivas (2015)]: <http://www.dna.caltech.edu/Papers/NiranjanSrinivas_2015_thesis.pdf>
[Song et al. (2016)]: <http://dx.doi.org/10.1021/acssynbio.6b00144>

