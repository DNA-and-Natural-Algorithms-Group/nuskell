# Nuskell -- A verifying and optimizing compiler for toehold-mediated strand displacement circuits.

**Nuskell** is a compiler written in python2.7. It translates formal chemical
reaction networks (CRNs) into a domain-level specification using a translation
scheme.  Full documentation will become available at [Nuskell].

### Examples

If you want to implement a formal CRN using a particular translation-scheme:

```
  $ echo "A + B <=> X + Y; X -> A" | nuskell --ts scheme.ts [options]
```

If you want to verify the equivalence between a formal CRN and an implementation CRN:

```
  $ echo "A + B <=> X + Y; X -> A" | nuskell --compare implementation.crn --verify bisimulation 
```

## Installation
```
  python setup.py install
```

### local installation
```
  python setup.py install --user
```
  
## Version
0.2.0

### Development / Unittests
```
  python setup.py test
```

### Build the documentation
```
  sphinx-build -b html docs ~/your/html/sourcedir/
```

### License
MIT

[//]: References
[nuskell]: <https://dna.caltech.edu/nuskell>
[peppersuite]: <https://dna.caltech.edu/peppersuite>

