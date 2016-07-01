# Nuskell -- A Comiler for DNA strand displacement computation

**Nuskell** is a compiler written in python2.7. The Nuskell compiler translates 
formal chemical reaction networks into a domain-level specification using
a translation scheme written in the nuskell programming language.
Full documentation will be available at [Nuskell]

### Examples

  $ ./compile translation_scheme.nus [options] < formal.crn > dom_design.pil

## Installation
  python setup.py install

### local installation
  python setup.py install --prefix=/your/destination/nuskell-x.y.r
  
Do not forget to set your environment variables when using local installations:
  
  export PATH=/your/destination/nuskell-x.y.r/bin:$PATH
  export PYTHONPATH=/your/destination/nuskell-x.y.r/lib/python2.7/site-packages/:$PYTHONPATH
  
## Version
0.2.0

### Devolment / Unittests
  python setup.py test

### Build the documentation
  sphinx-build -b html docs ~/your/html/sourcedir/

### Todos

 - Write Tests
 - Sphinx Documentation
 - Link to [peppersuite]
 - Release

### License
MIT

[//]: References
[nuskell]: <https://dna.caltech.edu/nuskell>
[peppersuite]: <https://dna.caltech.edu/peppersuite>

