#### 'nuskell': A Comiler for DNA strand displacement computation ####

Full documentation will be available at: [http://www.dna.caltech.edu/nuskell/](http://www.dna.caltech.edu/nuskell/)

usage:

  ./compile translation_scheme.nus [options] < formal.crn > dom_design.pil

nuskell is a compiler written in python2.7. The nuskell compiler translates 
formal chemical reaction networks into a domain-level specification using
a translation scheme written in the nuskell programming language.

#### Build the documentation ####

sphinx-build -b html docs ~/your/html/sourcedir/

##### Installation #####

python setup.py install

#### local installation ####
python setup.py install --prefix=/your/destination/nuskell-x.y.r

**Do not forget to set your environment variables when using local installations**

export PATH=/your/destination/nuskell-x.y.r/bin:$PATH
export PYTHONPATH=/your/destination/nuskell-x.y.r/lib/python2.7/site-packages/:$PYTHONPATH

##### Devolment / Unittests #####

python setup.py test

