This section will get you started using the enumerator API as quickly as possible. For more details on how the code is layed out, see the Architecture section below. 

To create an ``Enumerator`` object: 

*	Load an enumerator from an input file, using one of the input functions in the `input` module::
	
		from input import input_pil
		enum = input_pil('my_example_system.pil')

*	Alternatively, just write some complexes in Kernel notation::

		from input import enum_from_kernel
		input_string = """
		a(b(+ c d) e) f g
		f(g(+))
		"""
		enum = enum_from_kernel(input_string)

*	Or manually create some domains, strands, and complexes::

		from enumerator import Enumerator
		from utils import Complex, Strand, Domain

		domains = [Domain("a"), Domain("a", is_complement=True)]
		strands = [Strand("s1", domains[:])]
		complexes = [Complex("c1", strands[:], [[(0,1),(0,0)]] )]
		enum = Enumerator(domains, strands, complexes)

To enumerate reactions::

	enum.enumerate()

To get the resulting... ::

	# complexes
	enum.complexes

	# reactions
	enum.reactions

	# reaction rate constants
	[r.rate() for r in enum.reactions]

	# resting states
	enum.resting_states

To condense::

	from condense import condense_resting_states
	results = condense_resting_states(enum, compute_rates = True)
	
	# condensed reactions
	results['reactions']

To write output::

	from output import output_pil

	# un-condensed version
	output_pil(enum, 'my_example_system-enum.pil')

	# condensed version
	output_pil(enum, 'my_example_system-enum-condensed.pil', output_condensed = True)

.. note:: If you do _not_ call ``enum.enumerate()`` before running one of the output functions, you'll get some errors. If you just want to write output without enumerating, call ``enum.dry_run()`` first. 

