
Examples
========

Here we will give a demonstration of how to use some of `smact`'s features. For a full set of
work-through examples in Jupyter notebook form check out
`the examples section of our GitHub repo <https://github.com/WMD-group/SMACT/tree/master/examples>`_.
For workflows that have been used in real examples and in published work, visit our
`separate repository <https://github.com/WMD-group/smact_workflows>`_.

===========================
Element and species classes
===========================

The element and species classes are at the heart of :mod:`smact`'s functionality. Elements are the
elements of the periodic table. Species are elements, with some additional information; the
oxidation state and the coordination environment (if known). So for example the element iron
can have many oxidation states and those oxidation states can have many coordination
environments.

.. code:: python

    import smact

    iron = smact.Element('Fe')
    print("The element %s has %i oxidation states. They are %s." %
    (iron.symbol, len(iron.oxidation_states), iron.oxidation_states))

    The element Fe has 8 oxidation states. They are [-2, -1, 1, 2, 3, 4, 5, 6].

When an element has an oxidation state and coordination environment then it has additional
features. For example the Shannon radius [1]_ of the element, this is often useful for calculating
radius ratio rules [2]_, or for training neural networks [3]_ .

.. code:: python

    iron_square_planar = smact.Species('Fe', 2, '4_n')
    print('Square planar iron has a Shannon radius of %s Angstrom' % iron_square_planar.shannon_radius)

    Square planar iron has a Shannon radius of 0.77 Angstrom

=============
List building
=============

Often when using :mod:`smact` the aim will to be to search over combinations of a set of elements. This
is most efficiently achieved by setting up a dictionary of the elements that you want to search
over. The easiest way to achieve this in :mod:`smact` is to first create a list of the symbols of the elements
that you want to include, then to build a dictionary of the corresponding element objects.

The list can be built by hand, or if you want to cover a given range there is a helper function.

.. code:: python

    import smact

    elements = smact.ordered_elements(13, 27)
    print(elements)

    ['Al','Si','P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co']

For doing searches across combinations of elements it is then quickest to load the element objects into
a dictionary and search by key. This avoids having to repopulate the element class at each iteration of
the search.

.. code:: python

    element_list = smact.element_dictionary(elements)
    print(element_list)

    {'Al': <smact.Element at 0x10ecc5890>,
     'Ar': <smact.Element at 0x10ecc5cd0>,
     'Ca': <smact.Element at 0x10ecc5a10>,
     'Cl': <smact.Element at 0x10ecc5d90>,
     'Co': <smact.Element at 0x10ecc5f90>,
     'Cr': <smact.Element at 0x10ecc5ed0>,
     'Fe': <smact.Element at 0x10ecc5f50>,
     'K': <smact.Element at 0x10ecc5e90>,
     'Mn': <smact.Element at 0x10ecc5f10>,
     'P': <smact.Element at 0x10ecc5990>,
     'S': <smact.Element at 0x10ecc5e10>,
     'Sc': <smact.Element at 0x10ecc5150>,
     'Si': <smact.Element at 0x10e8bf190>,
     'Ti': <smact.Element at 0x10ecc5dd0>,
     'V': <smact.Element at 0x10ecc5e50>}

====================
Neutral combinations
====================

One of the most basic tests for establishing sensible combinations of elements is that they should form charge neutral
combinations. This is a straightforward combinatorial problem of comparing oxidation states and allowed stoichiometries.

:math:`\Sigma_i Q_in_i = 0`

where :math:`i` are the elements in the compound and :math:`Q` are the charges. We have a special function, ``smact_filter``,
which does this checking for a list of elements. The ``smact_filter`` also ensures that all elements specified to be anions
have electronegitivities greater than all elements specified to be cations.

As input ``smact_filter`` takes:

* ``els`` : a tuple of the elements to search over (required)
* ``threshold``: the upper limit of the stoichiometric ratios (default = 8)
* ``species_unique``: whether or not we want to consider elements in different oxidation states as unique in our results (default is False).

We can look for neutral combos.

.. code:: python

    import smact.screening

    elements = ['Ti', 'Al', 'O']
    space = smact.element_dictionary(elements)
    # We just want the element items from the dictionary
    eles = [e[1] for e in space.items()]
    # We set a threshold for the stoichiometry of 4
    allowed_combinations = smact.screening.smact_filter(eles, threshold=4)
    print(allowed_combinations)

    [(('Ti', 'Al', 'O'), (1, 3, 3)),
     (('Ti', 'Al', 'O'), (2, 3, 4)),
     (('Ti', 'Al', 'O'), (3, 1, 4)),
     (('Ti', 'Al', 'O'), (1, 4, 4)),
     (('Ti', 'Al', 'O'), (3, 1, 2)),
     (('Ti', 'Al', 'O'), (3, 2, 4)),
     (('Ti', 'Al', 'O'), (1, 2, 3)),
     (('Ti', 'Al', 'O'), (1, 3, 4)),
     (('Ti', 'Al', 'O'), (2, 4, 3)),
     (('Ti', 'Al', 'O'), (2, 1, 3)),
     (('Ti', 'Al', 'O'), (4, 2, 3)),
     (('Ti', 'Al', 'O'), (1, 3, 2)),
     (('Ti', 'Al', 'O'), (1, 2, 4)),
     (('Ti', 'Al', 'O'), (1, 1, 2)),
     (('Ti', 'Al', 'O'), (1, 2, 2)),
     (('Ti', 'Al', 'O'), (1, 1, 4)),
     (('Ti', 'Al', 'O'), (3, 1, 3)),
     (('Ti', 'Al', 'O'), (2, 1, 4)),
     (('Ti', 'Al', 'O'), (1, 1, 1)),
     (('Ti', 'Al', 'O'), (2, 2, 3)),
     (('Ti', 'Al', 'O'), (4, 1, 3)),
     (('Ti', 'Al', 'O'), (1, 1, 3)),
     (('Ti', 'Al', 'O'), (1, 4, 3)),
     (('Ti', 'Al', 'O'), (2, 1, 2))]

There is `an example <https://github.com/WMD-group/SMACT/blob/master/examples/Counting/Generate_compositions_lists.ipynb>`_
of how this function can be combined with multiprocessing to rapidly explore large subsets of chemical space.

==========================
Compound electronegativity
==========================

One property that is often used in high-throughput screening where band alignment is important is the
compound electronegativity. Ginley and Butler showed how the simple geometric mean of the
electronegitivities of a compound could be used to predict flat band potentials [4]_. :mod:`smact` has a built
in function to calculate this property for a given composition.

.. code:: python

    import smact.properties

    compound_electronegs = [smact.properties.compound_electroneg(elements = a[0], stoichs = a[1]) for \\
    a in allowed_combinations]

    print(compound_electronegs)

    [4.319343517137848,
     4.729831837874991,
     4.462035251666306,
     4.337155845378665,
     5.0575817742802025,
     4.777171739263751,
     4.427325394494835,
     5.34030430325585,
     4.583732423414276,
     4.980129115226567,
     4.652147502981397,
     5.284089129411956,
     4.726884428924315,
     4.373001170931816,
     4.808336266651247,
     5.041995471272069,
     4.587722671269271,
     5.437592861777965,
     5.010966817423813,
     4.964781503487637,
     4.768922515748819,
     4.409142747625072,
     5.74200359520417,
     4.677126472294396]

===============================
Interfacing to machine learning
===============================

When preparing to do machine learning, we have to convert the compositions that we have into
something that can be fed into an algorithm. Many of the properties provided in :mod:`smact` are suitable for this,
one can take properties like electronegativity, mass, electron affinity etc etc (for the full list see
:ref:`smact_module`).

One useful representation that is often used in machine learning is the one-hot-vector formulation. A similar
construction to this can be used to encode a chemical formula. A vector of length of the periodic table is
set up and each element set to be a number corresponding to the stoichiometric ratio of that element in the compound.
For example we could convert :math:`Ba(OH)_2`

.. code:: python

   ml_vector = smact.screening.ml_rep_generator(['Ba', 'H', 'O'], stoichs=[1, 2, 2])

There is also `an example <https://github.com/WMD-group/SMACT/blob/master/examples/Counting/Generate_compositions_lists.ipynb>`_
demonstrating the conversion of charge neutral compositions produced by `smact` to a list of formulas using Pymatgen,
or to a Pandas dataframe, both of which could then be used as input for a machine learning algorithm.
For a full machine learning example that uses `smact`, there is a repository `here <https://github.com/WMD-group/Solar_oxides_data>`_ 
which demonstrates a search for solar energy materials from the four-component (quaternary) oxide materials space.

.. [1]  "Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides".
         Acta Crystallogr A. 32: 751–767, 1976

.. [2]  "Crystal Structure and Chemical Constitution" Trans. Faraday Soc. 25, 253-283, 1929.

.. [3] "Deep neural networks for accurate predictions of crystal stability" Nat. Comms. 9, 3800, 2018.

.. [4] "Prediction of Flatband Potentials at Semiconductor‐Electrolyte Interfaces from Atomic Electronegativities"
       J. Electrochem. Soc. 125, 228-32, 1975.
