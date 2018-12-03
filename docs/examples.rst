
Examples
========

Here we will give a demonstration of how to use some of SMACT's features. For a full set of 
work-through tutorials in Jupyter notebook form check out `the tutorials section of our GitHub
repo <https://github.com/WMD-group/SMACT/tree/master/examples>`_

===========================
Element and species classes
===========================

The element and species classes are at the heart of SMACT's functionality. Elements are the
elements of the periodic table. Species are elements, with some additional information; the
oxidation state and the coordination environment (if known). So for example the element iron
can have many oxidation states and those oxidation states can have many coordination 
environments.

.. code:: python

    import smact

    iron = smact.Element('Fe')
    print("The element %s has %i oxidation states. They are %s." % 
    (iron.symbol, len(iron.oxidation_states), iron.oxidation_states))

When an element has an oxidation state and coordination environment then it has additional
features. For example the Shannon radius [1]_ of the element, this is often useful for calculating
radius ratio rules [2]_, or for training neural networks [3]_ .

=============
List building
=============

Often when using SMACT the aim will to be to search over combinations of a set of elements. This
is most efficiently achieved by setting up a dictionary of the elements that you want to search
over. The easiest way to achive this in SMACT is to first create a list of the symbols of the elements
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

.. math:: \Sum_i Q_in_i = 0

where :math:`i` are the elements in the compound and :math:`Q` are the charges. We have a special function, ``smact_test``,
which does this checking for a list of elements. As input ``smact_test`` takes:

* ``els`` : a tuple of the elements to search over (required)
* ``threshold``: the upper limit of the stoichiometric ratios (default = 8)
* ``include``: a list of elements that **must** be in every compound (optional)

We can look for neutral combos.

.. code:: python

    import smact.screening

    elements = ['Ti', 'Al', 'O']
    space = smact.element_dictionary(elements)
    # We just want the element items from the dictornary
    eles = [e[1] for e in space.items()]
    # With an 
    


.. [1]  "Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides". 
         Acta Crystallogr A. 32: 751–767, 1976

.. [2]  "Crystal Structure and Chemical Constitution" Trans. Faraday Soc. 25, 253-283, 1929.

.. [3] "Deep neural networks for accurate predictions of crystal stability" Nat. Comms. 9, 3800, 2018.

