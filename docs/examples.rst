
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
over 

.. [1]  "Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides". 
         Acta Crystallogr A. 32: 751–767, 1976

.. [2]  "Crystal Structure and Chemical Constitution" Trans. Faraday Soc. 25, 253-283, 1929.

.. [3] "Deep neural networks for accurate predictions of crystal stability" Nat. Comms. 9, 3800, 2018.

