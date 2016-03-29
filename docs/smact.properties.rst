smact.properties module
=======================

A collection of tools for estimating useful properties.

The "electronegativity of a compound" computed with
:func:`compound_electroneg` is the rescaled geometric mean of
electronegativity used in Nethercot's recipe for estimating the
photoelectric threshold: [1]_

.. math:: \Phi^{AB} = 2.86(\chi_{A}\chi_{B})^{1/2} + E_{g} / 2.

In other words, the computed group
:math:`2.86(\chi_{A}\chi_{B})^{1/2}`
is the mid-gap energy and the VBM/CBM positions can be estimated by
subtracting/adding half of the band gap :math:`E_g`.

.. [1] Nethercot, A. H. (1974). *Phys. Rev. Lett.*, **33**, 1088–1091. http://dx.doi.org/10.1103/PhysRevLett.33.1088

.. automodule:: smact.properties
    :members:
    :undoc-members:
    :show-inheritance:
