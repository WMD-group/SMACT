smact.properties module
=======================

A collection of tools for estimating physical properties
based on chemical composition.

The "electronegativity of a compound" computed with
:func:`compound_electroneg` is the rescaled geometric mean of
electronegativity used in Nethercot's recipe for estimating the
photoelectric threshold: [1]_

.. math:: \Phi^{AB} = 2.86(\chi_{A}\chi_{B})^{1/2} + E_{g} / 2.

In other words, the computed group
:math:`2.86(\chi_{A}\chi_{B})^{1/2}`
is the mid-gap energy. The valence band maximum/conduction band minimum positions
can be estimated by subtracting/adding half of the band gap :math:`E_g`.
This is an extension Mulliken's electronegativity scale in which
:math:`\chi_{A} = (I_{A} + E_{A})/2` (where :math:`I` and :math:`E`
are respectively the ionisation potential and electron affinity.) [2]_

.. [1] Nethercot, A. H., *Prediction of Fermi energies and photoelectric thresholds based on electronegativity concepts* Phys. Rev. Lett. **33**, 1088–1091 (1974). http://dx.doi.org/10.1103/PhysRevLett.33.1088

.. [2] Mulliken, R. S., *A new electroaffinity scale; together with data on valence states and on valence ionization potentials and electron affinities* J. Chem. Phys. **2**, 782 (1934). http://dx.doi.org/10.1063/1.1749394

.. automodule:: smact.properties
    :members:
    :undoc-members:
    :show-inheritance:
