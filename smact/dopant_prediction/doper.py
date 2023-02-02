### This Jupyter notebook creates ntype ptype possiblie dopants for input species
### using SMACT structure prediction
### Working with Kieth from SCIML and Anthony

###Doper ver 2

# Now 'Doper' can generate possible n-type p-type dopants for multicomponent materials (i.e. Ternary, Quaternary etc).
# Can plot the result of doping search within a single step
# """ex) test= Doper(('Cu1+','Zn2+','Ge4+','S2-'))
#         test.get_dopants(num_dopants = 10, plot_heatmap = True)"""


from typing import List, Tuple

from pymatgen.util import plotting

import smact
from smact.structure_prediction import mutation, utilities


class Doper:
    """
    A class to search for n & p type dopants
    Methods: get_dopants, plot_dopants

    Attributes:
        original_species: A tuple which describes the constituent species of a material. For example:

            >>> test= Doper(('Cu1+','Zn2+','Ge4+','S2-'))
            >>> test.original_species
                ('Cu1+','Zn2+','Ge4+','S2-')

    """

    def __init__(self, original_species: Tuple[str, ...]):
        """
        Intialise the `Doper` class with a tuple of species

        Args:
            original_species: See :class:`~.Doper`.

        """
        self.original_species = original_species

    def _get_cation_dopants(
        self, element_objects: List[smact.Element], cations: List[str]
    ):
        poss_n_type_cat = []
        poss_p_type_cat = []

        for element in element_objects:
            # [-2, -1, 0, +1, +2]
            oxi_state = element.oxidation_states
            el_symbol = element.symbol
            for state in oxi_state:
                for cation in cations:
                    _, charge = utilities.parse_spec(cation)
                    if state > charge:
                        poss_n_type_cat.append(
                            utilities.unparse_spec((el_symbol, state))
                        )
                    elif state < charge and state > 0:
                        poss_p_type_cat.append(
                            utilities.unparse_spec((el_symbol, state))
                        )

        return poss_n_type_cat, poss_p_type_cat

    def _get_anion_dopants(
        self, element_objects: List[smact.Element], anions: List[str]
    ):
        poss_n_type_an = []
        poss_p_type_an = []

        for element in element_objects:
            oxi_state = element.oxidation_states
            el_symbol = element.symbol
            for state in oxi_state:
                for anion in anions:
                    _, charge = utilities.parse_spec(anion)
                    if state > charge and state < 0:
                        poss_n_type_an.append(
                            utilities.unparse_spec((el_symbol, state))
                        )
                    elif state < charge:
                        poss_p_type_an.append(
                            utilities.unparse_spec((el_symbol, state))
                        )
        return poss_n_type_an, poss_p_type_an

    def _plot_dopants(self, results: dict):
        """
        Uses pymatgen plotting utilities to plot the results of doping search
        """
        for key, val in results.items():
            dict_results = {utilities.parse_spec(x)[0]: y for x, y in val}
            plotting.periodic_table_heatmap(
                elemental_data=dict_results,
                cmap="rainbow",
                blank_color="gainsboro",
                edge_color="white",
            )

    def get_dopants(
        self,
        num_dopants: int = 5,
        plot_heatmap: bool = False,
    ) -> dict:
        """
        Args:
            num_dopants (int): The number of suggestions to return for n- and p-type dopants.
            plot_heatmap (bool): If True, the results of the doping search are plotted as heatmaps

        Returns:
            (dict): Dopant suggestions, given as a dictionary with keys
            "n_type_cation", "p_type_cation", "n_type_anion", "p_type_anion".

        Examples:
            >>> test = Doper(('Ti4+','O2-'))
            >>> print(test.get_dopants(num_dopants=2))
                {'n-type cation substitutions': [('Ta5+', 8.790371775858281e-05),
                ('Nb5+', 7.830035204694342e-05)],
                'p-type cation substitutions': [('Na1+', 0.00010060400812977031),
                ('Zn2+', 8.56373996146833e-05)],
                'n-type anion substitutions': [('F1-', 0.01508116810515677),
                ('Cl1-', 0.004737202729901607)],
                'p-type anion substitutions': [('N3-', 0.0014663800608945628),
                ('C4-', 9.31310255126729e-08)]}
        """

        cations = []
        anions = []
        try:
            for ion in self.original_species:
                _, charge = utilities.parse_spec(ion)
                if charge > 0:
                    cations.append(ion)
                elif charge < 0:
                    anions.append(ion)
        except Exception as e:
            print(e, "charge is not defined")

        CM = mutation.CationMutator.from_json()

        # call all elements
        element_objects = list(smact.element_dictionary().values())

        poss_n_type_cat, poss_p_type_cat = self._get_cation_dopants(
            element_objects, cations
        )
        poss_n_type_an, poss_p_type_an = self._get_anion_dopants(
            element_objects, anions
        )

        n_type_cat, p_type_cat, n_type_an, p_type_an = [], [], [], []
        for cation in cations:
            for n_specie, p_specie in zip(poss_n_type_cat, poss_p_type_cat):
                n_type_cat.append((n_specie, CM.sub_prob(cation, n_specie)))
                p_type_cat.append((p_specie, CM.sub_prob(cation, p_specie)))

        for anion in anions:
            for n_specie, p_specie in zip(poss_n_type_an, poss_p_type_an):
                n_type_an.append((n_specie, CM.sub_prob(anion, n_specie)))
                p_type_an.append((p_specie, CM.sub_prob(anion, p_specie)))

        # [('B3+', 0.003), ('C4+', 0.001), (), (), ...] : list(tuple(str, float))
        # sort by probability
        n_type_cat.sort(key=lambda x: x[1], reverse=True)
        p_type_cat.sort(key=lambda x: x[1], reverse=True)
        n_type_an.sort(key=lambda x: x[1], reverse=True)
        p_type_an.sort(key=lambda x: x[1], reverse=True)

        results = {
            "n-type cation substitutions": n_type_cat[:num_dopants],
            "p-type cation substitutions": p_type_cat[:num_dopants],
            "n-type anion substitutions": n_type_an[:num_dopants],
            "p-type anion substitutions": p_type_an[:num_dopants],
        }

        # plot heatmap

        if plot_heatmap:
            self._plot_dopants(results)

        # return the top (num_dopants) results for each case
        return results
