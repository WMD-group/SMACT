### This Jupyter notebook creates ntype ptype possiblie dopants for input species
### using SMACT structure prediction 
### Working with Kieth from SCIML and Anthony 

from smact.structure_prediction import mutation, utilities
import smact
import re
from pymatgen.util import plotting

class Doper:
    '''
    A class to search for n & p type dopants
    Methods: get_dopants, plot_dopants
    '''
    def __init__(self, original_species: tuple[str], 
                  num_dopants=5
                  #match_oxi_sign=False):
                ):
      self.original_species = original_species
      self.num_dopants = num_dopants
      #self.match_oxi_sign = match_oxi_sign
    def get_dopants(self) -> dict:
        '''
    Note currently limited to binaries
    Args:
        ex) get_dopants(('Ti4+','O2-'))
        
        original_species (tuple(str)) = ('Cd2+', 'O2-')
        num_dopants (int) = The number of suggestions to return for n- and p-type dopants.
        match_oxi_sign (bool) = ? shoud do some tests
    
    Returns:
        (dict): Dopant suggestions, given as a dictionary with keys 
        "n_type_cation", "p_type_cation", "n_type_anion", "p_type_anion".
    '''
    
        cat, an = self.original_species
    # utilities.parse_spec('Cd2+') -> ('Cd', 2):(str, int)
    # utilities.parse_spec('O2-') -> ('O', -2):(str, int)
        higher_charge = {'element': utilities.parse_spec(cat)[0],
                        'charge': utilities.parse_spec(cat)[1],}
        lower_charge = {'element': utilities.parse_spec(an)[0],
                        'charge': utilities.parse_spec(an)[1],}  
    
        CM = mutation.CationMutator.from_json()
    
        #call all elements
        element_objects = list(smact.element_dictionary().values())
        poss_n_type_cat = []
        poss_p_type_cat = []
        poss_n_type_an = []
        poss_p_type_an = []
        
        # cation loop to identify species
        for element in element_objects:
            #[-2, -1, 0, +1, +2]
            oxi_state = element.oxidation_states
            els_symbol = element.symbol
            for state in oxi_state:
                if state > higher_charge['charge']:
                    poss_n_type_cat.append((els_symbol, state))
                elif state < higher_charge['charge'] and state > 0:
                    poss_p_type_cat.append((els_symbol, state))
        
        #anion loop to identify species
        for element in element_objects:
            oxi_state = element.oxidation_states
            els_symbol = element.symbol
            for state in oxi_state:
                if state > lower_charge['charge'] and state < 0:
                    poss_n_type_an.append((els_symbol, state))
                elif state < lower_charge['charge']:
                    poss_p_type_an.append((els_symbol, state))
        
        # cation loop to calculate substitution probability 
        parsed_poss_n_type_cat = [utilities.unparse_spec(specie) for specie in poss_n_type_cat]
        parsed_poss_p_type_cat = [utilities.unparse_spec(specie) for specie in poss_p_type_cat]
        # anion loop to calculate substitution probability 
        parsed_poss_n_type_an = [utilities.unparse_spec(specie) for specie in poss_n_type_an]
        parsed_poss_p_type_an = [utilities.unparse_spec(specie) for specie in poss_p_type_an]
        
        #calculate substitutional probability for all elememts stored in the list
        n_type_cat = [(specie, CM.sub_prob(cat, specie))
                 for specie in parsed_poss_n_type_cat]
        p_type_cat = [(specie, CM.sub_prob(cat, specie))
                 for specie in parsed_poss_p_type_cat]
        n_type_an = [(specie, CM.sub_prob(an, specie))
                 for specie in parsed_poss_n_type_an]
        p_type_an = [(specie, CM.sub_prob(an, specie))
                 for specie in parsed_poss_p_type_an]
        
        #[('B3+', 0.003), ('C4+', 0.001), (), (), ...] : list(tuple(str, float))
        #sort by probability
        n_type_cat.sort(key=lambda x: x[1], reverse=True)
        p_type_cat.sort(key=lambda x: x[1], reverse=True)
        n_type_an.sort(key=lambda x: x[1], reverse=True)
        p_type_an.sort(key=lambda x: x[1], reverse=True)
        
        self.results = {"n-type cation substitutions": n_type_cat[:self.num_dopants], 
                "p-type cation substitutions": p_type_cat[:self.num_dopants],
                "n-type anion substitutions": n_type_an[:self.num_dopants],
                "p-type anion substitutions": p_type_an[:self.num_dopants]}
        
        # return the top (num_dopants) results for each case 
        return self.results
    
    
    def plot_dopants(self):
        '''
        Uses pymatgen plotting utilities to plot the results of doping search
        '''
        for val in self.results.values():
            dict_results = dict((utilities.parse_spec(x)[0], y) for x, y in val)
            plotting.periodic_table_heatmap(elemental_data=dict_results, cmap='rainbow',
                                                blank_color='gainsboro', edge_color='white')
    


