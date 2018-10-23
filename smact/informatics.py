###############################################################################
# Copyright Daniel Davies (2018)                                              #
#                                                                             #
# This file is part of SMACT: informatics.py is free software: you can        #
# redistribute it and/or modify it under the terms of the GNU General Public  #
# License as published by the Free Software Foundation, either version 3 of   #
# the License, or (at your option) any later version.  This program is        #
# distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;   #
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A       #
# PARTICULAR PURPOSE.  See the GNU General Public License for more details.   #
# You should have received a copy of the GNU General Public License along     #
# with this program.  If not, see <http://www.gnu.org/licenses/>.             #
#                                                                             #
###############################################################################

###############################################################################
# Collection of functions for manipulating lists of Pymatgen structures       #
# and compositions. The outputs from many of these functions can be used      #
# directly in the oxidationstates module.                                     #
# NB:                                                                         #
# This module currently depends on Pymatgen, although this is not currently a #
# dependency of SMACT as a whole. See http://pymatgen.org/.                   #
###############################################################################

###  Imports
# General
import os, re, json
from tqdm import tqdm
from collections import Counter
import itertools

# pymatgen
from pymatgen import Structure, Specie, MPRester
from pymatgen.analysis.structure_prediction.substitutor import Substitutor
from pymatgen.io.cif import CifWriter
from smact import ordered_elements, Element, neutral_ratios, metals
from smact.screening import pauling_test

### Functions
def get_struc_list(json_name):
    """Import pymatgen Structure objects from a json as a list.
    TODO DWD: Move to more appropriate place. This is not oxidation state specific.
    Args:
        json_name (string): Path to json file containing a list of dicts in which one key is
        'structure'. The value of this entry should be a Pymatgen Structure object as a dict.
    Returns:
        struc_list (list): list of dicts containing 'id' and 'structure'.
    """
    with open(json_name, 'r') as f:
        saved_strucs = json.load(f)

    struc_list = []
    for i, entry in enumerate(tqdm(saved_strucs)):
        struc_list.append({'structure': Structure.from_dict(entry['structure']),
        'id': entry['id'] })
    return(struc_list)


def mp_filter(criteria, api_key=None):
    """Get a list of Materials Project task ids that match a set of given criteria.
    The criteria should be in the format used for a MP query.
    TODO DWD: Move to more appropriate place. This is not oxidation state specific.
    Args:
        criteria (dict): Criteria that can be used in an MPRester query
        api_key (str): Materials Project API key (from your MP dashboard)
    """
    if not api_key:
        print('You need to supply an api key.')
    else:
        m = MPRester(os.environ.get("MP_API_KEY"))
        properties = ['task_id']
        struc_filter = m.query(criteria,properties)
        id_list = [i['task_id'] for i in struc_filter]
        return id_list


def get_unique_species(structures, ordering='ptable', reverse=False,
                       cation_only = True, metal_only = True):
    """Given a set of pymatgen structures, in the form of dictionaries where the
    Structure is keyed as 'structure', returns a list of all the different
    Species present in that set.
    Args:
        structures (list): Dictionaries containing pymatgen Structures.
        ordering('string'): How to order the Species:
            ptable: order by periodic table position.
            Can be set to None.
        reverse (bool): Whether to reverse the ordering (descending order).
        cation_only (bool): Whether to only consider species in positive oxidation
            states.
        metal_only (bool): Whether to only consider metal elements.
    Returns:
        species_list (list): Unique species that are exhibited in the structures.

    """
    # Initially comb through the structures for all unique species
    species_list = []
    for i in structures:
        for sp in i['structure'].composition:
            species_list.append((sp))
    species_list = list(set(species_list))

    ordered_el = ordered_elements(1,103)
    # Turn into tuples for easy sorting
    species_list = [(i.symbol, i.oxi_state) for i in species_list]
    if ordering == 'ptable':
        species_list.sort(key = lambda x: (ordered_el.index(x[0]),x[1]), reverse=reverse)
        print("Species ordered by periodic table position.")
    else:
        print('Did not reorder the list of species...')

    # Turn back into Species objects
    species_list = [Specie(i[0], i[1]) for i in species_list]

    if metal_only:
        print('Metals only: ON')
        species_list = [i for i in species_list if (i.symbol in metals)]

    if cation_only:
        print('Cations only: ON')
        species_list = [i for i in species_list if (i.oxi_state > 0)]

    print("First species: {0}  last species: {1}".format(species_list[0], species_list[-1]))
    return species_list


def species_totals(structures, count_elements=False, anions=[],
                   edit_structures_dicts=True, return_species_list=False):
    """Given a set of pymatgen structures in the form of dictionaries where
    the Structure is keyed as 'structure', returns the number
    of compounds that features each Species.
    Args:
        structures (list): dictionaries containing pymatgen Structures.
        count_elements (bool): switch to counting elements not species.
        anions (list): Pymatgen.Species anions of interestself.
        edit_structure_dicts (bool): Modify the dicts in the structures list
        to add a 'most_eneg_anion' key.
    Returns:
        totals (dict): Totals of each species in structure list.
        or an_containing (dict): Totals of each species separated by anion.
        species_list (optional): List of species for structures as generated by
        get_unique_species.
    """
    # Simple method if simply counting all species or elements
    if not anions:
        totals = []
        if count_elements:
            for i in structures:
                comp = [j.symbol for j in i['structure'].composition]
                totals.append(comp)
            totals = [i for sublist in totals for i in sublist]
            totals = dict(Counter(totals))
        else:
            for i in structures:
                comp = [j for j in i['structure'].composition]
                totals.append(comp)
            totals = [i for sublist in totals for i in sublist]
            totals = dict(Counter(totals))
    # Method used if collecting count per anion 
    else:
        totals = {}
        for anion in tqdm(anions):
            an_containing = []
            for i in structures:
                if anion in i['structure'].composition:
                    # Check whether anion is most electronegative element
                    an_eneg = Element(anion.symbol).pauling_eneg
                    all_enegs = [Element(sp.symbol).pauling_eneg for \
                    sp in i['structure'].composition]
                    if all(eneg <= an_eneg for eneg in all_enegs):
                        comp = [j for j in i['structure'].composition]
                        an_containing.append(comp)
                        if edit_structures_dicts:
                            i['most_eneg_anion'] = anion

            an_containing = [i for sublist in an_containing for i in sublist]
            an_containing = dict(Counter(an_containing))
            an_containing.pop(anion)
            totals[anion] = an_containing

    # Return objects based on whether species list required
    if return_species_list:
        return(totals,get_unique_species(structures))
    else:
        return(totals)

def ternary_smact_combos(position1, position2, position3, threshold = 8):
    """ Combinatorially generate Pymatgen Species compositions using SMACT when up to three different
        lists are needed to draw species from (e.g. Ternary metal halides.)
    Args:
        position(n) (list of species): Species to be considered iteratively for each
                                     position.
        threshold (int): Max stoichiometry threshold.
    Returns:
        species_comps (list): Compositions as tuples of Pymatgen Species objects.
        """

    initial_comps_list = []
    for sp1, sp2, an in tqdm(itertools.product(position1, position2, position3)):
        e1, oxst1 = sp1.symbol, int(sp1.oxi_state)
        eneg1 = Element(e1).pauling_eneg
        e2, oxst2 = sp2.symbol, int(sp2.oxi_state)
        eneg2 = Element(e2).pauling_eneg
        e3, oxst3 = an.symbol, int(an.oxi_state)
        eneg3 = Element(e3).pauling_eneg

        symbols = [e1,e2,e3]
        ox_states = [oxst1, oxst2, oxst3]
        cn_e, cn_r = neutral_ratios(ox_states, threshold=threshold)

        if cn_e:
            enegs = [eneg1,eneg2,eneg3]
            eneg_ok = pauling_test(ox_states, enegs, symbols=symbols, repeat_cations=False)
            if eneg_ok:
                for ratio in cn_r:
                    comp = (symbols, ox_states, list(ratio))
                    initial_comps_list.append(comp)
    print('Number of compositions before reduction:  {}'.format(len(initial_comps_list)))

    # Create a list of pymatgen species for each comp
    print('Converting to Pymatgen Species...')
    species_comps = []
    for i in tqdm(initial_comps_list):
        comp = {}
        for sym,ox,ratio in zip(i[0],i[1],i[2]):
            comp[Specie(sym,ox)] = ratio
        comp_list = [[key]*val for key,val in comp.items()]
        comp_list = [item for sublist in comp_list for item in sublist]
        species_comps.append(comp_list)

    # Sort and ditch duplicates
    print('Ditching duplicates (sorry to have got your hopes up with the big numbers)...')
    for i in species_comps:
        i.sort()
        i.sort(key=lambda x: x.oxi_state, reverse=True)
    species_comps = list(set([tuple(i) for i in species_comps]))
    print('Total number of new compounds unique compositions: {0}'.format(len(species_comps)))
    return species_comps

def predict_structure(species, struc_list, check_dir=False, threshold = 0.00001):
    """ Predicted structures for set of pymatgen species using the Pymatgen structure predictor
    and save as cif file.
    TODO: This will be superceded by our own implementation of the structure prediction algorithm
    in future versions of SMACT.
    Args:
        species (list): Pymatgen Species for which structure should be predicted.
        struc_list (list): Pymatgen Structure objects to consider as parent structures in the
        substitution algorithm.
        check_dir (bool): check if directory already exists and only carry out
        prediction and write new files if it doesn't.
        threshold (float): Log-probability threshold for the Pymatgen structure predictor.
    Returns:
        Saves cif files of possible structures in new directory along with a summary .txt file
        containing info including probabilities.
    """
    sub = Substitutor(threshold = 0.00001)
    print('{}  ........'.format(species))
    dirname = ''.join([str(i) for i in species])
    path_exists = True if os.path.exists('./SP_results/{0}'.format(dirname)) else False

    if (check_dir and path_exists):
        print('Already exists')
    else:
        print('{} not already there'.format(dirname))
        suggested_strucs = sub.pred_from_structures(target_species=species, structures_list=struc_list,
                                                 remove_existing = False, remove_duplicates = True)
        suggested_strucs = sorted(suggested_strucs,
        key=lambda k: k.other_parameters['proba'], reverse = True)

        # Save the structures as cifs
        if not path_exists:
            os.makedirs('./SP_results/{0}'.format(dirname))

        for i, d in enumerate(suggested_strucs):
            cw = CifWriter(d.final_structure)
            cw.write_file("SP_results/{0}/{1}_BasedOn_{2}.cif".format(dirname, i, d.history[0]['source']))

        # Save the summary as a text file
        with open ('SP_results/{0}/{0}_summary.txt'.format(dirname), 'w') as f:
            f.write('Formula,  Probability,  Based on,  Substitutions \n')
            for i, struc in enumerate(suggested_strucs):
                f.write(' {0}, {1:12},   {2:.5f},         {3:9},    {4} \n'.format(i,struc.composition.reduced_formula,
                                                                               struc.other_parameters['proba'],
                                                        struc.history[0]['source'],
                                                        re.sub('Specie', '', str(struc.history[1]['species_map']))))
    print('Done.')

def add_probabilities(strucs):
    """ Add probabilities from summary text files for a list of dicts containing structures.
        Dicts must contain Structure, based_on (str).
        Args:
            strucs (list): Dicts containing Pymatgen Structures (keyed by 'struc')
            and a string of the parent structure formula (keyed by 'based_on').
        Returns:
            strucs (list): As supplied to function but with the additional probability info.
    """
    for i in tqdm(strucs):
        ions = ''.join([str(j) for j in i['struc'].composition])
        with open('SP_results/{}/{}_summary.txt'.format(ions, ions), 'r') as f:
            for lines in f:
                line = lines.split(',')
                if line[3].strip() == i['based_on']:
                    i['probability'] = float(line[2].strip())
    return strucs
