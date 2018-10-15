###############################################################################
# Copyright Daniel Davies (2018)                                              #
#                                                                             #
# This file is part of SMACT: oxidationstates.py is free software: you can    #
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
# Collection of functions for the statistical analysis of oxidation states    #
# tidied neatly into one module.                                              #
# NB:                                                                         #
# This module currently depends on Pymatgen, although this is not currently a #
# dependency of SMACT. See http://pymatgen.org/.                              #
###############################################################################

###  Imports
# General
import os, re, json
from tqdm import tqdm
from collections import Counter
import itertools

# Maths and plotting
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# pymatgen
from pymatgen import Structure, Specie, MPRester
from pymatgen.analysis.structure_prediction.substitutor import Substitutor
from pymatgen.io.cif import CifWriter
from smact import ordered_elements, Element, neutral_ratios
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


def get_unique_species(structures):
    """Given a set of pymatgen structures, in the form of dictionaries where the
    Structure is keyed as 'structure', returns a list of all the different
    Species present in that set.
    Args:
        structures (list): Dictionaries containing pymatgen Structures
    Returns:
        species_list (list): Unique species that are exhibited in the structures.
    """
    species_list = []
    for i in structures:
        for sp in i['structure'].composition:
            species_list.append((sp))
    species_list = list(set(species_list))
    return(species_list)


def sort_species(species_list, ordering='ptable', reverse=False):
    """Given a list of pymatgen Species, will order them according to a given
    rule and return the ordered list.
    Args:
        species_list (list): Pymatgen species objects
        ordering('string'): How to order the Species:
            ptable: order by periodic table position
        reverse (bool): Whether to reverse the ordering (descending order)
    Returns:
        species_list (list): ordered Pymatgen species objects
    """
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
    print("First species: {0}  last species: {1}".format(species_list[0], species_list[-1]))
    return species_list


def species_totals(structures, count_elements=False):
    """Given a set of pymatgen structures in the form of dictionaries where
    the Structure is keyed as 'structure', returns the number
    of compounds that features each Species.
    Args:
        structures (list): dictionaries containing pymatgen Structures
        count_elements (bool): switch to counting elements not species
    Returns:
        totals (dict): Totals of each species.
    """
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

    return(totals)

def species_totals_for_anion(anion, structures, edit_structure_dicts=True):
    """Finds the number of instances of each species in a list of pymatgen
    Structure objects when a given anion is the most electronegative one
    present. Also adds most electronegative anion to the dictionary.
    Args:
        anion (Pymatgen.Species): Anion of interest
        structures (list): Dictionaries containing pymatgen Structures
        edit_structure_dicts (bool): Modify the dicts in the structures list
        to include a 'most_eneg_anion' key.
    Returns:
        an_containing (dict): Totals of each species.
    """
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
                if edit_structure_dicts:
                    i['most_eneg_anion'] = anion

    an_containing = [i for sublist in an_containing for i in sublist]
    an_containing = dict(Counter(an_containing))
    an_containing.pop(anion)
    return(an_containing)

def generate_scores(spec_list, anion_count_dict, spec_totals, scoring, el_totals=None):
    """Given a list of Species of interest, and a dictionary of the
    occurrence of each species with each anion as the most electronegative one
    present, can return various 'scores' based on this info.
    Args:
        spec_list (list): Pymatgen species we are interested in
            e.g. [Specie Li+, Specie Na+, Specie Mg2+...].
        anion_count_dict (dict): Occurrence of each species with each anion as the most
        electronegiative element present (as generated by species_totals_for_anion)
            e.g. {Specie:F-: {Specie Rb+: 164, Specie K+...}...}.
        spec_totals (dict): Total number of compounds containing each Species (as
        generated by species_totals)
            e.g. {Specie Cs+: 1026, Specie Ag+: 521...}.
        scoring (str): Can take values:
            spec_fraction - Totals in dictionary / total compounds containing
            that species.
            an_norm - As in spec_frac but also divided by total compounds
            containing that anion.
            overall_score - As in an_norm but multiplied by
            (compounds containing species / compounds containing element)
            spec_prob - Calculate the probability of the compound based on
            available oxidation states.
        el_totals (dict): Total number of compounds containing each element.
    Returns:
        list_vals (dict) of list: Scores for each species, dict is keyed by anion.
        WARNING: lists of scores are in the same order as the supplied spec_list.
     """

    if scoring == 'spec_fraction':
        list_vals = {}
        for an, subdict in anion_count_dict.items():
            an_vals = []
            for i in spec_list:
                if i in anion_count_dict[an]:
                    an_vals.append(anion_count_dict[an][i]/spec_totals[i])
                else:
                    an_vals.append(0.)
            list_vals[an] = an_vals

    elif scoring == 'an_norm':
        list_vals = {}
        for an, subdict in anion_count_dict.items():
            an_vals = []
            for i in spec_list:
                if i in anion_count_dict[an]:
                    an_vals.append(anion_count_dict[an][i]/(spec_totals[i]*\
                    spec_totals[an])*1000)
                else:
                    an_vals.append(0.)
            list_vals[an] = an_vals

    elif scoring == 'overall_score':
        list_vals = {}
        for an, subdict in anion_count_dict.items():
            an_vals = []
            for i in spec_list:
                if i in anion_count_dict[an]:
                    an_vals.append(anion_count_dict[an][i]/(el_totals[i.symbol]*\
                    spec_totals[an])*2000)
                else:
                    an_vals.append(0.)
            list_vals[an] = an_vals

    elif scoring == 'spec_prob':
        list_vals = {}
        for an, subdict in anion_count_dict.items():
            an_vals = []
            for i in spec_list:
                if i in anion_count_dict[an]:
                    element_an_counts = [count for sp,count in subdict.items() if i.symbol==sp.symbol]
                    element_an_tot = sum(element_an_counts)
                    an_vals.append(anion_count_dict[an][i]/(element_an_tot))
                else:
                    an_vals.append(0.)
            list_vals[an] = an_vals
    return list_vals

    return(list_vals)


def plot_all_scores(list_scores, spec_list, break_points, figure_filename, raw_totals=None):
    """ Plots correlation plots for all the species considered.
    Args:
        list_scores (dict): Lists of scores for the species in spec_list
        keyed by anion (as generated by generate_scores).
        spec_list (list): Pymatgen species in same order as corresponding
        scores in each list within the dict list_scores (as used by generate_scores).
        break_points (list): ints, positions in lists to start each new plot.
        figure_filename (str): Name of .png file to save plot.
        raw_totals (list) of dicts: Two dicts: First dict is element totals, second
        dict is species totals. If included, will display these values on the plot.i
    Returns:
        Saves and shows plots. Data should be chunked into 14 subplots (3 x 4 + 2) and this
        is not currently flexible.
         """
    arr = np.array([list_scores[key] for key in list_scores])
    arr = arr.transpose()

    # y-axis
    y_labels = []
    for i in spec_list:
        if raw_totals:
            spec_tot = raw_totals[1][i]
            el_tot = raw_totals[0][i.symbol]
            y_lab = r"$_{{{0}/{1}}}$ {2}$^{{+{3}}}$".format(spec_tot,el_tot,i.symbol,int(i.oxi_state))
        else:
            y_lab = r"{0}$^{{+{1}}}$".format(i.symbol,int(i.oxi_state))
        y_labels.append(y_lab)
    # x-axis clean up
    x_labels = [key.symbol for key, val in list_scores.items()]

    # Chunk up data
    chunks = {}
    spec_labels = {}
    for i in range(0, len(break_points)-1):
        chunk = arr[break_points[i]: break_points[i+1]]
        label_chunk = y_labels[break_points[i]: break_points[i+1]]
        chunks['{}'.format(i)] = chunk
        spec_labels['{}'.format(i)] = label_chunk

    # Plotting
    f, axarr = plt.subplots(5, 3)
    plot = plt.setp(axarr, xticks=np.arange(0,8), xticklabels=x_labels)

    axarrs = [axarr[0,0],axarr[0,1], axarr[0,2],
          axarr[1,0], axarr[1,1],axarr[1,2],
          axarr[2,0],axarr[2,1], axarr[2,2],
          axarr[3,0],axarr[3,1], axarr[3,2],
          axarr[4,0],axarr[4,1]   ]

    '''axarrs = [ axarr[0,0], axarr[0,1],
            axarr[1,0], axarr[1,1],
            axarr[2,0], axarr[2,1],
            axarr[3,0], axarr[3,1],
            axarr[4,0], axarr[4,1],
            axarr[5,0], axarr[5,1],
            axarr[6,0], axarr[6,1]
            ]
'''
    for ax, i in zip(axarrs, list(chunks.keys())):
        ax1 =ax.imshow(chunks[i], interpolation='None',
        cmap = mpl.cm.get_cmap('inferno_r'), vmin = 0, vmax =1)
        plt.sca(ax)
        plt.yticks([i for i,n in enumerate(chunks[i])], spec_labels[i])

    # Tidy up
    plt.tight_layout()
    f.subplots_adjust(right = 0.85)
    cbar_ax = f.add_axes([0.9,0.05,0.04,0.3])
    f.colorbar(ax1, cax=cbar_ax)
    plt.savefig(figure_filename, dpi = 300)
    plt.show()


def plot_metal(metal, list_scores, spec_list, show_legend = False):
    """ Plot distribution of species for individual metals.
    Args:
        metal (str): metal element of interest to plot.
        list_scores (dict): Lists of scores for the species in spec_list
        keyed by anion (as generated by genrate_scores).
        spec_list (list): Pymatgen species in same order as corresponding
        scores in each list within the dict list_scores.
        show_legend (bool): display legend on plot.
    Returns:
        Shows and saves plot.
     """
    # Set initial very daft x and y range values to be adjusted below
    min_x, max_x = 20, -20
    max_y = 0

    overall_list = []
    for anion in list_scores.keys():
        an_dict = {}
        for spec,score in zip(spec_list,list_scores[anion]):
            if spec.symbol == metal:
                an_dict[spec.oxi_state] = score
        sorted_list = sorted(an_dict.items())
        if an_dict:
            x,y = zip(*sorted_list)
            x = list(x)
            y = list(y)
            overall_list.append([x,y,anion])
            if min(x) < min_x:
                min_x = min(x)
            if max(x) > max_x:
                max_x = max(x)
            if max(y) > max_y:
                max_y = max(y)
        else:
            overall_list.append([[1],[0],anion])

    # Plotting
    labels = ["$F^-$","$O^{2-}$","$Cl^-$","$Br^-$","$I^-$","$S^{2-}$","$Se^{2-}$","$Te^{2-}$"]

    # Aesthetics and bar positions
    width = 1./10.
    pos = 0
    colours = ['#E51200', '#DF5400', '#DA9300', '#D4CE00',
               '#97CF00', '#57C900', '#1BC400', '#00BF1D']

    for col,anion in enumerate(overall_list):
        pos += width
        plt.bar(np.array(anion[0]) + pos, anion[1], width, label=labels[col], color = colours[col])

    if show_legend:
        plt.legend(prop={'size':24})

    plt.xticks(np.arange(min_x, max_x+1)+0.5, np.arange(min_x,max_x+1, dtype=np.int_))
    if min_x < 0:
        min_x = 0
    plt.xlim(min_x, max_x+1)
    plt.ylim(0, 1.19)
    plt.xlabel('Oxidation state')
    plt.ylabel('Species fraction')

    at_no = int(Element(metal).number)
    mass = Element(metal).mass
    plt.text(np.mean([max_x,min_x]) + 0.5, 1.0, "$^{{{0}}}$\n  {1}\n$_{{{2:.2f}}}$".format(at_no,metal, mass),
            bbox={'facecolor':'gray', 'alpha':0.3, 'pad':20}, size=28)
    plt.tight_layout()
    plt.savefig('OxidationState_score_{0}'.format(metal, dpi=300))
    plt.show()

def assign_prob(structures, list_scores, species_list, scoring = 'overall_score', verbose=False,
        edit_struc_dict = True):
    """ Assigns probability values to structures based on the list of score values.
    Args:
        structures (list): Dictionaries containing pymatgen Structures, keyed under 'structure'.
        list_scores (dict): Lists of scores for the species in spec_list keyed by anion
        (as produced by generate_scores).
        species_list (list): Pymatgen species in same order as corresponding lists in list_scores.
        scoring (str): Can be either:
                        overall_score - Mean species-anion score for each species of interest
                        in the composition.
                        limiting_score - As above but minimum species-anion score.
                        probability - Product of scores.
                        probability_simple - Product of scores for different species only (set(comp))
        verbose (bool): Explicitly print any compounds that were skipped over due to the elements
        they contain.
        edit_struc_dict (bool): Add the probability to the dicts in the structures list.
    Returns:
        probabilities_list (list): Score for each structure in structures.
    """
    scores_dict = {}
    for key in list_scores.keys():
        an = {}
        for spec, val in zip(species_list, list_scores[key]):
            an[spec] = val
        scores_dict[key] = an

    probabilities_list = []
    for struc in tqdm(structures):
        scores = []
        comp = list(struc['structure'].species)
        try:
            scores = [scores_dict[struc['most_eneg_anion']][sp] for sp in comp \
            if sp in species_list]
            if scoring == 'overall_score':
                overall_score = np.mean(scores)
            elif scoring == 'limiting_score':
                overall_score = min(scores)
            elif scoring == 'probability':
                overall_score = np.prod(scores)
            elif scoring == 'probability_simple':
                scores = [scores_dict[struc['most_eneg_anion']][sp] for sp in list(set(comp)) \
                if sp in species_list]
                overall_score = np.prod(scores)
        except:
            if verbose:
                print('Could not get score for: {}'.format(comp))
            overall_score = 0
        if edit_struc_dict = True:
            struc['probability'] = overall_score
        probabilities_list.append(overall_score)

    return probabilities_list

def assign_prob_new(compositions, list_scores, species_list, scoring, verbose=False):
    """ Assign probabilities to novel compositions based on the list of score values.
    Args:
        compositions (list): list of pymatgen species (possibly as generated by smact)
        list_scores (dict): Lists of scores for the species in spec_list
        spec_list (list): Pymatgen species in same order as corresponding
        verbose (bool): Explicitly print any compounds that were skipped over
        scoring (str): Can be either:
                        overall_score: mean species-anion score for each species of interest
                        in the composition
                        limiting_score: as above but minimum species-anion score
                        probability: product of scores
                        probability_simple: product of scores for different species only (set(comp))
        """
    scores_dict = {}
    for key in list_scores.keys():
        an = {}
        for spec, val in zip(species_list, list_scores[key]):
            an[spec] = val
        scores_dict[key] = an

    probabilities_list = []
    for comp in tqdm(compositions):
        scores = []
        # Assumes the last species is the most electronegative
        most_eneg = comp[-1]
        try:
            scores = [scores_dict[most_eneg][sp] for sp in comp \
            if sp in species_list]
            if scoring == 'overall_score':
                overall_score = np.mean(scores)
            elif scoring == 'limiting_score':
                overall_score = min(scores)
            elif scoring == 'probability':
                overall_score = np.prod(scores)
            elif scoring == 'probability_simple':
                scores = [scores_dict[most_eneg][sp] for sp in list(set(comp)) \
                if sp in species_list]
                overall_score = np.prod(scores)
        except:
            if verbose:
                print('Could not get score for: {}'.format(comp))
            overall_score = 0
        probabilities_list.append(overall_score)
    return probabilities_list

def plot_scores_hist(scores, bins = 100, plot_title='plot', xlim = None, ylim = None):
    """ Plot histogram of a list of scores.        """
    print('Number of individual values: {}'.format(len(scores)))
    print('Number of zero values: {}'.format(scores.count(0.0)))
    print('Maximum score: {}'.format(max(scores)))

    plt.hist(scores, bins, color='orange', alpha=0.5)
    if xlim:
        plt.xlim(xlim[0],xlim[1])
    plt.title(plot_title)
    plt.xlabel('Score')
    if ylim:
        plt.ylim(ylim[0],ylim[1])
    plt.tight_layout()
    plt.savefig('{}.png'.format(plot_title), dpi=300)
    plt.show()

def ternary_smact_combos(position1, position2, position3, threshold = 8):
    """ Get Pymatgen species compositions using SMACT when up to three different
        lists are needed to draw species from. E.g. Ternary metal halides...
    Args:
        position(n) (list of species): Species to be considered iteratively for each
                                     position.
        threshold (int): Max stoichiometry threshold
        """

    initial_comps_list = []
    for sp1, sp2, an in tqdm(itertools.product(position1, position2, position3)):
        m1, oxst1 = sp1.symbol, int(sp1.oxi_state)
        eneg1 = Element(m1).pauling_eneg
        m2, oxst2 = sp2.symbol, int(sp2.oxi_state)
        eneg2 = Element(m2).pauling_eneg
        hal, oxst3 = an.symbol, int(an.oxi_state)
        eneg3 = Element(hal).pauling_eneg

        symbols = [m1,m2,hal]
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

def predict_structure(species, struc_list, check_dir=False):
    """ Save cif files of predicted structures for set of pymatgen species.
    Args:
        species (list): pymatgen species to predict structures for
        struc_list (list): pymatgen structures to substitute species into
        check_dir (bool): check if directory already exists and only carry out
        prediction if it doesn't.
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
    """
    for i in tqdm(strucs):
        ions = ''.join([str(j) for j in i['struc'].composition])
        with open('SP_results/{}/{}_summary.txt'.format(ions, ions), 'r') as f:
            for lines in f:
                line = lines.split(',')
                if line[3].strip() == i['based_on']:
                    i['probability'] = float(line[2].strip())
    return strucs


oxistate_prob_table = {('F', -1): [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0392156862745098, 0.35294117647058826, 0.6078431372549019, 0.21739130434782608, 0.45652173913043476, 0.17391304347826086, 0.15217391304347827, 0.13157894736842105, 0.5789473684210527, 0.19736842105263158, 0.09210526315789473, 0.0, 0.0, 0.47959183673469385, 0.41836734693877553, 0.10204081632653061, 0.0, 0.0, 0.0, 0.0, 0.23333333333333334, 0.7666666666666667, 0.0, 0.0, 0.6428571428571429, 0.30952380952380953, 0.047619047619047616, 0.0, 0.7288135593220338, 0.13559322033898305, 0.13559322033898305, 0.08080808080808081, 0.8686868686868687, 0.050505050505050504, 1.0, 0.0, 0.0, 1.0, 0.043478260869565216, 0.0, 0.9565217391304348, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.08333333333333333, 0.19444444444444445, 0.7222222222222222, 0.0, 0.2972972972972973, 0.08108108108108109, 0.02702702702702703, 0.5945945945945946, 0.13333333333333333, 0.13333333333333333, 0.13333333333333333, 0.4666666666666667, 0.13333333333333333, 0.05, 0.4, 0.55, 0.7142857142857143, 0.047619047619047616, 0.23809523809523808, 0.4, 0.4666666666666667, 0.13333333333333333, 1.0, 0.045454545454545456, 0.0, 0.9545454545454546, 0.4430379746835443, 0.02531645569620253, 0.5316455696202531, 0.35036496350364965, 0.021897810218978103, 0.6131386861313869, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.9375, 0.0625, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.47368421052631576, 0.5263157894736842, 0.0, 1.0, 0.0, 0.0, 0.3684210526315789, 0.631578947368421, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.26666666666666666, 0.7333333333333333, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.125, 0.0, 0.0, 0.875, 0.0, 0.0, 0.08333333333333333, 0.041666666666666664, 0.25, 0.625, 0.2, 0.5, 0.1, 0.2, 0.25, 0.75, 0.8444444444444444, 0.15555555555555556, 0.8888888888888888, 0.1111111111111111, 0.0, 0.0, 0.7142857142857143, 0.2857142857142857, 0.0, 1.0, 0.0, 0.034482758620689655, 0.41379310344827586, 0.1724137931034483, 0.3793103448275862], ('O', -2): [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.014814814814814815, 0.9851851851851852, 0.0025188916876574307, 0.0982367758186398, 0.8992443324937027, 0.008368200836820083, 0.1297071129707113, 0.26778242677824265, 0.5941422594142259, 0.0960960960960961, 0.40540540540540543, 0.04804804804804805, 0.10510510510510511, 0.34534534534534533, 0.0, 0.5506216696269982, 0.22380106571936056, 0.19005328596802842, 0.019538188277087035, 0.007104795737122558, 0.008880994671403197, 0.0019193857965451055, 0.2399232245681382, 0.7428023032629558, 0.015355086372360844, 0.0117096018735363, 0.7634660421545667, 0.16627634660421545, 0.0585480093676815, 0.0033003300330033004, 0.8679867986798679, 0.10231023102310231, 0.026402640264026403, 0.145748987854251, 0.8205128205128205, 0.033738191632928474, 1.0, 0.004464285714285714, 0.0, 0.9955357142857143, 0.005763688760806916, 0.002881844380403458, 0.9769452449567724, 1.0, 1.0, 0.0, 0.004149377593360996, 0.995850622406639, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0029154518950437317, 0.043731778425655975, 0.061224489795918366, 0.892128279883382, 0.0, 0.02872531418312388, 0.059245960502693, 0.10771992818671454, 0.8043087971274686, 0.029940119760479042, 0.0718562874251497, 0.3413173652694611, 0.49700598802395207, 0.059880239520958084, 0.0, 0.6949152542372882, 0.3050847457627119, 0.881578947368421, 0.07894736842105263, 0.039473684210526314, 0.9866666666666667, 0.0044444444444444444, 0.008888888888888889, 1.0, 0.0, 0.017857142857142856, 0.9821428571428571, 0.2875, 0.00625, 0.70625, 0.321285140562249, 0.0, 0.5863453815261044, 1.0, 1.0, 0.002680965147453083, 0.02680965147453083, 0.9705093833780161, 0.00819672131147541, 0.7049180327868853, 0.28688524590163933, 0.028409090909090908, 0.9715909090909091, 0.02158273381294964, 0.9784172661870504, 0.011560693641618497, 0.9884393063583815, 0.25, 0.75, 0.006329113924050633, 0.9936708860759493, 0.0, 0.0, 0.8709677419354839, 0.12903225806451613, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0136986301369863, 0.9863013698630136, 0.08333333333333333, 0.9166666666666666, 1.0, 0.0, 1.0, 0.0, 0.0, 0.01282051282051282, 0.02564102564102564, 0.9615384615384616, 0.0, 0.0, 0.00980392156862745, 0.029411764705882353, 0.9607843137254902, 0.0, 0.0, 0.05357142857142857, 0.26785714285714285, 0.19047619047619047, 0.4880952380952381, 0.046511627906976744, 0.4883720930232558, 0.3953488372093023, 0.06976744186046512, 0.23417721518987342, 0.7658227848101266, 0.7772020725388601, 0.22279792746113988, 0.8956043956043956, 0.1043956043956044, 0.002702702702702703, 0.0, 0.9297297297297298, 0.05945945945945946, 0.025, 0.975, 0.0, 0.010344827586206896, 0.07586206896551724, 0.07241379310344828, 0.8413793103448276], ('Cl', -1): [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.045454545454545456, 0.3181818181818182, 0.6363636363636364, 0.2727272727272727, 0.45454545454545453, 0.2727272727272727, 0.5, 0.2857142857142857, 0.14285714285714285, 0.07142857142857142, 0.5909090909090909, 0.4090909090909091, 0.0, 0.0, 0.0, 0.0, 0.9285714285714286, 0.0, 0.07142857142857142, 0.0, 0.0, 0.0, 0.0, 0.4375, 0.5625, 0.0, 0.0, 0.8148148148148148, 0.18518518518518517, 0.0, 0.0, 1.0, 0.0, 0.0, 0.4418604651162791, 0.5116279069767442, 0.046511627906976744, 1.0, 0.043478260869565216, 0.08695652173913043, 0.8695652173913043, 0.3333333333333333, 0.0, 0.5, 1.0, 1.0, 0.1111111111111111, 0.0, 0.8888888888888888, 0.03333333333333333, 0.23333333333333334, 0.3, 0.43333333333333335, 0.0, 0.0, 0.3333333333333333, 0.16666666666666666, 0.5, 0.21621621621621623, 0.1891891891891892, 0.21621621621621623, 0.2702702702702703, 0.10810810810810811, 0.1, 0.6, 0.3, 0.0, 0.0, 0.0, 0.75, 0.25, 0.8387096774193549, 0.0, 0.16129032258064516, 0.9473684210526315, 0.05263157894736842, 0.0, 1.0, 0.4166666666666667, 0.08333333333333333, 0.5, 0.4166666666666667, 0.0, 0.5833333333333334, 0.3181818181818182, 0.0, 0.5681818181818182, 1.0, 1.0, 0.0, 0.05, 0.95, 0.0, 0.875, 0.125, 0.0, 1.0, 0.0, 1.0, 0.3333333333333333, 0.6666666666666666, 0.4444444444444444, 0.5555555555555556, 0.16666666666666666, 0.8333333333333334, 0.2, 0.0, 0.8, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.5, 0.5, 1.0, 0.2857142857142857, 0.7142857142857143, 0.0, 0.0, 0.08695652173913043, 0.2608695652173913, 0.6521739130434783, 0.1111111111111111, 0.2222222222222222, 0.2777777777777778, 0.16666666666666666, 0.2222222222222222, 0.03225806451612903, 0.3225806451612903, 0.3225806451612903, 0.22580645161290322, 0.03225806451612903, 0.06451612903225806, 0.6666666666666666, 0.16666666666666666, 0.16666666666666666, 0.0, 0.020833333333333332, 0.9791666666666666, 0.875, 0.125, 0.8125, 0.1875, 0.0, 0.04, 0.96, 0.0, 0.0, 1.0, 0.0, 0.26666666666666666, 0.4, 0.26666666666666666, 0.06666666666666667], ('Br', -1): [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.5, 0.5, 0.5, 0.3333333333333333, 0.16666666666666666, 0.8333333333333334, 0.16666666666666666, 0.0, 0.0, 0.42857142857142855, 0.5714285714285714, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4444444444444444, 0.5555555555555556, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.7391304347826086, 0.17391304347826086, 0.08695652173913043, 1.0, 0.0, 0.1, 0.9, 0.5555555555555556, 0.0, 0.3333333333333333, 1.0, 1.0, 0.0, 0.0, 1.0, 0.1111111111111111, 0.2222222222222222, 0.3333333333333333, 0.3333333333333333, 0.0, 0.058823529411764705, 0.35294117647058826, 0.11764705882352941, 0.47058823529411764, 0.5, 0.4166666666666667, 0.08333333333333333, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5384615384615384, 0.15384615384615385, 0.3076923076923077, 1.0, 0.0, 0.0, 1.0, 0.43478260869565216, 0.2608695652173913, 0.30434782608695654, 0.625, 0.0, 0.375, 0.7777777777777778, 0.0, 0.0, 1.0, 1.0, 0.05, 0.2, 0.75, 0.0, 1.0, 0.0, 0.125, 0.875, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.2, 0.8, 0.14285714285714285, 0.14285714285714285, 0.7142857142857143, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.8, 0.2, 1.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.6, 0.4166666666666667, 0.08333333333333333, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.0, 0.5454545454545454, 0.09090909090909091, 0.18181818181818182, 0.09090909090909091, 0.09090909090909091, 1.0, 0.0, 0.0, 0.0, 0.023809523809523808, 0.9761904761904762, 0.9047619047619048, 0.09523809523809523, 1.0, 0.0, 0.1, 0.1, 0.75, 0.0, 0.14285714285714285, 0.8571428571428571, 0.0, 0.25, 0.625, 0.125, 0.0], ('I', -1): [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.5, 0.5, 0.5, 0.3333333333333333, 0.16666666666666666, 1.0, 0.0, 0.0, 0.0, 0.7777777777777778, 0.2222222222222222, 0.0, 0.0, 0.0, 0.14285714285714285, 0.8571428571428571, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.2, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.9473684210526315, 0.0, 0.05263157894736842, 1.0, 0.0, 0.2, 0.8, 0.5, 0.0, 0.16666666666666666, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.3076923076923077, 0.5384615384615384, 0.15384615384615385, 0.0, 0.1111111111111111, 0.2777777777777778, 0.1111111111111111, 0.5, 0.5, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3333333333333333, 0.2222222222222222, 0.4444444444444444, 1.0, 0.0, 0.0, 1.0, 0.5, 0.05555555555555555, 0.4444444444444444, 0.8, 0.0, 0.2, 0.7647058823529411, 0.0, 0.0, 1.0, 1.0, 0.043478260869565216, 0.34782608695652173, 0.6086956521739131, 0.5, 0.5, 0.0, 0.16666666666666666, 0.8333333333333334, 0.6666666666666666, 0.3333333333333333, 0.3333333333333333, 0.6666666666666666, 1.0, 0.0, 0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.3333333333333333, 0.6666666666666666, 0.0, 1.0, 1.0, 0.6666666666666666, 0.3333333333333333, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.2222222222222222, 0.1111111111111111, 0.6666666666666666, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.02857142857142857, 0.9714285714285714, 0.9642857142857143, 0.03571428571428571, 1.0, 0.0, 0.125, 0.0625, 0.8125, 0.0, 0.0, 1.0, 0.0, 0.2857142857142857, 0.7142857142857143, 0.0, 0.0], ('S', -2): [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.05263157894736842, 0.9473684210526315, 0.03125, 0.28125, 0.6875, 0.10204081632653061, 0.32653061224489793, 0.3673469387755102, 0.20408163265306123, 0.15789473684210525, 0.8070175438596491, 0.03508771929824561, 0.0, 0.0, 0.0, 0.9423076923076923, 0.038461538461538464, 0.019230769230769232, 0.0, 0.0, 0.0, 0.0, 0.47058823529411764, 0.5098039215686274, 0.0196078431372549, 0.0, 0.6333333333333333, 0.26666666666666666, 0.1, 0.0, 0.7407407407407407, 0.25925925925925924, 0.0, 0.8060344827586207, 0.1939655172413793, 0.0, 1.0, 0.0, 0.05555555555555555, 0.9444444444444444, 0.030927835051546393, 0.041237113402061855, 0.9175257731958762, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.125, 0.041666666666666664, 0.8333333333333334, 0.0, 0.10638297872340426, 0.1702127659574468, 0.3829787234042553, 0.3404255319148936, 0.10526315789473684, 0.42105263157894735, 0.3157894736842105, 0.0, 0.15789473684210525, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.8, 0.2, 1.0, 0.0, 0.0, 0.9919354838709677, 0.008064516129032258, 0.0, 1.0, 0.02666666666666667, 0.04, 0.9333333333333333, 0.3229166666666667, 0.041666666666666664, 0.6354166666666666, 0.8271604938271605, 0.0, 0.12345679012345678, 1.0, 1.0, 0.0, 0.04081632653061224, 0.9591836734693877, 0.023809523809523808, 0.9285714285714286, 0.047619047619047616, 0.041666666666666664, 0.9583333333333334, 0.058823529411764705, 0.9411764705882353, 0.029411764705882353, 0.9705882352941176, 0.9166666666666666, 0.08333333333333333, 0.05, 0.95, 0.0, 0.05555555555555555, 0.9444444444444444, 0.0, 0.045454545454545456, 0.9545454545454546, 0.05263157894736842, 0.9473684210526315, 1.0, 0.05555555555555555, 0.9444444444444444, 0.3333333333333333, 0.6666666666666666, 1.0, 0.07692307692307693, 0.9230769230769231, 0.022222222222222223, 0.08888888888888889, 0.022222222222222223, 0.3333333333333333, 0.5333333333333333, 0.0, 0.0, 0.2857142857142857, 0.0, 0.7142857142857143, 0.0, 0.6, 0.4, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 1.0, 0.9887640449438202, 0.011235955056179775, 1.0, 0.0, 0.0, 0.04918032786885246, 0.9508196721311475, 0.0, 0.08333333333333333, 0.9166666666666666, 0.027777777777777776, 0.1388888888888889, 0.7777777777777778, 0.027777777777777776, 0.027777777777777776], ('Se', -2): [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.14285714285714285, 0.8571428571428571, 0.13333333333333333, 0.4, 0.4666666666666667, 0.1111111111111111, 0.5, 0.2222222222222222, 0.16666666666666666, 0.07692307692307693, 0.9230769230769231, 0.0, 0.0, 0.0, 0.0, 0.9285714285714286, 0.07142857142857142, 0.0, 0.0, 0.0, 0.0, 0.0, 0.47058823529411764, 0.5294117647058824, 0.0, 0.0, 0.631578947368421, 0.2631578947368421, 0.10526315789473684, 0.05263157894736842, 0.7894736842105263, 0.15789473684210525, 0.0, 0.8391608391608392, 0.16083916083916083, 0.0, 1.0, 0.0, 0.1111111111111111, 0.8888888888888888, 0.07017543859649122, 0.08771929824561403, 0.8421052631578947, 1.0, 1.0, 0.0, 0.0, 1.0, 0.08333333333333333, 0.0, 0.0, 0.9166666666666666, 0.034482758620689655, 0.0, 0.20689655172413793, 0.3103448275862069, 0.4482758620689655, 0.3, 0.2, 0.4, 0.0, 0.1, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.6666666666666666, 0.3333333333333333, 1.0, 0.0, 0.0, 0.9850746268656716, 0.014925373134328358, 0.0, 1.0, 0.022727272727272728, 0.06818181818181818, 0.9090909090909091, 0.25862068965517243, 0.017241379310344827, 0.7241379310344828, 0.7586206896551724, 0.0, 0.13793103448275862, 1.0, 1.0, 0.0, 0.05555555555555555, 0.9444444444444444, 0.043478260869565216, 0.9130434782608695, 0.043478260869565216, 0.06666666666666667, 0.9333333333333333, 0.06666666666666667, 0.9333333333333333, 0.05, 0.95, 1.0, 0.0, 0.08333333333333333, 0.9166666666666666, 0.0, 0.058823529411764705, 0.9411764705882353, 0.0, 0.06666666666666667, 0.9333333333333333, 0.07692307692307693, 0.9230769230769231, 1.0, 0.1111111111111111, 0.8888888888888888, 0.3333333333333333, 0.6666666666666666, 1.0, 0.0, 1.0, 0.058823529411764705, 0.14705882352941177, 0.029411764705882353, 0.38235294117647056, 0.38235294117647056, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.9787234042553191, 0.02127659574468085, 1.0, 0.0, 0.0, 0.0967741935483871, 0.9032258064516129, 0.0, 0.0, 1.0, 0.05555555555555555, 0.05555555555555555, 0.8888888888888888, 0.0, 0.0], ('Te', -2): [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.125, 0.125, 0.75, 0.25, 0.5, 0.25, 0.5, 0.25, 0.25, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.8823529411764706, 0.11764705882352941, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75, 0.25, 0.0, 0.0, 0.8333333333333334, 0.0, 0.16666666666666666, 0.0, 0.8571428571428571, 0.14285714285714285, 0.0, 0.6764705882352942, 0.3235294117647059, 0.0, 1.0, 0.0, 0.1111111111111111, 0.8888888888888888, 0.7333333333333333, 0.06666666666666667, 0.2, 1.0, 1.0, 0.0, 0.0, 1.0, 0.09090909090909091, 0.36363636363636365, 0.09090909090909091, 0.45454545454545453, 0.0, 0.25, 0.3333333333333333, 0.25, 0.16666666666666666, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.045454545454545456, 0.09090909090909091, 0.8636363636363636, 0.36363636363636365, 0.0, 0.6363636363636364, 0.5263157894736842, 0.0, 0.0, 1.0, 1.0, 0.0, 0.16666666666666666, 0.8333333333333334, 0.2, 0.8, 0.0, 0.2857142857142857, 0.7142857142857143, 0.14285714285714285, 0.8571428571428571, 0.16666666666666666, 0.8333333333333334, 1.0, 0.0, 0.5, 0.5, 0.0, 0.2, 0.8, 0.0, 0.16666666666666666, 0.8333333333333334, 0.16666666666666666, 0.8333333333333334, 1.0, 0.3333333333333333, 0.6666666666666666, 1.0, 0.0, 1.0, 0.0, 1.0, 0.11764705882352941, 0.29411764705882354, 0.17647058823529413, 0.11764705882352941, 0.29411764705882354, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.08333333333333333, 0.9166666666666666, 0.0, 0.0, 1.0, 0.1111111111111111, 0.1111111111111111, 0.7777777777777778, 0.0, 0.0]}
species_list = [('Li', 1.0), ('Be', 2.0), ('Na', 1.0), ('Mg', 2.0), ('Al', 3.0), ('K', 1.0), ('Ca', 2.0), ('Sc', 1.0), ('Sc', 2.0), ('Sc', 3.0), ('Ti', 2.0), ('Ti', 3.0), ('Ti', 4.0), ('V', 2.0), ('V', 3.0), ('V', 4.0), ('V', 5.0), ('Cr', 2.0), ('Cr', 3.0), ('Cr', 4.0), ('Cr', 5.0), ('Cr', 6.0), ('Mn', 1.0), ('Mn', 2.0), ('Mn', 3.0), ('Mn', 4.0), ('Mn', 5.0), ('Mn', 6.0), ('Mn', 7.0), ('Fe', 1.0), ('Fe', 2.0), ('Fe', 3.0), ('Fe', 4.0), ('Co', 1.0), ('Co', 2.0), ('Co', 3.0), ('Co', 4.0), ('Ni', 1.0), ('Ni', 2.0), ('Ni', 3.0), ('Ni', 4.0), ('Cu', 1.0), ('Cu', 2.0), ('Cu', 3.0), ('Zn', 2.0), ('Ga', 1.0), ('Ga', 2.0), ('Ga', 3.0), ('Ge', 2.0), ('Ge', 3.0), ('Ge', 4.0), ('Rb', 1.0), ('Sr', 2.0), ('Y', 1.0), ('Y', 2.0), ('Y', 3.0), ('Zr', 1.0), ('Zr', 2.0), ('Zr', 3.0), ('Zr', 4.0), ('Nb', 1.0), ('Nb', 2.0), ('Nb', 3.0), ('Nb', 4.0), ('Nb', 5.0), ('Mo', 2.0), ('Mo', 3.0), ('Mo', 4.0), ('Mo', 5.0), ('Mo', 6.0), ('Ru', 2.0), ('Ru', 3.0), ('Ru', 4.0), ('Ru', 5.0), ('Ru', 6.0), ('Rh', 1.0), ('Rh', 3.0), ('Rh', 4.0), ('Pd', 2.0), ('Pd', 3.0), ('Pd', 4.0), ('Ag', 1.0), ('Ag', 2.0), ('Ag', 3.0), ('Cd', 2.0), ('In', 1.0), ('In', 2.0), ('In', 3.0), ('Sn', 2.0), ('Sn', 3.0), ('Sn', 4.0), ('Sb', 3.0), ('Sb', 4.0), ('Sb', 5.0), ('Cs', 1.0), ('Ba', 2.0), ('La', 1.0), ('La', 2.0), ('La', 3.0), ('Ce', 2.0), ('Ce', 3.0), ('Ce', 4.0), ('Pr', 2.0), ('Pr', 3.0), ('Nd', 2.0), ('Nd', 3.0), ('Sm', 2.0), ('Sm', 3.0), ('Eu', 2.0), ('Eu', 3.0), ('Gd', 2.0), ('Gd', 3.0), ('Tb', 1.0), ('Tb', 2.0), ('Tb', 3.0), ('Tb', 4.0), ('Dy', 2.0), ('Dy', 3.0), ('Ho', 2.0), ('Ho', 3.0), ('Er', 3.0), ('Tm', 2.0), ('Tm', 3.0), ('Yb', 2.0), ('Yb', 3.0), ('Lu', 3.0), ('Hf', 2.0), ('Hf', 4.0), ('Ta', 1.0), ('Ta', 2.0), ('Ta', 3.0), ('Ta', 4.0), ('Ta', 5.0), ('W', 2.0), ('W', 3.0), ('W', 4.0), ('W', 5.0), ('W', 6.0), ('Re', 2.0), ('Re', 3.0), ('Re', 4.0), ('Re', 5.0), ('Re', 6.0), ('Re', 7.0), ('Ir', 3.0), ('Ir', 4.0), ('Ir', 5.0), ('Ir', 6.0), ('Hg', 1.0), ('Hg', 2.0), ('Tl', 1.0), ('Tl', 3.0), ('Pb', 2.0), ('Pb', 4.0), ('Bi', 1.0), ('Bi', 2.0), ('Bi', 3.0), ('Bi', 5.0), ('Th', 3.0), ('Th', 4.0), ('U', 2.0), ('U', 3.0), ('U', 4.0), ('U', 5.0), ('U', 6.0)]
