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
# Currently this module is quite specific to a current piece of work and      #
# some functions need to be made more general.                                #
# This module also depends on Pymatgen, although this is not currently a      #
# dependency of SMACT. See http://pymatgen.org/.                              #
###############################################################################

###  Imports
import json
from tqdm import tqdm
from collections import Counter
import os
import numpy as np
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
from pymatgen import Structure, Specie, MPRester
from smact import ordered_elements, Element, neutral_ratios
from smact.screening import pauling_test

def get_struc_list(cifpath, json_name):
    """Import pymatgen Structure objects from a json.
    Args:
        cifpath (string): Filepath to json file
        json_name (string): Name of json file
    """
    with open('{0}{1}'.format(cifpath,json_name)) as f:
        saved_strucs = json.load(f)

    struc_list = []
    for i, entry in enumerate(tqdm(saved_strucs)):
        struc_list.append({'structure': Structure.from_dict(entry['structure']),
        'id': entry['id'] })
    return(struc_list)


def mp_filter(criteria, return_elements=False):
    """Get a list of materials project mpids that match a set of given criteria.
    The criteria should be in the format used for a MP query.
    Args:
        criteria (dict): Criteria that can be used in an MPRester query
        return_elements (bool): Also return the elements for that MP entry
    """
    m = MPRester(os.environ.get("MP_API_KEY"))
    if return_elements:
        properties = ['task_id', 'elements']
        struc_filter = m.query(criteria,properties)
        return struc_filter
    else:
        properties = ['task_id']
        struc_filter = m.query(criteria,properties)
        id_list = [i['task_id'] for i in struc_filter]
        return id_list


def species_count(structures):
    """Given a set of pymatgen structures, in the form of dictionaries where the
    Structure is keyed as 'structure', returns a list of all the different
    Species present in that set.
    Args:
        structures (list): Dictionaries containing pymatgen Structures
    """
    species_list = []
    for i in structures:
        for sp in i['structure'].composition:
            species_list.append((sp))
    species_list = list(set(species_list))
    return(species_list)


def sort_species(species_list, ordering):
    """Given a list of pymatgen Species, will order them according to a given
    rule and return the ordered list.
    ordering is a string that can take values: 'ptable': order by periodic table position.
    Args:
        species_list (list): Pymatgen species objects
    """
    ordered_el = ordered_elements(1,103)
    # Turn into tuples for easy sorting
    species_list = [(i.symbol, i.oxi_state) for i in species_list]
    if ordering == 'ptable':
        species_list.sort(key = lambda x: (ordered_el.index(x[0]),x[1]))
        print("Species ordered by periodic table position.")
    else:
        print('Did not reorder the list of species...')

    # Turn back into Species objects
    species_list = [Specie(i[0], i[1]) for i in species_list]
    print("First species: {0}  last species: {1}".format(species_list[0], species_list[-1]))
    return species_list


def species_totals(structures, count_elements=False):
    """Given a set of pymatgen structures in the form of dictionaries where
    the Structure is keyed as 'structure', returns a dictionary of the number
    of compounds that features each Species.
    Args:
        structures (list): dictionaries containing pymatgen Structures
        count_elements (bool): switch to counting elements not species
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

def find_instances(anion, structures):
    """Finds the number of instances of each species in a list of pymatgen
    Structure objects when a given anion is the most electronegative one
    present. Also adds most electronegative anion to the dictionary.
    Args:
        anion (Pymatgen Species): Anion of interest
        structures (list): Dictionaries containing pymatgen Structures
    """
    an_containing = []
    for i in structures:
        if anion in i['structure'].composition:
            # Check whether anion most electronegative element
            an_eneg = Element(anion.symbol).pauling_eneg
            all_enegs = [Element(sp.symbol).pauling_eneg for \
            sp in i['structure'].composition]
            if all(eneg <= an_eneg for eneg in all_enegs):
                comp = [j for j in i['structure'].composition]
                i['most_eneg_anion'] = anion
                an_containing.append(comp)

    an_containing = [i for sublist in an_containing for i in sublist]
    an_containing = dict(Counter(an_containing))
    an_containing.pop(anion)
    return(an_containing)

def generate_scores(spec_list, anion_count_dict, spec_totals, scoring, el_totals=None):
    """Given a list of Species of interest in order and a dictionary of the
    occurrence of each species with each anion as the most electronegative one
    present, can return various 'scores' based on this info.
    Args:
        spec_list (list): Pymatgen species we are interested in (in order)
        anion_count_dict (dict): Totals as described above
        spec_totals (dict): Total number of compounds containing each Species
        scoring (str): Can take values:
            spec_fraction  Totals in dictionary / total compounds containing
            that species.
            an_norm  as in spec_frac but also divided by total compounds
            containing that anion
            overall_score  as in an_norm but multiplied by
            (compounds containing species / compounds containing element)
        el_totals (dict): Total number of compounds containing each element
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
    return(list_vals)


def plot_all_scores(list_scores, spec_list, break_points, figure_filename, raw_totals=None):
    """ Plots correlation plots for all the species considered.
    Args:
        list_scores (dict): Lists of scores for the species in spec_list
        keyed by anion.
        spec_list (list): Pymatgen species in same order as corresponding
        scores in each list within the dict list_scores
        break_points (list): ints, positions in lists to start each new plot
        figure_filename (str): Name of .png file to save plot
        raw_totals (list of dicts): First dict is element totals, second
        dict is species totals. If included, will display these values on the plot.
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
    #HACK f, axarr = plt.subplots(5, 3)
    f, axarr = plt.subplots(7,2)
    plot = plt.setp(axarr, xticks=np.arange(0,8), xticklabels=x_labels)

    #HACK axarrs = [axarr[0,0],axarr[0,1], axarr[0,2],
    #HACK      axarr[1,0], axarr[1,1],axarr[1,2],
    #HACK      axarr[2,0],axarr[2,1], axarr[2,2],
    #HACK      axarr[3,0],axarr[3,1], axarr[3,2],
    #HACK      axarr[4,0]   ]

    axarrs = [ axarr[0,0], axarr[0,1],
            axarr[1,0], axarr[1,1],
            axarr[2,0], axarr[2,1],
            axarr[3,0], axarr[3,1],
            axarr[4,0], axarr[4,1],
            axarr[5,0], axarr[5,1],
            axarr[6,0]
            ]

    for ax, i in zip(axarrs, list(chunks.keys())):
        ax1 =ax.imshow(chunks[i], interpolation='None',
        cmap = mpl.cm.get_cmap('inferno_r'), vmin = 0, vmax =1)
        plt.sca(ax)
        plt.yticks([i for i,n in enumerate(chunks[i])], spec_labels[i])

    # Tidy up
    plt.tight_layout()
    f.subplots_adjust(right = 0.85)
    cbar_ax = f.add_axes([0.9,0.05,0.05,0.2])
    f.colorbar(ax1, cax=cbar_ax)
    plt.savefig(figure_filename, dpi = 300)
    plt.show()


def plot_metal(metal, list_scores, spec_list, show_legend = False):
    """ Plot distribution of species for individual metals.
    Args:
        metal (str): metal element of interest to plot
        list_scores (dict): Lists of scores for the species in spec_list
        keyed by anion.
        spec_list (list): Pymatgen species in same order as corresponding
        scores in each list within the dict list_scores
        show_legend (bool): display legend on plot
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
        plt.legend(prop={'size':18})

    plt.xticks(np.arange(min_x, max_x+1)+0.5, np.arange(min_x,max_x+1))
    if min_x < 0:
        min_x = 0
    plt.xlim(min_x, max_x+1)
    plt.ylim(0, max_y + 0.2)
    plt.xlabel('Oxidation state')
    plt.ylabel('Score')

    at_no = int(Element(metal).number)
    mass = Element(metal).mass
    plt.text(np.mean([max_x,min_x]) + 0.5, max_y, "$^{{{0}}}$\n  {1}\n$_{{{2:.2f}}}$".format(at_no,metal, mass),
            bbox={'facecolor':'gray', 'alpha':0.3, 'pad':10}, size=18)
    plt.tight_layout()
    plt.savefig('OxidationState_score_{0}'.format(metal, dpi=300))
    plt.show()

def assign_prob(structures, list_scores, species_list, scoring = 'overall_score', verbose=False):
    """ Assigns probability values to structures based on the list of score values.
    Args:
        structures (list): Dictionaries containing pymatgen Structures
        list_scores (dict): Lists of scores for the species in spec_list
        spec_list (list): Pymatgen species in same order as corresponding
        scoring (str): Can be either:
                        overall_score: mean species-anion score for each species of interest
                        in the composition
                        limiting_score: as above but minimum species-anion score
        verbose (bool): Explicitly print any compounds that were skipped over

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
        except:
            if verbose:
                print('Could not get score for: {}'.format(comp))
            overall_score = 0
        struc['probability'] = overall_score
        probabilities_list.append(overall_score)

    return probabilities_list

def assign_prob_new(compositions, list_scores, species_list, verbose=False):
    """ Assign probabilities to novel compositions based on the list of score values.
    Args:
        compositions (list): list of pymatgen species (possibly as generated by smact)
        list_scores (dict): Lists of scores for the species in spec_list
        spec_list (list): Pymatgen species in same order as corresponding
        verbose (bool): Explicitly print any compounds that were skipped over
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
            overall_score = np.mean(scores)
        except:
            if verbose:
                print('Could not get score for: {}'.format(comp))
            overall_score = 0
        probabilities_list.append(overall_score)
    return probabilities_list

def plot_scores_hist(scores, bins = 100, plot_title='plot', xlim = None):
    """ Plot histogram of a list of scores.        """
    print('Number of individual values: {}'.format(len(scores)))
    print('Number of zero values: {}'.format(scores.count(0.0)))
    print('Maximum score: {}'.format(max(scores)))

    plt.hist(scores, bins, color='orange', alpha=0.5)
    if xlim:
        plt.xlim(xlim[0],xlim[1])
    plt.title(plot_title)
    plt.xlabel('Score')
    plt.savefig('{}.png'.format(plot_title), dpi=300)
    plt.show()

def ternary_smact_combos(position1, position2, position3, threshold = 8):
    """ Get Pymatgen species compositions using SMACT when up to three different
        lists are needed to draw species from. E.g. Ternary metal halides...
    Args:
        positionn (list of species): Species to be considered iteratively for each
                                     position.
        """

    print('Generating combinations...')
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
