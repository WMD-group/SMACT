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
import matplotlib as mpl
import matplotlib.pyplot as plt
from pymatgen import Structure, Specie, MPRester
from smact import ordered_elements, Element

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


def mp_filter(criteria):
    """Get a list of materials project mpids that match a set of given criteria.
    The criteria should be in the format used for a MP query.
    Args:
        criteria (dict): Criteria that can be used in an MPRester query
    """
    m = MPRester(os.environ.get("MP_API_KEY"))
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


def species_totals(structures):
    """Given a set of pymatgen structures in the form of dictionaries where
    the Structure is keyed as 'structure', returns a dictionary of the number
    of compounds that features each Species.
    Args:
        structures (list): Dictionaries containing pymatgen Structures
    """
    totals = []
    for i in structures:
        comp = [j for j in i['structure'].composition]
        totals.append(comp)
    totals = [i for sublist in totals for i in sublist]
    totals = dict(Counter(totals))
    return(totals)

def find_instances(anion, structures):
    """Finds the number of instances of each species in a list of pymatgen
    Structure objects when a given anion is the most electronegative one
    present.
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
                an_containing.append(comp)

    an_containing = [i for sublist in an_containing for i in sublist]
    an_containing = dict(Counter(an_containing))
    an_containing.pop(anion)
    return(an_containing)

def generate_scores(spec_list, anion_count_dict, spec_totals, scoring):
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
                    an_vals.append(anion_count_dict[an][i]/spec_totals[i]*\
                    spec_totals[an])
                else:
                    an_vals.append(0.)
            list_vals[an] = an_vals
    return(list_vals)


def plot_all_scores(list_scores, spec_list, break_points, figure_filename):
    """ Plots correlation plots for all the species considered.
    Args:
        list_scores (dict): Lists of scores for the species in spec_list
        keyed by anion.
        spec_list (list): Pymatgen species in same order as corresponding
        scores in each list within the dict list_scores
        break_points (list): ints, positions in lists to start each new plot
        figure_filename (str): Name of .png file to save plot
         """
    arr = np.array([list_scores[key] for key in list_scores])
    arr = arr.transpose()

    # y-axis
    y_labels = []
    for i in spec_list:
        y_lab = '{0}(+{1})'.format(i.symbol,i.oxi_state)
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
          axarr[4,0]   ]

    for ax, i in zip(axarrs, list(chunks.keys())):
        ax1 =ax.imshow(chunks[i], interpolation='None',
        cmap = mpl.cm.get_cmap('inferno_r'), vmin = 0, vmax =1)
        plt.sca(ax)
        plt.yticks([i for i,n in enumerate(chunks[i])], spec_labels[i])

    # Tidy up
    plt.tight_layout()
    f.subplots_adjust(right = 0.85)
    cbar_ax = f.add_axes([0.9,0.05,0.05,0.7])
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
