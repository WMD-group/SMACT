from smact.lattice import *
#from smact.builder import *
#from ase.io import *
from smact.properties.compound_electroneg import *
from smact.properties.compound_electroneg_pauling import *
from smact.properties.Band_gap_simple import *
import copy
from smact.data import get_ionic, get_pauling, get_covalent
import os.path

# Get correct directory if calling from elsewhere
# (if statement needed in case of empty response)
this_directory = os.path.dirname(__file__)
if this_directory:
    this_directory = this_directory + '/'
smact_directory = this_directory + '../smact/'

# Generate a dictionary elements, form the dataset oxidationstates.data
# Dictionary contains elements and their oxidation states
# Reduce the regions of the periodic table to visit, by using search_space
search_space = {'Li','Be','Na','Mg','K','Ca','Rb','Sr','Cs','Ba','Al','Si','Ga','Ge','As','In','Sn','Sb','Te','Tl','Pb','Bi','Po','At','S','O','Se','F','Cl','Br','Zn','Cu','I'}
# Get the list of possible constituent elements

# Generate list of binary compositions which satisfy charge neutrality
neutral_species = []
for i, ele_a in enumerate(search_space):
    for ox_a in core.Element(ele_a).oxidation_states:
	for j, ele_b in enumerate(search_space):
	    if j > i:
	   	for ox_b in core.Element(ele_b).oxidation_states:
		    neutral_exists, neutral_ratios = core.charge_neutrality([ox_a,ox_b],threshold = 4)
		    if neutral_exists:
			combo = [core.Element(ele_a).symbol, core.Element(ele_b).symbol]
			neutral_species.append([combo,neutral_ratios])
print len(neutral_species)

# Calculate Band Gap and Mulliken electronegativity for each composition and write to file
# using compound_electroneg function from compound_electroneg.py
out = open("output.txt","w")
print "no. ", "band gap ", "Mulliken e-neg ", "internuclear dist"
out.write("no. " +", " + "band gap " +", " + "Mulliken e-neg " +", " + "internuclear dist" + "\n")

element_compositions = copy.deepcopy(neutral_species)

for combo in neutral_species:
    dist = core.Element(combo[0][0]).covalent_radius + core.Element(combo[0][1]).covalent_radius
    bandgap = band_gap_simple(False, combo[0][0], combo[0][1], dist)

    print combo[0][0], combo[0][1], bandgap
