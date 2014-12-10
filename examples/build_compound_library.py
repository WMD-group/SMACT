from smact.lattice import *
#from smact.builder import *
#from ase.io import *
import copy

import os # get correct path for datafiles when called from another directory
this_directory = os.path.dirname(__file__)
if this_directory:
    this_directory = this_directory + '/'
smact_directory = this_directory + '../smact/'

# Generate a dictionary elements, form the dataset oxidationstates.data
# Dictionary contains elements and their oxidation states
# Reduce the regions of the periodic table to visit, by using search_space
search_space = {'Ti','O-','Sn','F-','C','Sr','Mg','Cu','Li','S','Si','Ge'}
# Get the list of possible constituant elements
elements = {}
f = open(smact_directory + 'data/oxidationstates.data','r')
lines = f.readlines()
f.close()
for line in lines:
    row = line.split()
    for item in search_space:
        if row[0].startswith(item): 
    	    key = inp[0]
    	    elements[key] = inp[1]

# Generate list of compositions which satisfy charge neutrality
perovskite = Lattice(["A","B","C"],[1,1,3],[[1,2],[2,3,4],[-1,-2]])
perovskite_compositions = possible_compositions(perovskite, elements)


