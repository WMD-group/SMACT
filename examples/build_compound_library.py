from smact.lattice import *
import smact.builder as builder
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
search_space = ['Ti','O','Sn','F','Sr','Mg','Cu','Li','S','Si','Ge']
# Generate list of compositions which satisfy charge neutrality
#perovskite = Lattice(["A","B","C"],[1,1,3],[[1,2],[2,3,4],[-1,-2]])
perovskite = builder.cubic_perovskite(['H','H','H'])   # H is just a place holder, which must be a legitimate atom type
perovskite_compositions = possible_compositions(perovskite, search_space)


