from smact_lattice import *
from smact_builder import *
from ase.io import *


# Generate a dictionary elements, from the dataset oxidationstates.data
# Dictionary contains elements and their oxidation states
# Reduce the regions of the periodic table to visit, by using search_space
search_space = {'Pb','Ti','O-','Cs','Sn','F-','Cl','I-','Ca','Sr'}

elements = {}
f = open('test/oxidationstates.data','r')
lines = f.readlines()
f.close()
for line in lines:
    inp = line.split()
    for item in search_space:
        if inp[0].startswith(item): 
    	    key = inp[0]
    	    elements[key] = inp[1]


perovskite = Lattice(["A","B","C"],[1,1,3],[[1,2],[2,3,4],[-1,-2]])
perovskite_compositions = possible_compositions(perovskite, elements)

i = 0
for composition in perovskite_compositions:
    print composition
#    system = cubic_perovskite(composition)
#    write('%s.cif'%i, system, format='cif')
    i = i + 1
