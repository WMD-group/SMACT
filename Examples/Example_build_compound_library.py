from smact_lattice import *
#from smact_builder import *
#from ase.io import *
from compound_electroneg import *
from Band_gap_full import *
from Band_gap_simple import *
import copy

# Generate a dictionary elements, form the dataset oxidationstates.data
# Dictionary contains elements and their oxidation states
# Reduce the regions of the periodic table to visit, by using search_space
search_space = {'Ti','O-','Sn','F-','C','Sr','Mg','Cu','Li','S','Si','Ge'}
# Get the list of possible constituant elements
elements = {}
f = open('oxidationstates.data','r')
lines = f.readlines()
f.close()
for line in lines:
    inp = line.split()
    for item in search_space:
        if inp[0].startswith(item): 
    	    key = inp[0]
    	    elements[key] = inp[1]

# Generate list of compositions which satisfy charge neutrality
perovskite = Lattice(["A","B","C"],[1,1,3],[[1,2],[2,3,4],[-1,-2]])
perovskite_compositions = possible_compositions(perovskite, elements)


# Convert elements + oxidation states to element symbols only
for i in perovskite_compositions:
	for j in range(len(i)):
		i[j]= ''.join(c if (c not in map(str,range(0,10))) and (c != '-' ) else "" for c in i[j])
		
# Calculate Mulliken electronegativity for each composition and write to file
# using compound_electroneg function from compound_electroneg.py
out = open("output.txt","w")
element_compositions = copy.deepcopy(perovskite_compositions)
counter = 1
for combination, index in zip(perovskite_compositions, element_compositions):
	result = compound_electroneg(False,combination,[1,1,3])
	out.write(str(counter) + ". " + str(index) + ", " + str(result) + "\n")

	counter = counter + 1
out.close()

#    system = cubic_perovskite(composition)
#    write('%s.cif'%i, system, format='cif')
#    i = i + 1
