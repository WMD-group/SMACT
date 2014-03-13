from smact_lattice import *
#from smact_builder import *
#from ase.io import *
from compound_electroneg import *
from compound_electroneg_pauling import *
from Band_gap_simple import *
import copy
from smact_data import get_ionic, get_pauling, get_covalent

# Generate a dictionary elements, form the dataset oxidationstates.data
# Dictionary contains elements and their oxidation states
# Reduce the regions of the periodic table to visit, by using search_space
search_space = {'Li','Be','Na','Mg','K','Ca','Rb','Sr','Cs','Ba','Al','Si','Ga','Ge','As','In','Sn','Sb','Te','Tl','Pb','Bi','Po','At','S','O','Se','F','Cl','Br','Zn','Cu','I'}
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
simple = Lattice(["A","B"],[1,1],[[1,2],[-1,-2]])
simple_compositions = possible_compositions(simple, elements)


# Convert elements + oxidation states to element symbols only
for i in simple_compositions:
	for j in range(len(i)):
		i[j]= ''.join(c if (c not in map(str,range(0,10))) and (c != '-' ) else "" for c in i[j])
		

# Calculate Band Gap and Mulliken electronegativity for each composition and write to file
# using compound_electroneg function from compound_electroneg.py
out = open("output.txt","w")
print "no. ", "band gap ", "Mulliken e-neg ", "internuclear dist"
out.write("no. " +", " + "band gap " +", " + "Mulliken e-neg " +", " + "internuclear dist" + "\n")

element_compositions = copy.deepcopy(simple_compositions)

counter = 1
for combination, index in zip(simple_compositions, element_compositions):
	if float(abs(get_pauling(combination[1])-get_pauling(combination[0]))) >= 1:
		dist = float(get_ionic(combination[1])+get_ionic(combination[0]))
	else:
		dist = float(get_covalent(combination[1])+get_covalent(combination[0]))
	bandgap = band_gap_simple(False, combination[1], combination[0], dist)
	eneg = compound_electroneg_pauling(False,combination,[1,1])
	
	print str(index), "--", str(bandgap), "--", str(eneg), "--", str(dist)
	out.write(str(index) + ", " + str(bandgap) + ", " + str(eneg) + ", " + str(dist) + "\n")

#	result = compound_electroneg(False,combination,[1,1,3])
#	out.write(str(counter) + ". " + str(index) + ", " + str(result) + "\n")
#	print str(counter), "--", str(index), "--", str(result)
#	counter = counter + 1
out.close()

#    system = cubic_perovskite(composition)
#    write('%s.cif'%i, system, format='cif')
#    i = i + 1
