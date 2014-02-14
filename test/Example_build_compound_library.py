from smact_lattice import *
from CellBuilder import *
from ase.io import *

elements ={'Cu' : 2, 'Pb' : 2, 'Ti' : 4, 'I' : -1, 'O' : -2, 'Nb' : 4, 'Cl' : -1, 'Sn' : 1, 'S' : -2}

perovskite = Lattice(["A","B","C"],[1,1,3],[[1,2],[2,3,4],[-1,-2]])
perovskite_compositions = possible_compositions(perovskite)

i = 0
for composition in perovskite_compositions:
    print composition
#    system = cubic_perovskite(composition)
#    write('%s.cif'%i, system, format='cif')
    i = i + 1
