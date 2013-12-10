#!/usr/bin/env python
from numpy import product,sqrt
"""get_eig function not yet working"""
from smact_data import get_covalent
from scipy.constants import hbar, m_e
"""Get element names"""
An =raw_input("Enter Anion Symbol: ") 
Cat =raw_input("Enter Cation Symbol: ")
eig_an = 10.27 
eig_cat = 4.14

"""calculations"""
V2 = (2.16*(hbar)**2)/m_e*(((get_covalent(An)**-9))+((get_covalent(Cat))**-9)**2)

V3 =  (eig_cat - eig_an)/2 

Ev = ((eig_cat + eig_an)/2) - (sqrt((V2**2)+(V3**2)))

print An, get_covalent(An)
print Cat, get_covalent(Cat)
print "V2:", V2, "V3:", V3
print "Valence Band Maximum = ", Ev
