#!/usr/bin/env python
from numpy import product,sqrt
"""get_eig function not yet working"""
from smact_data import get_eig,get_covalent
from scipy.constants import hbar m_e

"""Get element names"""
An = 
Cat =

"""calculations"""
V2 = (2.16*(hbar)**2)/m_e*((get_covalent(An))+(get_covalent(Cat))**2)

V3 =  (get_eig(Cat) + get_eig(An))/2 

Ev = V3 + sqrt((V2**2)+(V3**2))

print "Valence Band Minimum = ", Ev
