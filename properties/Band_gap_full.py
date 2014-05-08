#!/usr/bin/env python

from numpy import product, sqrt
from smact_data import get_eig, get_eig_s, get_covalent
from scipy.constants import m_e, physical_constants

# Set constants
Eta_sss = -1.40
Eta_sps = 1.84
Eta_pps = 3.24
Eta_ppp = -0.81

hbarsq_over_me = 7.62

# Set Anion and Cation and 
Cat = raw_input("Enter Cation Symbol:")
An = raw_input("Enter Anion Symbol:")

d = float(raw_input("Enter distacne: "))

# Calculate V1 values
V1_Cat = (get_eig(Cat) - get_eig_s(Cat))/4
V1_An = (get_eig(An) - get_eig_s(An))/4

# Calculate V_sss, V_pps and V_ppp
V_sss = Eta_sss * hbarsq_over_me / (d**2)
V_pps = Eta_pps * hbarsq_over_me / (d**2)
V_ppp = Eta_ppp * hbarsq_over_me / (d**2)

# Calculate E_ss and E_xx
E_ss = V_sss
E_xx = (V_pps/3) + (2*V_ppp/3)

# Calculate Bracket 1 (B1) and Bracket 2 (B2) and B3
B1 = sqrt(((get_eig_s(Cat) - get_eig_s(An))/2)**2 + (4*E_ss)**2)
B2 = sqrt(((get_eig(Cat) - get_eig(An))/2)**2 + (4*E_xx)**2)
B3 = -2*(V1_Cat + V1_An)
# Calculate Band gap
Band_gap = -B1 + B2 - B3

print "Bracket 1: ", B1
print "Bracket 2: ", B2
print "2(v1cat + v1an ): ", B3
print "Cation p-eig: ", get_eig(Cat)
print "Cation s-eig: ", get_eig_s(Cat)
print "Anion p-eig: ", get_eig(An)
print "Anion s-eig: ", get_eig_s(An)

print "---", Band_gap, "---"
