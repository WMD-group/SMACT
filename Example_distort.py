import ase
import smact_builder as builder
import smact_distorter as distort
import numpy as np

# Build the input
test_case = builder.cubic_perovskite(['Ba','Ti','O'],[1,2,2])


# Get arrays of site labels and positions
positions = test_case.get_scaled_positions()
atomic_labels = test_case.get_chemical_symbols()

#Build a sub-lattice you wish to disorder [test case do the O sub-lattice]
sub_lattice = []
i = 0
for atom in atomic_labels:
    if atom == "Ba":
        sub_lattice.append(positions[i])
    i = i + 1
i = 0

# Enumerate the ineuivalent sites
inequivalent_sites = distort.get_inequivalent_sites(sub_lattice,test_case)

print "SUB LATTICE"
print "----- ----- ----- ----- ----- -----"
print sub_lattice
print "INEQUIVALENT SITES"
print "----- ----- ----- ----- ----- -----"
print inequivalent_sites
