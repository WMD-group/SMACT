# Example script using distorter, generate all possible 
# (symmetry inequivalent) substitutions of Sr on Ba  sites; 
# single and double substitutions.

import ase
import smact.builder as builder
import smact.distorter as distort
import numpy as np

# Build the input
test_case = builder.cubic_perovskite(['Ba','Ti','O'],[2,2,2])


print "------------------------------"
print "Original coordinates: ", test_case
print "------------------------------"


# Do the single substitution first, it is trivial as all Ba 
#  sites are equivalent we will choose the first Ba
subs_site = [0.0, 0.0, 0.0]
single_substitution = distort.make_substitution(test_case,subs_site,"Sr")
print "Single: ", single_substitution


#Build a sub-lattice you wish to disorder [test case do the Ba sub-lattice]
sub_lattice = distort.build_sub_lattice(single_substitution,"Ba")


# Enumerate the inequivalent sites
inequivalent_sites = distort.get_inequivalent_sites(sub_lattice,single_substitution)


# Replace Ba at inequivalent sites with Sr
for inequivalent_site in inequivalent_sites:
    print "------------------------------"
    print " Substituted coordinates" 
    #print test_case,inequivalent_site
    distorted = distort.make_substitution(single_substitution,inequivalent_site,"Sr")
    print distorted
    

print "------------------------------"
print "------------------------------"
