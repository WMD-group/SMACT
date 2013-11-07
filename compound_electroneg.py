#!/usr/bin/env python
from numpy import product
from smact_data import get_mulliken

"""Get element names"""
elements=raw_input("Enter elements (space separated): ").split(" ")
stoichs=raw_input("Enter stoichiometries (space separated): ").split(" ")
elementlist=list(elements)
stoichslist=list(stoichs)

"""Confirt stoichslist from string to float"""
stoichslist=map(float, stoichslist)

#print "List of elements=", elementlist
#print "Relative ratios=", stoichslist

"""Get mulliken values for each element"""
for i in range(0,len(elementlist)):
	elementlist[i]=get_mulliken(elementlist[i])

print "Electronegativities of elements=", elementlist

"""Raise each electronegativity to it's appropriate power"""
for i in range(0,len(elementlist)):
	elementlist[i]=[elementlist[i]**stoichslist[i]]

#print "Electronegativities raised to powers=", elementlist

"""Calculate final answer"""
prod = product(elementlist)
#print "Product=", prod
compelectroneg = (prod)**(1.0/(sum(stoichslist)))
print "Geometric mean = Compound 'electronegativity'=", compelectroneg
