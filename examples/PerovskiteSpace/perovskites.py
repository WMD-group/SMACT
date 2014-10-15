#!/usr/bin/env python 
# Using the ase spacegroup module this can build the structure, from 
# the composition, as defined in the smact_lattice module. 
#TO DO:
# Calculating the lattice parameters from the elemnts involved
################################################################################
# Copyright Tim Gauntlett   (2014)  				                           #
#                                                                              #
# This file is part of SMACT: parameters.py is free software: you can          #
# redistribute it and/or modify it under the terms of the GNU General Public   #
# License as published by the Free Software Foundation, either version 3 of    #
# the License, or (at your option) any later version.                          #
# This program is distributed in the hope that it will be useful, but WITHOUT  #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for    #
# more details.                                                                #
# You should have received a copy of the GNU General Public License along with #
# this program.  If not, see <http://www.gnu.org/licenses/>.                   #
#                                                                              #
################################################################################

import smact.core as core
import numpy as np
import csv

def E2_1(shannon_radius): #Cubic Pervoskite
    limiting_factors=[2*sum(shannon_radius[1:])]
    a = max(limiting_factors)
    b = a
    c = a
    space = a * np.sqrt(3) - 2 * shannon_radius[1]
    alpha = 90
    beta = 90
    gamma = 90
    return a,b,c,space,alpha,beta,gamma
e_lattices = [E2_1]

i=0
a_sym =[]
a_cood =[]
a_shannon =[]
b_sym =[]
b_cood =[]
b_shannon =[]
c_sym =['O']
c_cood =[2]
c_shannon =[1.21]



with open('data/shannon_radii.csv','rU') as f:
    reader = csv.reader(f)         
    reader.next() # Skip the first row which is titles      
    for row in reader:
        # print row                  
        if i <= 1000:  #########              ##     <<    <  ###<<< Safety valve  ######
            if row[1] == '2' and row[2] == '6_n'or'8_n':                
                a_sym.append(row[0])         
                #a_cood.append(row[2])
                a_shannon.append(row[3])
            if row[1] == '4' and '6_n':
                b_sym.append(row[0])
                b_cood.append(row[2])
                b_shannon.append(row[3])
        i+=1

for j, bsymbol in enumerate(b_sym):
    for i, symbol in enumerate(a_sym):
        print a_shannon[i],b_shannon[j]
        if a_shannon[i] > b_shannon[j]:  #### define how much greater.
            print symbol,bsymbol
            print a_sym[i],
            print "%.2f"%float(a_shannon[i]),
            print b_sym[j],
            print "%.2f"%float(b_shannon[j]),
            a = 2 * (float(b_shannon[0])) + 1.21
            space = a*(float(3**0.5)) - 2 * float(b_shannon[j])
            print "%.2f"%a,
            print "%.2f"%space
    

    
    
    
    
    
    
    
    
    
    
    
