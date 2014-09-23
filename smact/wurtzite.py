#!/usr/bin/env python 
# Using the ase spacegroup module this can build the structure, from 
# the composition, as defined in the smact_lattice module. 
#TO DO:
# Calculating the lattice parameters from the elemnts involved
################################################################################
# Copyright Tim Gauntlett, Adam J. Jackson   (2014)                            #
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
import os
import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm



# get correct path for datafiles when called from another directory
smact_directory = os.path.dirname(__file__)
# Append a trailing slash to make coherent directory name - this would select the
# root directory in the case of no prefix, so we need to check
if smact_directory:
    smact_directory = smact_directory + '/'


#B4    
def wurtzite(shannon_radius):
    #print shannon_radius
    ##limiting_factors=[2*(max(shannon_radius)*np.sqrt(2)), 4*(shannon_radius[0] + shannon_radius[1])**(1./3.)]
    a = 2*0.817*(shannon_radius[0]+shannon_radius[1])  # 0.817 is sin(109.6/2)
    b = a
    c = (shannon_radius[0]+shannon_radius[1])*(2+2*0.334)  # 0.334 is sin(109.5-90)
    alpha = 90
    beta = 90
    gamma = 120
    inner_space = a * (6**0.5) - (4*shannon_radius[0])
    return a,b,c,alpha,beta,gamma,inner_space

##get info from tables
def radius_covalent(element):
    with open(smact_directory + 'data/Covalent_radii.csv','rU') as f:
        reader = csv.reader(f)         
        reader.next() # Skip the first row which is titles   
        r_covalent=0 
        for row in reader:
            symbols= list(row[0])
            for char in symbols:
                if char == '(':
                    a=symbols.index("(")
                    symbols = symbols[:a]
            symbols="".join(symbols)
            if symbols == element:
                r_covalent = float(row[3])/ 100 #############################################
    return r_covalent
        
def formal_charge(element):
    with open('data/oxidation_states.txt','r') as f:
        data = f.readlines()
        oxidation_states = False
        for line in data:
            if not line.startswith('#'):
                l = line.split()
                if (l[0] == element):
                    if len(l)>1:
                        oxidation_states = l[1:]#############################################
                        for location, value in enumerate(oxidation_states):
                            oxidation_states[location] = int(value)
                    else:
                        oxidation_states = []
    return oxidation_states

def Radius_shannon(element,oxidation):
    with open('data/shannon_radii.csv','rU') as f:
        reader = csv.reader(f)
        r_shannon=False
        for row in reader:
            if element == row[0] and float(oxidation)==float(row[1]) and row[2]=='4_n': #only for Wurtzite
                r_shannon=row[3]
    return r_shannon                 
    
            
#------------------------------------------------------------------------------------------
output=[]

# Processing Lattices
#------------------------------------------------------------------------------------------
lattice_charge=-1

element=[]
all_elements=[]
elements=core.ordered_elements(1,100)     ## The safety valve!
holes=np.zeros(shape=(100,100))
matrix_hhi=np.zeros(shape=(100,100))
with open('data/ordered_periodic.txt','rU') as f:
    reader = csv.reader(f)
    for row in reader:
        all_elements.append(row[0])

tetra_cood=[]
for elm in elements:
    with open('data/shannon_radii.csv','rU') as f:
        reader = csv.reader(f)
        r_shannon=False
        for row in reader:
            if row[2]=="4_n" and row[0]==elm:
                tetra_cood.append(row[0])
elements=[]                
#print tetra_cood
for el in tetra_cood:
    if el not in elements:
        elements.append(el)          
#print elements
for a_element in elements:
    
    #print a_element
    a_formal_charge = formal_charge(a_element)
    if a_formal_charge:#i.e. a_formal_charge has a value
        for a_ox in a_formal_charge:
            a_shan=Radius_shannon(a_element,a_ox)
            a_cov=radius_covalent(a_element)
            if a_shan == False:  ####Covelent radius
                #print a_ox,radius_covalent(a_element)
                for b_element in elements:
                    if a_element != b_element:
                        #print " ",a_element, b_element
                        b_formal_charge = formal_charge(b_element)  ## check again
                        if b_formal_charge: #i.e. b_formal_charge has a value
                            for b_ox in b_formal_charge:
                                if a_ox+b_ox==lattice_charge:
                                    b_cov=radius_covalent(b_element)
                                    #print a_cov,b_cov
                                    g=[]
                                    g.append(float(a_cov))
                                    g.append(float(b_cov))
                                    g.sort(reverse=True)#Needed for calculating cavity
                                    #print g
                                    method="covalent radii"
                                    a,b,c,alpha,beta,gamma,inner_space = wurtzite(g)
                                    hhi_a=core.Element(a_element).HHI_p
                                    hhi_b=core.Element(b_element).HHI_p
                                    hhi_rp=(float(hhi_a)*float(hhi_b))**0.5
                                    holes[all_elements.index(a_element),all_elements.index(b_element)]=inner_space
                                    matrix_hhi[all_elements.index(a_element),all_elements.index(b_element)]=hhi_rp
                                    '''print a_element,b_element,method
                                    print "%.2f"%a,"%.2f"%b,"%.2f"%c,alpha,beta,gamma
                                    print "%.2f"%inner_space'''
                                    #output.append([a_element,b_element,hhi_a,hhi_b,method,a,b,c,inner_space])
           
            else:
                for b_element in elements:
                    if a_element != b_element:
                        #print " ",a_element, b_element
                        b_formal_charge = formal_charge(b_element)  ## check again
                        if b_formal_charge:
                            for b_ox in b_formal_charge:
                                if Radius_shannon(b_element,b_ox) == False:
                                    if a_ox+b_ox==lattice_charge:                       ## Charge of lattice
                                        b_cov=radius_covalent(b_element)
                                        #print a_cov,b_cov
                                        g=[]
                                        g.append(float(a_cov))
                                        g.append(float(b_cov))
                                        g.sort(reverse=True)#Needed for calculating cavity
                                        #print g
                                        method="covalent radii"
                                        a,b,c,alpha,beta,gamma,inner_space = wurtzite(g)
                                        hhi_a=core.Element(a_element).HHI_p
                                        hhi_b=core.Element(b_element).HHI_p
                                        hhi_rp=(float(hhi_a)*float(hhi_b))**0.5
                                        holes[all_elements.index(a_element),all_elements.index(b_element)]=inner_space
                                        matrix_hhi[all_elements.index(a_element),all_elements.index(b_element)]=hhi_rp
                                        '''print a_element,b_element,method
                                        print "%.2f"%a,"%.2f"%b,"%.2f"%c,alpha,beta,gamma
                                        print "%.2f"%inner_space'''
                                        #output.append([a_element,b_element,hhi_a,hhi_b,method,a,b,c,inner_space])
                                    
                                    
                                else:   ### Shannon radius
                                    b_shan=Radius_shannon(b_element,b_ox)
                                    if a_ox+b_ox==lattice_charge:
                                        #print a_shan,b_shan
                                        g=[]
                                        g.append(float(a_shan))
                                        g.append(float(b_shan))
                                        g.sort(reverse=True)#Needed for calculating cavity
                                        #print g
                                        method="shannon radii"
                                        a,b,c,alpha,beta,gamma,inner_space = wurtzite(g)
                                        hhi_a=core.Element(a_element).HHI_p
                                        hhi_b=core.Element(b_element).HHI_p
                                        hhi_rp=(float(hhi_a)*float(hhi_b))**0.5
                                        holes[all_elements.index(a_element),all_elements.index(b_element)]=inner_space
                                        matrix_hhi[all_elements.index(a_element),all_elements.index(b_element)]=hhi_rp
                                        '''print a_element,b_element,method
                                        print "%.2f"%a,"%.2f"%b,"%.2f"%c,alpha,beta,gamma
                                        print "%.2f"%inner_space'''
                                        #output.append([a_element,b_element,hhi_a,hhi_b,method,a,b,c,inner_space])
'''
with open('wurtzite.csv', 'w') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for structure in output:
        spamwriter.writerow(structure)
'''
#np.savetxt('matrix_neg_1.dat',holes)

'''
f=open('matrix.dat',"r")
lines = f.readlines()
f.close
print lines
holes = np.zeros(shape=(100,100))
i = 0
for line in lines:
    inp = line.split()
    print inp
    holes[i,:] = inp[:]
    i = i + 1
 '''  


plt.imshow(matrix_hhi, cmap=cm.BuPu, interpolation='nearest')
plt.xlabel('Elements(A) by proton numbers',fontsize=22)
plt.ylabel('Elements(B) by proton numbers',fontsize=22)
plt.colorbar()
plt.savefig('wurtzite_hhi.png')
plt.show()
plt.imshow(holes, cmap=cm.BuGn, interpolation='nearest')
plt.xlabel('Elements(A) by proton numbers',fontsize=22)
plt.ylabel('Elements(B) by proton numbers',fontsize=22)
plt.colorbar()
plt.savefig('wurtzite_holes.png')
plt.show()


