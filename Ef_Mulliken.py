#!/usr/bin/env python
from smact_data import get_pauling

Cat_list = ['Si','Ge','In','Ga','In','Al','Ga','In','Zn','Zn','Cd','Hg','Zn','Cd','Hg','Zn','Cd','Hg','K','Na','K','Cs','K','Cs','Li','Na','K','Rb','Cs']
An_list = ['Si','Ge','P','As','As','Sb','Sb','Sb','O','S','S','S','Se','Se','Se','Te','Te','Te','F','Cl','Cl','Cl','Br','Br','I','I','I','I','I']
print len(Cat_list)
print len(An_list)

i=0
for i in range(0,len(Cat_list)):
	Ef = 2.86*(get_pauling(Cat_list[i])*get_pauling(An_list[i]))**0.5
	print Cat_list[i],An_list[i],Ef


#def compound_pauling(Anion,Cation):
#Cat = raw_input("enter cation: ")
#An = raw_input("enter anion: ")
#Ef = 2.86*(get_pauling(Cat)*get_pauling(An))**0.5

