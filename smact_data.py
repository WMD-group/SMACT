#!/usr/bin/env python

def get_mulliken(element):
	"""Gets Mulliken electroneg from the IE and EA"""
	from chemlab.db import chemlabdb
	A = chemlabdb.ChemlabDB()
	ionpot=A.get('data','ionpotdict')[element]
	eaff=A.get('data','eaffdict')[element]
	mulliken=(ionpot+eaff)/2.0
	return mulliken

def get_pauling(element):

	from chemlab.db import chemlabdb 

	A = chemlabdb.ChemlabDB() 
	pauling=A.get('data','paulingenegdict')[element]
	return pauling

def get_covalent(element):

	from chemlab.db import chemlabdb
	
	A= chemlabdb.ChemlabDB()
	covalent=A.get('data','covalentdict')[element]
	return covalent
"""Need to link to database of eigenvalues (solid_properties.txt) """
def get_eig(element):

	from **** import ****
	A = ****
	eig= ****
	return eig
