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

