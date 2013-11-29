#!/usr/bin/env python

################################################################################
# Copyright Daniel Davies, Adam J. Jackson (2013)                              #
#                                                                              #
#  This file is part of SMACT: smact_data.py is free software: you can         #
#  redistribute it and/or modify it under the terms of the GNU General Public  #
#  License as published by the Free Software Foundation, either version 3 of   #
#  the License, or (at your option) any later version.                         #
#  This program is distributed in the hope that it will be useful, but WITHOUT #
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       #
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for   #
#  more details.                                                               #
#  You should have received a copy of the GNU General Public License along with#
#  this program.  If not, see <http://www.gnu.org/licenses/>.                  #
################################################################################

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

# AJJ: commenting this out while incomplete to avoid breaking things
#"""Need to link to database of eigenvalues (solid_properties.txt) """
#def get_eig(element):
#
#	from **** import ****
#	A = ****
#	eig= ****
#	return eig
