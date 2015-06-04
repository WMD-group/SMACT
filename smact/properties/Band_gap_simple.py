#!/usr/bin/env python

################################################################################
# Copyright Daniel Davies (2013)                              				   #
#                                                                              #
# This file is part of SMACT: Band_gap_simple.py is free software: you can 	   #
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
import sys
from numpy import sqrt
from smact.data import get_eig, get_eig_s
import smact.core as core

def band_gap_simple(verbose=False,Anion=None,Cation=None,Distance=None):

	"""Estimates the band gap from elemental data.

	Args:
		Anion: Element symbol of the dominant anion in the system.
		Cation: Element symbol of the the dominant cation in the system.
		Distance: Nuclear separation between anion and cation i.e. sum of ionic
				radii.
				verbose: An optional True/False flag. If True, additional information is
				printed to the standard output. [Defult: False]

				Returns: Band gap as a float in eV. 

	"""
	# Set constants
	hbarsq_over_m = 7.62
	
	# Get anion and cation
	if Anion:
		An = Anion
	if Cation:
		Cat = Cation
	if not Anion:
		An = raw_input("Enter Anion Symbol:")
	if not Cation:
		Cat = raw_input("Enter Cation Symbol:")
	if Distance:
		d = float(Distance)
	if not Distance:
		d = float(raw_input("Enter internuclear separation (Angstroms): "))
	
	# Calculate values of equation components 
	V1_Cat = (core.Element(Cat).eig - core.Element(Cat).eig_s)/4
	V1_An = (core.Element(An).eig - core.Element(An).eig_s)/4
	V1_bar = (V1_An + V1_Cat)/2
	V2 = 2.16 * hbarsq_over_m / (d**2)
	V3 = (core.Element(Cat).eig - core.Element(An).eig)/2
	alpha_m = (1.11*V1_bar)/sqrt(V2**2 + V3**2)

	# Calculate Band gap
	Band_gap = (3.60/3)*(sqrt(V2**2 + V3**2))*(1-alpha_m)
	if verbose:
		print "V1_bar = ", V1_bar
		print "V2 = ", V2
		print "alpha_m = ", alpha_m
		print "V3 = ", V3
#	print "--- band_gap = ", Band_gap, "---"
	return Band_gap
	
if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(
		description ="""Compound band gap estimates from elemental data.""")
	parser.add_argument("-a", "--Anion", type=str, help ="""Element symbol for anion.""")
	parser.add_argument("-c", "--Cation", type=str, help ="""Element symbol for cation.""")
	parser.add_argument("-d", "--Distance", type=float, help ="""Internuclear separation.""")
	parser.add_argument("-v", "--Verbose", help ="""More Verbose output.""",action="store_true")
	parser.add_argument("-q", "--quiet", help="Quiet output", action="store_true")
	args = parser.parse_args()
	if args.quiet:
		verbose_flag = False
	else:
		verbose_flag = True
	print band_gap_simple(verbose=verbose_flag,Anion=args.Anion,Cation=args.Cation,Distance=args.Distance)
	
	
