#!/usr/bin/env python 
# Surfaces module for SMACT. Currently cuts a series of common surfaces
# TODO:
# 	Systematically determine all possible inequivalent cuts for a surface
# 	Develop surface matching algorithm
################################################################################
# Copyright Keith T Butler    (2013)                                           #
#                                                                              #
# This file is part of SMACT: surface.py is free software: you can             #
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

from ase.utils.geometry import cut

def cut100(crystal,nlayers=6):
    """Cleave the 100 surface"""
    crystal100 = cut(crystal,[0,0,1],[0,1,0],nlayers=nlayers)
    return crystal100

def cut010(crystal,nlayers=6):
    """Cleave the 010 surface"""
    crystal010 = cut(crystal,[1,0,0],[0,0,1],nlayers=nlayers)
    return crystal010

def cut001(crystal,nlayers=6):
    """Cleave the 001 surface"""
    crystal001 = cut(crystal,[0,1,0],[0,0,1],nlayers=nlayers)
    return crystal001

def cut110(crystal,nlayers=6):
    """Cleave the 110 surface"""
    crystal110 = cut(crystal,[0,0,1],[1,1,0],nlayers=nlayers)
    return crystal110

def cut111(crystal,nlayers=6):
    """Cleave the 111 surface"""
    crystal111 = cut(crystal,[1,0,1],[1,1,0],nlayers=nlayers)
    return crystal111
