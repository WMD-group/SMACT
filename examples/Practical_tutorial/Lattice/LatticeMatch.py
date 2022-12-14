import glob
import itertools
import math
import re

# import argparse
from optparse import OptionParser

import ase.build.surface as surface

# import smact.core as core
import ase.io as io
import numpy as np


# We need a class "pair" which contains the information about a matching interface pair
class Pair:
    """Class providing standard nformation on interface matching pairs."""

    def __init__(
        self,
        material1,
        material2,
        surface1,
        surface2,
        multiplicity1,
        multiplicity2,
        strains,
    ):
        """
        Attributes:
                Pair.material1 : name of first material
                Pair.material2 : name of second material
                Pair.surface1 : miller index of surface 1 Format: [1,1,1]
                Pair.surface2 : miller index of surface 2 Format: [1,1,1]
                Pair.multiplicity1 : multiplicity of u,v in surface1 Format: [1,1]
                Pair.multiplicity2 : multiplicity of u,v in surface2 Format: [1,1]
                Pair.strains : stains in u,v and gamma Format: [0.1,0.1,0.1]
        """
        self.material1 = material1
        self.material2 = material2
        self.surface1 = surface1
        self.surface2 = surface2
        self.multiplicity1 = multiplicity1
        self.multiplicity2 = multiplicity2
        self.strains = strains


# parser = argparse.ArgumentParser()
parser = OptionParser()
parser.add_option(
    "-a",
    "--matera",
    action="store",
    type="string",
    dest="mater1",
    default="material.cif",
    help="The first material. Default material.cif",
)
parser.add_option(
    "-b",
    "--materb",
    action="store",
    type="string",
    dest="mater2",
    default="material-b.cif",
    help="The second material. Default material-b.cif",
)
parser.add_option(
    "-s",
    "--strain",
    action="store",
    type="float",
    dest="strain",
    default="0.04",
    help="Maximum strain between the interfaces. Default 0.04",
)
parser.add_option(
    "-l",
    "--limit",
    action="store",
    type="int",
    dest="limit",
    default=5,
    help="Maximum number of supercell expansions for each interface. Default 5",
)
parser.add_option("-v", action="store_true", dest="verbose")
# parser.add_option("-v", "--verbose", dest="verbose", type='bool', default=True,  help="increase output verbosity")
# parser.add_argument("--print_strains", help="increase output verbosity of strian output")
(options, args) = parser.parse_args()

# Define some basic algebra


def dotproduct(v1, v2):
    return sum((a * b) for a, b in zip(v1, v2))


def length(v):
    return math.sqrt(dotproduct(v, v))


def angle(v1, v2):
    return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))


# Define the tests as required by Fig 2 of Zur and McGill
def test1(a, b):
    if np.dot(a, b) < 0:
        return False
    else:
        return True


def test2(a, b):
    if np.linalg.norm(a) > np.linalg.norm(b):
        return False
    else:
        return True


def test3(a, b):
    if np.linalg.norm(b) > np.linalg.norm(b + a):
        return False
    else:
        return True


def test4(a, b):
    if np.linalg.norm(b) > np.linalg.norm(b - a):
        return False
    else:
        return True


# Switching of vectors as prescribed by Fig 2 of Zur and McGill
def cond1(a, b):
    if np.dot(a, b) < 0:
        b = -b
    return a, b


def cond2(a, b):
    if np.linalg.norm(a) > np.linalg.norm(b):
        c = b
        b = a
        a = c
        cond1(a, b)
    return a, b


def cond3(a, b):
    if np.linalg.norm(b) > np.linalg.norm(b + a):
        b = b + a
    return a, b


def cond4(a, b):
    if np.linalg.norm(b) > np.linalg.norm(b - a):
        b = b - a
    return a, b


# §
def reduce_vectors(va, vb):
    """Reduce the surface vectors to their minimal cell, as outlined in figure 2 of
    J. Appl. Phys. 55, 380 (1984)
    Args:
    va,vb: lists of real numbers/integers. Minimum length 2
    Returns:
    a,b: arrays of dimension (2x1). The reduced vectors

    """
    a = np.asarray(va[0:2])
    b = np.asarray(vb[0:2])
    test_truth = 0  # A test that all 4 necessary conditions are met
    while test_truth < 4:
        a, b = cond1(a, b)
        a, b = cond2(a, b)
        a, b = cond3(a, b)
        a, b = cond4(a, b)
        truths = [test1(a, b), test2(a, b), test3(a, b), test4(a, b)]
        test_truth = sum(truths)

    return a, b


def surface_vectors(lattice, miller):
    """Given the xtal as defined with a (3x3) cell uses the ase surface module to cut the required surface.
    Args:
            lattice : ase Atoms object
            miller  : miller indices of the surface, a tuple of integers, length 3.
    Returns:
            vectors[0/1] : the surface vectors (u,v), list of real numbers.
    """
    surf = surface(lattice, miller, layers=1)
    vectors = surf.cell[0:2]
    return vectors[0], vectors[1]


def surface_ratios(surface_a, surface_b, threshold=0.05, limit=5):
    """Given two surfaces a and b defined by vectors (u,v) tests to see if there is a ratio of r1/r2 which gives a lattice vector with mismatch less than the threshold.
    Args:
            surfaces_ : the surface vectors and angle. A (3) tuple of real numbers.
            threshold : the limit for lattice mismatch.
            limit : the maximum number to multiply lattices by to obtain match.
    Returns:
            exists : a bool of wheter the match was found
            multiplicity : a list (1x2) with the integer values to multiply by.
    """

    epitaxy = False
    u_satisfied = False
    v_satisfied = False
    angles_match = False

    super_a = (0, 0)
    super_b = (0, 0)
    strains = [0, 0, 0]
    if options.verbose:
        print(" Entering matching ratios routine")
        print("------   -------  -------  -------")
        print("Surface a: ", surface_a)
        print("Surface b: ", surface_b)
        print("------   -------  -------  -------")
    #   Ensure that the vector angles match (up to an arbitrary 180 degree rotation)
    if surface_a[2] - surface_b[2] <= 0.2:
        angles_match = True
    strains[2] = surface_a[2] - surface_b[2]
    if options.verbose:
        print("Angles are within tolerance")
    #   Run the test for the first lattice vector set
    if angles_match:
        for i in range(1, limit + 1):
            if u_satisfied:
                break
            for j in range(1, limit + 1):
                r1 = float(i * surface_a[0])
                r2 = float(j * surface_b[0])
                strain = 2 * abs(r1 - r2) / (r1 + r2)
                if strain < threshold:
                    a = (i, j)
                    u_satisfied = True
                    strains[0] = strain
                    break

    #   Run the test for the second lattice vector set
    if u_satisfied:  # Dont bother with v if u not satisfied
        for i in range(1, limit + 1):
            if v_satisfied:
                break
            for j in range(1, limit + 1):
                r1 = float(i * surface_a[1])
                r2 = float(j * surface_b[1])
                strain = 2 * abs(r1 - r2) / (r1 + r2)
                if strain < threshold:
                    b = (i, j)
                    v_satisfied = True
                    strains[1] = strain
                    break
    if angles_match and u_satisfied and v_satisfied:
        super_a = (a[0], b[0])
        super_b = (a[1], b[1])
    epitaxy = True
    return epitaxy, super_a, super_b, strains


material1 = re.sub(r"\.cif$", "", options.mater1)
material2 = re.sub(r"\.cif$", "", options.mater2)

# Code output header, for God's sake clean this up!
# God's will be done.
"""
print "###-----------------------------------###"
print "###                                   ###"
print "       Scanning for srufaces of "
print "         ",material1, material2
print "###                                   ###"
print "###-----------------------------------###"
"""

# Scan through the various miller indices
indices_a = list(itertools.product([0, 1], repeat=3))
indices_a = [(1, 0, 0), (1, 1, 0)]  # Just consider these two for MAPI
xtalA = io.read(options.mater1)
xtalB = io.read(options.mater2)
# mathced_pairs is a list containing the Pairs class for all unique interfaces
matched_pairs = []
print(" ")
print(
    "-------------------------------------------------------------------------------------------------"
)
print(
    "-------------------------------------------------------------------------------------------------"
)
print(
    "(miller 1) (miller 2) (mult 1) (mult 2) [Strain in vector u, Strain in vector v, angle mis-match]"
)
# Set the  values for material A
for index_a in indices_a:
    if index_a != (0, 0, 0):
        vec1, vec2 = surface_vectors(xtalA, index_a)
        r_vec1, r_vec2 = reduce_vectors(vec1, vec2)
        surface_vector_1 = (length(r_vec1), length(r_vec2), angle(r_vec1, r_vec2))
    # Set the  values for material B
    indices_b = list(itertools.product([0, 1], repeat=3))
    for index_b in indices_b:
        if index_b != (0, 0, 0):
            vec1, vec2 = surface_vectors(xtalB, index_b)
            r_vec1, r_vec2 = reduce_vectors(vec1, vec2)
            surface_vector_2 = (length(r_vec1), length(r_vec2), angle(r_vec1, r_vec2))
            epitaxy, a, b, strains = surface_ratios(
                surface_vector_1,
                surface_vector_2,
                threshold=options.strain,
                limit=options.limit,
            )
            # print index_a, index_b, surface_vector_1, surface_vector_2
            if epitaxy:
                if options.verbose:
                    print("------   ------   ------   ------   ------")
                    print(
                        "Surface A  multiplicity: ",
                        a,
                        "semiconductor multiplicity: ",
                        b,
                    )
                    print("Strain vector: ", strains)
                    surface_super_cell_a = np.asarray(
                        surface_vector_1[0:2]
                    ) * np.asarray(a)
                    print("Surface super-cell vector: ", surface_super_cell_a)
                    surface_super_cell_b = np.asarray(
                        surface_vector_2[0:2]
                    ) * np.asarray(b)
                    print("Surface super-cell vector: ", surface_super_cell_b)
                    print("------   ------   ------   ------   ------")
                else:
                    new = Pair(material1, material2, index_a, index_b, a, b, strains)
                    isnew = True
                    if len(matched_pairs) == 0 and new.strains[2] == 0.0:
                        matched_pairs.append(new)
                    for old in matched_pairs:
                        if (
                            new.multiplicity1 == old.multiplicity1
                            and new.multiplicity2 == old.multiplicity2
                            and new.strains == old.strains
                        ):
                            isnew = False
                        else:
                            continue
                if len(matched_pairs) > 0 and isnew and new.strains[2] == 0.0:
                    matched_pairs.append(new)
                # print("%s ; %s ; %s; %s ; %s; %s" % ('Found',a,index_a,b,index_b,strains))
print(material1, material2, len(matched_pairs))
for pair in matched_pairs:
    print(
        pair.surface1,
        pair.surface2,
        pair.multiplicity1,
        pair.multiplicity2,
        pair.strains,
    )

print(
    "#################################################################################################"
)
# print "###-----------------------------------###"
# print "###                                   ###"
# print "###-----------------------------------###"
