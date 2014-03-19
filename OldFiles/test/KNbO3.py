import ase
import smact_builder as builder
import numpy as np
from ase.io import *

# Build the input
test_case = builder.cubic_perovskite(['K','Nb','O'],[4.07, 4.07, 4.07, 90, 90, 90])

write('KNbO3_Cubic.cif',test_case)
