import ase
from ase.io import *
import smact.builder as build
import smact.surface as surface

PbZrO = build.cubic_perovskite(['Pb','Zr','O'])

S100 = surface.cut100(PbZrO)
S110 = surface.cut110(PbZrO)
S111 = surface.cut111(PbZrO)

write('100.cif',S100)
write('110.cif',S110)
write('111.cif',S111)
write('Bulk.cif',PbZrO)
