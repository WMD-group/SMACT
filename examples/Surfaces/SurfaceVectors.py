
# coding: utf-8

# In[55]:

import smact.core as core
import ase.io.vasp as io
import ase.lattice.general_surface as surface
import math


# In[56]:

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))
def length(v):
  return math.sqrt(dotproduct(v, v))
def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))


# In[57]:

def surface_vectors(lattice,miller):
    surf = surface.surface(lattice,miller,layers=1)
    vectors = surf.cell[0:2]
    return (length(vectors[0]),length(vectors[1]),np.degrees(angle(vectors[0],vectors[1])))


# In[65]:

xtal = io.read_vasp('CONTCAR')
vec = surface_vectors(xtal,(1,1,1))
vec


# In[ ]:



