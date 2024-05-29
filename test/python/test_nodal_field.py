import sys
sys.path.append('../../build/python/')


import uguca as ug
import numpy as np


length = 1.
nb_elements = 4

mesh = ug.SimpleMesh(Lx=length, Nx=nb_elements)

coords = mesh.getLocalCoords()

test_field = ug.NodalField(mesh, {ug.x, ug.y})
test_field.setAllValuesTo(2, 0)
test_field.storage()[1, :] = 1. 
print(test_field.storage())


print(np.linalg.norm(test_field.storage(), axis=0))