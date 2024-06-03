import sys
sys.path.append('../../build/python/')


import uguca as ug
import numpy as np


length = 1.
nb_elements = 4

mesh = ug.mesh.SimpleMesh(Lx=length, Nx=nb_elements)

coords = mesh.getLocalCoords()

test_field = ug.mesh.NodalField(mesh, {ug.x, ug.y})
test_field.setAllValuesTo(2, 0)
test_field.array()[1, :] = 1. 
print(test_field.array())


print(np.linalg.norm(test_field.array(), axis=0))
