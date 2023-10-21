import puguca
import numpy as np

lengths = [0.25, 0.25]
nb_elements = [128, 128]

mesh = puguca.SimpleMesh(Lx=lengths[0], Nx=nb_elements[0],
                         Lz=lengths[1], Nz=nb_elements[1])

coords_x = mesh.getLocalCoords(0)
coords_z = mesh.getLocalCoords(2)
print(coords_x)
print(coords_z)

law = puguca.BarrasLaw(mesh, 3.5e6, 2e-5)

top_mat = puguca.Material(7e9, 0.33, 2000)
top_mat.readPrecomputedKernels();

interface = puguca.DefRigInterface(mesh, top_mat, law)

loads = interface.getLoad()
loads.component(1)[:] = 1e6

total_duration = 4.5e-4
time_factor = 0.3
time_step = time_factor*interface.getStableTimeStep()
interface.setTimeStep(time_step)
interface.init(False)

crack_length = 0.05
indexes = np.where((np.abs(coords_x - lengths[0]/2.) < crack_length/2.) &
                   (np.abs(coords_z - lengths[1]/2.) < crack_length/2.))[0]
print(indexes)
tau_max = law.getTauMax()
tau_max[indexes] = 0.


interface.initDump('fracture_3d_example', '.')
interface.registerDumpFields('cohesion_0,cohesion_1,cohesion_2,top_disp_0,top_disp_1,top_disp_2')
interface.dump(0, 0)


nb_time_steps = int(total_duration/time_step)

for s in range(nb_time_steps):
    interface.advanceTimeStep(True)

    if s%60 == 0:
        print(s)
        interface.dump(s, s*time_step)
