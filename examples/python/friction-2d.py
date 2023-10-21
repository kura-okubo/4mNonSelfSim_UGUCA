import puguca
import numpy as np
import sys

length = 1.
nb_elements = 1024

mu_s = 0.8
mu_k = 0.2
mus_ampl = 0.15

mesh = puguca.SimpleMesh(Lx=length, Nx=nb_elements)

coords_x = mesh.getLocalCoords(0)

law = puguca.LinearCoulombFrictionLaw(mesh, mu_s, mu_k, 1e-5)

top_mat = puguca.Material(7e9, 0.33, 2000)
top_mat.readPrecomputedKernels();

interface = puguca.UnimatShearInterface(mesh, top_mat, law)

total_duration = 2.3e-4
time_factor = 0.4
time_step = time_factor*interface.getStableTimeStep()
interface.setTimeStep(time_step)

mus = law.getMuS()
mus_min = sys.float_info.max
mus[:] = mu_s + mus_ampl*(1.0*np.sin(2*coords_x*np.pi) +
                          0.5*np.cos(12*coords_x*np.pi) +
                          0.7*np.sin(18*coords_x*np.pi))

mus_min = min(mus_min, np.min(mus))


normal_load = -5e6
shear_load_rate = 1e10

loads = interface.getLoad()
loads.component(0)[:] = mus_min*np.abs(normal_load)
loads.component(1)[:] = normal_load

interface.init(False)

interface.initDump('friction_2d_example', '.')
interface.registerDumpFields('cohesion_0,top_disp_0,top_velo_0,mu_s')
interface.dump(0, 0)


nb_time_steps = int(total_duration/time_step)

for s in range(nb_time_steps):
    loads.component(0)[:] = mus_min*np.abs(normal_load) + s*time_step*shear_load_rate 
    interface.advanceTimeStep(True)

    if s%10 == 0:
        print(s)
        interface.dump(s, s*time_step)
