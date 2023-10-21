import puguca
import numpy as np

length = 1.
nb_elements = 1024

mesh = puguca.SimpleMesh(Lx=length, Nx=nb_elements)

coords_x = mesh.getLocalCoords(0)

law = puguca.BarrasLaw(mesh, 3.5e6, 2e-5)

top_mat = puguca.Material(7e9, 0.33, 2000)
top_mat.readPrecomputedKernels();

bot_mat = puguca.Material(5e9, 0.25, 1200)
bot_mat.readPrecomputedKernels();

interface = puguca.BimatInterface(mesh, top_mat, bot_mat, law)

loads = interface.getLoad()
loads.component(0)[:] = 2e6
loads.component(1)[:] = 1e6

total_duration = 2.6e-4
time_factor = 0.4
time_step = time_factor*interface.getStableTimeStep()
interface.setTimeStep(time_step)
interface.init(False)

crack_length = 0.05
indexes = np.where(np.abs(coords_x - length/2.) < crack_length/2.)[0]

tau_max = law.getTauMax()
tau_max[indexes] = 0.


interface.initDump('fracture_2d_example-py', '.')
interface.registerDumpFields('cohesion_0,cohesion_1,top_disp_0,top_disp_1,bot_disp_0,bot_disp_1,tau_max')
interface.dump(0, 0)


nb_time_steps = int(total_duration/time_step)

for s in range(nb_time_steps):
    interface.advanceTimeStep(True)

    if s%100 == 0:
        print(s)
        interface.dump(s, s*time_step)
