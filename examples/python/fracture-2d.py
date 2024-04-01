import sys
sys.path.append('../../build/python/')
import uguca as ug
import numpy as np

length = 1.
nb_elements = 1024

mesh = ug.SimpleMesh(Lx=length, Nx=nb_elements)

coords = mesh.getLocalCoords()

law = ug.BarrasLaw(mesh, 3.5e6, 2e-5)

top_mat = ug.Material(7e9, 0.33, 2000)
top_mat.readPrecomputedKernels();

bot_mat = ug.Material(5e9, 0.25, 1200)
bot_mat.readPrecomputedKernels()

interface = ug.BimatInterface(mesh, {ug.x, ug.y},  top_mat, bot_mat, law)

loads = interface.getLoad()

loads[0, :] = 2e6
loads[1, :] = 1e6

total_duration = 2.6e-4
time_factor = 0.4
time_step = time_factor*interface.getStableTimeStep()
interface.setTimeStep(time_step)
interface.init(False)

crack_length = 0.05
indexes = np.where(np.abs(coords[0, :] - length/2.) < crack_length/2.)[0]

tau_max = law.getTauMax()
tau_max[indexes] = 0.


interface.initDump('fracture_2d_example', '.')
interface.registerDumpFields('cohesion,top_disp,bot_disp,tau_max')
interface.dump(0, 0)

nb_time_steps = int(total_duration/time_step)

for s in range(nb_time_steps):
    interface.advanceTimeStep(True)

    if s%10 == 0:
        print(s)
        interface.dump(s, s*time_step)
