import uguca as ug
import numpy as np

class CustomLaw(ug.InterfaceLaw):
    def __init__(self, mesh, tau_max, delta_c):
        super().__init__(mesh)
        self.tau_max = np.full(mesh.getNbLocalNodesAlloc(), tau_max)
        self.delta_c = np.full(mesh.getNbLocalNodesAlloc(), delta_c)
        self.gap = ug.NodalField(mesh)


    def computeCohesiveForces(self, cohesion, predicting):
        cohesion = self.getCohesion()
        interface = self.getInterface()

        # find current gap
        interface.computeGap(self.gap, predicting)
        gaps_array = np.array((self.gap.component(0), self.gap.component(1)))
        gap_norm = np.linalg.norm(gaps_array, axis=0)

        # find forces needed to close normal gap
        cohesion_force = cohesion.component(1)
        interface.closingNormalGapForce(cohesion_force, predicting)

        # find force needed to maintain shear gap
        interface.maintainShearGapForce(cohesion)
        tau_shear = np.abs(cohesion.component(0))

        alpha_field = np.zeros(self.getMesh().getNbLocalNodesAlloc())
        strength = self.tau_max[:] * np.maximum(0., 1.-gap_norm/self.delta_c)

        cohesion_force[:] = np.minimum(cohesion_force, strength)
        alpha_field = np.minimum(1, np.abs(strength/tau_shear))

        cohesion.component(0)[:] *= alpha_field[:]


length = 1.
nb_elements = 1024

mesh = ug.SimpleMesh(Lx=length, Nx=nb_elements)
coords_x = mesh.getLocalCoords(0)

law = CustomLaw(mesh, tau_max=3.5e6, delta_c=2e-5)

top_mat = ug.Material(7e9, 0.33, 2000)
top_mat.readPrecomputedKernels()

bot_mat = ug.Material(5e9, 0.25, 1200)
bot_mat.readPrecomputedKernels()


interface = ug.BimatInterface(mesh, top_mat, bot_mat, law)

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

law.tau_max[indexes] = 0.


interface.initDump('fracture_2d_example', '.')
interface.registerDumpFields('cohesion_0,cohesion_1,top_disp_0,top_disp_1,bot_disp_0,bot_disp_1')
interface.dump(0, 0)


nb_time_steps = int(total_duration/time_step)

for s in range(nb_time_steps):
    interface.advanceTimeStep(True)

    if s%10 == 0:
        print(s)
        interface.dump(s, s*time_step)
