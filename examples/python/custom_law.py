import uguca as ug
import numpy as np

class CustomLaw(ug.InterfaceLaw):
    def __init__(self, mesh):
        super().__init__(mesh)
        self.tau_max = np.zeros(mesh.getNbLocalNodesAlloc(), dtype=float)
        self.delta_c = np.zeros(mesh.getNbLocalNodesAlloc(), dtype=float)
        self.gap_norm = np.zeros(mesh.getNbLocalNodesAlloc(), dtype=float)
        
    def computeCohesiveForces(self, predicting):
        
        cohesion = self.getCohesion()
        print(cohesion.component(1))

        interface = self.getInterface()
        load = interface.getLoad()
        print(load.component(1))
        
        gap = ug.NodalField(super().getMesh())
        
        
        py_cohesion = np.zeros(super().getMesh().getNbLocalNodesAlloc())
        interface.closingNormalGapForce(py_cohesion, predicting)
        print(py_cohesion)

        

        
length = 1.
nb_elements = 1024

mesh = ug.SimpleMesh(Lx=length, Nx=nb_elements)
law = CustomLaw(mesh)

top_mat = ug.Material(7e9, 0.33, 2000)
top_mat.readPrecomputedKernels();

bot_mat = ug.Material(5e9, 0.25, 1200)
bot_mat.readPrecomputedKernels();


interface = ug.BimatInterface(mesh, top_mat, bot_mat, law)


loads = interface.getLoad()
loads.component(0)[:] = 2e6
loads.component(1)[:] = 1e6

total_duration = 2.6e-4
time_factor = 0.4
time_step = time_factor*interface.getStableTimeStep()
interface.setTimeStep(time_step)
interface.init(False)
