"""
half_space submodule
"""

from __future__ import annotations
import puguca.mesh

__all__ = ["HalfSpace", "Material"]

class HalfSpace:
    def __init__(
        self,
        material: Material,
        mesh: puguca.mesh.FFTableMesh,
        side_factor: int,
        components: ...,
        std: ...,
        std: ...,
        name: str = "half_space",
    ) -> None: ...
    def getDisp(self, predicting: bool = False) -> puguca.mesh.FFTableNodalField: ...
    def getVelo(self, predicting: bool = False) -> puguca.mesh.NodalField: ...

class Material:
    def __init__(
        self, E: float, nu: float, rho: float, pstress: bool = False
    ) -> None: ...
    def readPrecomputedKernels(
        self,
        path: str = "/home/florez/work/projects/continuously-weakening-friction/uguca/kernels/",
    ) -> None: ...
