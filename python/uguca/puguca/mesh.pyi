"""
mesh submodule
"""
from __future__ import annotations
import numpy
import typing
__all__ = ['BaseMesh', 'DistributedFFTableMesh', 'FFTableMesh', 'FFTableNodalField', 'NodalField', 'SimpleMesh']
class BaseMesh:
    def getLocalCoords(self) -> numpy.ndarray:
        ...
class DistributedFFTableMesh(FFTableMesh):
    pass
class FFTableMesh(BaseMesh):
    pass
class FFTableNodalField(NodalField):
    def __init__(self, name: str = 'unnamed') -> None:
        ...
    def array(self) -> numpy.ndarray:
        ...
class NodalField:
    def __init__(self, mesh: BaseMesh, components: ..., std: ..., std: ..., name: str = 'unnamed') -> None:
        ...
    def array(self) -> numpy.ndarray:
        ...
    def setAllValuesTo(self, arg0: float, arg1: int) -> None:
        ...
class SimpleMesh(DistributedFFTableMesh):
    @typing.overload
    def __init__(self, Lx: float, Nx: int) -> None:
        ...
    @typing.overload
    def __init__(self, Lx: float, Nx: int, Lz: float, Nz: int) -> None:
        ...
    def getNbLocalNodesAlloc(self) -> int:
        ...
