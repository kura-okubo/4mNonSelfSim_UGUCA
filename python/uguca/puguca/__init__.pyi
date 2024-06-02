"""
python module of uguca
"""
from __future__ import annotations
import typing
from . import half_space
from . import interface
from . import io
from . import mesh
__all__ = ['SolverMethod', 'SpatialDirection', 'adaptive', 'dynamic', 'half_space', 'interface', 'io', 'mesh', 'quasi_dynamic', 'spatial_dir_count', 'static', 'x', 'y', 'z']
class SolverMethod:
    """
    Members:
    
      static
    
      quasi_dynamic
    
      dynamic
    
      adaptive
    """
    __members__: typing.ClassVar[dict[str, SolverMethod]]  # value = {'static': <SolverMethod.static: 0>, 'quasi_dynamic': <SolverMethod.quasi_dynamic: 1>, 'dynamic': <SolverMethod.dynamic: 2>, 'adaptive': <SolverMethod.adaptive: 3>}
    adaptive: typing.ClassVar[SolverMethod]  # value = <SolverMethod.adaptive: 3>
    dynamic: typing.ClassVar[SolverMethod]  # value = <SolverMethod.dynamic: 2>
    quasi_dynamic: typing.ClassVar[SolverMethod]  # value = <SolverMethod.quasi_dynamic: 1>
    static: typing.ClassVar[SolverMethod]  # value = <SolverMethod.static: 0>
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class SpatialDirection:
    """
    Members:
    
      x
    
      y
    
      z
    
      spatial_dir_count
    """
    __members__: typing.ClassVar[dict[str, SpatialDirection]]  # value = {'x': <SpatialDirection.x: 0>, 'y': <SpatialDirection.y: 1>, 'z': <SpatialDirection.z: 2>, 'spatial_dir_count': <SpatialDirection.spatial_dir_count: 3>}
    spatial_dir_count: typing.ClassVar[SpatialDirection]  # value = <SpatialDirection.spatial_dir_count: 3>
    x: typing.ClassVar[SpatialDirection]  # value = <SpatialDirection.x: 0>
    y: typing.ClassVar[SpatialDirection]  # value = <SpatialDirection.y: 1>
    z: typing.ClassVar[SpatialDirection]  # value = <SpatialDirection.z: 2>
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
adaptive: SolverMethod  # value = <SolverMethod.adaptive: 3>
dynamic: SolverMethod  # value = <SolverMethod.dynamic: 2>
quasi_dynamic: SolverMethod  # value = <SolverMethod.quasi_dynamic: 1>
spatial_dir_count: SpatialDirection  # value = <SpatialDirection.spatial_dir_count: 3>
static: SolverMethod  # value = <SolverMethod.static: 0>
x: SpatialDirection  # value = <SpatialDirection.x: 0>
y: SpatialDirection  # value = <SpatialDirection.y: 1>
z: SpatialDirection  # value = <SpatialDirection.z: 2>
