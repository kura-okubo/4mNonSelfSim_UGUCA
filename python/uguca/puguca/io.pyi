"""
io submodule
"""
from __future__ import annotations
import typing
__all__ = ['BaseIO', 'Dumper', 'Format', 'Restart']
class BaseIO:
    def dump(self, arg0: int, arg1: float) -> None:
        ...
    def initIO(self, arg0: str, arg1: str, arg2: Format) -> None:
        ...
    def load(self, arg0: int) -> None:
        ...
class Dumper(BaseIO):
    pass
class Format:
    """
    Members:
    
      ASCII
    
      CSV
    
      Binary
    """
    ASCII: typing.ClassVar[Format]  # value = <Format.ASCII: 0>
    Binary: typing.ClassVar[Format]  # value = <Format.Binary: 2>
    CSV: typing.ClassVar[Format]  # value = <Format.CSV: 1>
    __members__: typing.ClassVar[dict[str, Format]]  # value = {'ASCII': <Format.ASCII: 0>, 'CSV': <Format.CSV: 1>, 'Binary': <Format.Binary: 2>}
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
class Restart(BaseIO):
    def __init__(self, bname: str, path: str, format: Format = ...) -> None:
        ...
