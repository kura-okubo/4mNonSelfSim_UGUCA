# uguca python interface

Uses pybind11 to create a python interface for the uguca library.

Uses `pybind11-stubgen` to generate static type hints for the python interface.

## Installation

1. Compile the uguca library.
2. `cd build/python`
3. `pip install .` (or `pip install -e . --config-settings editable_mode=compat` for editable mode)
