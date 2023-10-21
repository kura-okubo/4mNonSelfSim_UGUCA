__copyright__ = (
    "Copyright (©) 2021 ETH Zurich (David S. Kammer)"
    "Laboratory (CMBM - Computational Mechanics of Building Materials)"
)

from . import puguca as _puguca

private_keys = set(dir(_puguca)) - set(dir())

for k in private_keys:
    globals()[k] = getattr(_puguca, k)
