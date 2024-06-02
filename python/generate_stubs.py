#!/bin/python3
"""Generate stubs for the uguca module.

Dependencies:
    - pybind11-stubgen
    - black
"""

import os
import subprocess

PATH_TO_COMPILED_EXTENSION = "build/python/uguca"
EXTENSION_NAME = "puguca"
PACKAGE_NAME = "uguca"

# Generate stubs for the uguca module
print("Generating stubs for the uguca module")

# change directory to root of the project where the .git folder is located
current_dir = os.getcwd()
while not os.path.exists(".git"):
    os.chdir("..")
    if os.getcwd() == "/":
        raise RuntimeError("Error: Could not find the root of the project")

os.chdir(PATH_TO_COMPILED_EXTENSION)

subprocess.run(
    ["pybind11-stubgen", "-o", os.path.join(current_dir, PACKAGE_NAME), EXTENSION_NAME],
    check=True,
)

# change back to the original directory
os.chdir(current_dir)
