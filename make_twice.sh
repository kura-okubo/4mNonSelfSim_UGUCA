#!bin/sh


export PATH=$PATH:/opt/openmpi-5.0.2/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/openmpi-5.0.2/bin

cmake ..
make -j8
cmake ..
make -j8
#ctest
