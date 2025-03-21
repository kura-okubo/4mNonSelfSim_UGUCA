#!bin/sh

rm ./CMakeCache.txt
rm -r ./example_4mNonSelfSim
rm -r ./src
rm -r ./test

cmake -DCMAKE_C_COMPILER=`which mpicc` -DCMAKE_CXX_COMPILER=`which mpicxx` ..
make -j8
cmake -DCMAKE_C_COMPILER=`which mpicc` -DCMAKE_CXX_COMPILER=`which mpicxx` ..
make -j8
