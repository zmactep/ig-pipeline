#!/bin/bash

mkdir build
mkdir lib
cd ./build
cmake ../src/ -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Debug
make
make install
