#!/bin/bash

rm -rf ./build
rm -rf ./lib
rm -rf ./bin
mkdir build
mkdir bin
cd ./build
cmake ../src/ -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Debug
make
make install
