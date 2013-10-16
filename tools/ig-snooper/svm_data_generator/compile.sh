#!/bin/bash

echo "Building libs"
cd ../../common_lib
./clean_all.sh
./compile.sh
echo "Done"

echo "Building svm_data_generator"
cd ../ig-snooper/svm_data_generator
mkdir ./bin
mkdir ./build
cd build
cmake ../src/ -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Debug
make
make install
