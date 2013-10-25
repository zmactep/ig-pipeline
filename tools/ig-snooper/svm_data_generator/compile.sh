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
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS="-std=c++11 -stdlib=libc++ -lc++abi" ../src/ 
make 
make install
