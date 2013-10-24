#!/bin/bash

cd ./third_party
./clean_all.sh
./compile.sh
cd ..
rm -rf ./build
rm -rf ./lib
rm -rf ./bin
mkdir build
mkdir lib
mkdir bin
cd ./build
cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS="-std=c++11 -stdlib=libc++" ../src/
make
make install
