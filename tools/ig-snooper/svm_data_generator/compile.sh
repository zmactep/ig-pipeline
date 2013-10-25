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

OS=$(uname -rs)
if [[ "$OS" == "Darwin 13.0.0" ]]
then
  echo "You are using OS X Mavericks"
  cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS="-std=c++11 -stdlib=libc++" ../src/
elif [[ "$OS" == *ARCH* ]]
then 
  echo "You are using ARCH"
  cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS="-std=c++11 -stdlib=libc++ -lc++abi" ../src/
else
  echo "You are using OS X Mountain Lion or other"
  cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS="-std=c++11" ../src/
fi  

make 
make install
