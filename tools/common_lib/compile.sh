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
