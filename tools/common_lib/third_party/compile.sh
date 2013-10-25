echo "Cleaning lib folder"
rm -rf ./lib
mkdir ./lib

echo "Building gtest..."
cd ./gtest-1.7.0
rm -rf ./build
mkdir ./build
cd ./build

OS=$(uname -rs)
if [[ "$OS" == "Darwin 13.0.0" ]]
then
  echo "You are using OS X Mavericks"
  cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS="-std=c++11 -stdlib=libc++" ..
elif [[ "$OS" == *ARCH* ]]
then 
  echo "You are using ARCH"
  cmake -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS="-std=c++11 -stdlib=libc++ -lc++abi" ..
else
  echo "You are using OS X Mountain Lion or other"
  cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_FLAGS="-std=c++11" ..
fi  

make
cp *.a ../../lib 
cd ..
rm -rf ./build
cd ..
echo "Done"

echo "Building libsvm"
cd ./libsvm-3.17
make
cd ./../libsvm-string-3.17
make
echo "Done"
