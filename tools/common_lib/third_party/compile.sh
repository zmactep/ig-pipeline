echo "Cleaning lib folder"
rm -rf ./lib
mkdir ./lib

echo "Building gtest..."
cd ./gtest-1.6.0
rm -rf ./build
mkdir ./build
cd ./build
cmake -DCMAKE_CXX_FLAGS="-std=c++11" ..
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
