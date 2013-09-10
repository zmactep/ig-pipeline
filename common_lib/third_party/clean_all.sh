echo "Removing all libs and build files.."
rm -rf ./lib
rm -rf ./gtest-1.6.0/build
cd ./libsvm-3.17
make clean
cd ./../libsvm-string-3.17
make clean
echo "Done."