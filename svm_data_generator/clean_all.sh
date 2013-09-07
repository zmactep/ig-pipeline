#!/bin/bash

echo "Cleaning common_lib"
cd ../common_lib/
./clean_all.sh

cd ../svm_data_generator
echo "Cleaning build and bin directories"
rm -rf ./build
rm -rf ./bin
echo "Done"
