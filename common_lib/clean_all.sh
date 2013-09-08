#!/bin/bash

echo "Cleaning bin, build and lib directories"

cd ./third_party
./clean_all.sh
cd ..

rm -rf ./build
rm -rf ./lib
rm -rf ./bin
echo "Done."
