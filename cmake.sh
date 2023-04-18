#!bin/bash
rm -r build/*
mkdir build
cd src
make clean
cd -
cd build
cmake ..
make
cd -
