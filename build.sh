#!/bin/bash

source $HOME/setup_lcg_cuda.sh

cd acts/build

cmake ..
sed -i 's/c++17/c++20/g' build.ninja
ninja -j64
