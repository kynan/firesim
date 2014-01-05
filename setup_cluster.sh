#!/bin/bash

#module load cmake
#module load boost
#module load irrlicht
#module load gcc-4.2.0

[ -z "$1" ] && echo "usage: $0 <build_type_1> [<build_type_2> ...]" && exit

export PLATFORM=cluster

for arg in $@; do
  rm -rf cluster"$arg"
  mkdir cluster"$arg"
  cmake -E chdir cluster"$arg" cmake -DPLATFORM=cluster -DCMAKE_BUILD_TYPE=$arg -DBOOST_ROOT=/central/boost ..
#make -C cluster"$arg"
done

#mkdir Release
#cmake -E chdir Release cmake -DCMAKE_BUILD_TYPE=Release -DBOOST_ROOT=/usr/local/boost ..
#mkdir Debug
#cmake -E chdir Debug cmake -DCMAKE_BUILD_TYPE=Debug -DBOOST_ROOT=/usr/local/boost ..
#mkdir ReleaseNSmago
#cmake -E chdir ReleaseNSmago cmake -DCMAKE_BUILD_TYPE=ReleaseNSmago -DBOOST_ROOT=/usr/local/boost ..
#mkdir DebugNSmago
#cmake -E chdir DebugNSmago cmake -DCMAKE_BUILD_TYPE=DebugNSmago -DBOOST_ROOT=/usr/local/boost ..
