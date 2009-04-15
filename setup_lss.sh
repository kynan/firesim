#!/bin/bash

[ -z "$1" ] && echo "usage: $0 <build_type_1> [<build_type_2> ...]" && exit

for arg in $@; do
  if [ -d "$arg" ]; then
    read -p "Directory $arg already exists. Purge and reinitialize? " && [ "$REPLY" != "yes" ] && continue
    rm -rf "$arg"
  fi
  mkdir "$arg"
  cmake -E chdir $arg cmake -DCMAKE_BUILD_TYPE=$arg -DBOOST_ROOT=/usr/local/boost ..
done

#mkdir Release
#cmake -E chdir Release cmake -DCMAKE_BUILD_TYPE=Release -DBOOST_ROOT=/usr/local/boost ..
#mkdir Debug
#cmake -E chdir Debug cmake -DCMAKE_BUILD_TYPE=Debug -DBOOST_ROOT=/usr/local/boost ..
#mkdir ReleaseNSmago
#cmake -E chdir ReleaseNSmago cmake -DCMAKE_BUILD_TYPE=ReleaseNSmago -DBOOST_ROOT=/usr/local/boost ..
#mkdir DebugNSmago
#cmake -E chdir DebugNSmago cmake -DCMAKE_BUILD_TYPE=DebugNSmago -DBOOST_ROOT=/usr/local/boost ..
