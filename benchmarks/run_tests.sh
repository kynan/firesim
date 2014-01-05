#!/bin/bash

[ -z "$PLATFORM" ] && echo "Need to set PLATFORM environment variable" && exit

BINS=( "../bin/$PLATFORM/Release/lbm_reference" "../bin/$PLATFORM/ReleaseNSmago/lbm_reference" )
TESTS=`ls *.cfg`

for bin in $BINS; do
  for tst in $TESTS; do
    $bin $tst
  done
done
