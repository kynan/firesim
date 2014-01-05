#!/bin/bash

for f in dynlight irrlicht_030; do
  for sprites in 1 2 4 8; do
    sed -e s/SPRITES/$sprites/ dynlight > dynlight_${sprites}.cfg
  done
done

for size in 2 3 4 5 6 7 8 9 10 11 12; do
 sed -e s/SIZE/$size/ ldc > ldc_${size}0.cfg
done

for f in ldc_120_*; do
 for rad in 05 10 15 20 25 30 35 40 45 50; do
   sed -e s/RADIUS/$rad/ $f > "$f"_"$rad".cfg
 done
done
