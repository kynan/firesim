#1/bin/bash

rm smago.dat nosmago.dat spheres.dat

#for i in `ls ldc_???.smago`; do
#  awk '/Domain/ { sz = $4 } /MLUP/ { mlup = $3; fmlup = $4 } END { print sz, mlup, fmlup }' $i >> mlup.smago.dat
#done
#for i in `ls ldc_???.nosmago`; do
#  awk '/Domain/ { sz = $4 } /MLUP/ { mlup = $3 } END { print sz, mlup }' $i>> mlup.nosmago.dat
#done
#for i in `ls ldc_???.smago`; do
#  awk '/Domain/ { sz = $4 } /MLUP/ { mlup = $4 } END { print sz, mlup }' $i >> fmlup.smago.dat
#done
#for i in `ls ldc_???.nosmago`; do
#  awk '/Domain/ { sz = $4 } /MLUP/ { mlup = $4 } END { print sz, mlup }' $i >> fmlup.nosmago.dat
#done

echo '"Domain size" "MLUP/s" "FMLUP/s" "Fluid" "Boundary total" "Boundary no-slip" "Boundary acceleration" "Runtime total [ms]" "Runtime fluid [ms]"' >> smago.dat
for i in ldc_???.smago; do
  awk '/Domain/ { sz = $4 } /MLUP/ { mlup = $3; fmlup = $4 } /share/ { f = $5; b = $6; n = $7; a = $8 } /\[s\]/ { r = $4; rf = $5} END { print sz, mlup, fmlup, f, b, n, a, r*10, rf*10 }' $i >> smago.dat
done

echo '"Domain size" "MLUP/s" "FMLUP/s" "Fluid" "Boundary total" "Boundary no-slip" "Boundary acceleration" "Runtime total [ms]" "Runtime fluid [ms]"' >> nosmago.dat
for i in ldc_???.nosmago; do
  awk '/Domain/ { sz = $4 } /MLUP/ { mlup = $3; fmlup = $4 } /share/ { f = $5; b = $6; n = $7; a = $8 } /\[s\]/ { r = $4; rf = $5} END { print sz, mlup, fmlup, f, b, n, a, r*10, rf*10 }' $i >> nosmago.dat
done

echo '"Number of cells" "Runtime staircase [ms]" "Runtime curved [ms]" "Runtime moving [ms]" "Share staircase" "Share curved" "Share moving" "MLUP/s staircase" "MLUP/s curved" "MLUP/s moving"' >> spheres.dat
for i in 05 10 15 20 25 30 35 40 45 50; do
  cat ldc_120_sphere_moving_$i.nosmago ldc_120_sphere_staircase_$i.nosmago ldc_120_sphere_stationary_$i.nosmago | awk '/\[s\]/ { if ($12 > 0) rc = $12; if ($13 > 0) rs = $13; if ($14 > 0) rm = $14 } /cells/ { if ($14 > 0) c = $14 } /share/ { if ($12 > 0) sc = $12; if ($13 > 0) ss = $13; if ($14 > 0) sm = $14 } END { print c, rs*10, rc*10, rm*10, ss, sc, sm, c/(rs*10000), c/(rc*10000), c/(rm*10000) }' >> spheres.dat
done