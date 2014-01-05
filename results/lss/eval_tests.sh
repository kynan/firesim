#1/bin/bash

rm *.dat

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

echo '"Domain size" "MLUP/s" "FMLUP/s" "Fluid" "Boundary total" "Boundary no-slip" "Boundary acceleration"' >> smago.dat
for i in `ls ldc_???.smago`; do
  awk '/Domain/ { sz = $4 } /MLUP/ { mlup = $3; fmlup = $4 } /share/ { f = $5; b = $6; n = $7; a = $8 } END { print sz, mlup, fmlup, f, b, n, a }' $i >> smago.dat
done

echo '"Domain size" "MLUP/s" "FMLUP/s" "Fluid" "Boundary total" "Boundary no-slip" "Boundary acceleration"' >> nosmago.dat
for i in `ls ldc_???.nosmago`; do
  awk '/Domain/ { sz = $4 } /MLUP/ { mlup = $3; fmlup = $4 } /share/ { f = $5; b = $6; n = $7; a = $8 } END { print sz, mlup, fmlup, f, b, n, a }' $i >> nosmago.dat
done

for i in *_staircase_