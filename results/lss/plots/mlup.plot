set output 'mlup.ps'
set title "MLUP/s with and without Smagorinsky turbulence correction"
set xlabel "Domain size"
set ylabel "MLUP/s"
set key invert reverse Left outside
set key autotitle columnheader
set yrange [0:2]
plot 'nosmago.dat' u 1:2 w l ti "No Smagorinsky turbulence correction", 'smago.dat' u 1:2 w l ti "Smagorinsky turbulence correction"