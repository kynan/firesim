set output 'fmlup.ps'
set title "FMLUP/s with and without Smagorinsky turbulence correction"
set xlabel "Domain size"
set ylabel "FMLUP/s"
set key invert reverse Left outside
set key autotitle columnheader
set yrange [0:2]
plot 'nosmago.dat' u 1:3 w l ti "No Smagorinsky turbulence correction", 'smago.dat' u 1:3 w l ti "Smagorinsky turbulence correction"