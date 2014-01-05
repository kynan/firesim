#set title "MLUP/s with and without Smagorinsky turbulence correction"
set xlabel "Domain size"
set ylabel "MLUP/s"
set key reverse Left inside b l
set yrange [0:2.5]

plot 'nosmago.dat' u 1:2 w l lw 3 ti "Smagorinsky off, all cells", \
     'nosmago.dat' u 1:3 w l lw 3 ti "Smagorinsky off, fluid cells", \
     'smago.dat' u 1:2 w l lw 3 ti "Smagorinsky on, all cells", \
     'smago.dat' u 1:3 w l lw 3 ti "Smagorinsky on, fluid cells"