#set title "Runtime in ms of a LBM timestep with and without Smagorinsky turbulence correction"
set xlabel "Domain size"
set ylabel "Runtime [ms]"
set key reverse Left inside t l

plot 'smago.dat'   u 1:8 w l lw 3 ti "Smagorinsky on, complete", \
     'smago.dat'   u 1:9 w l lw 3 ti "Smagorinsky on, collide-stream", \
     'nosmago.dat' u 1:8 w l lw 3 ti "Smagorinsky off, complete", \
     'nosmago.dat' u 1:9 w l lw 3 ti "Smagorinsky off, collide-stream"