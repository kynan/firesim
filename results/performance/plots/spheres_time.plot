#set title "Runtime of obstacle computation"
set xlabel "Obstacle boundary cells"
set ylabel "Runtime [ms]"
set key invert reverse Left inside t l

plot 'spheres.dat' u 1:2 w l lw 3 ti "staircase approximation", \
                '' u 1:3 w l lw 3 ti "curved boundary", \
                '' u 1:4 w l lw 3 ti "moving boundary"