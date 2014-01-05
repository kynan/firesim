#set title "MLUP/s of obstacle computation"
set xlabel "Obstacle boundary cells"
set ylabel "MLUP/s"
set key reverse Left inside t l
set yrange [0:7.5]

plot 'spheres.dat' u 1:8  w l lw 3 ti "staircase approximation", \
                '' u 1:9  w l lw 3 ti "curved boundary", \
                '' u 1:10 w l lw 3 ti "moving boundary"