set title "Runtime of particle emission, LBM and particle update computations"
set ylabel "Runtime in sec"
set xlabel "Number of particles"
set key invert reverse Left outside
set key autotitle columnheader
#set yrange[0:100]
set auto x
#set style data lines
set output 'particle_time_030.ps'
plot 'particle_030.dat' using 1:3 ti column(3) w filledcurve x1, \
                        '' using 1:3:($3+$4) ti column(4) w filledcu, \
                        '' using 1:($3+$4):($3+$4+$5) ti column(5) w filledcu
set output 'particle_time_060.ps'
plot 'particle_060.dat' using 1:3 ti column(3) w filledcurve x1, \
                        '' using 1:3:($3+$4) ti column(4) w filledcu, \
                        '' using 1:($3+$4):($3+$4+$5) ti column(5) w filledcu
set output 'particle_time_120.ps'
plot 'particle_120.dat' using 1:3 ti column(3) w filledcurve x1, \
                        '' using 1:3:($3+$4) ti column(4) w filledcu, \
                        '' using 1:($3+$4):($3+$4+$5) ti column(5) w filledcu