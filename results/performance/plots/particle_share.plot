#set title "Share of particle emission, LBM and particle update computations"
set xlabel "% of computation time"
set ylabel "Number of particles"
set key invert reverse Left outside
set key autotitle columnheader
set yrange[0:100]
set auto x
#set style data lines

set output 'particle_share_030.eps'
plot 'particle_030.dat' using 1:($3*100/$2) ti column(3) w filledcurve x1, \
                        '' using 1:($3*100/$2):(($3+$4)*100/$2) ti column(4) w filledcu, \
                        '' using 1:(($3+$4)*100/$2):(($3+$4+$5)*100/$2) ti column(5) w filledcu

set output 'particle_share_060.eps'
plot 'particle_060.dat' using 1:($3*100/$2) ti column(3) w filledcurve x1, \
                        '' using 1:($3*100/$2):(($3+$4)*100/$2) ti column(4) w filledcu, \
                        '' using 1:(($3+$4)*100/$2):(($3+$4+$5)*100/$2) ti column(5) w filledcu

set output 'particle_share_120.eps'
plot 'particle_120.dat' using 1:($3*100/$2) ti column(3) w filledcurve x1, \
                        '' using 1:($3*100/$2):(($3+$4)*100/$2) ti column(4) w filledcu, \
                        '' using 1:(($3+$4)*100/$2):(($3+$4+$5)*100/$2) ti column(5) w filledcu