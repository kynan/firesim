set title "Runtime of particle emission, LBM, particle update and rendering computations"
set ylabel "Runtime in sec"
set xlabel "Number of particles"
set key invert reverse Left outside
set key autotitle columnheader
#set yrange[0:100]
set auto x
#set style data lines
set output 'irrlicht_time_1.ps'
plot 'irrlicht_030_01.dat' using 1:4 ti column(4) w filledcurve x1, \
                        '' using 1:4:($4+$5) ti column(5) w filledcu, \
                        '' using 1:($4+$5):($4+$5+$6) ti column(6) w filledcu, \
                        '' using 1:($4+$5+$6):($4+$5+$6+$7) ti column(7) w filledcu
set output 'irrlicht_time_2.ps'
plot 'irrlicht_030_02.dat' using 1:4 ti column(4) w filledcurve x1, \
                        '' using 1:4:($4+$5) ti column(5) w filledcu, \
                        '' using 1:($4+$5):($4+$5+$6) ti column(6) w filledcu, \
                        '' using 1:($4+$5+$6):($4+$5+$6+$7) ti column(7) w filledcu
set output 'irrlicht_time_4.ps'
plot 'irrlicht_030_04.dat' using 1:4 ti column(4) w filledcurve x1, \
                        '' using 1:4:($4+$5) ti column(5) w filledcu, \
                        '' using 1:($4+$5):($4+$5+$6) ti column(6) w filledcu, \
                        '' using 1:($4+$5+$6):($4+$5+$6+$7) ti column(7) w filledcu
set output 'irrlicht_time_8.ps'
plot 'irrlicht_030_08.dat' using 1:4 ti column(4) w filledcurve x1, \
                        '' using 1:4:($4+$5) ti column(5) w filledcu, \
                        '' using 1:($4+$5):($4+$5+$6) ti column(6) w filledcu, \
                        '' using 1:($4+$5+$6):($4+$5+$6+$7) ti column(7) w filledcu
