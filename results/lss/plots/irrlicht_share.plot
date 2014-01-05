set title "Share of particle emission, LBM, particle update and rendering computations"
set xlabel "% of computation time"
set ylabel "Number of particles"
set key invert reverse Left outside
set key autotitle columnheader
set yrange[0:100]
set auto x
#set style data lines
set output 'irrlicht_share_1.ps'
plot 'irrlicht_030_01.dat' using 1:($4*100/$3) ti column(4) w filledcurve x1, \
                        '' using 1:($4*100/$3):(($4+$5)*100/$3) ti column(5) w filledcu, \
                        '' using 1:(($4+$5)*100/$3):(($4+$5+$6)*100/$3) ti column(6) w filledcu, \
                        '' using 1:(($4+$5+$6)*100/$3):(($4+$5+$6+$7)*100/$3) ti column(7) w filledcu
set output 'irrlicht_share_2.ps'
plot 'irrlicht_030_02.dat' using 1:($4*100/$3) ti column(4) w filledcurve x1, \
                        '' using 1:($4*100/$3):(($4+$5)*100/$3) ti column(5) w filledcu, \
                        '' using 1:(($4+$5)*100/$3):(($4+$5+$6)*100/$3) ti column(6) w filledcu, \
                        '' using 1:(($4+$5+$6)*100/$3):(($4+$5+$6+$7)*100/$3) ti column(7) w filledcu
set output 'irrlicht_share_4.ps'
plot 'irrlicht_030_04.dat' using 1:($4*100/$3) ti column(4) w filledcurve x1, \
                        '' using 1:($4*100/$3):(($4+$5)*100/$3) ti column(5) w filledcu, \
                        '' using 1:(($4+$5)*100/$3):(($4+$5+$6)*100/$3) ti column(6) w filledcu, \
                        '' using 1:(($4+$5+$6)*100/$3):(($4+$5+$6+$7)*100/$3) ti column(7) w filledcu
set output 'irrlicht_share_8.ps'
plot 'irrlicht_030_08.dat' using 1:($4*100/$3) ti column(4) w filledcurve x1, \
                        '' using 1:($4*100/$3):(($4+$5)*100/$3) ti column(5) w filledcu, \
                        '' using 1:(($4+$5)*100/$3):(($4+$5+$6)*100/$3) ti column(6) w filledcu, \
                        '' using 1:(($4+$5+$6)*100/$3):(($4+$5+$6+$7)*100/$3) ti column(7) w filledcu
