#set title "Share of particle emission, LBM, particle update and rendering computations"
set ylabel "% of computation time"
set xlabel "Number of particles"
set key horizontal Left reverse invert outside c b
set tics out
set yrange[0:100]

set output 'irrlicht_share_1.eps'
plot 'irrlicht_030_01.dat' using 1:($4*100/$3) ti "Particle emission" w filledcurve x1, \
                        '' using 1:($4*100/$3):(($4+$5)*100/$3) ti "Rendering" w filledcu, \
                        '' using 1:(($4+$5)*100/$3):(($4+$5+$6)*100/$3) ti "LBM" w filledcu, \
                        '' using 1:(($4+$5+$6)*100/$3):(($4+$5+$6+$7)*100/$3) ti "Particle update" w filledcu

set output 'irrlicht_share_2.eps'
plot 'irrlicht_030_02.dat' using 1:($4*100/$3) ti "Particle emission" w filledcurve x1, \
                        '' using 1:($4*100/$3):(($4+$5)*100/$3) ti "Rendering" w filledcu, \
                        '' using 1:(($4+$5)*100/$3):(($4+$5+$6)*100/$3) ti "LBM" w filledcu, \
                        '' using 1:(($4+$5+$6)*100/$3):(($4+$5+$6+$7)*100/$3) ti "Particle update" w filledcu

set output 'irrlicht_share_4.eps'
plot 'irrlicht_030_04.dat' using 1:($4*100/$3) ti "Particle emission" w filledcurve x1, \
                        '' using 1:($4*100/$3):(($4+$5)*100/$3) ti "Rendering" w filledcu, \
                        '' using 1:(($4+$5)*100/$3):(($4+$5+$6)*100/$3) ti "LBM" w filledcu, \
                        '' using 1:(($4+$5+$6)*100/$3):(($4+$5+$6+$7)*100/$3) ti "Particle update" w filledcu

set output 'irrlicht_share_8.eps'
plot 'irrlicht_030_08.dat' using 1:($4*100/$3) ti "Particle emission" w filledcurve x1, \
                        '' using 1:($4*100/$3):(($4+$5)*100/$3) ti "Rendering" w filledcu, \
                        '' using 1:(($4+$5)*100/$3):(($4+$5+$6)*100/$3) ti "LBM" w filledcu, \
                        '' using 1:(($4+$5+$6)*100/$3):(($4+$5+$6+$7)*100/$3) ti "Particle update" w filledcu
