#set title "Runtime of particle emission, LBM, particle update and rendering computations"
set ylabel "Runtime [ms]"
set xlabel "Number of particles"
set key reverse invert Left inside t l
set tics out

# set output 'irrlicht_time_1.eps'
# plot 'irrlicht_030_01.dat' using 1:4 ti column(4) w filledcurve x1, \
#                         '' using 1:4:($4+$5) ti column(5) w filledcu, \
#                         '' using 1:($4+$5):($4+$5+$6) ti column(6) w filledcu, \
#                         '' using 1:($4+$5+$6):($4+$5+$6+$7) ti column(7) w filledcu
#
# set output 'irrlicht_time_2.eps'
# plot 'irrlicht_030_02.dat' using 1:4 ti column(4) w filledcurve x1, \
#                         '' using 1:4:($4+$5) ti column(5) w filledcu, \
#                         '' using 1:($4+$5):($4+$5+$6) ti column(6) w filledcu, \
#                         '' using 1:($4+$5+$6):($4+$5+$6+$7) ti column(7) w filledcu
#
# set output 'irrlicht_time_4.eps'
# plot 'irrlicht_030_04.dat' using 1:4 ti column(4) w filledcurve x1, \
#                         '' using 1:4:($4+$5) ti column(5) w filledcu, \
#                         '' using 1:($4+$5):($4+$5+$6) ti column(6) w filledcu, \
#                         '' using 1:($4+$5+$6):($4+$5+$6+$7) ti column(7) w filledcu
#
# set output 'irrlicht_time_8.eps'
# plot 'irrlicht_030_08.dat' using 1:4 ti column(4) w filledcurve x1, \
#                         '' using 1:4:($4+$5) ti column(5) w filledcu, \
#                         '' using 1:($4+$5):($4+$5+$6) ti column(6) w filledcu, \
#                         '' using 1:($4+$5+$6):($4+$5+$6+$7) ti column(7) w filledcu

set output 'irrlicht_time_1.eps'
plot "<awk '{print $1, $4*1000, ($4+$5)*1000, ($4+$5+$6)*1000, ($4+$5+$6+$7)*1000}' irrlicht_030_01.dat" \
using 1:2 w filledcurve x1 ti "Particle emission", \
                        '' using 1:2:3 w filledcu ti "Rendering", \
                        '' using 1:3:4 w filledcu ti "LBM", \
                        '' using 1:4:5 w filledcu ti "Particle update"

set output 'irrlicht_time_2.eps'
plot "<awk '{print $1, $4*1000, ($4+$5)*1000, ($4+$5+$6)*1000, ($4+$5+$6+$7)*1000}' irrlicht_030_02.dat" \
using 1:2 w filledcurve x1 ti "Particle emission", \
                        '' using 1:2:3 w filledcu ti "Rendering", \
                        '' using 1:3:4 w filledcu ti "LBM", \
                        '' using 1:4:5 w filledcu ti "Particle update"

set output 'irrlicht_time_4.eps'
plot "<awk '{print $1, $4*1000, ($4+$5)*1000, ($4+$5+$6)*1000, ($4+$5+$6+$7)*1000}' irrlicht_030_04.dat" \
using 1:2 w filledcurve x1 ti "Particle emission", \
                        '' using 1:2:3 w filledcu ti "Rendering", \
                        '' using 1:3:4 w filledcu ti "LBM", \
                        '' using 1:4:5 w filledcu ti "Particle update"

set output 'irrlicht_time_8.eps'
plot "<awk '{print $1, $4*1000, ($4+$5)*1000, ($4+$5+$6)*1000, ($4+$5+$6+$7)*1000}' irrlicht_030_08.dat" \
using 1:2 w filledcurve x1 ti "Particle emission", \
                        '' using 1:2:3 w filledcu ti "Rendering", \
                        '' using 1:3:4 w filledcu ti "LBM", \
                        '' using 1:4:5 w filledcu ti "Particle update"
