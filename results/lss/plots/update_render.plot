set title "Runtime of particle update vs. particle number"
set ylabel "Runtime in sec"
set xlabel "Number of particles"
set key invert reverse Left outside
set key autotitle columnheader
#set yrange[0:100]
set auto x
#set style data lines
set output 'update_render.ps'
plot 'irrlicht_030_01.dat' using 1:7 ti column(7) w l, \
                        '' using 1:5 ti column(5) w l