#set title "Runtime of particle update vs. particle number"
set ylabel "Runtime [ms]"
set xlabel "Number of particles"
set key invert reverse Left inside t l
set tics out

plot 'particle_030.dat'    u 1:($4*1000) w l lw 3 ti "LBM 30^3 lattice", \
     'particle_060.dat'    u 1:($4*1000) w l lw 3 ti "LBM 60^3 lattice", \
                        '' u 1:($5*1000) w l lw 3 ti "Update no sprites", \
     'irrlicht_030_01.dat' u 1:($7*1000) w l lw 3 ti "Update 1 sprite", \
                        '' u 1:($5*1000) w l lw 3 ti "Rendering 1 sprite", \
     'irrlicht_030_02.dat' u 1:($7*1000) w l lw 3 lc 9 ti "Update 2 sprites", \
                        '' u 1:($5*1000) w l lw 3 ti "Rendering 2 sprites" #, \
#     'irrlicht_030_04.dat' u 1:($7*1000) w l lw 3 ti "Update 4 sprites", \
#                        '' u 1:($5*1000) w l lw 3 ti "Rendering 4 sprites", \
#     'irrlicht_030_08.dat' u 1:($7*1000) w l lw 3 ti "Update 8 sprites", \
#                        '' u 1:($5*1000) w l lw 3 ti "Rendering 8 sprites"
