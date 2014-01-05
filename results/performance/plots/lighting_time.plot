#set title "Rendering time for dynamic lighting turned on and off"
set ylabel "Rendering time [ms]"
set xlabel "Number of sprites"
set key invert reverse Left inside t l
set tics out

plot 'lighting_00.dat' using 1:($5*1000) w l lw 3 ti "No dynamic lighting", \
     'lighting_01.dat' using 1:($5*1000) w l lw 3 ti "Dynamic lighting (1 Object)", \
     'lighting_02.dat' using 1:($5*1000) w l lw 3 ti "Dynamic lighting (2 Objects)", \
     'lighting_04.dat' using 1:($5*1000) w l lw 3 ti "Dynamic lighting (4 Objects)", \
     'lighting_08.dat' using 1:($5*1000) w l lw 3 ti "Dynamic lighting (8 Objects)"