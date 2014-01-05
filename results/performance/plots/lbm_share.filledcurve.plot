#set title "Share of collide-stream and boundary computations for LBM step"
set xlabel "% of computation time"
set ylabel "Domain size"
set key invert reverse Left outside
set key autotitle columnheader
set yrange[0:100]
set auto x
#set style data lines

plot 'nosmago.dat' using 1:($4*100) ti column(4) w filledcurve x1, '' using 1:($4*100):(($4+$6)*100) ti column(6) w filledcu, '' using 1:(($4+$6)*100):(($4+$6+$7)*100) ti column(7) w filledcu
set output 'lbm_share.filledcurve.smago.ps'
plot 'smago.dat' using 1:($4*100) ti column(4) w filledcurve x1, '' using 1:($4*100):(($4+$6)*100) ti column(6) w filledcu, '' using 1:(($4+$6)*100):(($4+$6+$7)*100) ti column(7) w filledcu
