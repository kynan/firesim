#set title "Share of collide-stream and boundary computations"
set xlabel "% of computation time"
set ylabel "Domain size"
set key invert reverse Left outside
set key autotitle columnheader
set yrange [0:100]
set auto x
unset xtics
set xtics nomirror
set style data histogram
set style histogram rowstacked
set style fill solid border -1
set boxwidth 0.75

#plot 'share.nosmago.dat' using ($2*100):xtic(1), '' using ($4*100), '' using ($5*100)
plot newhistogram "Without Smagorinsky turbulence correction", 'nosmago.dat' using ($4*100):xtic(1) ti column(4), '' using ($6*100) ti column(6), '' using ($7*100) ti column(7), \
     newhistogram "With Smagorinsky turbulence correction", 'smago.dat' using ($4*100):xtic(1) ti column(4), '' using ($6*100) ti column(6), '' using ($7*100) ti column(7)