#
#
set title "Karthikeyan et al., Fig 5 Redlich-Kister fit"
set xlabel "x in LixCoO2"
# set xlabel "DeltaG_f (cal/gmol)"
set ylabel "Open circuit voltage vs Li(S)"
set xtic auto
set ytic auto
# set key top right
# set key box
# set key left bottom Left box 3
# set key 0.00023,-14.0
# set logscale
# unset logscale; set logscale x
set pointsize 2.5
set xrange [ 0.45 : 1.0 ]
set yrange [ 3.6 : 4.3 ]
#set datafile separator ","
plot "LiCoO2_Thermo_Direct.out" using 1:3 with lines pointtype 7 \
	title "Open circuit vs extent (Direct)" , \
     "LiCoO2_Thermo.out" using 2:3 with linespoints ls 5 \
        title "Open circuit vs extent (Cantera) " , \
     "LiCoO2_dualfoil.out" using 2:3 with linespoints ls 3 \
        title "dualfoil fit "

#plot "temp.csv" using 8:7 with linespoints ls 5 
# plot "obj_func.dat" using 1:2 notitle
#     "FOFT_new-coord.dat" using 2:9 notitle, \
#     "FOFT_new-inj-rate.dat" using 2:4 notitle, \
#     "FOFT_new-inj-rate.dat" using 2:9 notitle
# plot "junk1" using 2:4 notitle
# plot '< join FOFT.dat FOFT_new-coord.dat' using 1:N+1 notitle
