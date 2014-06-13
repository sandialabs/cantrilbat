#
#
set title "Cantera-DAKOTA Coupling - Zinc Electrode Rest Potentials at 25 C vs Hg/HgO"
set xlabel "DG0_Zincate"
# set xlabel "DeltaG_f (cal/gmol)"
set ylabel "Open circuit voltage"
set xtic auto
set ytic auto
# set key top right
# set key box
# set key left bottom Left box 3
# set key 0.00023,-14.0
# set logscale
# unset logscale; set logscale x
set pointsize 2.5
set xrange [ 0.0 : 1.0 ]
set datafile separator ","
plot "MCMB_globalResults_0_0.csv" using 10:6 with points pointtype 7 \
	title "Open circuit vs extent"
#plot "temp.csv" using 8:7 with linespoints ls 5 
# plot "obj_func.dat" using 1:2 notitle
#     "FOFT_new-coord.dat" using 2:9 notitle, \
#     "FOFT_new-inj-rate.dat" using 2:4 notitle, \
#     "FOFT_new-inj-rate.dat" using 2:9 notitle
# plot "junk1" using 2:4 notitle
# plot '< join FOFT.dat FOFT_new-coord.dat' using 1:N+1 notitle
