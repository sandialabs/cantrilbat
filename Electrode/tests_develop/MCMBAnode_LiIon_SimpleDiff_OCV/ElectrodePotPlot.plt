#
#   Plot of the bottom line results
#   Should mimic Fig. 1 Isaacson
#
#set terminal png
set terminal jpeg 
set title "Cantera-DAKOTA Coupling - Zinc Electrode Rest Potentials at 25 C vs Hg/HgO"
set xlabel "Molality KOH at 0.1 Zincate"
set ylabel "E_Zn (mV) vs Hg/HgO"
set xtic auto
set ytic auto
# set key top right
# set key box
# set key left bottom Left box 3
#set key 0.00023,-14.0
#set logscale x
# unset logscale; set logscale x
set yrange [-1.20 : -1.44 ] reverse
set pointsize 2.5
#
# linestyle
set style line 5 lt rgb "cyan" lw 3 pt 6
#
#
plot "num_data.txt" using 3:1 with linespoints ls 5  \
	title "KOH Dependence" , \
     "exp_data.txt" using 3:1 with points pt 4 title "KOH data"

# plot "obj_func.dat" using 1:2 notitle
#     "FOFT_new-coord.dat" using 2:9 notitle, \
#     "FOFT_new-inj-rate.dat" using 2:4 notitle, \
#     "FOFT_new-inj-rate.dat" using 2:9 notitle
# plot "junk1" using 2:4 notitle
# plot '< join FOFT.dat FOFT_new-coord.dat' using 1:N+1 notitle
