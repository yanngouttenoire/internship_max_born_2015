set terminal postscript eps enhanced color solid font 'Helvetica,10'
set output 'comparison[6-12].eps'
set multiplot  layout 1, 1
set xlabel 'Asymptotic energy (au)' offset 0,4
set ylabel 'Probability (log)'
set xrange [0:1]
set key on outside left bmargin box title sprintf("nField = 1000, nVYPerp = 1, nVZPrimPerp = 1000, hydrogen potential, \n binsWidth = 0.001, ErrorMax = 1e-10, dtMin = 1e-20, \n linear field, fieldAmplMax = 0.06au, waveLenght = 1064nm, ")
plot 'spectrum6.dat' index 0 using 1:2 w l lc rgb 'royalblue' title 'weightMinThreshold = 0, weightTooSmallNbr = 0 %, duration = 1 h 0 min', \
'spectrum9.dat' index 0 using 1:2 w l lc rgb 'dark-violet' title 'weightMinThreshold = 5E-5, weightTooSmallNbr = 31%, duration = 0 h 47 min', \
'spectrum8.dat' index 0 using 1:2 w l lc rgb 'orange-red' title 'weightMinThreshold = 1E-4, weightTooSmallNbr = 47%, duration = 0 h 39 min', \
'spectrum10.dat' index 0 using 1:2 w l lc rgb 'brown' title 'weightMinThreshold = 2E-4, weightTooSmallNbr = 65%, duration = 0 h 31 min', \
'spectrum11.dat' index 0 using 1:2 w l lc rgb 'sandybrown' title 'weightMinThreshold = 3E-4, weightTooSmallNbr = 77%, duration = 0 h 22 min', \
'spectrum12.dat' index 0 using 1:2 w l lc rgb 'gold' title 'weightMinThreshold = 4E-4, weightTooSmallNbr = 85.4%, duration = 0 h 15 min', \
'spectrum7.dat' index 0 using 1:2 w l lc rgb 'green' title 'weightMinThreshold = 5E-4, weightTooSmallNbr = 92%, duration = 0 h 9 min'
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
