set terminal postscript eps enhanced color solid font 'Helvetica,10'
set output 'comparison[29-36].eps'
set multiplot  layout 1, 1
set xlabel 'Asymptotic energy (au)' offset 0,4
set ylabel 'Probability (log)'
set xrange [0:1]

set key on outside center bmargin box title sprintf("hydrogen potential, weightMinThreshold = 0, \n binsWidth = 0.001, ErrorMax = 1e-10, dtMin = 1e-20, \n linear field, fieldAmplMax = 0.06au, waveLenght = 1064nm, ")
set key left
plot 'spectrum29.dat' using 1:2 w l lc rgb 'violet' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 100', \
'spectrum30.dat' using 1:2 w l lc rgb 'royalblue' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 500', \
'spectrum6.dat' using 1:2 w l lc rgb 'cyan' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 1000', \
'spectrum31.dat' using 1:2 w l lc rgb 'green' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 2000', \
'spectrum32.dat' using 1:2 w l lc rgb 'orange' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 5000', \
'spectrum33.dat' using 1:2 w l lc rgb 'black' title 'nField = 100, nVYPerp = 1, nVZPrimPerp = 1000', \
'spectrum34.dat' using 1:2 w l lc rgb 'black' title 'nField = 500, nVYPerp = 1, nVZPrimPerp = 1000', \
'spectrum6.dat' using 1:2 w l lc rgb 'black' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 1000', \
'spectrum35.dat' using 1:2 w l lc rgb 'black' title 'nField = 2000, nVYPerp = 1, nVZPrimPerp = 1000', \
'spectrum36.dat' using 1:2 w l lc rgb 'black' title 'nField = 5000, nVYPerp = 1, nVZPrimPerp = 1000'
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
