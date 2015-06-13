set terminal postscript eps enhanced color solid font 'Helvetica,10'
set output 'comparison[37-40].eps'
set multiplot  layout 1, 1
set xlabel 'Asymptotic energy (au)' offset 0,4
set ylabel 'Probability (log)'
set xrange [0:1]

set key on outside center bmargin box title sprintf("hydrogen potential, weightMinThreshold = 0, \n binsWidth = 0.001, ErrorMax = 1e-10, dtMin = 1e-20, \n ellipticity=0.1, fieldAmplMax = 0.06au, waveLenght = 1064nm, ")
set key left
plot 'spectrum37.dat' using 1:2 w l lc rgb 'violet' title 'nField = 100 000, nVYPerp = 1, nVZPrimPerp = 1, initial velocity of the eletron is null', \
'spectrum38.dat' using 1:2 w l lc rgb 'royalblue' title 'nField = 500 000, nVYPerp = 1, nVZPrimPerp = 1, initial velocity of the eletron is null', \
'spectrum39.dat' using 1:2 w l lc rgb 'cyan' title 'nField = 1 000 000, nVYPerp = 1, nVZPrimPerp = 1, initial velocity of the eletron is null', \
'spectrum40.dat' using 1:2 w l lc rgb 'green' title 'nField = 2 000 000, nVYPerp = 1, nVZPrimPerp = 1, initial velocity of the eletron is null'
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
