set terminal postscript eps enhanced color solid font 'Helvetica,10'
set output 'comparison[54-57].eps'
set multiplot  layout 1, 1
set xlabel 'Asymptotic energy (au)' offset 0,4
set ylabel 'Probability (log)'
set xrange [0:1]

set key on outside center bmargin box title sprintf("molecule potential, charges = 0.4 0.2 0.2 0.2 , bondLength = 2 1 3, \n weightMinThreshold = 0, binsWidth = 0.005, ErrorMax = 1e-10, dtMin = 1e-20, \n ellipticity=0.1, fieldAmplMax = 0.06au, waveLenght = 1064nm, ")
set key left
plot 'spectrum54.dat' using 1:2 w l lc rgb 'violet' title 'nField = 1000, nVYPerp = 100, nVZPrimPerp = 1', \
'spectrum55.dat' using 1:2 w l lc rgb 'royalblue' title 'nField = 1000, nVYPerp = 500, nVZPrimPerp = 1', \
'spectrum56.dat' using 1:2 w l lc rgb 'cyan' title 'nField = 1000, nVYPerp = 1000, nVZPrimPerp = 1', \
'spectrum57.dat' using 1:2 w l lc rgb 'green' title 'nField = 1000, nVYPerp = 2000, nVZPrimPerp = 1'
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
