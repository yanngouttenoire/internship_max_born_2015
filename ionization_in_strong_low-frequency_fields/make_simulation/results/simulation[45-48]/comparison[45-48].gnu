set terminal postscript eps enhanced color solid font 'Helvetica,10'
set output 'comparison[45-48].eps'
set multiplot  layout 1, 1
set xlabel 'Asymptotic energy (au)' offset 0,4
set ylabel 'Probability (log)'
set xrange [0:1]
set key on outside center bmargin Left reverse box title sprintf("molecule potential, charges = 0.4 0.2 0.2 0.2 , bondLength = 2 1 3 , \n weightThreshold = 0, binsWidth = 0.005, ErrorMax = 1e-10, dtMin = 1e-20, \n ellipticity = 0, fieldAmplMax = 0.06au, waveLenght = 1064nm")
set label "Vy positive" center at 0.2,2.25 tc rgb "red"
set label "Vy negative" center at 0.2,1.75 tc rgb "red"
plot 'spectrum45.dat' using 1:2 w l lc rgb 'violet' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 100', \
'spectrum46.dat' using 1:2 w l lc rgb 'royalblue' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 500', \
'spectrum47.dat' using 1:2 w l lc rgb 'green' title 'nField = 100, nVYPerp = 1, nVZPrimPerp = 1000', \
'spectrum48.dat' using 1:2 w l lc rgb 'red' title 'nField = 500, nVYPerp = 1, nVZPrimPerp = 1000'
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
