set terminal postscript eps enhanced color solid font 'Helvetica,10'
set output 'spectrum[58-59].eps'
set multiplot  layout 1, 1
set xlabel 'Asymptotic energy (au)' offset 0,4
set ylabel 'Probability (log)'
set xrange [0.2:0.4]
set yrange [1.5:2.5]
set key on outside center bmargin Left reverse box title sprintf("nField = 100, nVYPerp = 1000, nVZPrimPerp = 100, \nmolecule potential, charges = 0.4 0.2 0.2 0.2 , bondLength = 2 1 3 , \n weightThreshold = 0, unexpectedStopNbr = 4.4%, binsWidth = 0.001, ErrorMax = 1e-10, dtMin = 1e-20, \n ellipticity = 0.1, fieldAmplMax = 0.06au, waveLenght = 1064nm, duration = 18 h 3 min 41 s")

plot 'spectrum58.dat' using 1:2  lc rgb 'violet' title 'nField = 100, nVYPerp = 1000, nVZPrimPerp = 100, duration=18h4min', \
'spectrum59.dat' using 1:2  lc rgb 'royalblue' title 'nField = 100, nVYPerp = 100, nVZPrimPerp = 100, duration=2h15min', \
'spectrum62.dat' using 1:2  lc rgb 'green' title 'nField = 100, nVYPerp = 100, nVZPrimPerp = 100, binsWidth = 0.005'
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
