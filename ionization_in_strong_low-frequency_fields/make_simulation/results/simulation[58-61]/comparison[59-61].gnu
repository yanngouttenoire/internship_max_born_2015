set terminal postscript eps enhanced color solid font 'Helvetica,10'
set output 'spectrum[59-61].eps'
set multiplot  layout 1, 1
set xlabel 'Asymptotic energy (au)' offset 0,4
set ylabel 'Probability (log)'
set xrange [0:1]

set key on outside center bmargin Left reverse box title sprintf("nField = 100, nVYPerp = 1000, nVZPrimPerp = 100, \nmolecule potential, charges = 0.4 0.2 0.2 0.2 , bondLength = 2 1 3 , \n weightThreshold = 0, unexpectedStopNbr = 4.4%, binsWidth = 0.001, ErrorMax = 1e-10, dtMin = 1e-20, \n ellipticity = 0.1, fieldAmplMax = 0.06au, waveLenght = 1064nm, duration = 2 h 15 min")

plot 'spectrum59.dat' using 1:2 w l lc rgb 'violet' title 'ellipticity=0.1', \
'spectrum60.dat' using 1:2 w l lc rgb 'royalblue' title 'ellipticity=0.2', \
'spectrum61.dat' using 1:2 w l lc rgb 'cyan' title 'ellipticity=0.3'
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
