set terminal postscript eps enhanced color solid font 'Helvetica,10'
set output 'comparison[49-53].eps'
set multiplot  layout 1, 1
set xlabel 'Asymptotic energy (au)' offset 0,4
set ylabel 'Probability (log)'
set xrange [0:1]
set key on outside center bmargin Left reverse box title sprintf("nField = 1000, nVYPerp = 1, nVZPrimPerp = 100, \n molecular potential, charges = 0.4 0.2 0.2 0.2 , bondLength = 2 1 3 , \n weightThreshold = 0, binsWidth = 0.005, ErrorMax = 1e-10, dtMin = 1e-20, \n elliptic field, fieldAmplMax = 0.06au, waveLenght = 1064nm")
set label "Vy positive" center at 0.4,1.75 tc rgb "gold"
set label "Vy negative" center at 0.4,-1 tc rgb "gold"
plot 'spectrum45.dat' using 1:2 w l lc rgb 'violet' title 'ellipticity=0', \
'spectrum49.dat' using 1:2 w l lc rgb 'royalblue' title 'ellipticity=0.1', \
'spectrum50.dat' using 1:2 w l lc rgb 'cyan' title 'ellipticity=0.2', \
'spectrum51.dat' using 1:2 w l lc rgb 'green' title 'ellipticity=0.3', \
'spectrum52.dat' using 1:2 w l lc rgb 'orange' title 'ellipticity=0.4', \
'spectrum53.dat' using 1:2 w l lc rgb 'red' title 'ellipticity=0.5'
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
