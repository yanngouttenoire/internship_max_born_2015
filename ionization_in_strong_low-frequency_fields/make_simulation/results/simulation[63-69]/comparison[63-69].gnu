set terminal postscript eps enhanced color solid font 'Helvetica,10'
set output 'comparison[63-69].zoom.eps'
set multiplot  layout 1, 1
set xlabel 'Asymptotic energy (au)' offset 0,4
set ylabel 'Probability (log)'
set xrange [0:0.7]
set yrange [2:3.1]
set key on outside center bmargin box title sprintf("nField = 100, nVYPerp = 1000, nVZPrimPerp = 100, charges = 0.4 0.2 0.2 0.2 , bondLength = 2 1 3 , \n weightThreshold = 0, binsWidth = 0.002, \n ErrorMax = 1e-10, dtMin = 1e-20, ellipticity = 0.1, fieldAmplMax = 0.06au, waveLenght = 1064nm, \n duration = 60 h")
set key left
plot 'spectrum63.dat' using 1:2 w l lc rgb 'violet' title 'chiral molecule W', \
'spectrum64.dat' using 1:2 w l lc rgb 'royalblue' title 'chiral molecule X1', \
'spectrum65.dat' using 1:2 w l lc rgb 'cyan' title 'chiral molecule X2', \
'spectrum66.dat' using 1:2 w l lc rgb 'green' title 'chiral molecule Z1', \
'spectrum67.dat' using 1:2 w l lc rgb 'yellow' title 'chiral molecule Z2', \
'spectrum68.dat' using 1:2 w l lc rgb 'orange' title 'chiral molecule Y1', \
'spectrum69.dat' using 1:2 w l lc rgb 'red' title 'chiral molecule Y2'
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
