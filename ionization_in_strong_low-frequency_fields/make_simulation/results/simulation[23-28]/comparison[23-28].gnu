set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'comparison[23-28].eps'
set multiplot  layout 1, 1
set xlabel 'Asymptotic energy (au)' offset 0,4
set ylabel 'Probability (log)'
set xrange [0:2]
set key on outside center bmargin box title sprintf("nField = 1000, nVYPerp = 1, nVZPrimPerp = 1000, \n hydrogen potential, weightMinThreshold = 0, \n binsWidth = 0.001, ErrorMax = 1e-10, dtMin = 1e-20, \n fieldAmplMax = 0.06au, waveLenght = 1064nm, duration = 1h29min")
plot 'spectrum23.dat' index 0 using 1:2 w l lc rgb 'violet' title 'ellipticity=0.0', \
'spectrum24.dat' index 0 using 1:2 w l lc rgb 'royalblue' title 'ellipticity=0.1', \
'spectrum25.dat' index 0 using 1:2 w l lc rgb 'cyan' title 'ellipticity=0.2', \
'spectrum26.dat' index 0 using 1:2 w l lc rgb 'green' title 'ellipticity=0.3', \
'spectrum27.dat' index 0 using 1:2 w l lc rgb 'orange' title 'ellipticity=0.4', \
'spectrum28.dat' index 0 using 1:2 w l lc rgb 'red' title 'ellipticity=0.5'
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
