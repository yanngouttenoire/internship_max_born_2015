set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'spectrum.eps'
set multiplot  layout 1, 1
set xlabel 'Asymptotic energy (au)'
set ylabel 'Probability (log)'
set xrange [0:1]
set key on outside left bmargin box title sprintf("nField = 1000, nVYPerp = 1, nVZPrimPerp = 1000, weightMinThreshold = 0.0004, weightTooSmallNbr = 85.4326%, \n spectraPointNbr = 96433, unexpectedStopNbr = 4.9241%, binsWidth = 0.001, ErrorMax = 1e-10, dtMin = 1e-20, \n linear field, fieldAmplMax = 0.0608629au, waveLenght = 1064nm, duration = 0 h 14 min 55 s , ")
plot 'data.dat' index 0 using 1:2 w l lc rgb 'violet' title 'Photo-electron spectrum with vY positive', \
'data.dat' index 1 using 1:2 w l lc rgb 'violet' title 'Photo-electron spectrum with vY negative' 
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
