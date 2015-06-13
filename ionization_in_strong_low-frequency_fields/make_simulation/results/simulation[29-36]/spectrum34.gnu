set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'spectrum.eps'
set multiplot  layout 1, 1
set key width -15
set xlabel 'Asymptotic energy (au)'
set ylabel 'Probability (log)'
set xrange [0:1]
set key on outside center bmargin Left reverse box title sprintf("nField = 500, nVYPerp = 1, nVZPrimPerp = 1000, weightThreshold = 0, weightTooSmallNbr = 0%, \n spectraPointNbr = 464782, unexpectedStopNbr = 7%, binsWidth = 0.001, ErrorMax = 1e-10, dtMin = 1e-20, \n ellipticity = 0, fieldAmplMax = 0.0608629au, waveLenght = 1064nm, duration = 1 h 51 min 7 s , ")
plot 'data.dat' index 0 using 1:2 w l lc rgb 'violet' title 'Photo-electron spectrum with vY positive', \
'data.dat' index 1 using 1:2 w l lc rgb 'violet' title 'Photo-electron spectrum with vY negative' 
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
