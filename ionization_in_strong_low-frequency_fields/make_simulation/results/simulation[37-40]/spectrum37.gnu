set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'spectrum.eps'
set multiplot  layout 1, 1
set key width -15
set xlabel 'Asymptotic energy (au)'
set ylabel 'Probability (log)'
set xrange [0:1]
set key on outside center bmargin Left reverse box title sprintf("nField = 1000, nVYPerp = 100, nVZPrimPerp = 1, charges = 0.4 0.2 0.2 0.2 , bondLength = 2 1 3 , \n weightThreshold = 0, weightTooSmallNbr = 0%, spectraPointNbr = 100000, unexpectedStopNbr = 0%, binsWidth = 0.005, \n ErrorMax = 1e-10, dtMin = 1e-20, ellipticity = 0.1, fieldAmplMax = 0.0608629au, waveLenght = 1064nm, \n duration = 0 h 19 min 40 s , ")
plot 'data.dat' index 0 using 1:2 w l lc rgb 'violet' title 'Photo-electron spectrum with vY positive', \
'data.dat' index 1 using 1:2 w l lc rgb 'violet' title 'Photo-electron spectrum with vY negative' 
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
