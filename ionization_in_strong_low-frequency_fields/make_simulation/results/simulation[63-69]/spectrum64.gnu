set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'spectrum.eps'
set multiplot  layout 1, 1
set key width -15
set xlabel 'Asymptotic energy (au)'
set ylabel 'Probability (log)'
set xrange [0:1]
set key on outside center bmargin Left reverse box title sprintf("nField = 100, nVYPerp = 1000, nVZPrimPerp = 100, charges = 0.4 0.2 0.2 0.2 , bondLength = 2 3 -1 , \n weightThreshold = 0, weightTooSmallNbr = 0%, spectraPointNbr = 9501080, unexpectedStopNbr = 4.9%, binsWidth = 0.002, \n ErrorMax = 1e-10, dtMin = 1e-20, ellipticity = 0.1, fieldAmplMax = 0.0608629au, waveLenght = 1064nm, \n duration = 56 h 9 min 41 s , ")
plot 'data.dat' index 0 using 1:2 w l lc rgb 'violet' title 'Photo-electron spectrum with vY positive', \
'data.dat' index 1 using 1:2 w l lc rgb 'violet' title 'Photo-electron spectrum with vY negative' 
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
