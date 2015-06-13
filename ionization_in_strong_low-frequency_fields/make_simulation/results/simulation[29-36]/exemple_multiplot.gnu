set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'comparison[29-36].eps'
set multiplot  layout 2, 2
set rmargin 0
set xlabel 'Asymptotic energy (au)'
set ylabel 'Probability (log)'
set xrange [0:1]

set key on outside center bmargin box title sprintf("weightMinThreshold = 0, weightTooSmallNbr = 0%, \n binsWidth = 0.001, ErrorMax = 1e-10, dtMin = 1e-20, \n linear field, fieldAmplMax = 0.0608629au, waveLenght = 1064nm, ")
set key left
plot 'spectrum29.dat' index 0 using 1:2 w l lc rgb 'violet' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 100, duration = 29min', \
'spectrum30.dat' index 0 using 1:2 w l lc rgb 'royalblue' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 500, duration =  1h52min', \
'spectrum6.dat' index 0 using 1:2 w l lc rgb 'cyan' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 1000, duration =  1h0min', \
'spectrum31.dat' index 0 using 1:2 w l lc rgb 'green' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 2000, duration =  5h14min', \
'spectrum32.dat' index 0 using 1:2 w l lc rgb 'orange' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 5000, duration =  8h39min'
unset ylabel
set lmargin 0
set rmargin 0
plot 'spectrum33.dat' index 0 using 1:2 w l lc rgb 'violet' title 'nField = 100, nVYPerp = 1, nVZPrimPerp = 1000, duration = 29min', \
'spectrum34.dat' index 0 using 1:2 w l lc rgb 'royalblue' title 'nField = 500, nVYPerp = 1, nVZPrimPerp = 1000, duration =  1h51min', \
'spectrum6.dat' index 0 using 1:2 w l lc rgb 'cyan' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 1000, duration =  1h0min', \
'spectrum35.dat' index 0 using 1:2 w l lc rgb 'green' title 'nField = 2000, nVYPerp = 1, nVZPrimPerp = 1000, duration =  5h15min', \
'spectrum36.dat' index 0 using 1:2 w l lc rgb 'orange' title 'nField = 5000, nVYPerp = 1, nVZPrimPerp = 1000, duration =  8h40min'
plot 'spectrum29.dat' index 0 using 1:2 w l lc rgb 'violet' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 100', \
'spectrum30.dat' index 0 using 1:2 w l lc rgb 'royalblue' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 500', \
'spectrum6.dat' index 0 using 1:2 w l lc rgb 'cyan' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 1000', \
'spectrum31.dat' index 0 using 1:2 w l lc rgb 'green' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 2000', \
'spectrum32.dat' index 0 using 1:2 w l lc rgb 'orange' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 5000', \
'spectrum33.dat' index 0 using 1:2 w l lc rgb 'black' title 'nField = 100, nVYPerp = 1, nVZPrimPerp = 1000', \
'spectrum34.dat' index 0 using 1:2 w l lc rgb 'black' title 'nField = 500, nVYPerp = 1, nVZPrimPerp = 1000', \
'spectrum6.dat' index 0 using 1:2 w l lc rgb 'black' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 1000', \
'spectrum35.dat' index 0 using 1:2 w l lc rgb 'black' title 'nField = 2000, nVYPerp = 1, nVZPrimPerp = 1000', \
'spectrum36.dat' index 0 using 1:2 w l lc rgb 'black' title 'nField = 5000, nVYPerp = 1, nVZPrimPerp = 1000'
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
