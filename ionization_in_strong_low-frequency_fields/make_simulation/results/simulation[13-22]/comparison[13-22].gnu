set terminal postscript eps enhanced color solid font 'Helvetica,10'
set output 'comparison[13-22].eps'
set multiplot  layout 1, 1
set xlabel 'Asymptotic energy (au)' offset 0,4
set ylabel 'Probability (log)'
set xrange [0:1]
set key on outside center bmargin box title sprintf("hydrogen potential, weightMinThreshold = 0, \n binsWidth = 0.001, ErrorMax = 1e-10, dtMin = 1e-20, \n linear field, fieldAmplMax = 0.06au, waveLenght = 1064nm, ")
plot 'spectrum13.dat' index 0 using 1:2 w l lc rgb 'royalblue' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 100, duration = 29min', \
'spectrum14.dat' index 0 using 1:2 w l lc rgb 'skyblue' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 500, duration =  2h2min', \
'spectrum6.dat' index 0 using 1:2 w l lc rgb 'skyblue' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 1000, duration =  1h0min', \
'spectrum15.dat' index 0 using 1:2 w l lc rgb 'orange-red' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 2000, duration =  6h24min', \
'spectrum16.dat' index 0 using 1:2 w l lc rgb 'brown' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 5000, duration =  12h12min', \
'spectrum17.dat' index 0 using 1:2 w l lc rgb 'sandybrown' title 'nField = 1000, nVYPerp = 1, nVZPrimPerp = 10000, duration = 17h20min', \
'spectrum18.dat' index 0 using 1:2 w l lc rgb 'black' title 'nField = 100, nVYPerp = 1, nVZPrimPerp = 1000, duration = 29min', \
'spectrum19.dat' index 0 using 1:2 w l lc rgb 'black' title 'nField = 500, nVYPerp = 1, nVZPrimPerp = 1000, duration =  2h2min', \
'spectrum20.dat' index 0 using 1:2 w l lc rgb 'black' title 'nField = 2000, nVYPerp = 1, nVZPrimPerp = 1000, duration =  6h27min', \
'spectrum21.dat' index 0 using 1:2 w l lc rgb 'black' title 'nField = 5000, nVYPerp = 1, nVZPrimPerp = 1000, duration =  12h10min', \
'spectrum22.dat' index 0 using 1:2 w l lc rgb 'black' title 'nField = 10000, nVYPerp = 1, nVZPrimPerp = 1000, duration = 17h24min'
unset multiplot
pause -1
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'courbe.eps'
replot
