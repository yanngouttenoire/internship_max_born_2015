#Line style for axes, grey color
set style line 100 linecolor rgb '#808080'
#Line style for grid, dashed, grey color
set style line 101 linecolor rgb '#808080' linetype 0 

#Line styles for curves
#RED
set style line 1 lc rgb '#A00000' pt 6 ps 1 lt 1 lw 2
#GREEN-BLUE (cyan) (Complementary Color)
set style line 2 lc rgb '#00A0A0' pt 6 ps 1 lt 1 lw 2
#GREEN
set style line 3 lc rgb '#00A000' pt 6 ps 1 lt 1 lw 2
#RED-BLUE (purple) (Complementary Color)
set style line 6 lc rgb '#A000A0' pt 6 ps 1 lt 1 lw 2
#BLUE
set style line 4 lc rgb '#0000A0' pt 6 ps 1 lt 1 lw 2
#RED-GREEN (gold) (Complementary Color)
set style line 5 lc rgb '#A0A000' pt 6 ps 1 lt 1 lw 2

set terminal epscairo font 'Gill Sans,9' rounded fontscale 0.4

#Remove border on bottom and right, these borders are useless and make it harder to see plotted lines near the border
#Also, put it in grey, no need for so much emphasis on a border 
set border 6 back linestyle 100
set grid linestyle 101
#We also put the key box in grey
set key box ls 100

set xtics scale 0
set format x '' 
#set xrange [0:20]
set x2tics 
set x2label 'Asymptotic energy (eV)'
set ylabel 'Probability (linear scale)'

