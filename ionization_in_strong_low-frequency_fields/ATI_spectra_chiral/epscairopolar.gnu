set terminal epscairo font 'Gill Sans,9' rounded fontscale 0.4
unset border
set polar
set angles degrees

set style line 100 linecolor rgb '#808080'
set style line 101 linecolor rgb '#808080' linetype 0  lw 2
set style line 1 lc rgb '#A00000' pt 6 ps 1 lt 1 lw 4
set style line 2 lc rgb '#00A0A0' pt 6 ps 1 lt 1 lw 4
set style line 3 lc rgb '#00A000' pt 6 ps 1 lt 1 lw 4
set style line 6 lc rgb '#A000A0' pt 6 ps 1 lt 1 lw 4
set style line 4 lc rgb '#0000A0' pt 6 ps 1 lt 1 lw 4
set style line 5 lc rgb '#A0A000' pt 6 ps 1 lt 1 lw 4

set grid polar 30 ls 101

set size square
set border linestyle 101
unset border
set xtics scale 0 0,.5,10
set ytics scale 0 0,.5,10
set format x ""
set format y ""
set key bmargin center
set format r ""
set rtics scale 0
unset raxis

amplitude=0.55
set_label_angle(x, text,tag) = sprintf("set label %f '%s' at graph (0.5+amplitude*cos(%f)), (0.5+amplitude*sin(%f)) center", tag, text, x, x)

#here all labels are created
eval set_label_angle(0, "0°",100)
#eval set_label_angle(60, "60°",101)
eval set_label_angle(90, "90°",102)
#eval set_label_angle(120, "120°",103)
eval set_label_angle(180, "180°",104)
#eval set_label_angle(240, "240°",105)
eval set_label_angle(270, "270°",106)
#eval set_label_angle(300, "300°",107)


