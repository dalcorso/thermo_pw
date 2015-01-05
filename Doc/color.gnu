set xrange [-1.0:8]
set yrange [-1:6]
set arrow from 0.0,2.0 to 1.0,2.0 lw 3 lc rgb "red"
set arrow from 2.0,2.0 to 3.0,2.0 lw 3 lc rgb "green"
set arrow from 4.0,2.0 to 5.0,2.0 lw 3 lc rgb "blue"
set arrow from 6.0,2.0 to 7.0,2.0 lw 3 lc rgb "cyan"

set arrow from 0.0,1.5 to 1.0,1.5 lw 3 lc rgb "magenta"
set arrow from 2.0,1.5 to 3.0,1.5 lw 3 lc rgb "gold"
set arrow from 4.0,1.5 to 5.0,1.5 lw 3 lc rgb "pink"
set arrow from 6.0,1.5 to 7.0,1.5 lw 3 lc rgb "black"

set arrow from 0.0,1.0 to 1.0,1.0 lw 3 lc rgb "olive"
set arrow from 2.0,1.0 to 3.0,1.0 lw 3 lc rgb "brown"
set arrow from 4.0,1.0 to 5.0,1.0 lw 3 lc rgb "light-blue"
set arrow from 6.0,1.0 to 7.0,1.0 lw 3 lc rgb "orange"

set label "1" at -0.2,2.0 center
set label "2" at  1.8,2.0 center
set label "3" at  3.8,2.0 center
set label "4" at  5.8,2.0 center
set label "5" at -0.2,1.5 center
set label "6" at  1.8,1.5 center
set label "7" at  3.8,1.5 center
set label "8" at  5.8,1.5 center
set label "9" at -0.2,1.0 center
set label "10" at  1.7,1.0 center
set label "11" at  3.7,1.0 center
set label "12" at  5.7,1.0 center
plot x+10

