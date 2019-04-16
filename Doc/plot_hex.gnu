set encoding iso_8859_15
set terminal postscript enhanced solid color "AvantGarde-Book" 20
set output "figure_hex.ps"

set nokey
set size 1.0,1.0
a=1.0
sqrt3=sqrt(3.0)
set xrange[-2.3*a:2.3*a]
set yrange[-2.3*a:2.3*a]
set arrow from -2.0*a,0.0 to 2.0*a lw 3 lc rgb "blue"
set arrow from 0.0,-2.0*a to 0.0,2.0*a lw 3 lc rgb "blue"

set arrow from a,0.0 to a*0.5,a*sqrt3*0.5 lw 3  lc rgb "black" nohead
set arrow from -a,0.0 to -a*0.5,a*sqrt3*0.5 lw 3  lc rgb "black" nohead
set arrow from -a*0.5,a*sqrt3*0.5 to a*0.5,a*sqrt3*0.5 lw 3 lc rgb "black" nohead
set arrow from a,0.0 to a*0.5,-a*sqrt3*0.5 lw 3  lc rgb "black" nohead
set arrow from -a,0.0 to -a*0.5,-a*sqrt3*0.5 lw 3  lc rgb "black" nohead
set arrow from -a*0.5,-a*sqrt3*0.5 to a*0.5,-a*sqrt3*0.5 lw 3 lc rgb "black" nohead

unset xtics
unset ytics
unset border

b=a*0.8
set arrow from -b,-b*sqrt3 to b,b*sqrt3 lw 3 lc rgb "red" nohead
set arrow from -b,b*sqrt3 to b,-b*sqrt3 lw 3 lc rgb "olive" nohead

set arrow from -b*sqrt3,-b to b*sqrt3,b lw 3 lc rgb "olive" nohead
set arrow from -b*sqrt3,b to b*sqrt3,-b lw 3 lc rgb "red" nohead

set label "[2,1,0]" at 1.5,0.8 font "Times-Roman,32"
set label "[1,1,0]" at 0.83,1.5 font "Times-Roman,32"
set label "[0,1,0]" at -1.5,1.5 font "Times-Roman,32"
set label "[1,-1,0]" at 1.5,-0.8 font "Times-Roman,32"

set label "x" at 1.9,-0.20 font "Times-Roman,32"
set label "y" at 0.1,1.9 font "Times-Roman,32"
plot x+10
