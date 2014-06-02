set term postscript eps enhanced color

set xlabel 'Temperature'
set ylabel 'G_2'
set xrange [25:80]

set output 'coherence.eps'
set title "Second order coherence function for\n0.1mm apertures at 21mm from the optic axis"
plot 'coherence.dat' w lp pt 7 ps 0.5 lc 3
