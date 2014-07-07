set term postscript eps enhanced color

g(x)=a+x*(b+x*c)
xmin=1.5;xmax=3.0

set xlabel 'External angle [degrees]'
set xtics 0.1
set ylabel 'G_2'
set xrange [xmin:xmax]

set output 'coherence_angle.eps'
set title "Second order coherence function"
plot 'coherence.dat' u 1:3 w lp pt 7 ps 0.5 lc 3

set output 'coherence_norm.eps'
set title "Second order coherence function"
plot 'coherence.dat' u 1:3 w lp pt 7 ps 0.5 lc 3

