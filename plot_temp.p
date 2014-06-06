set term postscript eps enhanced color

g(x)=a+x*(b+x*c)
xmin=25;xmax=80

set xlabel 'Temperature'
set ylabel 'G_2'
set xrange [xmin:xmax]

set output 'coherence.eps'
set title "Second order coherence function"
plot 'coherence.dat' w lp pt 7 ps 0.5 lc 3


set xlabel 'Phase-mismatch [PI]'
set xrange [-5:5]
set xtics 1
a=-165279;#@2.28degrees
a=-16.4346;#@2.3degrees
b=0.259301;c=0.00117307
set output 'coherence_phase.eps'
plot 'coherence.dat' u (g($1)):2 w lp pt 7 ps 0.5 lc 3