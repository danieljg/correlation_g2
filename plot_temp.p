set term postscript eps enhanced color

g(x)=a+x*(b+x*c)
xmin=25;xmax=80

set xlabel 'Temperature'
set ylabel 'G_2'
set xrange [xmin:xmax]

set title "Second order coherence function"

#signal_apertures=(1  1  1   0.1 0.1 0.1 5  5  5   5)
#idler_apertures=( 5  10 15  5   10  15  5  10 15  25)

signal_ap="1"
set output './singles_convergence/temp_signal_'.signal_ap.'.eps'
plot for [idler_ap in "5 10 15"] 'singles_convergence/temp/signal_'.signal_ap.'_idler_'.idler_ap.'.dat' w lp lw 2 ps 0.5 ti 'signal aperture ='.signal_ap.', idler aperture ='.idler_ap

signal_ap="0.1"
set output 'singles_convergence/temp_signal_'.signal_ap.'.eps'
plot for [idler_ap in "5 10 15"] 'singles_convergence/temp/signal_'.signal_ap.'_idler_'.idler_ap.'.dat' w lp lw 2 ps 0.5 ti 'signal aperture ='.signal_ap.', idler aperture ='.idler_ap

signal_ap="5"
set output 'singles_convergence/temp_signal_'.signal_ap.'.eps'
plot for [idler_ap in "5 10 15 25"] 'singles_convergence/temp/signal_'.signal_ap.'_idler_'.idler_ap.'.dat' w lp lw 2 ps 0.5 ti 'signal aperture ='.signal_ap.', idler aperture ='.idler_ap


set xlabel 'Phase-mismatch [PI]'
set xrange [-5:5]
set xtics 1

signal_ap="1"
set output 'singles_convergence/phase_signal_'.signal_ap.'.eps'
plot for [idler_ap in "5 10 15"] 'singles_convergence/temp/signal_'.signal_ap.'_idler_'.idler_ap.'.dat' u 3:2 w lp lw 2 ps 0.5 ti 'signal aperture ='.signal_ap.', idler aperture ='.idler_ap

signal_ap="0.1"
set output 'singles_convergence/phase_signal_'.signal_ap.'.eps'
plot for [idler_ap in "5 10 15"] 'singles_convergence/temp/signal_'.signal_ap.'_idler_'.idler_ap.'.dat' u 3:2 w lp lw 2 ps 0.5 ti 'signal aperture ='.signal_ap.', idler aperture ='.idler_ap

signal_ap="5"
set output 'singles_convergence/phase_signal_'.signal_ap.'.eps'
plot for [idler_ap in "5 10 15 25"] 'singles_convergence/temp/signal_'.signal_ap.'_idler_'.idler_ap.'.dat' u 3:2 w lp lw 2 ps 0.5 ti 'signal aperture ='.signal_ap.', idler aperture ='.idler_ap

