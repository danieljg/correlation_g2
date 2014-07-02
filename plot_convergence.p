set term postscript eps enhanced color

xmin=0;xmax=20

set xlabel 'Idler aperture [mm]'
set ylabel 'G_2'
set xrange [xmin:xmax]

# apertures to try
#signal_ap=(0.1 1 3 5)
#idler_ap=( 0.1 0.5 1 2 3 4 5 6 7 8 9 10 12.5 15)
# bandwidths to try
#idler_bw=( 20 50 100 200 )
# configuration of detectors
#signal_angle=2.4
#idler_angle=2.4
# working temp
#temp=51.5

do for [k in " 0.1 1 3 5"]{ 
set output 'singles/convergence/signal_ap_'.k.'.eps'
set title "Convergence test for idler aperture diameters and idler filter bandwidths\n T=51.5C, external signal and idler angle=2.4degrees, signal aperture =".k."mm"
plot for [l in "20 50 100 200] './singles_convergence/'.l.'nm/signal_ap_'.k.'.dat' u 1:3  w lp ti 'idler bandwidth = '.l.'nm'
}
