#!/bin/bash
echo compiling
make coherence_temp
# apertures to try
signal_ap=(0.1 1 3 5)
idler_ap=( 5 7.5 10 15)
# bandwidths to try
idler_bw=( 20 50 100 200 )
# configuration of detectors
signal_angle=2.4
idler_angle=2.4
# working temp
temp=51.5
# cleaning up previous data, if any
rm singles_convergence/temp/*.dat
echo running
for ((k=0;k<4;k++))
 do for ((j=0;j<4;k++))
  do for ((l=0;l<4;l++))
   time ./coherence_convergence ${signal_ap[$k]} ${idler_ap[$j]} ${signal_angle} ${idler_angle} ${idler_bw[$l]} ${temp}
   cat coherence.dat >> singles_convergence/temp/signal_ap_${signal_ap[$k]}_idler_bw_${idler_bw[$l]}.dat
  done
 done
done
gnuplot plot_convergence.p

