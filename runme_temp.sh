#!/bin/bash
echo compiling
make coherence_temp
signal_apertures=(1  1  1   0.1 0.1 0.1 5  5  5   5)
idler_apertures=( 5  10 15  5   10  15  5  10 15  25)
echo running
for ((k=0;k<10;k++))
do
 time ./coherence_temp ${signal_apertures[$k]} ${idler_apertures[$k]} 2.4 2.4
 cp coherence.dat singles_convergence/temp/signal_${signal_apertures[$k]}_idler_${idler_apertures[$k]}.dat
done
gnuplot plot_temp.p

