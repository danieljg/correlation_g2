#!/bin/bash
echo compiling
make coherence_angle
signal_apertures=(1 0.1 0.1 5)
idler_apertures=(1 0.1 5 5)
echo running
for ((k=0;k<=3;k++))
do
 time ./coherence_angle ${signal_apertures[$k]} ${idler_apertures[$k]}
 gnuplot plot_angle.p
 cp coherence_angle.eps graphs/angle/signal_${signal_apertures[$k]}_idler_${idler_apertures[$k]}.eps
 cp coherence.dat archive/angle/signal_${signal_apertures[$k]}_idler_${idler_apertures[$k]}.dat
done
