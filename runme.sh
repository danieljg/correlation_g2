#!/bin/bash
echo compiling
make coherence_temp
signal_apertures=(1 0.1 0.1 5)
idler_apertures=(1 0.1 5 5)
echo running
for ((k=0;k<=3;k++))
do
 time ./coherence_temp ${signal_apertures[$k]} ${idler_apertures[$k]} 21.0 21.0
 gnuplot plot.p
 cp coherence.eps graphs/temp/signal_${signal_apertures[$k]}_idler_${idler_apertures[$k]}.eps
 cp coherence.dat archive/temp/signal_${signal_apertures[$k]}_idler_${idler_apertures[$k]}.dat
done

