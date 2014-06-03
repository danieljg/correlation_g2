#!/bin/bash
echo compiling
make
signal_apertures=(1 3 5)
idler_apertures=(1 3 5)
echo running
for ((k=0;k<=2;k++))
do
 time ./coherence ${signal_apertures[$k]} ${idler_apertures[$k]} 21.0 21.0
 gnuplot plot.p
 cp coherence.eps graphs/signal_${signal_apertures[$k]}_idler_${idler_apertures[$k]}.eps
 cp coherence.dat archive/signal_${signal_apertures[$k]}_idler_${idler_apertures[$k]}.dat
done

