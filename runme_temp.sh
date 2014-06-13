#!/bin/bash
echo compiling
make coherence_temp
signal_apertures=(1 0.1 0.1 5)
idler_apertures=(1 0.1 5 5)
echo running
for ((k=0;k<=3;k++))
do
 time ./coherence_temp ${signal_apertures[$k]} ${idler_apertures[$k]} 2.4 2.4
 gnuplot plot_temp.p
 cp coherence.eps graphs/temp/signal_${signal_apertures[$k]}_idler_${idler_apertures[$k]}.eps
 cp coherence_phase.eps graphs/temp/phase_signal_${signal_apertures[$k]}_idler_${idler_apertures[$k]}.eps
 cp coherence.dat archive/temp/signal_${signal_apertures[$k]}_idler_${idler_apertures[$k]}.dat
done
