#!/bin/bash

echo compiling
cd build/
make coherence_convergence
mv coherence_convergence ../
cd ..

# apertures to try
signal_ap=(0.1 1 3 5)
idler_ap=( 0.1 0.5 1 2 3 4 5 6 7 8 9 10 12.5 15)
# bandwidths to try
idler_bw=( 20 50 100 200 )
# configuration of detectors
signal_angle=2.4
idler_angle=2.4
# working temp
temp=51.5

# cleaning up previous data, if any
for((l=0;l<4;l++))
do
mkdir singles_convergence/${idler_bw[$l]}nm
rm singles_convergence/${idler_bw[$l]}nm/*.dat
done

echo running
for((k=0;k<4;k++))
do
for((j=0;j<13;j++))
do
for((l=0;l<4;l++))
do
time ./coherence_convergence ${signal_ap[$k]} ${idler_ap[$j]} ${signal_angle} ${idler_angle} ${idler_bw[$l]} ${temp}
cat coherence.dat >> singles_convergence/${idler_bw[$l]}nm/signal_ap_${signal_ap[$k]}.dat
done
done
done

gnuplot plotting_scripts/plot_convergence.p
rm coherence.dat

