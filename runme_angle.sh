#!/bin/bash

echo compiling
cd build/
make coherence_angle
mv coherence_angle ../
cd ..

# apertures to try
signal_ap=(0.3 1.0 3 5)
idler_ap=0.3 #( 0.3 1.0 3 5)
# configuration of detectors
signal_angle=2.4
idler_angle=2.4
# working temp
temp=51.5

# cleaning up previous data, if any
mkdir data/angle
rm data/angle/*.dat

echo running
for((k=0;k<5;k++))
do
time ./coherence_angle ${signal_ap[$k]} ${idler_ap} ${idler_angle} ${temp}
cat coherence.dat >> data/angle/signal_${signal_ap[$k]}_idler_${idler_ap}.dat
done

gnuplot plotting_scripts/plot_angle.p
rm coherence.dat

