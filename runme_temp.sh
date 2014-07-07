#!/bin/bash

echo compiling
cd build/
make coherence_temp
mv coherence_temp ../
cd ..

# apertures to try
signal_ap=(0.3 1.0 5)
idler_ap=0.3 #( 0.3 1.0 5)
# configuration of detectors
signal_angle=2.4
idler_angle=2.4

# cleaning up previous data, if any
mkdir data/temp
rm data/temp/*.dat

echo running
for((k=0;k<5;k++))
do
time ./coherence_temp ${signal_ap[$k]} ${idler_ap} ${signal_angle} ${idler_angle}
cat coherence.dat >> data/temp/signal_${signal_ap[$k]}_idler_${idler_ap}.dat
done

gnuplot plotting_scripts/plot_temp.p
rm coherence.dat

