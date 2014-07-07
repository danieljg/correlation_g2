#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -o outputfile
#PBS -e errorfile
#PBS -q default
#PBS -V
#PBS -N ERMjob_pcs
echo compiling
cd build/
make coherence_angle
mv coherence_angle ../
cd ..

# apertures to try
signal_ap=0.3 #(0.3 1.0 5)
idler_ap=( 0.3 1.0 3.0 5)
# configuration of detectors
idler_angle=2.4
# working temp
temp=51.5

# cleaning up previous data, if any
mkdir data/angle
rm data/angle/*.dat

echo running
for((k=0;k<5;k++))
do
time ./coherence_angle ${signal_ap} ${idler_ap[$k]} ${idler_angle} ${temp}
cat coherence.dat >> data/angle/signal_${signal_ap}_idler_${idler_ap[$k]}.dat
done

gnuplot plotting_scripts/plot_angle.p
rm coherence.dat

