#!/bin/sh
#PBS-lwalltime=1:00:00
#PBS-lselect=1:ncpus=1:mpiprocs=1:ompthreads=1:mem=23500mb
module load netcdf
ulimit -s unlimited
cd /scratch/rw408/ForcingData/CFSR/CFSR2ROMS

./process_dlw > process_dlw.out
./process_dsw > process_dsw.out




