#!/bin/sh
#PBS-lwalltime=1:00:00
#PBS-lselect=1:ncpus=1:mpiprocs=1:ompthreads=1:mem=23500mb
module load netcdf
ulimit -s unlimited
cd /scratch/rw408/ForcingData/CFSR/CFSR2ROMS

./process_dlw > process_dlw.out
./process_dsw > process_dsw.out
./process_cloud > process_cloud.out
./process_rain > process_rain.out
./process_spe_hum > process_spe_hum.out
./process_tair > process_tair.out
./process_wind > process_wind.out
./process_msl > process_msl.out




