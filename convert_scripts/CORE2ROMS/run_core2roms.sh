#!/bin/sh
#PBS-lwalltime=72:00:00
#PBS-lselect=1:ncpus=1:mpiprocs=11:ompthreads=1:mem=23000mb
module load intel-suite mpi
module load netcdf

cd /scratch/rw408/ForcingData/CORE2/CORE2ROMS/

./core2roms_monthly

./core2roms_daily

./core2roms_6hr




