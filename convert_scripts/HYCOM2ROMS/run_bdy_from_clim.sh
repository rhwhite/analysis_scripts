#!/bin/sh
#PBS-lwalltime=72:00:00
#PBS-lselect=1:ncpus=1:mpiprocs=1:ompthreads=1:mem=800mb
cd /scratch/rw408/ForcingData/HYCOM/HYCOM2ROMS/

mpiexec bry_from_climatology



