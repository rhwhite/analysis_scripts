#!/bin/sh
#PBS-lwalltime=72:00:00
#PBS-lselect=1:ncpus=1:mpiprocs=1:ompthreads=1:mem=800mb
cd /scratch/rw408/HYCOM/HYCOM2ROMS/

mpiexec smooth_u

mpiexec smooth_v

mpiexec smooth_temp

mpiexec smooth_ssh

mpiexec smooth_salt





