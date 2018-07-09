#!/bin/bash

cd /scratch/rw408/ForcingData/SODA/DATA/SODA_ACT_V4

ncrcat SODA_ACT_V4.3_clim_199707.nc \
       SODA_ACT_V4.3_clim_199708.nc \
       SODA_ACT_V4.3_clim_199709.nc \
       SODA_ACT_V4.3_clim_199710.nc \
       SODA_ACT_V4.3_clim_199711.nc \
       SODA_ACT_V4.3_clim_199712.nc \
       SODA_ACT_V4.3_clim_199801.nc \
       SODA_ACT_V4.3_clim_199802.nc \
ACT_V4.3_SODAclim_Jul1997_6mo.nc
echo 'Climfile created' 

