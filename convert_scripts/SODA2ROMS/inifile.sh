#!/bin/bash
module load nco
cd /scratch/rw408/ForcingData/SODA/DATA/SODA_ACT_V5

ncks -d clim_time,0 ACT_V5_SODAclim_Jul1997_11yr.nc SODAclim_temp.nc
ncrename -v clim_time,ocean_time SODAclim_temp.nc
ncrename -d clim_time,ocean_time SODAclim_temp.nc
ncap -s ocean_time=ocean_time+15.5 SODAclim_temp.nc SODAini_temp.nc
ncap -s ocean_time=ocean_time*24.*60.*60. SODAini_temp.nc ACT_V5_SODAini_Jul1997.nc
rm SODAclim_temp.nc
rm SODAini_temp.nc
echo 'inifile created'

