#!/bin/bash
ncks -d clim_time,0 SODA_FSC_clim_15y_coarse_V2.nc SODA_FSC_ini_V2.nc
ncrename -v clim_time,ocean_time SODA_FSC_ini_V2.nc
ncrename -d clim_time,ocean_time SODA_FSC_ini_V2.nc
ncap2 -s ocean_time=ocean_time+15.5 SODA_FSC_ini_V2.nc SODA_FSC_ini_coarse_V2.nc
rm SODA_FSC_ini_V2.nc
echo 'inifile created'
