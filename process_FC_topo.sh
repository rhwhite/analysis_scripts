#!/bin/sh
# Script to calculate annual averages, climatological averages and running
# averages on raw monthly data from CESM (coupled) 
#
# Author: Rachel White, rachel.white@cantab.net
#
# Created: June 2018

# get top level, only Atlantic, specific variables and concatenate
dir="/glade/scratch/rachelwh/archive/"
exp="CAM4POP_NoMThosing"
printf -v yrstart "%04d" 1
printf -v yrend "%04d" 200
echo $yrstart

cd ${dir}/${exp}/ocn/hist/
mkdir raw

mv ${exp}.pop.h.*.nc raw/
cd raw
mkdir Extra
mv Extra/* .

mv ${exp}.pop.h.????-??.nc Extra
for iyear in `seq -f "%04g" $startyear $endyear`; do
    # record average across all monthly files for that year    
    mv Extra/${exp}.pop.h.????-??.nc .
done

ncrcat -O -d z_t,0,0 -d transport_reg,1,1 -v MOC,HMXL,RHO,SALT,TEMP,BSF,N_HEAT,N_SALT,TAUX,TAUY,TAREA raw/${exp}.pop.h.0* cat_limited_${exp}.pop.h.${yrstart}-${yrend}.nc

# Average over transport_ref dimension to get rid of it
ncwa -O -a transport_reg cat_limited_${exp}.pop.h.${yrstart}-${yrend}.nc cat_limited_${exp}.pop.h.${yrstart}-${yrend}.nc

# Shifttime to get correct months (otherwise cdo reads end of Jan as 1-Feb)
cdo yearmean -shifttime,-1mo cat_limited_${exp}.pop.h.${yrstart}-${yrend}.nc annmean_limited_${exp}.pop.h.${yrstart}-${yrend}.nc

# Need to shifttime BEFORE selecting month, which means putting it AFTER in the command
cdo yearmean -selmon,1,2,12 -shifttime,-1mo cat_limited_${exp}.pop.h.${yrstart}-${yrend}.nc annmean_limited_${exp}.pop.h.${yrstart}-${yrend}.nc 
cdo yearmean -selmon,3,4,5 -shifttime,-1mo cat_limited_${exp}.pop.h.${yrstart}-${yrend}.nc annmean_limited_${exp}.pop.h.${yrstart}-${yrend}.nc
cdo yearmean -selmon,6,7,8 -shifttime,-1mo cat_limited_${exp}.pop.h.${yrstart}-${yrend}.nc annmean_limited_${exp}.pop.h.${yrstart}-${yrend}.nc
cdo yearmean -selmon,9,10,11 -shifttime,-1mo cat_limited_${exp}.pop.h.${yrstart}-${yrend}.nc annmean_limited_${exp}.pop.h.${yrstart}-${yrend}.nc

# Calculate climatologies
ncra -d time,100,200 annmean_limited_CAM4POP_f19g16C_noMT.pop.h.${yrstart}-${yrend}.nc ANNclim_CAM4POP_f19g16C_noMT.pop.h.0001-0300.nc
ncra -d time,100,200 DJF_limited_CAM4POP_f19g16C_noMT.pop.h.${yrstart}-${yrend}.nc DJFclim_CAM4POP_f19g16C_noMT.pop.h.0001-0300.nc
ncra -d time,100,200 MAM_limited_CAM4POP_f19g16C_noMT.pop.h.${yrstart}-${yrend}.nc MAMclim_CAM4POP_f19g16C_noMT.pop.h.0001-0300.nc
ncra -d time,100,200 JJA_limited_CAM4POP_f19g16C_noMT.pop.h.${yrstart}-${yrend}.nc JJAclim_CAM4POP_f19g16C_noMT.pop.h.0001-0300.nc
ncra -d time,100,200 SON_limited_CAM4POP_f19g16C_noMT.pop.h.${yrstart}-${yrend}.nc SONclim_CAM4POP_f19g16C_noMT.pop.h.0001-0300.nc


# Regrid
export NCL_dir=${dir}/${exp}/ocn/hist/
export NCL_file='clim_CAM4POP_f19g16C_noMT.pop.h.0001-0300.nc'

# Need to change directory so NCL can find the weights files
#cd /home/disk/eos4/rachel/git/NCL/Regrid/
#ncl regrid_POP_seas.ncl



# Now do atmospheric variables
#cd ${dir}/${exp}/atm/hist/

#ncrcat -v PRECT,TS,TREFHT,TMQ,TAUX,TAUY,SHFLX,LHFLX,RHREFHT,QFLX,PSL,PS,PRECC,PRECL,PBLH,ICEFRAC,FSUTOA,FSNTOAC,FSNTOA,FSNTC,FSNT,FSNSC,FSNS,FSDSC,FSDS,FLUTC,FLUT,FLNTC,FLNT,FLNSC,FLNS,FLDSC,FLDS,CLDTOT,CLDMED,CLDLOW,CLDHGH raw/${exp}.cam2.h0.*.nc cat_${exp}.cam2.h0.${yrstart}-${yrend}.nc

#cdo yearmean -shifttime,-1mo cat_${exp}.cam2.h0.${yrstart}-${yrend}.nc annmean_${exp}.cam2.h0.${yrstart}-${yrend}.nc


