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

yrstart=1
yrend=200

climstart=101
climend=200

version="cam"

printf -v styrstart "%04d" $yrstart
printf -v styrend "%04d" $yrend

printf -v stclimstart "%04d" $climstart
printf -v stclimend "%04d" $climend

echo $styrstart

# Create directory structure and move files
cd ${dir}/${exp}/ocn/hist/
mkdir raw

mv ${exp}.pop.h.*.nc raw/
cd raw
mkdir toavg

# move only required files into toavg
for iyear in `seq -f "%04g" $startyear $endyear`; do
    # record average across all monthly files for that year    
    mv ${exp}.pop.h.????-??.nc toavg
done

ncrcat -O -d z_t,0,0 -d transport_reg,1,1 -v MOC,HMXL,RHO,SALT,TEMP,BSF,N_HEAT,N_SALT,TAUX,TAUY,TAREA,ROFF_F toavg/${exp}.pop.h.0* ../cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc

# move everything back out of toavg directory
mv toavg/* .

# move out of raw directory
cd ../
# Average over transport_ref dimension to get rid of it
ncwa -O -a transport_reg cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc

# Calculate annual mean with Shifttime to get correct months (otherwise cdo reads end of Jan as 1-Feb)
cdo yearmean -shifttime,-1mo cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc annmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc

# Make monthly averages
# Need to shifttime BEFORE selecting month, which means putting it AFTER in the command
cdo yearmean -selmon,1,2,12 -shifttime,-1mo cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc DJFmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc 
cdo yearmean -selmon,3,4,5 -shifttime,-1mo cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc MAMmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc
cdo yearmean -selmon,6,7,8 -shifttime,-1mo cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc JJAmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc
cdo yearmean -selmon,9,10,11 -shifttime,-1mo cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc SONmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc

# Calculate climatologies
ncra -d time,${stclimstart},${stclimend} annmean_limited_CAM4POP_f19g16C_noMT.pop.h.${styrstart}-${styrend}.nc ANNclim_CAM4POP_f19g16C_noMT.pop.h.${stclimstart}-${stclimend}.nc
ncra -d time,${stclimstart},${stclimend} DJF_limited_CAM4POP_f19g16C_noMT.pop.h.${styrstart}-${styrend}.nc DJFclim_CAM4POP_f19g16C_noMT.pop.h.${stclimstart}-${stclimend}.nc
ncra -d time,${stclimstart},${stclimend} MAM_limited_CAM4POP_f19g16C_noMT.pop.h.${styrstart}-${styrend}.nc MAMclim_CAM4POP_f19g16C_noMT.pop.h.${stclimstart}-${stclimend}.nc
ncra -d time,${stclimstart},${stclimend} JJA_limited_CAM4POP_f19g16C_noMT.pop.h.${styrstart}-${styrend}.nc JJAclim_CAM4POP_f19g16C_noMT.pop.h.${stclimstart}-${stclimend}.nc
ncra -d time,${stclimstart},${stclimend} SON_limited_CAM4POP_f19g16C_noMT.pop.h.${styrstart}-${styrend}.nc SONclim_CAM4POP_f19g16C_noMT.pop.h.${stclimstart}-${stclimend}.nc


# Regrid
export NCL_dir=${dir}/${exp}/ocn/hist/
export NCL_file='clim_CAM4POP_f19g16C_noMT.pop.h.0001-0300.nc'

# Need to change directory so NCL can find the weights files
#cd /home/disk/eos4/rachel/git/NCL/Regrid/
#ncl regrid_POP_seas.ncl



# Now do atmospheric variables
cd ${dir}/${exp}/atm/hist/

# Create directory structure and move files
mkdir raw

mv ${exp}.${version}.h*.nc raw/
cd raw
mkdir toavg

# move only required files into toavg
for iyear in `seq -f "%04g" $startyear $endyear`; do
    # record average across all monthly files for that year    
    mv ${exp}.${version}.h0.????-??.nc toavg
done

ncrcat -v PRECT,TS,TREFHT,TMQ,TAUX,TAUY,SHFLX,LHFLX,RHREFHT,QFLX,PSL,PS,PRECC,PRECL,PBLH,ICEFRAC,FSUTOA,FSNTOAC,FSNTOA,FSNTC,FSNT,FSNSC,FSNS,FSDSC,FSDS,FLUTC,FLUT,FLNTC,FLNT,FLNSC,FLNS,FLDSC,FLDS,CLDTOT,CLDMED,CLDLOW,CLDHGH toavg/${exp}.cam*.h0.*.nc ../cat_${exp}.${version}.h0.${styrstart}-${styrend}.nc

# move everything back out of toavg directory
mv toavg/* .

# move out of raw directory
cd ../

cdo yearmean -shifttime,-1mo cat_${exp}.${version}.h0.${styrstart}-${styrend}.nc annmean_${exp}.${version}.h0.${styrstart}-${styrend}.nc


