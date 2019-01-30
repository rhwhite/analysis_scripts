#!/bin/sh
# Script to calculate annual averages, climatological averages and running
# averages on raw monthly data from CESM (coupled) 
#
# Author: Rachel White, rachel.white@cantab.net
#
# Created: June 2018

# get top level for TEMP, RHO, SALT 
# for MOC get only Atlantic
# get specific variables and concatenate

module load cdo

dir="/glade/scratch/rachelwh/archive/"
exp="CAM4POP_NoR_f19"

yrstart=1
yrend=300

climstart=251
climend=300

version="cam"

runocn=1
runatm=1

printf -v styrstart "%04d" $yrstart
printf -v styrend "%04d" $yrend

printf -v stclimstart "%04d" $climstart
printf -v stclimend "%04d" $climend

# If runocn = 1 then run the ocean analysis:
if [ $runocn -eq 1 ]; then
    # Create directory structure and move files
    cd ${dir}/${exp}/ocn/hist/
    mkdir raw

    mv ${exp}.pop.h.*.nc raw/
    cd raw
    mkdir toavg
    mv toavg/* .

    # move only required files into toavg
    for iyear in `seq -f "%04g" $yrstart $yrend`; do
        # record average across all monthly files for that year    
        mv limited_${exp}.pop.h.${iyear}-??.nc toavg
    done

    ncrcat -O -d z_t,0,0 -d transport_reg,1,1 -v MOC,HMXL,RHO,SALT,TEMP,BSF,N_HEAT,N_SALT,TAUX,TAUY,ROFF_F toavg/limited_${exp}.pop.h.0* ../cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc

    ncrcat -O -d z_t,0,0 -d transport_reg,0,0 -v MOC toavg/limited_${exp}.pop.h.0* ../cat_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc

    # move everything back out of toavg directory
    mv toavg/* .

    # move out of raw directory
    cd ../
    # Average over transport_ref dimension to get rid of it
    ncwa -O -a transport_reg cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc
    ncwa -O -a transport_reg cat_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc cat_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc

    # Calculate annual mean with Shifttime to get correct months (otherwise cdo reads end of Jan as 1-Feb)
    cdo yearmonmean -shifttime,-1mo cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc annmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc
    cdo yearmonmean -shifttime,-1mo cat_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc annmean_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc

    # Make monthly averages, using yearmonmean to take into account monthly weighting
    # Need to shifttime BEFORE selecting month, which means putting it AFTER in the command
    cdo yearmonmean -selmon,1,2,12 -shifttime,-1mo cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc DJFmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc 
    cdo yearmonmean -selmon,3,4,5 -shifttime,-1mo cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc MAMmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc
    cdo yearmonmean -selmon,6,7,8 -shifttime,-1mo cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc JJAmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc
    cdo yearmonmean -selmon,9,10,11 -shifttime,-1mo cat_limited_${exp}.pop.h.${styrstart}-${styrend}.nc SONmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc

    cdo yearmonmean -selmon,1,2,12 -shifttime,-1mo cat_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc DJFmean_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc
    cdo yearmonmean -selmon,3,4,5 -shifttime,-1mo cat_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc MAMmean_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc
    cdo yearmonmean -selmon,6,7,8 -shifttime,-1mo cat_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc JJAmean_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc
    cdo yearmonmean -selmon,9,10,11 -shifttime,-1mo cat_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc SONmean_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc

    # Calculate climatologies for seasons - each season weighted equally here
    ncra -d time,${stclimstart},${stclimend} annmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc ANNclim_${exp}.pop.h.${stclimstart}-${stclimend}.nc
    ncra -d time,${stclimstart},${stclimend} DJFmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc DJFclim_${exp}.pop.h.${stclimstart}-${stclimend}.nc
    ncra -d time,${stclimstart},${stclimend} MAMmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc MAMclim_${exp}.pop.h.${stclimstart}-${stclimend}.nc
    ncra -d time,${stclimstart},${stclimend} JJAmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc JJAclim_${exp}.pop.h.${stclimstart}-${stclimend}.nc
    ncra -d time,${stclimstart},${stclimend} SONmean_limited_${exp}.pop.h.${styrstart}-${styrend}.nc SONclim_${exp}.pop.h.${stclimstart}-${stclimend}.nc

    ncra -d time,${stclimstart},${stclimend} annmean_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc ANNclim_MOC_${exp}.pop.h.${stclimstart}-${stclimend}.nc
    ncra -d time,${stclimstart},${stclimend} DJFmean_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc DJFclim_MOC_${exp}.pop.h.${stclimstart}-${stclimend}.nc
    ncra -d time,${stclimstart},${stclimend} MAMmean_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc MAMclim_MOC_${exp}.pop.h.${stclimstart}-${stclimend}.nc
    ncra -d time,${stclimstart},${stclimend} JJAmean_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc JJAclim_MOC_${exp}.pop.h.${stclimstart}-${stclimend}.nc
    ncra -d time,${stclimstart},${stclimend} SONmean_MOC_${exp}.pop.h.${styrstart}-${styrend}.nc SONclim_MOC_${exp}.pop.h.${stclimstart}-${stclimend}.nc


    # Regrid
    export NCL_dir=${dir}/${exp}/ocn/hist/
    export NCL_file='clim_${exp}.pop.h.0001-0300.nc'

    # Need to change directory so NCL can find the weights files
    #cd /home/disk/eos4/rachel/git/NCL/Regrid/
    #ncl regrid_POP_seas.ncl
fi


# Now if runatm=1 then do atmospheric variables
if [ $runatm -eq 1 ]; then
    cd ${dir}/${exp}/atm/hist/

    # Create directory structure and move files
    mkdir raw

    mv ${exp}.${version}.h*.nc raw/
    cd raw
    mkdir toavg
    mv toavg/* .

    # move only required files into toavg
    for iyear in `seq -f "%04g" $yrstart $yrend`; do
        # record average across all monthly files for that year    
        mv ${exp}.${version}.h0.${iyear}-??.nc toavg
    done

    #ncrcat -v PRECT,TS,TREFHT,TMQ,TAUX,TAUY,SHFLX,LHFLX,RHREFHT,QFLX,PSL,PS,PRECC,PRECL,PBLH,ICEFRAC,FSUTOA,FSNTOAC,FSNTOA,FSNTC,FSNT,FSNSC,FSNS,FSDSC,FSDS,FLUTC,FLUT,FLNTC,FLNT,FLNSC,FLNS,FLDSC,FLDS,CLDTOT,CLDMED,CLDLOW,CLDHGH toavg/${exp}.cam*.h0.*.nc ../cat_${exp}.${version}.h0.${styrstart}-${styrend}.nc

    ncrcat -v TS,TREFHT,TMQ,TAUX,TAUY,SHFLX,LHFLX,QFLX,PSL,PS,PRECC,PRECL,PBLH,ICEFRAC,FSUTOA,FSNTOAC,FSNTOA,FSNTC,FSNT,FSNSC,FSNS,FSDSC,FSDS,FLUTC,FLUT,FLNTC,FLNT,FLNSC,FLNS,FLDSC,FLDS,CLDTOT,CLDMED,CLDLOW,CLDHGH toavg/${exp}.cam*.h0.*.nc ../cat_${exp}.${version}.h0.${styrstart}-${styrend}.nc

    # move everything back out of toavg directory
    mv toavg/* .

    # move out of raw directory
    cd ../

    cdo yearmonmean -shifttime,-1mo cat_${exp}.${version}.h0.${styrstart}-${styrend}.nc annmean_${exp}.${version}.h0.${styrstart}-${styrend}.nc

fi
