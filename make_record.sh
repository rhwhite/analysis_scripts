#!/bin/sh

names=("WACCM_f19_CTL")
months=("06" "07" "08")
startyr=2
endyr=41

for iname in ${names[@]}; do
    cd /home/disk/eos4/rachel/CESM_outfiles/HYAK/WACCM_f19_CTL/atm/hist/TakNak/
    #cd /home/disk/eos4/rachel/Obs/ERAI/Daily/
 
    for yr in `seq -f "%04g" $startyr $endyr`; do    
        for mon in ${months[@]}; do
            ncks -O --mk_rec_dmn time TakNakfluxes_daily_${yr}-${mon}_WACCM_f19_CTL.cam2.h2.nc TakNakfluxes_daily_${yr}-${mon}_WACCM_f19_CTL.cam2.h2.nc
        done
    done

    ncrcat TakNakfluxes_daily_????-06_WACCM_f19_CTL.cam2.h2.nc TakNakfluxes_daily_????-07_WACCM_f19_CTL.cam2.h2.nc TakNakfluxes_daily_????-08_WACCM_f19_CTL.cam2.h2.nc TakNakfluxes_dailyclim_JJA_${startyr}-${endyr}_WACCM_f19_CTL.cam2.h2.nc
done


