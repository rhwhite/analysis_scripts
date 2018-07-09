#!/bin/sh

Experiments=("CESMnoTf19") #"CESMnoT2f19" "CESMnoT4f19" "CESMtopof19" "CESMnotopof19")

AMJ=1
ON=0

for experi in ${Experiments[@]}; do
	cd /home/disk/rachel/CESM_outfiles/${experi}/atm/hist/
	if [ $AMJ -eq 1 ]; then
		echo ncea ClimMon04_${experi}_2-41_Xpl_cam2h0.nc ClimMon05_${experi}_2-41_Xpl_cam2h0.nc ClimMon06_${experi}_2-41_Xpl_cam2h0.nc AMJ2_40_mean_Xpl_${experi}.cam2.h0.nc
		ncea ClimMon04_${experi}_2-41_Xpl_cam2h0.nc ClimMon05_${experi}_2-41_Xpl_cam2h0.nc ClimMon06_${experi}_2-41_Xpl_cam2h0.nc AMJ2_40_mean_Xpl_${experi}.cam2.h0.nc

		ncea ClimMon04_${experi}_2-41_cam2h0.nc ClimMon05_${experi}_2-41_cam2h0.nc ClimMon06_${experi}_2-41_cam2h0.nc AMJ2_40_mean_${experi}.cam2.h0.nc
	fi

	if [ $ON -eq 1 ]; then
		ncea ClimMon10_${experi}_2-41_Xpl_cam2h0.nc ClimMon11_${experi}_2-41_Xpl_cam2h0.nc ON2_40_mean_Xpl_${experi}.cam2.h0.nc

                ncea ClimMon10_${experi}_2-41_cam2h0.nc ClimMon11_${experi}_2-41_cam2h0.nc ON2_40_mean_${experi}.cam2.h0.nc

	fi
done 

