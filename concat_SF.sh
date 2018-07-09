#!/bin/sh

startyr=2
endyr=71
dir=/home/disk/eos4/rachel/CESM_outfiles/
Experiments=("CESMtopof19") # "CESMnoT2f19" "CESMnoT4f19" "CESMtopof19")

SFZA=1
Surf=0

for experi in ${Experiments[@]}; do
        cd ${dir}${experi}/atm/hist/raw

	mkdir Temp
	for (( iyear=$startyr; iyear<=$endyr; iyear ++ )); do
		mv ${experi}.cam2.h0.`printf %04g $iyear`*_pl.nc Temp
                mv ${experi}.cam2.h0.`printf %04g $iyear`*.nc Temp
	done
	
	for month in `seq -f "%02g" 1 12`; do

		if [ $SFZA -eq 1 ]; then 
			#echo ncrcat -v SFza Temp/${experi}.cam2.h0.*${month}_pl.nc ../SFza_${experi}_${startyr}-${endyr}_${month}.nc
			ncrcat -v SFza Temp/${experi}.cam2.h0.*${month}_pl.nc ../SFza_${experi}_${startyr}-${endyr}_${month}.nc
		fi

		if [ $Surf -eq 1 ]; then
	
			#echo ncrcat -v TS,PRECT,TREFHT,U10 Temp/${experi}.cam2.h0.*${month}.nc ../Surf_${experi}_${startyr}-${endyr}_${month}.nc
			ncrcat -v TS,PRECT,TREFHT,U10 Temp/${experi}.cam2.h0.*${month}.nc ../Surf_${experi}_${startyr}-${endyr}_${month}.nc
		fi
	done
	mv Temp/* .
done 

