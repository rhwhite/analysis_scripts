#!/bin/sh

startdec=1
enddec=4
Experiments=("CAM4POP_f45g37" "CAM4POP_f45g37_10All_75S-0S_0" "CAM4POP_f45g37_10All_75S-45S_0" "CAM4POP_f45g37_10At25_55N_1" "CAM4POP_f45g37_10At25_55N_2" "CAM4POP_f45g37_10At25_55N_3" "CAM4POP_f45g37_10At25_55N_4")

for experi in ${Experiments[@]}; do
	echo $experi
	cd /home/disk/eos4/rachel/CESM_outfiles/HYAK/FluxEXP/$experi
	#mkdir Extra
	#mkdir ocn
	#mkdir atm
	#mv *POP_f45g37*.pop.* ocn
        #mv *POP_f45g37*.cam2.* atm

	cd ocn
	#mkdir Temp
	#mv *.nc Temp
        for ((idec=$enddec;idec>=$startdec;idec--)); do
		echo $idec
		mv Temp/CAM4POP_f45g37*pop.h.00$idec* .
		ncra CAM4POP_f45g37*pop.h.* ClimAvg_${experi}.pop.h.00${idec}0-00${enddec}9.nc

	done
done 

