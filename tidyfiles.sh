#!/bin/sh

names=("CAM4POP_f45g37" "CAM4POP_f45g37_10At25_55N_1" "CAM4POP_f45g37_10At25_55N_2" "CAM4POP_f45g37_10At25_55N_3" "CAM4POP_f45g37_10At25_55N_4" "CAM4POP_f45g37_10All_75S-45S_0" "CAM4POP_f45g37_10All_75S-0S_0")

for iname in ${names[@]}; do
	mv *Avg*_${iname}.cam2.h0.0*.nc ${iname}/atm/
        mv *Avg*_${iname}.pop.h.0*.nc ${iname}/ocn/
        mv ${iname}.pop.h.0*.nc ${iname}/ocn/

done


