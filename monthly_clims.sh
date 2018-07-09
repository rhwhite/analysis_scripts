#!/bin/sh

startyrs=(411)  #500
endyrs=(440)  #529
#startyrs=(32)
#endyrs=(71)

Experiments=("CAM4POP_f19g16C_noR") #"CESMtopof19" "CESMnoT2f19") #CAM4POP_B1850_NCAR")	#("CAM4POP_f19g16_NoTopo")

stmonth=1
endmonth=12

annual=1
months=1
std=0
surf=1
modlev=0
ocn=0
atm=1

expcount=0
for experi in ${Experiments[@]}; do
	startyr=${startyrs[$expcount]}
	endyr=${endyrs[$expcount]}
	expcount=$((expcount+1))
	echo $startyr $endyr
	if [ "$experi" = "CAM4POP_B1850_NCAR" ]; then
		experiname=b40.1850.track1.2deg.003
	else
		experiname=$experi
	fi

	if [ $atm -eq 1 ]; then
		cd ${experi}/atm/hist/raw/
		mkdir Temp
		for (( iyear=$startyr; iyear<=$endyr; iyear ++ )); do
			echo mv ${experiname}.cam2.h0.`printf %04g $iyear`*.nc Temp
			mv ${experiname}.cam2.h0.`printf %04g $iyear`*.nc Temp
		done

		pwd	

                cd ../	# get out of raw output directory

		pwd
                if [ $months -eq 1 ]; then

			for month in `seq -f "%02g" $stmonth $endmonth`; do

				if [ $std -eq 1 ]; then 
					echo ncra raw/Temp/${experiname}.cam2.h0.*-${month}_pl.nc Clim_${experiname}.cam.h0.${startyr}-${endyr}_${month}_pl.nc
					ncra raw/Temp/${experiname}.cam2.h0.*-${month}_pl.nc Clim_${experiname}.cam.h0.${startyr}-${endyr}_${month}_pl.nc
					ncrcat raw/Temp/${experiname}.cam2.h0.*-${month}_pl.nc Mon_${experiname}.cam.h0.${startyr}-${endyr}_${month}_pl.nc
				fi


                                if [ $modlev -eq 1 ]; then

                                        echo ncrcat -v VT,VTH2d,VTH3d,VQ raw/Temp/${experiname}.cam2.h0.*${month}.nc Mon_${experiname}_${startyr}-${endyr}_${month}.nc
					ncrcat -v VT,VTH2d,VTH3d,VQ raw/Temp/${experiname}.cam2.h0.*${month}.nc Clim_${experiname}_${startyr}-${endyr}_${month}.nc
                                        ncrcat -v VT,VTH2d,VTH3d,VQ raw/Temp/${experiname}.cam2.h0.*${month}.nc Mon_${experiname}_${startyr}-${endyr}_${month}.nc

                                fi
				if [ $surf -eq 1 ]; then
			
					echo ncrcat -v TS,PRECT,TREFHT,U10 raw/Temp/${experiname}.cam2.h0.*${month}.nc Surf_${experiname}_${startyr}-${endyr}_${month}.nc
					ncrcat -v U10,TS,TROP_Z,PRECT,TREFHT,TAUX,TAUY,SRFRAD,SHFLX,PS,PHIS,LHFLX,FSNTOA,FSNT,FSNS,FLNT,FLNS raw/Temp/${experiname}.cam2.h0.*${month}.nc Surf_${experiname}_${startyr}-${endyr}_${month}.nc
					ncra -v U10,TS,TROP_Z,PRECT,TREFHT,TAUX,TAUY,SRFRAD,SHFLX,PS,PHIS,LHFLX,FSNTOA,FSNT,FSNS,FLNT,FLNS raw/Temp/${experiname}.cam2.h0.*${month}.nc ClimSurf_${experiname}_${startyr}-${endyr}_${month}.nc
				fi
			done
		fi

		if [ $annual -eq 1 ]; then
			ncra raw/Temp/${experiname}.cam2.h0.*-*.nc ClimAnn_${experiname}.cam.h0.${startyr}-${endyr}.nc
		fi
		mv raw/Temp/* raw/
		cd ../../../	# get back out to original directory
	fi
        if [ $ocn -eq 1 ]; then
                cd ${experi}/ocn/hist/raw/
                mkdir Temp
		pwd
                for (( iyear=$startyr; iyear<=$endyr; iyear ++ )); do
                        mv ${experiname}.pop.h.`printf %04g $iyear`*.nc Temp
                done

                cd ../ 	# get out of raw directory
		if [ $months -eq 1 ]; then
			for month in `seq -f "%02g" $stmonth $endmonth`; do

				if [ $std -eq 1 ]; then
					echo ncra raw/Temp/${experiname}.pop.h.*-${month}.nc Clim_${experiname}.pop.h.${startyr}-${endyr}_${month}.nc
					ncra raw/Temp/${experiname}.pop.h.*-${month}.nc Clim_${experiname}.pop.h.${startyr}-${endyr}_${month}.nc
				fi

			done
		fi
		if [ $annual -eq 1 ]; then
			echo ncra raw/Temp/${experiname}.pop.h.*-*.nc ClimAnn_${experiname}.pop.h.${startyr}-${endyr}.nc
	                ncra raw/Temp/${experiname}.pop.h.*-*.nc ClimAnn_${experiname}.pop.h.${startyr}-${endyr}.nc
		fi
	        mv raw/Temp/* raw/
		cd ../../../
        fi
done 

