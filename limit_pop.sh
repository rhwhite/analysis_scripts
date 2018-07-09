#!/bin/sh

while getopts ":d:x:s:e:t:f:" opt; do
    case $opt in
    # Main directory of experiments
    d) dir="$OPTARG"
    ;;
    # Name of experiment
    x) exp="$OPTARG"
    ;;
    # start and end year for annual and running mean analysis
    s) startyear="$OPTARG"
    ;;
    e) endyear="$OPTARG"
    ;;
    # start and end year for climatological average
    t) startyearkeep="$OPTARG"
    ;;
    f) endyearkeep="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG"
    ;;
    esac
done

stmonth=1
endmonth=12

if [ $endyear -gt $startyearkeep ]; then 
    echo incorrect years, endyear is greater than startyearkeep and will delete files you want to keep!
    exit 1
fi
echo ${dir}/${exp}/ocn/hist/raw/
cd ${dir}/${exp}/ocn/hist/raw/

echo creating limited monthly means, and removing originals for spin-up period $startyear - $endyear
for yr in `seq -f "%04g" $startyear $endyear`; do
    echo $yr
    
    for month in `seq -f "%02g" $stmonth $endmonth`; do
        ncks -v N_SALT,N_HEAT,MOC,BSF,HBLT,HMXL,Q,PV,TEMP,SALT,RHO,RHO_VINT,SSH,SHF,SHF_QSW,SFWF,TAUX,TAUY,FW,ROFF_F,SENH_F,LWUP_F,LWDN_F,MELTH_F,IFRAC ${exp}.pop.h.${yr}-${month}.nc limited_${exp}.pop.h.${yr}-${month}.nc && rm ${exp}.pop.h.${yr}-${month}.nc

    done
done

echo creating limited monthly means, and KEEPING originals for final period $startyearkeep - $endyearkeep
for yr in `seq -f "%04g" $startyearkeep $endyearkeep`; do
    echo $yr

    for month in `seq -f "%02g" $stmonth $endmonth`; do
        ncks -v N_SALT,N_HEAT,MOC,BSF,HBLT,HMXL,Q,PV,TEMP,SALT,RHO,RHO_VINT,SSH,SHF,SHF_QSW,SFWF,TAUX,TAUY,FW,ROFF_F,SENH_F,LWUP_F,LWDN_F,MELTH_F,IFRAC ${exp}.pop.h.${yr}-${month}.nc limited_${exp}.pop.h.${yr}-${month}.nc

    done
done


