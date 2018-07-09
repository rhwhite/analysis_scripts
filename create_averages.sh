#!/bin/sh
# Script to calculate annual averages, climatological averages and running
# averages on raw monthly data from CESM (coupled) 
#
# Author: Rachel White, rachel.white@cantab.net
#
# Created: June 2016

# Set defaults: do nothing!
ifatm=0
ifocn=0

# If ocean data is monthly (rather than annually) set this to 1:
ifocnmon=0

# Set AnnAvg to 1 if you want to calculate annual averages for atm or ocean 
# data from monthly averages
AnnAvg=0

# Set ClimAvg to 1 if you want to calculate climatological averages from annual
# averages, for years startyearclim to endyearclim 
ClimAvg=0
# Set RunAvg to 1 if you want to calculate running averages using runmean years
# You can specify ocean and atm separately using ifatm and ifocn
RunAvg=0

# Get arguments
while getopts ":d:x:s:t:e:f:m:p:n:c:r:a:o:z:" opt; do
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
    t) startyearclim="$OPTARG"
    ;;
    f) endyearclim="$OPTARG"
    ;;
    # Number of years to include in the running mean
    m) runmean="$OPTARG"
    ;;
    # Number of years of spinup for ocean
    p) ocnspin="$OPTARG"
    ;;
    #) Annual averages
    n) AnnAvg="$OPTARG"
    ;;
    # Climatology
    c) ClimAvg="$OPTARG"
    ;;
    #) Running means
    r) RunAvg="$OPTARG"
    ;;
    #) Do atmosphere?
    a) ifatm="$OPTARG"
    ;;
    #) Do ocean?
    o) ifocn="$OPTARG"
    ;;
    # Is ocean data monthly averages?
    z) ifocnmon="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG"
    ;;
    esac
done

echo Running create_yearly_avg.sh for ${exp} in ${dir};
if [ $ifatm -eq 1 ]; then echo computing atmospheric values; fi 
if [ $ifocn -eq 1 ]; then echo computing oceanic values; fi

if [ $AnnAvg -eq 1 ]; then echo for annual averages; fi
if [ $RunAvg -eq 1 ]; then echo for running means of ${runmean} years; fi
if [ $ClimAvg -eq 1 ]; then echo climatology of $startyearclim to $endyearclim; fi

# get strings for climatological averages
printf -v startyearclimS "%03g" $startyearclim
printf -v endyearS "%03g" $endyearclim

iyear=$startyear

if [ $AnnAvg -eq 1 ]; then
    if [ $ifatm -eq 1 ]; then

       
        cd ${dir}/${exp}/atm/hist/
        echo pwd
        mkdir raw
        mv ${exp}.cam2.h0.*.nc raw/
        cd raw

        for iyear in `seq -f "%04g" $startyear $endyear`; do
            # record average across all monthly files for that year    
            ncra -O ${exp}.cam2.h0.${iyear}-??.nc AnnAvg_${exp}.cam2.h0.${iyear}.nc
        done
    fi

    if [ $ifocn -eq 1 ]; then
        if [ $ifocnmon -eq 1 ]; then
            cd ${dir}/${exp}/ocn/hist/
            mkdir raw/
            mv limited* raw/
            mv ${exp}.pop.h.*.nc raw/
            cd raw/
            for iyear in `seq -f "%04g" $startyear $endyear`; do
                # record average across all monthly files for that year    
                if [ $iyear -le $ocnspin ]; then
                    ncra -O limited_${exp}.pop.h.${iyear}-??.nc AnnAvg_${exp}.pop.h.${iyear}.nc
                else
                    ncra -O ${exp}.pop.h.${iyear}-??.nc AnnAvg_${exp}.pop.h.${iyear}.nc
                fi
            done
        fi
    fi
fi

if [ $ClimAvg -eq 1 ]; then
    if [ $ifatm -eq 1 ]; then
        cd ${dir}/${exp}/atm/hist/raw/

        mkdir Temp
        mv Temp/* .
        mv AnnAvg_${exp}.cam2.h0.*.nc Temp
        
        for (( iyear=$startyearclim; iyear<=$endyearclim; iyear ++ )); do
            mv Temp/AnnAvg_${exp}.cam2.h0.0`printf %03g $iyear`.nc .
        done

        ncra -O AnnAvg_${exp}.cam2.h0.0*.nc ClimAvg_${exp}.cam2.h0.0`printf %03g $startyearclimS`-0`printf %03g $endyearclim`.nc
        mv Temp/* .
    fi

    cd ${dir}${exp}/ocn/hist/raw/

    if [ $ifocn -eq 1 ]; then
        mkdir Temp
        mv Temp/* .
        if [ $ifocnmon -eq 1 ]; then
            mv AnnAvg_${exp}.pop.h.0*.nc Temp
            for (( iyear=$startyearclim; iyear<=$endyear; iyear ++ )); do
                mv Temp/AnnAvg_${exp}.pop.h.0`printf %03g ${iyear}`.nc .
            done

            ncra -O AnnAvg_${exp}.pop.h.0*.nc ClimAvg_${exp}.pop.h.0`printf %03g ${startyearclim}`-0`printf %03g ${endyearclim}`.nc
            mv Temp/* .
        else
            mv ${exp}.pop.h.*.nc Temp
            for (( iyear=$startyearclim; iyear<=$endyear; iyear ++ )); do
                mv Temp/${exp}.pop.h.0`printf %03g ${iyear}`.nc .
            done

            ncra -O ${exp}.pop.h.0*.nc ClimAvg_${exp}.pop.h.0`printf %03g ${startyearclim}`-0`printf %03g ${endyearclim}`.nc
            mv Temp/* .

        fi
    fi
fi

if [ $RunAvg -eq 1 ]; then
    if [ $ifatm -eq 1 ]; then
        cd ${dir}/${exp}/atm/hist/

        mkdir Temp
        mv Temp/* .
        mv AnnAvg_${exp}.cam2.h0.*.nc Temp
        for (( iyear=$startyear; iyear<$startyear+$runmean; iyear++ )); do
            mv Temp/AnnAvg_${exp}.cam2.h0.0`printf %03g ${iyear}`.nc .
        done
        (( avgyear = ($startyear + $startyear+$runmean-1)/2 ))	
        ncra -O AnnAvg_${exp}.cam2.h0.0*.nc tRunAvg${runmean}_${exp}.cam2.h0.0`printf %03g ${avgyear}`.nc

        for (( firstyear=$startyear; firstyear<$endyear-$runmean; firstyear++ )) do
            (( newyear = $firstyear+$runmean ))
            (( avgyear = ($firstyear+1+$newyear)/2 ))
            mv AnnAvg_${exp}.cam2.h0.0`printf %03g ${firstyear}`.nc Temp
            mv Temp/AnnAvg_${exp}.cam2.h0.0`printf %03g ${newyear}`.nc .
            ncra -O AnnAvg_${exp}.cam2.h0.0*.nc tRunAvg${runmean}_${exp}.cam2.h0.0`printf %03g ${avgyear}`.nc
        done

        ncrcat tRunAvg${runmean}_${exp}.cam2.h0.0*.nc RunAvg${runmean}_${exp}.cam2.h0.0`printf %03g ${startyear}`-0`printf %03g ${endyear}`.nc
        rm -f tRunAvg${runmean}_${exp}.cam2.h0.0*.nc

        mv Temp/* .
    fi


    if [ $ifocn -eq 1 ]; then
        cd ${dir}/${exp}/ocn/hist/
        mkdir Temp
        mv Temp/* .
        #echo ${exp}*.nc Temp

        if [ $ifocnmon -eq 1 ]; then
            mv AnnAvg_${exp}*.nc Temp
            for (( iyear=$startyear; iyear<$startyear+$runmean; iyear++ )); do
                mv Temp/AnnAvg_${exp}.pop.h.0`printf %03g ${iyear}`.nc .
            done
            (( avgyear = ($startyear + $startyear+$runmean-1)/2 ))
            ncra -O AnnAvg_${exp}.pop.h.0*.nc tRunAvg${runmean}_${exp}.pop.h.0`printf %03g ${avgyear}`.nc

            for (( firstyear=$startyear; firstyear<$endyear-$runmean; firstyear++ )) do
                (( newyear = $firstyear+$runmean ))
                (( avgyear = ($firstyear+1+$newyear)/2 ))
                mv AnnAvg_${exp}.pop.h.0`printf %03g ${firstyear}`.nc Temp
                mv Temp/AnnAvg_${exp}.pop.h.0`printf %03g ${newyear}`.nc .
                ncra -O AnnAvg_${exp}.pop.h.0*.nc tRunAvg${runmean}_${exp}.pop.h.0`printf %03g ${avgyear}`.nc
            done

            ncrcat tRunAvg${runmean}_${exp}.pop.h.0*.nc RunAvg${runmean}_${exp}.pop.h.0`printf %03g ${startyear}`-0`printf %03g ${endyear}`.nc
            rm -f tRunAvg${runmean}_${exp}.pop.h.0*.nc

        else
            mv ${exp}*.nc Temp
            for (( iyear=$startyear; iyear<$startyear+$runmean; iyear++ )); do
                mv Temp/${exp}.pop.h.0`printf %03g ${iyear}`.nc .
            done
            (( avgyear = ($startyear + $startyear+$runmean-1)/2 ))
            ncra -O ${exp}.pop.h.0*.nc tRunAvg${runmean}_${exp}.pop.h.0`printf %03g ${avgyear}`.nc

            for (( firstyear=$startyear; firstyear<$endyear-$runmean; firstyear++ )) do
                (( newyear = $firstyear+$runmean ))
                (( avgyear = ($firstyear+1+$newyear)/2 ))
                mv ${exp}.pop.h.0`printf %03g ${firstyear}`.nc Temp
                mv Temp/${exp}.pop.h.0`printf %03g ${newyear}`.nc .
                ncra -O ${exp}.pop.h.0*.nc tRunAvg${runmean}_${exp}.pop.h.0`printf %03g ${avgyear}`.nc
            done
            
            ncrcat tRunAvg${runmean}_${exp}.pop.h.0*.nc RunAvg${runmean}_${exp}.pop.h.0`printf %03g ${startyear}`-0`printf %03g ${endyear}`.nc
            rm -f tRunAvg${runmean}_${exp}.pop.h.0*.nc
        
        fi
            mv Temp/* .

    fi

fi


