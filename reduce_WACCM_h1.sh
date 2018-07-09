#!/bin/sh
# Script to calculate annual averages, climatological averages and running
# averages on raw monthly data from CESM (coupled) 
#
# Author: Rachel White, rachel.white@cantab.net
#
# Created: March 2018

dir=/home/disk/eos4/rachel/CESM_outfiles/HYAK/
exp=WACCM_f19_LGM
startyr=15
endyr=25


# Get arguments
while getopts ":d:x:s:e:" opt; do
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
    esac
done

#Loop through years
iyear=$startyear

echo ${dir}/${exp}/atm/hist/raw/
cd ${dir}/${exp}/atm/hist/raw/


for iyear in `seq -f "%04g" $startyear $endyear`; do
    cd ${dir}/${exp}/atm/hist/raw/
    # Remove 3D fields
    ncks -O -x -v QRS_TOT,T ${exp}.cam2.h1.${iyear}-01-01-21600.nc ${exp}.cam2.h1.${iyear}-01-01-21600.nc

done



