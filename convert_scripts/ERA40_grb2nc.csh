#!/bin/csh
module load intel-suite mpi
module load netcdf

cd /scratch/rw408/ERA40/1979/ap
echo 'cd /scratch/rw408/ERA40/1979/ap'
set year = 1979
set month = 1
set day = 1
set hour = 0


set dpm = (31 28 31 30 31 30 31 31 30 31 30 31)
set ii = 1
set six = 6
set zero = 0
echo 'all variables set'

while ($ii < 1461)
echo $ii
if ($month <10) then
   set smonth = $zero""$month
else 
   set smonth = $month
endif
if ($day <10) then
   set sday = $zero""$day
else
   set sday = $day
endif
if ($hour <12) then
   set shour = $zero""$hour
else
   set shour = $hour
endif
echo 'year' $year
echo 'month' $smonth
echo 'day' $sday
echo 'hour' $shour

ncl_convert2nc gpap$year$smonth$sday$shour""*.grb -u initial_time0_hours -U initial_time0_hours -itime
@ hour = $hour + 6
echo 'hour' $hour
if ($hour >18) then
echo 'hour great than 18' $hour
set hour = 0
@ day ++
endif
if ($day > $dpm[$month]) then
echo 'day greater than ' $dpm[$month] $day
set day = 1
@ month++
endif
if ($month >12) then
echo 'month great than 12' $month
set month = 1
@ year++
endif

@ ii++

end

ncl_convert2nc gpap????0229*.grb -u initial_time0_hours -U time initial_time0_hours -itime


