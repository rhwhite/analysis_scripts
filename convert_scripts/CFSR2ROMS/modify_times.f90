PROGRAM modify_times
! ______________________________________________
use netcdf
USE, INTRINSIC :: IEEE_ARITHMETIC
implicit none

interface
     FUNCTION jd(yy,mm,dd)
     implicit none
     INTEGER jd,yy,mm,dd
     end function jd
end interface

INTEGER iyear, itime, imonth
INTEGER year_st, year_end, month_st, month_end
INTEGER ncid_cfsr, status_cfsr
INTEGER LonDimID, LatDimID, TimeDimId,TimeVarId, FCTimeVarId
INTEGER nlon, nlat, ntime
INTEGER nvalue
INTEGER jdref,jdref_cfsr
INTEGER vars

REAL, ALLOCATABLE, DIMENSION(:,:) :: scr, scrp, emask
REAL, ALLOCATABLE, DIMENSION(:,:) :: error, work
INTEGER :: time, fctime
REAL :: tx
CHARACTER(len=160) pathin,filein
CHARACTER(len= 50):: ryear, rmonth, num_temp
REAL, PARAMETER :: undef = 2.E+35            ! Undefined land value
REAL, PARAMETER :: critx=0.01
REAL, PARAMETER :: cor=1.6
INTEGER, PARAMETER :: mxs=100
CHARACTER(len = 5), dimension(4) :: varnames

tx=0.9*undef

varnames(1) = 'dlw'
varnames(2) = 'dsw'
varnames(3) = 'rain'
varnames(4) = 'cloud'

OPEN(unit=10,file='modify_times.in')
READ(10,'(a)') pathin
READ(10,*) year_st, month_st 
READ(10,*) year_end, month_end
CLOSE(10)


DO vars = 1,4

iyear= year_st

DO imonth= month_st,12
 
   write(ryear,*) iyear
   ryear = adjustl(ryear)
   write(num_temp,*) imonth
   num_temp=ADJUSTL(num_temp)
   if (len_trim(num_temp).EQ.1) rmonth= '0'//TRIM(num_temp)
   if (len_trim(num_temp).EQ.2) rmonth= TRIM(num_temp)

   write(*,*) 'Treatment of year ',TRIM(ryear),' and month ',TRIM(rmonth),' in progress ...'
   
   filein=TRIM(pathin)//TRIM(varnames(vars))//'_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
 
   status_cfsr = nf90_open(TRIM(filein),nf90_write,ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   status_cfsr = nf90_inq_dimid(ncid_cfsr, "time", TimeDimID)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
 
   status_cfsr = nf90_Inquire_Dimension(ncid_cfsr,TimeDimID,len = ntime)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   
   
   jdref_cfsr=jd(iyear,imonth,1) 
   
   DO itime= 1,ntime
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_inq_varid(ncid_cfsr, "forecast_hour",FCTimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,FCTimeVarId,fctime,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

      time= time+1
      fctime = fctime + 1

      status_cfsr = nf90_put_var(ncid_cfsr,TimeVarId,time,(/itime/))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_put_var(ncid_cfsr,FCTimeVarId,fctime,(/itime/))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   

   ENDDO

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   write(*,*) 'Time modified by one hours for: ',filein
ENDDO

if(year_end.gt.year_st+1)then

DO iyear= year_st+1,year_end-1

   DO imonth= 1,12
 
      write(ryear,*) iyear
      ryear = adjustl(ryear)
      write(num_temp,*) imonth
      num_temp=ADJUSTL(num_temp)
      if (len_trim(num_temp).EQ.1) rmonth= '0'//TRIM(num_temp)
      if (len_trim(num_temp).EQ.2) rmonth= TRIM(num_temp)
      
      write(*,*) 'Treatment of year ',TRIM(ryear),' and month ',TRIM(rmonth),' in progress ...'
   
      filein=TRIM(pathin)//TRIM(varnames(vars))//'_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
      
      status_cfsr = nf90_open(TRIM(filein),nf90_write,ncid_cfsr)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      
      status_cfsr = nf90_inq_dimid(ncid_cfsr, "lon", LonDimID)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_inq_dimid(ncid_cfsr, "lat", LatDimID)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_inq_dimid(ncid_cfsr, "time", TimeDimID)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   
      status_cfsr = nf90_Inquire_Dimension(ncid_cfsr,LonDimID,len = nlon)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_Inquire_Dimension(ncid_cfsr,LatDimID,len = nlat)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_Inquire_Dimension(ncid_cfsr,TimeDimID,len = ntime)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      
   
      jdref_cfsr=jd(iyear,imonth,1) 
   
      DO itime= 1,ntime
         ! _________________________________________________________________
         ! Read in Time
         status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_inq_varid(ncid_cfsr, "forecast_hour",FCTimeVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         
         status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,FCTimeVarId,fctime,start=(/ itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

         time= time+1
         fctime = fctime + 1
         
         status_cfsr = nf90_put_var(ncid_cfsr,TimeVarId,time,(/itime/))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_put_var(ncid_cfsr,FCTimeVarId,fctime,(/itime/))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   

      ENDDO

      ! Close netcdf file
      status_cfsr = nf90_close(ncid_cfsr)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      write(*,*) 'Time modified by one hours for: ',filein
   ENDDO
ENDDO
end if

if (year_end.ne.year_st)then
iyear= year_end

DO imonth= 1,month_end
 
   write(ryear,*) iyear
   ryear = adjustl(ryear)
   write(num_temp,*) imonth
   num_temp=ADJUSTL(num_temp)
   if (len_trim(num_temp).EQ.1) rmonth= '0'//TRIM(num_temp)
   if (len_trim(num_temp).EQ.2) rmonth= TRIM(num_temp)

   write(*,*) 'Treatment of year ',TRIM(ryear),' and month ',TRIM(rmonth),' in progress ...'
   
      filein=TRIM(pathin)//TRIM(varnames(vars))//'_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
      
      status_cfsr = nf90_open(TRIM(filein),nf90_write,ncid_cfsr)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      
      status_cfsr = nf90_inq_dimid(ncid_cfsr, "lon", LonDimID)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_inq_dimid(ncid_cfsr, "lat", LatDimID)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_inq_dimid(ncid_cfsr, "time", TimeDimID)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   
      status_cfsr = nf90_Inquire_Dimension(ncid_cfsr,LonDimID,len = nlon)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_Inquire_Dimension(ncid_cfsr,LatDimID,len = nlat)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_Inquire_Dimension(ncid_cfsr,TimeDimID,len = ntime)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      
   
      jdref_cfsr=jd(iyear,imonth,1) 
   
      DO itime= 1,ntime
         ! _________________________________________________________________
         ! Read in Time
         status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_inq_varid(ncid_cfsr, "forecast_hour",FCTimeVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         
         status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,FCTimeVarId,fctime,start=(/ itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

         time= time+1
         fctime = fctime + 1
         
         status_cfsr = nf90_put_var(ncid_cfsr,TimeVarId,time,(/itime/))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_put_var(ncid_cfsr,FCTimeVarId,fctime,(/itime/))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   

      ENDDO

      ! Close netcdf file
      status_cfsr = nf90_close(ncid_cfsr)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      write(*,*) 'Time modified by one hours for: ',filein

   ENDDO
end if
ENDDO
! ______________________________________________
END PROGRAM modify_times



