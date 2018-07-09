PROGRAM process_6hr_rain
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

INTEGER iyear, itime, imonth, iyear2, imonth2, itime2
INTEGER year_st, year_end, month_st, month_end
INTEGER ncid_cfsr, status_cfsr, ncid_cfsr2
INTEGER LonDimID, LatDimID, TimeDimId,TimeVarId, VarId
INTEGER TimeDimId2,TimeVarId2, VarId2
INTEGER nlon, nlat, ntime, ntime2
INTEGER nvalue
INTEGER jdref,jdref_cfsr

REAL, ALLOCATABLE, DIMENSION(:,:) :: scr, scrp, emask
REAL, ALLOCATABLE, DIMENSION(:,:) :: error, work
DOUBLE PRECISION :: time
REAL :: tx
CHARACTER(len=160) pathin, pathout,filein,filein2,fileout, time_att_units
CHARACTER(len= 50):: ryear, rmonth, num_temp
CHARACTER(len= 50):: ryear2, rmonth2, num_temp2

REAL, PARAMETER :: undef = 2.E+35            ! Undefined land value
REAL, PARAMETER :: critx=0.01
REAL, PARAMETER :: cor=1.6
INTEGER, PARAMETER :: mxs=100

tx=0.9*undef

OPEN(unit=10,file='process.in')
READ(10,'(a)') pathin
READ(10,'(a)') pathout
READ(10,*) year_st, month_st 
READ(10,*) year_end, month_end
CLOSE(10)


jdref= jd(1970,1,1)
time_att_units="hours since 1970-01-01 00:00:00"

iyear= year_st

DO imonth= month_st,12
 
   write(ryear,*) iyear
   ryear = adjustl(ryear)
   write(num_temp,*) imonth
   num_temp=ADJUSTL(num_temp)
   if (len_trim(num_temp).EQ.1) rmonth= '0'//TRIM(num_temp)
   if (len_trim(num_temp).EQ.2) rmonth= TRIM(num_temp)
   
   if (imonth .eq. 12) then
      imonth2 = 1
      iyear2 = iyear +1
   else
      imonth2 = imonth + 1
      iyear2 = iyear
   end if
   write(ryear2,*) iyear2
   ryear2 = adjustl(ryear2)
   write(num_temp2,*) imonth2
   num_temp2=ADJUSTL(num_temp2)
   if (len_trim(num_temp2).EQ.1) rmonth2= '0'//TRIM(num_temp2)
   if (len_trim(num_temp2).EQ.2) rmonth2= TRIM(num_temp2)
   

   write(*,*) 'Treatment of year ',TRIM(ryear),' and month ',TRIM(rmonth),' in progress ...'
   
   filein=TRIM(pathout)//'/backup/'//'rain_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
   fileout =TRIM(pathout)//'rain_'//TRIM(ryear)//TRIM(rmonth)//'.nc'

   CALL system('cp ' // TRIM(filein) // ' ' // TRIM(fileout))

   filein2 =TRIM(pathout)//'rain_'//TRIM(ryear2)//TRIM(rmonth2)//'.nc'
   
   status_cfsr = nf90_open(TRIM(fileout),nf90_write,ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_open(TRIM(filein2),nf90_write,ncid_cfsr2)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)


   status_cfsr = nf90_inq_dimid(ncid_cfsr, "time", TimeDimID)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_Inquire_Dimension(ncid_cfsr,TimeDimID,len = ntime)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
        
   jdref_cfsr=jd(iyear,imonth,1) 
   ALLOCATE(scr(nlon,nlat))
   
   DO itime= 2,ntime
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

! Put time back in i-1 index
      status_cfsr = nf90_put_var(ncid_cfsr,TimeVarId,time,(/itime-1/))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)   
      ! _________________________________________________________________
      ! Read in DLW value
      ! Get DLW var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "PRATE_L1",VarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,VarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
! Put var back in i-1 index
      status_cfsr = nf90_put_var(ncid_cfsr,VarId,scr,start=(/1, 1,itime-1 /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)  

   ENDDO
! Get 1st time from next file
      status_cfsr = nf90_inq_varid(ncid_cfsr2, "time",TimeVarId2)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr2,TimeVarId2,time,start=(/ 1 /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
!Put this into original file
      status_cfsr = nf90_put_var(ncid_cfsr,TimeVarId,time,(/ntime/))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
!Repeat with Var
      status_cfsr = nf90_inq_varid(ncid_cfsr2, "PRATE_L1",VarId2)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr2,VarId2,scr,start=(/1, 1, 1 /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
!Put into file
      status_cfsr = nf90_put_var(ncid_cfsr,VarId,scr,start=(/1, 1,ntime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)  

   DEALLOCATE(scr)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_close(ncid_cfsr2)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   write(*,*) 'Rearanged times in: ','rain_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
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
   
   if (imonth .eq. 12) then
      imonth2 = 1
      iyear2 = iyear +1
   else
      imonth2 = imonth + 1
      iyear2 = iyear
   end if
   write(ryear2,*) iyear2
   ryear2 = adjustl(ryear2)
   write(num_temp2,*) imonth2
   num_temp2=ADJUSTL(num_temp2)
   if (len_trim(num_temp2).EQ.1) rmonth2= '0'//TRIM(num_temp2)
   if (len_trim(num_temp2).EQ.2) rmonth2= TRIM(num_temp2)
   

   write(*,*) 'Treatment of year ',TRIM(ryear),' and month ',TRIM(rmonth),' in progress ...'
   
   filein=TRIM(pathout)//'/backup/'//'rain_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
   fileout =TRIM(pathout)//'rain_'//TRIM(ryear)//TRIM(rmonth)//'.nc'

   CALL system('cp ' // TRIM(filein) // ' ' // TRIM(fileout))

   filein2 =TRIM(pathout)//'rain_'//TRIM(ryear2)//TRIM(rmonth2)//'.nc'
   
   status_cfsr = nf90_open(TRIM(fileout),nf90_write,ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_open(TRIM(filein2),nf90_write,ncid_cfsr2)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)


   status_cfsr = nf90_inq_dimid(ncid_cfsr, "time", TimeDimID)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_Inquire_Dimension(ncid_cfsr,TimeDimID,len = ntime)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_inq_dimid(ncid_cfsr2, "time", TimeDimID2)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_Inquire_Dimension(ncid_cfsr2,TimeDimID2,len = ntime2)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
        
   jdref_cfsr=jd(iyear,imonth,1) 
   ALLOCATE(scr(nlon,nlat))
   
   DO itime= 2,ntime
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

! Put time back in i-1 index
      status_cfsr = nf90_put_var(ncid_cfsr,TimeVarId,time,(/itime-1/))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)   
      ! _________________________________________________________________
      ! Read in DLW value
      ! Get DLW var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "PRATE_L1",VarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,VarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
! Put var back in i-1 index
      status_cfsr = nf90_put_var(ncid_cfsr,VarId,scr,start=(/1, 1,itime-1 /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)  

   ENDDO
! Get 1st time from next file
itime = 1
      status_cfsr = nf90_inq_varid(ncid_cfsr2, "time",TimeVarId2)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr2,TimeVarId2,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
!Put this into original file
      status_cfsr = nf90_put_var(ncid_cfsr,TimeVarId,time,(/ntime/))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
!Repeat with Var
      status_cfsr = nf90_inq_varid(ncid_cfsr2, "PRATE_L1",VarId2)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr2,VarId2,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
!Put into file
      status_cfsr = nf90_put_var(ncid_cfsr,VarId,scr,start=(/1, 1,ntime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)  

   DEALLOCATE(scr)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_close(ncid_cfsr2)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   write(*,*) 'Rearanged times in: ','rain_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
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
   
   if (imonth .eq. 12) then
      imonth2 = 1
      iyear2 = iyear +1
   else
      imonth2 = imonth + 1
      iyear2 = iyear
   end if
   write(ryear2,*) iyear2
   ryear2 = adjustl(ryear2)
   write(num_temp2,*) imonth2
   num_temp2=ADJUSTL(num_temp2)
   if (len_trim(num_temp2).EQ.1) rmonth2= '0'//TRIM(num_temp2)
   if (len_trim(num_temp2).EQ.2) rmonth2= TRIM(num_temp2)
   

   write(*,*) 'Treatment of year ',TRIM(ryear),' and month ',TRIM(rmonth),' in progress ...'
   
   filein=TRIM(pathout)//'/backup/'//'rain_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
   fileout =TRIM(pathout)//'rain_'//TRIM(ryear)//TRIM(rmonth)//'.nc'

   CALL system('cp ' // TRIM(filein) // ' ' // TRIM(fileout))

   filein2 =TRIM(pathout)//'rain_'//TRIM(ryear2)//TRIM(rmonth2)//'.nc'
   
   status_cfsr = nf90_open(TRIM(fileout),nf90_write,ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_open(TRIM(filein2),nf90_write,ncid_cfsr2)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)


   status_cfsr = nf90_inq_dimid(ncid_cfsr, "time", TimeDimID)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_Inquire_Dimension(ncid_cfsr,TimeDimID,len = ntime)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_inq_dimid(ncid_cfsr2, "time", TimeDimID2)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_Inquire_Dimension(ncid_cfsr2,TimeDimID2,len = ntime2)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
        
   jdref_cfsr=jd(iyear,imonth,1) 
   ALLOCATE(scr(nlon,nlat))
    
   DO itime= 2,ntime
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

! Put time back in i-1 index
      status_cfsr = nf90_put_var(ncid_cfsr,TimeVarId,time,(/itime-1/))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)   
      ! _________________________________________________________________
      ! Read in DLW value
      ! Get DLW var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "PRATE_L1",VarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,VarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
! Put var back in i-1 index
      status_cfsr = nf90_put_var(ncid_cfsr,VarId,scr,start=(/1, 1,itime-1 /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)  

   ENDDO
! Get 1st time from next file
itime = 1
      status_cfsr = nf90_inq_varid(ncid_cfsr2, "time",TimeVarId2)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr2,TimeVarId2,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
!Put this into original file
      status_cfsr = nf90_put_var(ncid_cfsr,TimeVarId,time,(/ntime/))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
!Repeat with Var
      status_cfsr = nf90_inq_varid(ncid_cfsr2, "PRATE_L1",VarId2)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr2,VarId2,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
!Put into file
      status_cfsr = nf90_put_var(ncid_cfsr,VarId,scr,start=(/1, 1,ntime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)  

   DEALLOCATE(scr)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_close(ncid_cfsr2)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   write(*,*) 'Rearanged times in: ','rain_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
ENDDO
end if
! ______________________________________________
END PROGRAM process_6hr_rain



