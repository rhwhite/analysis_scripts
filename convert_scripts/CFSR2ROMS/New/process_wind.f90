PROGRAM process_wind
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
INTEGER LonDimID, LatDimID, TimeDimId,TimeVarId, VarId, VarId2
INTEGER nlon, nlat, ntime
INTEGER nvalue
INTEGER jdref,jdref_cfsr

REAL, ALLOCATABLE, DIMENSION(:,:) :: scr, scrp, emask
REAL, ALLOCATABLE, DIMENSION(:,:) :: error, work
DOUBLE PRECISION :: time
REAL :: tx
CHARACTER(len=160) pathin, pathout,filein,fileout, time_att_units
CHARACTER(len= 50):: ryear, rmonth, num_temp
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

   write(*,*) 'Treatment of year ',TRIM(ryear),' and month ',TRIM(rmonth),' in progress ...'
   
   filein=TRIM(pathin)//'wnd10m.gdas.'//TRIM(ryear)//TRIM(rmonth)//'.grb2.nc'
   fileout=TRIM(pathout)//'wind10m_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
   
   CALL system('mv ' // TRIM(filein) // ' ' // TRIM(fileout))

   status_cfsr = nf90_open(TRIM(fileout),nf90_write,ncid_cfsr)
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
   
   
   ALLOCATE(scr(nlon,nlat))
   ALLOCATE(scrp(nlon,nlat))
   ALLOCATE(work(nlon,nlat))
   ALLOCATE(error(nlon,nlat))

   jdref_cfsr=jd(iyear,imonth,1) 
   
   DO itime= 1,ntime
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + time/24. - REAL(jdref)
      time= time*24.
      time= nint(time)
      status_cfsr = nf90_put_att(ncid_cfsr, TimeVarID,"units",TRIM(time_att_units))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_put_var(ncid_cfsr,TimeVarId,time,(/itime/))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)   
      ! _________________________________________________________________
      ! Read in value
      ! Get var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "U_GRD_L103",VarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,VarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Check for masked-out values 
      scrp = 0.
      scrp = scr
      if(maxval(scrp).gt.9999999.99) call write_err("missing values exist")

      status_cfsr = nf90_inq_varid(ncid_cfsr, "V_GRD_L103",VarId2)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,VarId2,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Check for masked-out values 
      scrp = 0.
      scrp = scr
      if(maxval(scrp).gt.9999999.99) call write_err("missing values exist")

   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   write(*,*) 'Check for missing values and time in hours added for: ','wind10m_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
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
      
      filein=TRIM(pathin)//'wnd10m.gdas.'//TRIM(ryear)//TRIM(rmonth)//'.grb2.nc'
      fileout=TRIM(pathout)//'wind10m_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
   
      CALL system('mv ' // TRIM(filein) // ' ' // TRIM(fileout))
      
      status_cfsr = nf90_open(TRIM(fileout),nf90_write,ncid_cfsr)
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
      
      
      ALLOCATE(scr(nlon,nlat))
      ALLOCATE(scrp(nlon,nlat))
      ALLOCATE(work(nlon,nlat))
      ALLOCATE(error(nlon,nlat))
      
      jdref_cfsr=jd(iyear,imonth,1) 
   
      DO itime= 1,ntime
         ! _________________________________________________________________
         ! Read in Time
         status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         time= REAL(jdref_cfsr) + time/24. - REAL(jdref)
         time= time*24.
         time= nint(time)
         status_cfsr = nf90_put_att(ncid_cfsr, TimeVarID,"units",TRIM(time_att_units))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_put_var(ncid_cfsr,TimeVarId,time,(/itime/))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)   
         ! _________________________________________________________________
         ! Read in value
         ! Get var id
         status_cfsr = nf90_inq_varid(ncid_cfsr, "U_GRD_L103",VarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,VarId,scr,start=(/1, 1,itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         ! Check for masked-out values 
         scrp = 0.
         scrp = scr
         if(maxval(scrp).gt.9999999.99) call write_err("missing values exist")
         
         status_cfsr = nf90_inq_varid(ncid_cfsr, "V_GRD_L103",VarId2)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,VarId2,scr,start=(/1, 1,itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         ! Check for masked-out values 
         scrp = 0.
         scrp = scr
         if(maxval(scrp).gt.9999999.99) call write_err("missing values exist")
      ENDDO
      
      DEALLOCATE(scr)
      DEALLOCATE(scrp)
      DEALLOCATE(work)
      DEALLOCATE(error)
      
      ! Close netcdf file
      status_cfsr = nf90_close(ncid_cfsr)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      write(*,*) 'Check for missing values and time in hours added for: ','wind10m_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
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
   
   filein=TRIM(pathin)//'wnd10m.gdas.'//TRIM(ryear)//TRIM(rmonth)//'.grb2.nc'
   fileout=TRIM(pathout)//'wind10m_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
   
   CALL system('mv ' // TRIM(filein) // ' ' // TRIM(fileout))

   status_cfsr = nf90_open(TRIM(fileout),nf90_write,ncid_cfsr)
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
   
   
   ALLOCATE(scr(nlon,nlat))
   ALLOCATE(scrp(nlon,nlat))
   ALLOCATE(work(nlon,nlat))
   ALLOCATE(error(nlon,nlat))
   
   jdref_cfsr=jd(iyear,imonth,1) 
   
   DO itime= 1,ntime
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + time/24. - REAL(jdref)
      time= time*24.
      time= nint(time)
      status_cfsr = nf90_put_att(ncid_cfsr, TimeVarID,"units",TRIM(time_att_units))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_put_var(ncid_cfsr,TimeVarId,time,(/itime/))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)   
      ! _________________________________________________________________
      ! Read in value
      ! Get var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "U_GRD_L103",VarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,VarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Check for masked-out values 
      scrp = 0.
      scrp = scr
      if(maxval(scrp).gt.9999999.99) call write_err("missing values exist")
      
      status_cfsr = nf90_inq_varid(ncid_cfsr, "V_GRD_L103",VarId2)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,VarId2,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Check for masked-out values 
      scrp = 0.
      scrp = scr
      if(maxval(scrp).gt.9999999.99) call write_err("missing values exist")
    ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   write(*,*) 'Check for missing values and time in hours added for: ','wind10m_'//TRIM(ryear)//TRIM(rmonth)//'.nc'

ENDDO
end if
! ______________________________________________
END PROGRAM process_wind



