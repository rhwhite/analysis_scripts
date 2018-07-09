PROGRAM smooth_cloud
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

INTEGER iyear, itime, imonth, irec
INTEGER year_st, year_end, month_st, month_end
INTEGER ncid_cfsr, status_cfsr, ncido, statuso
INTEGER LonDimID, LatDimID, TimeDimId,TimeVarId, CLOUDVarId, TimeOutId, CloudOutId, LatOutId, LonOutId, LonVarId, LatVarId
INTEGER nlon, nlat, ntime
INTEGER nvalue
INTEGER jdref,jdref_cfsr

REAL, ALLOCATABLE, DIMENSION(:,:) :: scr, scrp, emask
REAL, ALLOCATABLE, DIMENSION(:) :: lon,lat
REAL, ALLOCATABLE, DIMENSION(:,:) :: error, work
REAL :: time
REAL :: tx
CHARACTER(len=160) pathin, pathout,filein1, filein2,filein3,filein4,filein5,filein6,fileout, time_att_units, end_month, filecdl
CHARACTER(len= 50):: ryear, rmonth, num_temp
REAL, PARAMETER :: undef = 2.E+35            ! Undefined land value
REAL, PARAMETER :: critx=0.01
REAL, PARAMETER :: cor=1.6
INTEGER, PARAMETER :: mxs=100

tx=0.9*undef

OPEN(unit=10,file='smooth_cloud.in')
READ(10,'(a)') pathin
READ(10,'(a)') pathout
READ(10,*) year_st, month_st 
READ(10,*) year_end, month_end
CLOSE(10)


jdref= jd(1970,1,1)
time_att_units="hours since 1970-01-01 00:00:00"
filecdl='cloud.cdl'

iyear= year_st

DO imonth= month_st,12
 
   write(ryear,*) iyear
   ryear = adjustl(ryear)
   write(num_temp,*) imonth
   num_temp=ADJUSTL(num_temp)
   if (len_trim(num_temp).EQ.1) rmonth= '0'//TRIM(num_temp)
   if (len_trim(num_temp).EQ.2) rmonth= TRIM(num_temp)

   if((imonth.eq.1).or.(imonth.eq.3).or.(imonth.eq.5).or.(imonth.eq.7).or.(imonth.eq.8).or.(imonth.eq.10).or.(imonth.eq.12))then
      end_month='31.nc'
   elseif(imonth.eq.2)then
      if(MOD(iyear,4).eq.0)then 
         if(MOD(iyear,400).eq.0 )then
            end_month='29.nc'
         elseif(MOD(iyear,100).eq.0) then
            end_month='28.nc'
         else
            end_month='29.nc'
         endif
      else
         end_month='28.nc'
      endif
   else
      end_month='30.nc'
   endif

   write(*,*) 'Treament of year ',TRIM(ryear),' and month ',TRIM(rmonth),' in progress ...'
   
   filein1=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'01-'//TRIM(ryear)//TRIM(rmonth)//'05.nc'
   filein2=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'06-'//TRIM(ryear)//TRIM(rmonth)//'10.nc'
   filein3=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'11-'//TRIM(ryear)//TRIM(rmonth)//'15.nc'
   filein4=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'16-'//TRIM(ryear)//TRIM(rmonth)//'20.nc'
   filein5=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'21-'//TRIM(ryear)//TRIM(rmonth)//'25.nc'
   filein6=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'26-'//TRIM(ryear)//TRIM(rmonth)//TRIM(end_month)

   fileout=TRIM(pathout)//'cloud_'//TRIM(ryear)//TRIM(rmonth)//'.nc'

   
   CALL system('ncgen -o ' // TRIM(fileout) // ' '// TRIM(filecdl))

   statuso = nf90_open(TRIM(fileout),nf90_write,ncido)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   irec= 0
   statuso = nf90_inq_varid(ncido, "T_CDC_200_0_1_0",CloudOutId)
   if(statuso /= nf90_NoErr) call handle_err(statuso) 
   statuso = nf90_inq_varid(ncido, "time",TimeOutId)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_inq_varid(ncido, "lat",LatOutId)
   if(statuso /= nf90_NoErr) call handle_err(statuso) 
    statuso = nf90_inq_varid(ncido, "lon",LonOutId)
   if(statuso /= nf90_NoErr) call handle_err(statuso) 

   status_cfsr = nf90_open(TRIM(filein1),nf90_nowrite,ncid_cfsr)
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

   if (irec.eq.0)then
   ALLOCATE(lon(nlon))
   ALLOCATE(lat(nlat))
   status_cfsr = nf90_inq_varid(ncid_cfsr, "lon",LonVarId)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_get_var(ncid_cfsr,LonVarId,lon)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   statuso = nf90_put_var(ncido,LonOutId,lon)
   if(statuso /= nf90_NoErr) call handle_err(statuso)  
   status_cfsr = nf90_inq_varid(ncid_cfsr, "lat",LatVarId)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_get_var(ncid_cfsr,LatVarId,lat)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   statuso = nf90_put_var(ncido,LatOutId,lat)
   if(statuso /= nf90_NoErr) call handle_err(statuso)    
   DEALLOCATE(lon)
   DEALLOCATE(lat)  
   end if
   
   
   ALLOCATE(scr(nlon,nlat))
   ALLOCATE(scrp(nlon,nlat))
   ALLOCATE(work(nlon,nlat))
   ALLOCATE(error(nlon,nlat))

   jdref_cfsr=jd(iyear,imonth,1) 
   
   DO itime= 1,ntime
      irec = irec +1
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
      time= time*24.
      statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      statuso = nf90_put_var(ncido,TimeOutId,time,(/irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ! _________________________________________________________________
      ! Read in CLOUD value
      ! Get CLOUD var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Fill in masked-out values 
      scrp = 0.
      scrp = scr
      WHERE(scrp.gt.9999999.99) scrp = undef
      CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
      scr = scrp
      ! Smooth zonally
      CALL smooth_polar(scrp, scr, nlon, nlat)
      ! Write out in netCDF file
      statuso = nf90_put_var(ncido,CloudOutId,scr,(/1,1,irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   status_cfsr = nf90_open(TRIM(filein2),nf90_nowrite,ncid_cfsr)
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

   jdref_cfsr=jd(iyear,imonth,6) 
   
   DO itime= 1,ntime
      irec = irec +1
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
      time= time*24.
      statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      statuso = nf90_put_var(ncido,TimeOutId,time,(/irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ! _________________________________________________________________
      ! Read in CLOUD value
      ! Get CLOUD var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Fill in masked-out values 
      scrp = 0.
      scrp = scr
      WHERE(scrp.gt.9999999.99) scrp = undef
      CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
      scr = scrp
      ! Smooth zonally
      CALL smooth_polar(scrp, scr, nlon, nlat)
      ! Write out in netCDF file
      statuso = nf90_put_var(ncido,CloudOutId,scr,(/1,1,irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   status_cfsr = nf90_open(TRIM(filein3),nf90_nowrite,ncid_cfsr)
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

   jdref_cfsr=jd(iyear,imonth,11) 
   
   DO itime= 1,ntime
      irec = irec +1
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
      time= time*24.
      statuso = nf90_put_att(ncido,TimeOutId,"units",TRIM(time_att_units))
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      statuso = nf90_put_var(ncido,TimeOutId,time,(/irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ! _________________________________________________________________
      ! Read in CLOUD value
      ! Get CLOUD var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Fill in masked-out values 
      scrp = 0.
      scrp = scr
      WHERE(scrp.gt.9999999.99) scrp = undef
      CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
      scr = scrp
      ! Smooth zonally
      CALL smooth_polar(scrp, scr, nlon, nlat)
      ! Write out in netCDF file
      statuso = nf90_put_var(ncido,CloudOutId,scr,(/1,1,irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   status_cfsr = nf90_open(TRIM(filein4),nf90_nowrite,ncid_cfsr)
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

   jdref_cfsr=jd(iyear,imonth,16) 
   
   DO itime= 1,ntime
      irec = irec +1
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
      time= time*24.
      statuso = nf90_put_att(ncido,TimeOutId,"units",TRIM(time_att_units))
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      statuso = nf90_put_var(ncido,TimeOutId,time,(/irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ! _________________________________________________________________
      ! Read in CLOUD value
      ! Get CLOUD var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Fill in masked-out values 
      scrp = 0.
      scrp = scr
      WHERE(scrp.gt.9999999.99) scrp = undef
      CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
      scr = scrp
      ! Smooth zonally
      CALL smooth_polar(scrp, scr, nlon, nlat)
      ! Write out in netCDF file
      statuso = nf90_put_var(ncido,CloudOutId,scr,(/1,1,irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   status_cfsr = nf90_open(TRIM(filein5),nf90_nowrite,ncid_cfsr)
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

   jdref_cfsr=jd(iyear,imonth,21) 
   
   DO itime= 1,ntime
      irec = irec +1
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
      time= time*24.
      statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      statuso = nf90_put_var(ncido, TimeOutId,time,(/irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ! _________________________________________________________________
      ! Read in CLOUD value
      ! Get CLOUD var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Fill in masked-out values 
      scrp = 0.
      scrp = scr
      WHERE(scrp.gt.9999999.99) scrp = undef
      CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
      scr = scrp
      ! Smooth zonally
      CALL smooth_polar(scrp, scr, nlon, nlat)
      ! Write out in netCDF file
      statuso = nf90_put_var(ncido, CloudOutId,scr,(/1,1,irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   status_cfsr = nf90_open(TRIM(filein6),nf90_nowrite,ncid_cfsr)
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

   jdref_cfsr=jd(iyear,imonth,26) 
   
   DO itime= 1,ntime
      irec = irec +1
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
      time= time*24.
      statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      statuso = nf90_put_var(ncido,TimeOutId,time,(/irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ! _________________________________________________________________
      ! Read in CLOUD value
      ! Get CLOUD var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Fill in masked-out values 
      scrp = 0.
      scrp = scr
      WHERE(scrp.gt.9999999.99) scrp = undef
      CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
      scr = scrp
      ! Smooth zonally
      CALL smooth_polar(scrp, scr, nlon, nlat)
      ! Write out in netCDF file
      statuso = nf90_put_var(ncido,CloudOutId,scr,(/1,1,irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   ! Close netcdf file
   statuso = nf90_close(ncido)
   if(statuso /= nf90_NoErr) call handle_err(statuso)

   write(*,*) 'Laplacian interpolation and Polar smoothing done for file: ','cloud_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(year_end.gt.year_st+1)then

DO iyear= year_st+1,year_end-1

   DO imonth= 1,12
 
      write(ryear,*) iyear
      ryear = adjustl(ryear)
      write(num_temp,*) imonth
      num_temp=ADJUSTL(num_temp)
      if (len_trim(num_temp).EQ.1) rmonth= '0'//TRIM(num_temp)
      if (len_trim(num_temp).EQ.2) rmonth= TRIM(num_temp)

      if((imonth.eq.1).or.(imonth.eq.3).or.(imonth.eq.5).or.(imonth.eq.7).or.(imonth.eq.8).or.(imonth.eq.10).or.(imonth.eq.12))then
         end_month='31.nc'
      elseif(imonth.eq.2)then
         if(MOD(iyear,4).eq.0)then 
            if(MOD(iyear,400).eq.0 )then
               end_month='29.nc'
            elseif(MOD(iyear,100).eq.0)then
               end_month='28.nc'
            else
               end_month='29.nc'
            endif
         else
            end_month='28.nc'
         endif
      else
         end_month='30.nc'
      endif
      
      write(*,*) 'Treament of year ',TRIM(ryear),' and month ',TRIM(rmonth),' in progress ...'
   
      filein1=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'01-'//TRIM(ryear)//TRIM(rmonth)//'05.nc'
      filein2=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'06-'//TRIM(ryear)//TRIM(rmonth)//'10.nc'
      filein3=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'11-'//TRIM(ryear)//TRIM(rmonth)//'15.nc'
      filein4=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'16-'//TRIM(ryear)//TRIM(rmonth)//'20.nc'
      filein5=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'21-'//TRIM(ryear)//TRIM(rmonth)//'25.nc'
      filein6=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'26-'//TRIM(ryear)//TRIM(rmonth)//TRIM(end_month)
      
      fileout=TRIM(pathout)//'cloud_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
      
      CALL system('ncgen -o ' //TRIM(fileout)// ' '// TRIM(filecdl))
      
      statuso = nf90_open(TRIM(fileout),nf90_write,ncido)
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      irec= 0
      statuso = nf90_inq_varid(ncido, "T_CDC_200_0_1_0",CloudOutId)
      if(statuso /= nf90_NoErr) call handle_err(statuso) 
      statuso = nf90_inq_varid(ncido, "time",TimeOutId)
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      
      status_cfsr = nf90_open(TRIM(filein1),nf90_nowrite,ncid_cfsr)
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

      if (irec.eq.0)then
         ALLOCATE(lon(nlon))
         ALLOCATE(lat(nlat))
         status_cfsr = nf90_inq_varid(ncid_cfsr, "lon",LonVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,LonVarId,lon)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         statuso = nf90_put_var(ncido,LonOutId,lon)
         if(statuso /= nf90_NoErr) call handle_err(statuso)  
         status_cfsr = nf90_inq_varid(ncid_cfsr, "lat",LatVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,LatVarId,lat)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         statuso = nf90_put_var(ncido,LatOutId,lat)
         if(statuso /= nf90_NoErr) call handle_err(statuso)    
         DEALLOCATE(lon)
         DEALLOCATE(lat)  
      end if
      
      ALLOCATE(scr(nlon,nlat))
      ALLOCATE(scrp(nlon,nlat))
      ALLOCATE(work(nlon,nlat))
      ALLOCATE(error(nlon,nlat))
      
      jdref_cfsr=jd(iyear,imonth,1) 
      
      DO itime= 1,ntime
         irec = irec +1
         ! _________________________________________________________________
         ! Read in Time
         status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
         time= time*24.
         statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
         if(statuso /= nf90_NoErr) call handle_err(statuso)
         statuso = nf90_put_var(ncido, TimeOutId,time,(/irec/))
         if(statuso /= nf90_NoErr) call handle_err(statuso)   
         ! _________________________________________________________________
         ! Read in CLOUD value
         ! Get CLOUD var id
         status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         ! Fill in masked-out values 
         scrp = 0.
         scrp = scr
         WHERE(scrp.gt.9999999.99) scrp = undef
         CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
         scr = scrp
         ! Smooth zonally
         CALL smooth_polar(scrp, scr, nlon, nlat)
         ! Write out in netCDF file
         statuso = nf90_put_var(ncido,CloudOutId,scr,(/1,1,irec/))
         if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ENDDO
      
      DEALLOCATE(scr)
      DEALLOCATE(scrp)
      DEALLOCATE(work)
      DEALLOCATE(error)
      
      ! Close netcdf file
      status_cfsr = nf90_close(ncid_cfsr)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      
      status_cfsr = nf90_open(TRIM(filein2),nf90_nowrite,ncid_cfsr)
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
      
      jdref_cfsr=jd(iyear,imonth,6) 
      
      DO itime= 1,ntime
         irec = irec +1
         ! _________________________________________________________________
         ! Read in Time
         status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
         time= time*24.
         statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
         if(statuso /= nf90_NoErr) call handle_err(statuso)
         statuso = nf90_put_var(ncido,TimeOutId,time,(/irec/))
         if(statuso /= nf90_NoErr) call handle_err(statuso)   
         ! _________________________________________________________________
         ! Read in CLOUD value
         ! Get CLOUD var id
         status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime/))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         ! Fill in masked-out values 
         scrp = 0.
         scrp = scr
         WHERE(scrp.gt.9999999.99) scrp = undef
         CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
         scr = scrp
         ! Smooth zonally
         CALL smooth_polar(scrp, scr, nlon, nlat)
         ! Write out in netCDF file
         statuso = nf90_put_var(ncido, cloudOutId,scr,(/1,1,irec/))
         if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ENDDO
      
      DEALLOCATE(scr)
      DEALLOCATE(scrp)
      DEALLOCATE(work)
      DEALLOCATE(error)
      
      ! Close netcdf file
      status_cfsr = nf90_close(ncid_cfsr)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      
      status_cfsr = nf90_open(TRIM(filein3),nf90_nowrite,ncid_cfsr)
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
      
      jdref_cfsr=jd(iyear,imonth,11) 
      
      DO itime= 1,ntime
         irec = irec +1
         ! _________________________________________________________________
         ! Read in Time
         status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
         time= time*24.
         statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
         if(statuso /= nf90_NoErr) call handle_err(statuso)
         statuso = nf90_put_var(ncido,TimeOutId,time,(/irec/))
         if(statuso /= nf90_NoErr) call handle_err(statuso)   
         ! _________________________________________________________________
         ! Read in CLOUD value
         ! Get CLOUD var id
         status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         ! Fill in masked-out values 
         scrp = 0.
         scrp = scr
         WHERE(scrp.gt.9999999.99) scrp = undef
         CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
         scr = scrp
         ! Smooth zonally
         CALL smooth_polar(scrp, scr, nlon, nlat)
         ! Write out in netCDF file
         statuso = nf90_put_var(ncido, CloudOutId,scr,(/1,1,irec/))
         if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ENDDO
      
      DEALLOCATE(scr)
      DEALLOCATE(scrp)
      DEALLOCATE(work)
      DEALLOCATE(error)
      
      ! Close netcdf file
      status_cfsr = nf90_close(ncid_cfsr)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      
      status_cfsr = nf90_open(TRIM(filein4),nf90_nowrite,ncid_cfsr)
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
      
      jdref_cfsr=jd(iyear,imonth,16) 
      
      DO itime= 1,ntime
         irec = irec +1
         ! _________________________________________________________________
         ! Read in Time
         status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
         time= time*24.
         statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
         if(statuso /= nf90_NoErr) call handle_err(statuso)
         statuso = nf90_put_var(ncido, TimeOutId,time,(/irec/))
         if(statuso /= nf90_NoErr) call handle_err(statuso)   
         ! _________________________________________________________________
         ! Read in CLOUD value
         ! Get CLOUD var id
         status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         ! Fill in masked-out values 
         scrp = 0.
         scrp = scr
         WHERE(scrp.gt.9999999.99) scrp = undef
         CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
         scr = scrp
         ! Smooth zonally
         CALL smooth_polar(scrp, scr, nlon, nlat)
         ! Write out in netCDF file
         statuso = nf90_put_var(ncido, CloudOutId,scr,(/1,1,irec/))
         if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ENDDO
      
      DEALLOCATE(scr)
      DEALLOCATE(scrp)
      DEALLOCATE(work)
      DEALLOCATE(error)
      
      ! Close netcdf file
      status_cfsr = nf90_close(ncid_cfsr)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      
      status_cfsr = nf90_open(TRIM(filein5),nf90_nowrite,ncid_cfsr)
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
      
      jdref_cfsr=jd(iyear,imonth,21) 
      
      DO itime= 1,ntime
         irec = irec +1
         ! _________________________________________________________________
         ! Read in Time
         status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
         time= time*24.
         statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
         if(statuso /= nf90_NoErr) call handle_err(statuso)
         statuso = nf90_put_var(ncido, TimeOutId,time,(/irec/))
         if(statuso /= nf90_NoErr) call handle_err(statuso)   
         ! _________________________________________________________________
         ! Read in CLOUD value
         ! Get CLOUD var id
         status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         ! Fill in masked-out values 
         scrp = 0.
         scrp = scr
         WHERE(scrp.gt.9999999.99) scrp = undef
         CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
         scr = scrp
         ! Smooth zonally
         CALL smooth_polar(scrp, scr, nlon, nlat)
         ! Write out in netCDF file
         statuso = nf90_put_var(ncido, CloudOutId,scr,(/1,1,irec/))
         if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ENDDO
      
      DEALLOCATE(scr)
      DEALLOCATE(scrp)
      DEALLOCATE(work)
      DEALLOCATE(error)
      
      ! Close netcdf file
      status_cfsr = nf90_close(ncid_cfsr)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      
      
      status_cfsr = nf90_open(TRIM(filein6),nf90_nowrite,ncid_cfsr)
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
      
      jdref_cfsr=jd(iyear,imonth,26) 
      
      DO itime= 1,ntime
         irec = irec +1
         ! _________________________________________________________________
         ! Read in Time
         status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
         time= time*24.
         statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
         if(statuso /= nf90_NoErr) call handle_err(statuso)
         statuso = nf90_put_var(ncido, TimeOutId,time,(/irec/))
         if(statuso /= nf90_NoErr) call handle_err(statuso)   
         ! _________________________________________________________________
         ! Read in CLOUD value
         ! Get CLOUD var id
         status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
         if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
         ! Fill in masked-out values 
         scrp = 0.
         scrp = scr
         WHERE(scrp.gt.9999999.99) scrp = undef
         CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
         scr = scrp
         ! Smooth zonally
         CALL smooth_polar(scrp, scr, nlon, nlat)
         ! Write out in netCDF file
         statuso = nf90_put_var(ncido, CloudOutId,scr,(/1,1,irec/))
         if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ENDDO
      
      DEALLOCATE(scr)
      DEALLOCATE(scrp)
      DEALLOCATE(work)
      DEALLOCATE(error)
      
      ! Close netcdf file
      status_cfsr = nf90_close(ncid_cfsr)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      
      ! Close netcdf file
      statuso = nf90_close(ncido)
      if(statuso /= nf90_NoErr) call handle_err(statuso)

      write(*,*) 'Laplacian interpolation and Polar smoothing done for file: ','cloud_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
      
   ENDDO
ENDDO
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if (year_end.ne.year_st)then
iyear= year_end

DO imonth= 1,month_end
 
   write(ryear,*) iyear
   ryear = adjustl(ryear)
   write(num_temp,*) imonth
   num_temp=ADJUSTL(num_temp)
   if (len_trim(num_temp).EQ.1) rmonth= '0'//TRIM(num_temp)
   if (len_trim(num_temp).EQ.2) rmonth= TRIM(num_temp)

   if((imonth.eq.1).or.(imonth.eq.3).or.(imonth.eq.5).or.(imonth.eq.7).or.(imonth.eq.8).or.(imonth.eq.10).or.(imonth.eq.12))then
      end_month='31.nc'
   elseif(imonth.eq.2)then
      if(MOD(iyear,4).eq.0)then 
         if(MOD(iyear,400).eq.0 )then
            end_month='29.nc'
         elseif(MOD(iyear,100).eq.0)then 
            end_month='28.nc'
         else
            end_month='29.nc'
         endif
      else
         end_month='28.nc'
      endif
   else
      end_month='30.nc'
   endif

   write(*,*) 'Treament of year ',TRIM(ryear),' and month ',TRIM(rmonth),' in progress ...'
   
   filein1=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'01-'//TRIM(ryear)//TRIM(rmonth)//'05.nc'
   filein2=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'06-'//TRIM(ryear)//TRIM(rmonth)//'10.nc'
   filein3=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'11-'//TRIM(ryear)//TRIM(rmonth)//'15.nc'
   filein4=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'16-'//TRIM(ryear)//TRIM(rmonth)//'20.nc'
   filein5=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'21-'//TRIM(ryear)//TRIM(rmonth)//'25.nc'
   filein6=TRIM(pathin)//'flxf01.gdas.'//TRIM(ryear)//TRIM(rmonth)//'26-'//TRIM(ryear)//TRIM(rmonth)//TRIM(end_month)

   fileout=TRIM(pathout)//'cloud_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
   
   CALL system('ncgen -o ' //TRIM(fileout)// ' '// TRIM(filecdl))

   statuso = nf90_open(TRIM(fileout),nf90_write,ncido)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   irec= 0
   statuso = nf90_inq_varid(ncido, "T_CDC_200_0_1_0",CloudOutId)
   if(statuso /= nf90_NoErr) call handle_err(statuso) 
   statuso = nf90_inq_varid(ncido, "time",TimeOutId)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
 
   status_cfsr = nf90_open(TRIM(filein1),nf90_nowrite,ncid_cfsr)
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
   
   if (irec.eq.0)then
   ALLOCATE(lon(nlon))
   ALLOCATE(lat(nlat))
   status_cfsr = nf90_inq_varid(ncid_cfsr, "lon",LonVarId)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_get_var(ncid_cfsr,LonVarId,lon)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   statuso = nf90_put_var(ncido,LonOutId,lon)
   if(statuso /= nf90_NoErr) call handle_err(statuso)  
   status_cfsr = nf90_inq_varid(ncid_cfsr, "lat",LatVarId)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   status_cfsr = nf90_get_var(ncid_cfsr,LatVarId,lat)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
   statuso = nf90_put_var(ncido,LatOutId,lat)
   if(statuso /= nf90_NoErr) call handle_err(statuso)    
   DEALLOCATE(lon)
   DEALLOCATE(lat)  
   end if
   
   ALLOCATE(scr(nlon,nlat))
   ALLOCATE(scrp(nlon,nlat))
   ALLOCATE(work(nlon,nlat))
   ALLOCATE(error(nlon,nlat))

   jdref_cfsr=jd(iyear,imonth,1) 
   
   DO itime= 1,ntime
      irec = irec +1
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
      time= time*24.
      statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      statuso = nf90_put_var(ncido, TimeOutId,time,(/irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ! _________________________________________________________________
      ! Read in CLOUD value
      ! Get CLOUD var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime/))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Fill in masked-out values 
      scrp = 0.
      scrp = scr
      WHERE(scrp.gt.9999999.99) scrp = undef
      CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
      scr = scrp
      ! Smooth zonally
      CALL smooth_polar(scrp, scr, nlon, nlat)
      ! Write out in netCDF file
      statuso = nf90_put_var(ncido, CloudOutId,scr,(/1,1,irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   status_cfsr = nf90_open(TRIM(filein2),nf90_nowrite,ncid_cfsr)
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

   jdref_cfsr=jd(iyear,imonth,6) 
   
   DO itime= 1,ntime
      irec = irec +1
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
      time= time*24.
      statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      statuso = nf90_put_var(ncido,TimeOutId,time,(/irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ! _________________________________________________________________
      ! Read in CLOUD value
      ! Get CLOUD var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Fill in masked-out values 
      scrp = 0.
      scrp = scr
      WHERE(scrp.gt.9999999.99) scrp = undef
      CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
      scr = scrp
      ! Smooth zonally
      CALL smooth_polar(scrp, scr, nlon, nlat)
      ! Write out in netCDF file
      statuso = nf90_put_var(ncido, CloudOutId,scr,(/1,1,irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   status_cfsr = nf90_open(TRIM(filein3),nf90_nowrite,ncid_cfsr)
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

   jdref_cfsr=jd(iyear,imonth,11) 
   
   DO itime= 1,ntime
      irec = irec +1
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
      time= time*24.
      statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      statuso = nf90_put_var(ncido, TimeOutId,time,(/irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ! _________________________________________________________________
      ! Read in CLOUD value
      ! Get CLOUD var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Fill in masked-out values 
      scrp = 0.
      scrp = scr
      WHERE(scrp.gt.9999999.99) scrp = undef
      CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
      scr = scrp
      ! Smooth zonally
      CALL smooth_polar(scrp, scr, nlon, nlat)
      ! Write out in netCDF file
      statuso = nf90_put_var(ncido, CloudOutId,scr,(/1,1,irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   status_cfsr = nf90_open(TRIM(filein4),nf90_nowrite,ncid_cfsr)
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

   jdref_cfsr=jd(iyear,imonth,16) 
   
   DO itime= 1,ntime
      irec = irec +1
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
      time= time*24.
      statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      statuso = nf90_put_var(ncido, TimeOutId,time,(/irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ! _________________________________________________________________
      ! Read in CLOUD value
      ! Get CLOUD var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Fill in masked-out values 
      scrp = 0.
      scrp = scr
      WHERE(scrp.gt.9999999.99) scrp = undef
      CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
      scr = scrp
      ! Smooth zonally
      CALL smooth_polar(scrp, scr, nlon, nlat)
      ! Write out in netCDF file
      statuso = nf90_put_var(ncido, CloudOutId,scr,(/1,1,irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   status_cfsr = nf90_open(TRIM(filein5),nf90_nowrite,ncid_cfsr)
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

   jdref_cfsr=jd(iyear,imonth,21) 
   
   DO itime= 1,ntime
      irec = irec +1
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
      time= time*24.
      statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      statuso = nf90_put_var(ncido, TimeOutId,time,(/irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ! _________________________________________________________________
      ! Read in CLOUD value
      ! Get CLOUD var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Fill in masked-out values 
      scrp = 0.
      scrp = scr
      WHERE(scrp.gt.9999999.99) scrp = undef
      CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
      scr = scrp
      ! Smooth zonally
      CALL smooth_polar(scrp, scr, nlon, nlat)
      ! Write out in netCDF file
      statuso = nf90_put_var(ncido, CloudOutId,scr,(/1,1,irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)


   status_cfsr = nf90_open(TRIM(filein6),nf90_nowrite,ncid_cfsr)
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

   jdref_cfsr=jd(iyear,imonth,26) 
   
   DO itime= 1,ntime
      irec = irec +1
      ! _________________________________________________________________
      ! Read in Time
      status_cfsr = nf90_inq_varid(ncid_cfsr, "time",TimeVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,TimeVarId,time,start=(/ itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      time= REAL(jdref_cfsr) + REAL(time)/24. - REAL(jdref)
      time= time*24.
      statuso = nf90_put_att(ncido, TimeOutId,"units",TRIM(time_att_units))
      if(statuso /= nf90_NoErr) call handle_err(statuso)
      statuso = nf90_put_var(ncido, TimeOutId,time,(/irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
      ! _________________________________________________________________
      ! Read in CLOUD value
      ! Get CLOUD var id
      status_cfsr = nf90_inq_varid(ncid_cfsr, "T_CDC_200_0_1_0",CLOUDVarId)
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      status_cfsr = nf90_get_var(ncid_cfsr,CLOUDVarId,scr,start=(/1, 1,itime /))
      if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)
      ! Fill in masked-out values 
      scrp = 0.
      scrp = scr
      WHERE(scrp.gt.9999999.99) scrp = undef
      CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
      scr = scrp
      ! Smooth zonally
      CALL smooth_polar(scrp, scr, nlon, nlat)
      ! Write out in netCDF file
      statuso = nf90_put_var(ncido,CloudOutId,scr,(/1,1,irec/))
      if(statuso /= nf90_NoErr) call handle_err(statuso)   
   ENDDO

   DEALLOCATE(scr)
   DEALLOCATE(scrp)
   DEALLOCATE(work)
   DEALLOCATE(error)

   ! Close netcdf file
   status_cfsr = nf90_close(ncid_cfsr)
   if(status_cfsr /= nf90_NoErr) call handle_err(status_cfsr)

   ! Close netcdf file
   statuso = nf90_close(ncido)
   if(statuso /= nf90_NoErr) call handle_err(statuso)

   write(*,*) 'Laplacian interpolation and Polar smoothing done for file: ','cloud_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
ENDDO
end if
! ______________________________________________
END PROGRAM smooth_cloud



