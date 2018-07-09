PROGRAM smooth_rad_precip
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

INTEGER iyear, itime
INTEGER iday1, iday2, ndays, n6hr
INTEGER ncid_core, status_core
INTEGER ncid_core1
INTEGER LonDimID, LatDimID, MaskVarId, SWVarId, LWVarId, TimeDimID
INTEGER SnowVarId, RainVarId
INTEGER LonVarID, LatVarID, MLonVarID, MLatVarID
INTEGER nlon, nlat
INTEGER nvalue
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr, scrp, emask
REAL, ALLOCATABLE, DIMENSION(:,:) :: error, work
REAL, ALLOCATABLE, DIMENSION(:) :: core_lats1,core_lons1,core_lats2,core_lons2

REAL :: tx
CHARACTER(len=160) pathin, pathout,filein,fileout,maskfile,maskout
CHARACTER(len= 50):: ryear
REAL, PARAMETER :: undef = 2.E+35            ! Undefined land value
REAL, PARAMETER :: critx=0.01
REAL, PARAMETER :: cor=1.6
INTEGER, PARAMETER :: mxs=100

tx=0.9*undef

OPEN(unit=10,file='smooth_daily_fields.in')
READ(10,'(a)') pathin
READ(10,'(a)') pathout
CLOSE(10)

maskout=TRIM(pathout)//'land.sfc.gauss.nc'

! Get dimensions of CORE files
status_core = nf90_open(TRIM(maskout),nf90_write,ncid_core)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_inq_dimid(ncid_core, "lon", LonDimID)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_inq_dimid(ncid_core, "lat", LatDimID)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_Inquire_Dimension(ncid_core,LonDimID,len = nlon)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_Inquire_Dimension(ncid_core,LatDimID,len = nlat)
if(status_core /= nf90_NoErr) call handle_err(status_core)

ALLOCATE(scr(nlon,nlat))
ALLOCATE(scrp(nlon,nlat))
ALLOCATE(emask(nlon,nlat))
ALLOCATE(work(nlon,nlat))
ALLOCATE(error(nlon,nlat))
ALLOCATE(core_lats1(nlat))
ALLOCATE(core_lons1(nlon))
ALLOCATE(core_lats2(nlat))
ALLOCATE(core_lons2(nlon))

! Get land-sea mask
status_core = nf90_inq_varid(ncid_core, "land",MaskVarId)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_get_var(ncid_core,MaskVarId,emask)
if(status_core /= nf90_NoErr) call handle_err(status_core)

! Get lats and lons to check dimensions
status_core = nf90_inq_varid(ncid_core, "lat", MLatVarId)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_get_var(ncid_core,MLatVarId,core_lats1)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_inq_varid(ncid_core, "lon", MLonVarId)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_get_var(ncid_core,MLonVarId,core_lons1)
if(status_core /= nf90_NoErr) call handle_err(status_core)



! Do 6-hourly files

write(*,*) 'Treatment of daily rad file in progress ...'

filein=TRIM(pathin)// 'ncar_rad.1948-2009.23OCT2012.nc'
fileout=TRIM(pathout)// 'rad.1948-2009.nc'

CALL system('cp ' // TRIM(filein) // ' ' // TRIM(fileout))

status_core = nf90_open(TRIM(fileout),nf90_write,ncid_core1)
if(status_core /= nf90_NoErr) call handle_err(status_core)

! Get number of 6-hourly records to process
status_core = nf90_inq_dimid(ncid_core1, "TIME", TimeDimID)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_Inquire_Dimension(ncid_core1,TimeDimID,len = n6hr)

write(*,*) 'Number of daily fields ', n6hr

! Check lons and lats
status_core = nf90_inq_varid(ncid_core1, "LAT", LatVarId)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_get_var(ncid_core1,LatVarId,core_lats2)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_inq_varid(ncid_core1, "LON", LonVarId)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_get_var(ncid_core1,LonVarId,core_lons2)
if(status_core /= nf90_NoErr) call handle_err(status_core)

if (ALL(core_lons1.eq.core_lons2)) then
   write(*,*) 'Longitudes checked'
else
   write(*,*) 'Error on Longitudes, not equal with mask'
   stop "stopped"
endif

if (size(core_lats1) .ne. size(core_lats2)) then
  write(*,*) 'Error on Latitudes, not equal with mask'
  stop "stopped"
endif

if (core_lats1(1) * core_lats2(1) < 0) then
  write(*,*) 'Error on Latitudes, wrong way round!'
  stop "stopped"
endif


status_core = nf90_close(ncid_core)
if(status_core /= nf90_NoErr) call handle_err(status_core)


DO itime= 1,n6hr
! _________________________________________________________________
! Read in Q10 value
! Get Q10 var id
    status_core = nf90_inq_varid(ncid_core1, "LWDN_MOD",LWVarId)
    if(status_core /= nf90_NoErr) call handle_err(status_core)
    status_core = nf90_get_var(ncid_core1,LWVarId,scr,start=(/1, 1,itime /))
    if(status_core /= nf90_NoErr) call handle_err(status_core)
! Fill in masked-out values 
    scrp = 0.
    scrp = scr
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
    scr = scrp
! Smooth zonally
    CALL smooth_polar(scrp, scr, nlon, nlat)
! Write out in netCDF file
    status_core = nf90_put_var(ncid_core1,LWVarId,scr,(/1,1,itime/))
    if(status_core /= nf90_NoErr) call handle_err(status_core)   
! _________________________________________________________________
! Read in T10 value
! Get T10 var id
    status_core = nf90_inq_varid(ncid_core1, "SWDN_MOD",SWVarId)
    if(status_core /= nf90_NoErr) call handle_err(status_core)
    status_core = nf90_get_var(ncid_core1,SWVarId,scr,start=(/1, 1,itime /))
    if(status_core /= nf90_NoErr) call handle_err(status_core)
! Fill in masked-out values 
    scrp = 0.
    scrp = scr
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
    scr = scrp
! Smooth zonally
    CALL smooth_polar(scrp, scr, nlon, nlat)
! Write out in netCDF file
    status_core = nf90_put_var(ncid_core1,SWVarId,scr,(/1,1,itime/))
    if(status_core /= nf90_NoErr) call handle_err(status_core)   
 ENDDO

! Close netcdf files
    status_core = nf90_close(ncid_core1)
    if(status_core /= nf90_NoErr) call handle_err(status_core)

    write(*,*) 'Laplacian interpolation and Polar smoothing done for daily files'


! Do monthly files

write(*,*) 'Treatment of monthly precip file in progress ...'

filein=TRIM(pathin)// 'ncar_precip.1948-2009.23OCT2012.nc'
fileout=TRIM(pathout)// 'precip.1948-2009.nc'

CALL system('cp ' // TRIM(filein) // ' ' // TRIM(fileout))

status_core = nf90_open(TRIM(fileout),nf90_write,ncid_core1)
if(status_core /= nf90_NoErr) call handle_err(status_core)

! Get number of 6-hourly records to process
status_core = nf90_inq_dimid(ncid_core1, "TIME", TimeDimID)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_Inquire_Dimension(ncid_core1,TimeDimID,len = n6hr)

write(*,*) 'Number of monthly fields ', n6hr

! Check lons and lats
status_core = nf90_inq_varid(ncid_core1, "LAT", LatVarId)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_get_var(ncid_core1,LatVarId,core_lats2)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_inq_varid(ncid_core1, "LON", LonVarId)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_get_var(ncid_core1,LonVarId,core_lons2)
if(status_core /= nf90_NoErr) call handle_err(status_core)

if (ALL(core_lons1.eq.core_lons2)) then
   write(*,*) 'Longitudes checked'
else
   write(*,*) 'Error on Longitudes, not equal with mask'
   stop "stopped"
endif

if (size(core_lats1) .ne. size(core_lats2)) then
  write(*,*) 'Error on Latitudes, not equal with mask'
  stop "stopped"
endif

if (core_lats1(1) * core_lats2(1) < 0) then
  write(*,*) 'Error on Latitudes, wrong way round!'
  stop "stopped"
endif

DO itime= 1,n6hr
! _________________________________________________________________
! Read in Rain value
! Get Rain var id
    status_core = nf90_inq_varid(ncid_core1, "RAIN",RainVarId)
    if(status_core /= nf90_NoErr) call handle_err(status_core)
    status_core = nf90_get_var(ncid_core1,RainVarId,scr,start=(/1, 1,itime /))
    if(status_core /= nf90_NoErr) call handle_err(status_core)
! Fill in masked-out values 
    scrp = 0.
    scrp = scr
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
    scr = scrp
! Smooth zonally
    CALL smooth_polar(scrp, scr, nlon, nlat)
! Write out in netCDF file
    status_core = nf90_put_var(ncid_core1,RainVarId,scr,(/1,1,itime/))
    if(status_core /= nf90_NoErr) call handle_err(status_core)   
! _________________________________________________________________
! Read in Snow value
! Get Snow var id
    status_core = nf90_inq_varid(ncid_core1, "SNOW",SnowVarId)
    if(status_core /= nf90_NoErr) call handle_err(status_core)
    status_core = nf90_get_var(ncid_core1,SnowVarId,scr,start=(/1, 1,itime /))
    if(status_core /= nf90_NoErr) call handle_err(status_core)
! Fill in masked-out values 
    scrp = 0.
    scrp = scr
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
    scr = scrp
! Smooth zonally
    CALL smooth_polar(scrp, scr, nlon, nlat)
! Write out in netCDF file
    status_core = nf90_put_var(ncid_core1,SnowVarId,scr,(/1,1,itime/))
    if(status_core /= nf90_NoErr) call handle_err(status_core)   
 ENDDO

! Close netcdf files
    status_core = nf90_close(ncid_core1)
    if(status_core /= nf90_NoErr) call handle_err(status_core)

    write(*,*) 'Laplacian interpolation and Polar smoothing done for monthly files'



! ______________________________________________
  END PROGRAM smooth_rad_precip



