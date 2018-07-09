! Program to generate ROMS forcing files from ERAI netCDF files
! produced using cvt_gribfiles.bsh

PROGRAM cfsr2roms
use netcdf
implicit none

INTEGER nvars = 8

INTEGER i, j, k, l, imon, itime, ntimes, ia, ja, IFL
INTEGER irec
INTEGER iyear, iday, imonth,ihour, nb_records, indice
INTEGER greg(6)
INTEGER LonDimID, LatDimID, TimeDimId
INTEGER nlonw,nlatw,nlonmsl,nlatmsl,nlont2m,nlatt2m,nlonq2m,nlatq2m,nlonrain,nlatrain
INTEGER nloncloud,nlatcloud,nlondlw,nlatdlw,nlondsw,nlatdsw
INTEGER XiDimID, EtaDimID, Lp, Mp, Lp_dim, Mp_dim
INTEGER Lp_wind,Mp_wind,Lp_msl,Mp_msl,Lp_t2m,Mp_t2m,Lp_q2m,Mp_q2m,Lp_rain,Mp_rain
INTEGER Lp_cloud,Mp_cloud,Lp_dlw,Mp_dlw,Lp_dsw,Mp_dsw     

INTEGER LonVarIdwind, LatVarIdwind
INTEGER LonVarIdmsl, LatVarIdmsl
INTEGER LonVarIdq2m, LatVarIdq2m
INTEGER LonVarIdt2m, LatVarIdt2m
INTEGER LonVarIdcloud, LatVarIdcloud
INTEGER LonVarIdrain, LatVarIdrain
INTEGER LonVarIddsw, LatVarIddsw
INTEGER LonVarIddlw, LatVarIddlw


INTEGER LonRhoVarId, LatRhoVarId, AngleVarId
INTEGER U10VarId, V10VarId
INTEGER TimeVarId_wind,TimeVarId_msl,TimeVarId_t2m,TimeVarId_q2m,TimeVarId_rain,TimeVarId_cloud,TimeVarId_dlw,TimeVarId_dsw
INTEGER MSLVarId, T2MVarId, Q2MVarId, CLOUDVarId, RAINVarId, DLWVarId, DSWVarId
INTEGER year_st, month_st, day_st
INTEGER year_end, month_end, day_end
INTEGER year_ref, month_ref, day_ref

INTEGER XiRhoDimIdwind, EtaRhoDimIdwind, TimeOutDimIDwind
INTEGER XiRhoDimIdmsl, EtaRhoDimIdmsl, TimeOutDimIDmsl
INTEGER XiRhoDimIdt2m, EtaRhoDimIdt2m, TimeOutDimIDt2m
INTEGER XiRhoDimIdq2m, EtaRhoDimIdq2m, TimeOutDimIDq2m
INTEGER XiRhoDimIdcloud, EtaRhoDimIdcloud, TimeOutDimIDcloud
INTEGER XiRhoDimIdrain, EtaRhoDimIdrain, TimeOutDimIDrain
INTEGER XiRhoDimIddsw, EtaRhoDimIddsw, TimeOutDimIDdsw
INTEGER XiRhoDimIddlw, EtaRhoDimIddlw, TimeOutDimIDdlw


INTEGER TimeOutVarId, UwindVarId,VwindVarId,PairVarId
INTEGER TairVarId, QairVarId,cloudairVarId,rainairVarId, swradVarId, lwradVarId

INTEGER statuso, ncido
INTEGER ncidowind,ncidomsl,ncidot2m,ncidoq2m,ncidorain,ncidocloud,ncidodlw,ncidodsw
INTEGER status_cfsr, ncid_cfsr
INTEGER status_roms, ncid_roms
INTEGER status_cloud, ncid_cloud
INTEGER status_rain, ncid_rain
INTEGER status_wind, ncid_wind
INTEGER status_t2m, ncid_t2m
INTEGER status_q2m, ncid_q2m
INTEGER status_dlw, ncid_dlw
INTEGER status_dsw, ncid_dsw
INTEGER status_msl, ncid_msl


INTEGER i1_wind, i2_wind, j1_wind, j2_wind
INTEGER i1_msl, i2_msl, j1_msl, j2_msl
INTEGER i1_msltemp, i2_msltemp, j1_msltemp, j2_msltemp
INTEGER i1_t2m, i2_t2m, j1_t2m, j2_t2m
INTEGER i1_q2m, i2_q2m, j1_q2m, j2_q2m
INTEGER i1_rain, i2_rain, j1_rain, j2_rain
INTEGER i1_cloud, i2_cloud, j1_cloud, j2_cloud
INTEGER i1_dlw, i2_dlw, j1_dlw, j2_dlw
INTEGER i1_dsw, i2_dsw, j1_dsw, j2_dsw



INTEGER jdref, jdref_cfsr
INTEGER yes_no

REAL  rx, rxm, ry, rym
REAL  riminw, rimaxw, rjminw, rjmaxw
REAL  riminmsl, rimaxmsl, rjminmsl, rjmaxmsl
REAL  rimint2m, rimaxt2m, rjmint2m, rjmaxt2m
REAL  riminq2m, rimaxq2m, rjminq2m, rjmaxq2m
REAL  riminrain, rimaxrain, rjminrain, rjmaxrain
REAL  rimincloud, rimaxcloud, rjmincloud, rjmaxcloud
REAL  rimindlw, rimaxdlw, rjmindlw, rjmaxdlw
REAL  rimindsw, rimaxdsw, rjmindsw, rjmaxdsw

REAL, DIMENSION(:,:), allocatable :: ipos_wind
REAL, DIMENSION(:,:), allocatable :: jpos_wind
REAL, DIMENSION(:,:), allocatable :: ipos_msl
REAL, DIMENSION(:,:), allocatable :: jpos_msl
REAL, DIMENSION(:)  , allocatable :: ipos_msltemp
REAL, DIMENSION(:)  , allocatable :: jpos_msltemp
REAL, DIMENSION(:,:), allocatable :: ipos_t2m
REAL, DIMENSION(:,:), allocatable :: jpos_t2m
REAL, DIMENSION(:,:), allocatable :: ipos_q2m
REAL, DIMENSION(:,:), allocatable :: jpos_q2m
REAL, DIMENSION(:,:), allocatable :: ipos_rain
REAL, DIMENSION(:,:), allocatable :: jpos_rain
REAL, DIMENSION(:,:), allocatable :: ipos_dlw
REAL, DIMENSION(:,:), allocatable :: jpos_dlw
REAL, DIMENSION(:,:), allocatable :: ipos_dsw
REAL, DIMENSION(:,:), allocatable :: jpos_dsw
REAL, DIMENSION(:,:), allocatable :: ipos_cloud
REAL, DIMENSION(:,:), allocatable :: jpos_cloud

REAL, ALLOCATABLE, DIMENSION(:) :: lonw,latw
REAL, ALLOCATABLE, DIMENSION(:) :: lonwtemp,latwtemp
REAL, ALLOCATABLE, DIMENSION(:) :: IX_wind,IY_wind
REAL, ALLOCATABLE, DIMENSION(:,:) :: work_wind
REAL, ALLOCATABLE, DIMENSION(:) :: lonmsl,latmsl
REAL, ALLOCATABLE, DIMENSION(:) :: lonmsltemp,latmsltemp
REAL, ALLOCATABLE, DIMENSION(:) :: IX_msl,IY_msl
REAL, ALLOCATABLE, DIMENSION(:,:) :: work_msl
REAL, ALLOCATABLE, DIMENSION(:) :: lont2m,latt2m
REAL, ALLOCATABLE, DIMENSION(:) :: IX_t2m,IY_t2m
REAL, ALLOCATABLE, DIMENSION(:,:) :: work_t2m
REAL, ALLOCATABLE, DIMENSION(:) :: lonq2m,latq2m
REAL, ALLOCATABLE, DIMENSION(:) :: IX_q2m,IY_q2m
REAL, ALLOCATABLE, DIMENSION(:,:) :: work_q2m
REAL, ALLOCATABLE, DIMENSION(:) :: lonrain,latrain
REAL, ALLOCATABLE, DIMENSION(:) :: IX_rain,IY_rain
REAL, ALLOCATABLE, DIMENSION(:,:) :: work_rain
REAL, ALLOCATABLE, DIMENSION(:) :: loncloud,latcloud
REAL, ALLOCATABLE, DIMENSION(:) :: IX_cloud,IY_cloud
REAL, ALLOCATABLE, DIMENSION(:,:) :: work_cloud
REAL, ALLOCATABLE, DIMENSION(:) :: londlw,latdlw
REAL, ALLOCATABLE, DIMENSION(:) :: IX_dlw,IY_dlw
REAL, ALLOCATABLE, DIMENSION(:,:) :: work_dlw
REAL, ALLOCATABLE, DIMENSION(:) :: londsw,latdsw
REAL, ALLOCATABLE, DIMENSION(:) :: IX_dsw,IY_dsw
REAL, ALLOCATABLE, DIMENSION(:,:) :: work_dsw

REAL, ALLOCATABLE, DIMENSION(:,:) :: scr1_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr2_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: lon_rho
REAL, ALLOCATABLE, DIMENSION(:,:) :: lat_rho
REAL, ALLOCATABLE, DIMENSION(:,:) :: angle

REAL, ALLOCATABLE, DIMENSION(:,:) :: U10
REAL, ALLOCATABLE, DIMENSION(:,:) :: V10
REAL, ALLOCATABLE, DIMENSION(:,:) :: MSL
REAL, ALLOCATABLE, DIMENSION(:,:) :: Q2M
REAL, ALLOCATABLE, DIMENSION(:,:) :: T2M
REAL, ALLOCATABLE, DIMENSION(:,:) :: CLOUD
REAL, ALLOCATABLE, DIMENSION(:,:) :: RAIN
REAL, ALLOCATABLE, DIMENSION(:,:) :: DLW
REAL, ALLOCATABLE, DIMENSION(:,:) :: DSW

REAL, ALLOCATABLE, DIMENSION(:,:) :: Uwind
REAL, ALLOCATABLE, DIMENSION(:,:) :: Vwind
REAL, ALLOCATABLE, DIMENSION(:,:) :: Pair
REAL, ALLOCATABLE, DIMENSION(:,:) :: Pairtemp
REAL, ALLOCATABLE, DIMENSION(:,:) :: Tair
REAL, ALLOCATABLE, DIMENSION(:,:) :: Qair
REAL, ALLOCATABLE, DIMENSION(:,:) :: cloudair
REAL, ALLOCATABLE, DIMENSION(:,:) :: rainair
REAL, ALLOCATABLE, DIMENSION(:,:) :: swrad
REAL, ALLOCATABLE, DIMENSION(:,:) :: lwrad

DOUBLE PRECISION tday_wind, tday_msl, tday_t2m, tday_q2m, tday_rain, tday_cloud, tday_dlw, tday_dsw, tday
DOUBLE PRECISION date
REAL datetime

REAL, PARAMETER    :: undef = 2.E+35            ! Undefined land value
REAL pi, DTOR

CHARACTER(len=4) year_char
CHARACTER(len=2) month_char, day_char
CHARACTER(len=50)  ryear, num_temp, rmonth
CHARACTER(len=160) path_cfsr, maskfile, gridfile, fileout, time_units_att
CHARACTER(len=80) wind10m_file, msl_file, t2m_file, q2m_file, rain_file, cloud_file, dlw_file, dsw_file

interface
     FUNCTION jd(yy,mm,dd)
     implicit none
     INTEGER jd,yy,mm,dd
     end function jd
end interface

pi = ATAN(1.)*4.
DTOR = pi/180.

OPEN(unit=10,file='cfsr2roms.in')
READ(10,*)
READ(10,'(a)') path_cfsr
READ(10,*)
READ(10,'(a)') gridfile
READ(10,*)
READ(10,'(a)') fileout
READ(10,*)!interpolation on ROMS grid: yes=1/no=0
READ(10,*) yes_no
READ(10,*)!start of the simulation
READ(10,*) year_st, month_st, day_st
READ(10,*)!end of the simulation
READ(10,*) year_end, month_end, day_end
READ(10,*)! reference time of the simulation
READ(10,*) year_ref, month_ref, day_ref
CLOSE(10)

! Get reference Julian day for output
jdref = jd(year_ref,month_ref,day_ref)
! Get reference Julian day (1970-01-01 00:00:00Z) for input
jdref_cfsr = jd(1970,1,1)


write(num_temp,*) year_ref
year_char=ADJUSTL(num_temp)
write(num_temp,*) month_ref
num_temp=ADJUSTL(num_temp)
if (len_trim(num_temp).EQ.1) month_char= '0'//TRIM(num_temp)
if (len_trim(num_temp).EQ.2) month_char= TRIM(num_temp)
write(num_temp,*) day_ref
num_temp=ADJUSTL(num_temp)
if (len_trim(num_temp).EQ.1) day_char= '0'//TRIM(num_temp)
if (len_trim(num_temp).EQ.2) day_char= TRIM(num_temp)
time_units_att='seconds since '//TRIM(year_char)//'-'//TRIM(month_char)//'-'//TRIM(day_char)//' 00:00:00'

write(*,*) 'Grid file = ',TRIM(gridfile)
! Get info on ROMS grid
status_roms = nf90_open(TRIM(gridfile),nf90_nowrite,ncid_roms)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_inq_dimid(ncid_roms, "xi_rho", XiDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_inq_dimid(ncid_roms, "eta_rho", EtaDimID)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_Inquire_Dimension(ncid_roms,XiDimID,len = Lp)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_Inquire_Dimension(ncid_roms,EtaDimID,len = Mp)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
ALLOCATE(lon_rho(Lp,Mp))
ALLOCATE(lat_rho(Lp,Mp))
ALLOCATE(angle(Lp,Mp))
ALLOCATE(ipos_wind(Lp,Mp))
ALLOCATE(jpos_wind(Lp,Mp))
ALLOCATE(ipos_msl(Lp,Mp))
ALLOCATE(jpos_msl(Lp,Mp))
ALLOCATE(ipos_t2m(Lp,Mp))
ALLOCATE(jpos_t2m(Lp,Mp))
ALLOCATE(ipos_q2m(Lp,Mp))
ALLOCATE(jpos_q2m(Lp,Mp))
ALLOCATE(ipos_rain(Lp,Mp))
ALLOCATE(jpos_rain(Lp,Mp))
ALLOCATE(ipos_cloud(Lp,Mp))
ALLOCATE(jpos_cloud(Lp,Mp))
ALLOCATE(ipos_dlw(Lp,Mp))
ALLOCATE(jpos_dlw(Lp,Mp))
ALLOCATE(ipos_dsw(Lp,Mp))
ALLOCATE(jpos_dsw(Lp,Mp))
ALLOCATE(scr1_out(Lp,Mp))
ALLOCATE(scr2_out(Lp,Mp))
write(*,*) 'ROMS Lp, Mp = ',Lp,Mp
! Get lons and lats
status_roms = nf90_inq_varid(ncid_roms, "lon_rho",LonRhoVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_inq_varid(ncid_roms, "lat_rho",LatRhoVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_inq_varid(ncid_roms, "angle",AngleVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_get_var(ncid_roms,LonRhoVarId,lon_rho)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_get_var(ncid_roms,LatRhoVarId,lat_rho)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_get_var(ncid_roms,AngleVarId,angle)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_close(ncid_roms)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

nb_records = (jd(year_end,month_end,day_end)-jd(year_st,month_st,day_st)+1)*4

DO irec=1,nb_records

date= REAL(jd(year_st,month_st,day_st))+REAL(irec-1)*(6.d0/24.d0)+(1.d0/24.d0) + 0.00001
greg=0
CALL GR(greg,date)
iyear= greg(1)
imonth= greg(2)
iday= greg(3)
ihour= greg(4)

write(*,'(a24,i2,a1,i2,a1,i4,a4,i2,a1)')'Writting data for date: ',greg(3),'/',greg(2),'/',greg(1), ' at ',greg(4),'h'

if(ihour.eq.1)then
indice = 1
elseif(ihour.eq.7)then
indice = 2
elseif(ihour.eq.13)then
indice = 3
elseif(ihour.eq.19)then
indice = 4
end if

write(ryear,*) iyear
ryear = adjustl(ryear)
write(num_temp,*) imonth
num_temp=ADJUSTL(num_temp)
if (len_trim(num_temp).EQ.1) rmonth= '0'//TRIM(num_temp)
if (len_trim(num_temp).EQ.2) rmonth= TRIM(num_temp)
itime= (jd(iyear,imonth,iday)-jd(iyear,imonth,1))*4 + indice


! Get dimensions of CFSR files
wind10m_file = TRIM(path_cfsr)//'wind10m_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
status_wind = nf90_open(TRIM(wind10m_file),nf90_nowrite,ncid_wind)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
msl_file = TRIM(path_cfsr)//'msl_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
status_msl = nf90_open(TRIM(msl_file),nf90_nowrite,ncid_msl)
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
t2m_file = TRIM(path_cfsr)//'t2m_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
status_t2m = nf90_open(TRIM(t2m_file),nf90_nowrite,ncid_t2m)
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
q2m_file = TRIM(path_cfsr)//'q2m_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
status_q2m = nf90_open(TRIM(q2m_file),nf90_nowrite,ncid_q2m)
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
rain_file = TRIM(path_cfsr)//'rain_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
status_rain = nf90_open(TRIM(rain_file),nf90_nowrite,ncid_rain)
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
dlw_file = TRIM(path_cfsr)//'dlw_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
status_dlw = nf90_open(TRIM(dlw_file),nf90_nowrite,ncid_dlw)
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
dsw_file = TRIM(path_cfsr)//'dsw_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
status_dsw = nf90_open(TRIM(dsw_file),nf90_nowrite,ncid_dsw)
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)
cloud_file = TRIM(path_cfsr)//'cloud_'//TRIM(ryear)//TRIM(rmonth)//'.nc'
status_cloud = nf90_open(TRIM(cloud_file),nf90_nowrite,ncid_cloud)
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)

if(irec.eq.1)then
status_wind = nf90_inq_dimid(ncid_wind, "lon", LonDimID)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_wind = nf90_inq_dimid(ncid_wind, "lat", LatDimID)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_wind = nf90_Inquire_Dimension(ncid_wind,LonDimID,len = nlonw)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_wind = nf90_Inquire_Dimension(ncid_wind,LatDimID,len = nlatw)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_msl = nf90_inq_dimid(ncid_msl, "lon", LonDimID)
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
status_msl = nf90_inq_dimid(ncid_msl, "lat", LatDimID)
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
status_msl = nf90_Inquire_Dimension(ncid_msl,LonDimID,len = nlonmsl)
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
status_msl = nf90_Inquire_Dimension(ncid_msl,LatDimID,len = nlatmsl)
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
status_t2m = nf90_inq_dimid(ncid_t2m, "lon", LonDimID)
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
status_t2m = nf90_inq_dimid(ncid_t2m, "lat", LatDimID)
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
status_t2m = nf90_Inquire_Dimension(ncid_t2m,LonDimID,len = nlont2m)
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
status_t2m = nf90_Inquire_Dimension(ncid_t2m,LatDimID,len = nlatt2m)
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
status_q2m = nf90_inq_dimid(ncid_q2m, "lon", LonDimID)
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
status_q2m = nf90_inq_dimid(ncid_q2m, "lat", LatDimID)
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
status_q2m = nf90_Inquire_Dimension(ncid_q2m,LonDimID,len = nlonq2m)
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
status_q2m = nf90_Inquire_Dimension(ncid_q2m,LatDimID,len = nlatq2m)
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
status_rain = nf90_inq_dimid(ncid_rain, "lon", LonDimID)
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
status_rain = nf90_inq_dimid(ncid_rain, "lat", LatDimID)
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
status_rain = nf90_Inquire_Dimension(ncid_rain,LonDimID,len = nlonrain)
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
status_rain = nf90_Inquire_Dimension(ncid_rain,LatDimID,len = nlatrain)
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
status_dlw = nf90_inq_dimid(ncid_dlw, "lon", LonDimID)
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
status_dlw = nf90_inq_dimid(ncid_dlw, "lat", LatDimID)
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
status_dlw = nf90_Inquire_Dimension(ncid_dlw,LonDimID,len = nlondlw)
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
status_dlw = nf90_Inquire_Dimension(ncid_dlw,LatDimID,len = nlatdlw)
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
status_dsw = nf90_inq_dimid(ncid_dsw, "lon", LonDimID)
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)
status_dsw = nf90_inq_dimid(ncid_dsw, "lat", LatDimID)
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)
status_dsw = nf90_Inquire_Dimension(ncid_dsw,LonDimID,len = nlondsw)
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)
status_dsw = nf90_Inquire_Dimension(ncid_dsw,LatDimID,len = nlatdsw)
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)
status_cloud = nf90_inq_dimid(ncid_cloud, "lon", LonDimID)
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)
status_cloud = nf90_inq_dimid(ncid_cloud, "lat", LatDimID)
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)
status_cloud = nf90_Inquire_Dimension(ncid_cloud,LonDimID,len = nloncloud)
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)
status_cloud = nf90_Inquire_Dimension(ncid_cloud,LatDimID,len = nlatcloud)
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)

ALLOCATE(lonw(nlonw))
ALLOCATE(latw(nlatw))
ALLOCATE(IX_wind(nlonw))
ALLOCATE(IY_wind(nlatw))
ALLOCATE(work_wind(nlonw,nlatw))
ALLOCATE(U10(nlonw,nlatw))
ALLOCATE(V10(nlonw,nlatw))
ALLOCATE(lonmsl(nlonmsl))
ALLOCATE(latmsl(nlatmsl))
ALLOCATE(IX_msl(nlonmsl))
ALLOCATE(IY_msl(nlatmsl))
ALLOCATE(work_msl(nlonmsl,nlatmsl))
ALLOCATE(MSL(nlonmsl,nlatmsl))
ALLOCATE(lont2m(nlont2m))
ALLOCATE(latt2m(nlatt2m))
ALLOCATE(IX_t2m(nlont2m))
ALLOCATE(IY_t2m(nlatt2m))
ALLOCATE(work_t2m(nlont2m,nlatt2m))
ALLOCATE(T2M(nlont2m,nlatt2m))
ALLOCATE(lonq2m(nlonq2m))
ALLOCATE(latq2m(nlatq2m))
ALLOCATE(IX_q2m(nlonq2m))
ALLOCATE(IY_q2m(nlatq2m))
ALLOCATE(work_q2m(nlonq2m,nlatq2m))
ALLOCATE(Q2M(nlonq2m,nlatq2m))
ALLOCATE(lonrain(nlonrain))
ALLOCATE(latrain(nlatrain))
ALLOCATE(IX_rain(nlonrain))
ALLOCATE(IY_rain(nlatrain))
ALLOCATE(work_rain(nlonrain,nlatrain))
ALLOCATE(RAIN(nlonrain,nlatrain))
ALLOCATE(londlw(nlondlw))
ALLOCATE(latdlw(nlatdlw))
ALLOCATE(IX_dlw(nlondlw))
ALLOCATE(IY_dlw(nlatdlw))
ALLOCATE(work_dlw(nlondlw,nlatdlw))
ALLOCATE(DLW(nlondlw,nlatdlw))
ALLOCATE(londsw(nlondsw))
ALLOCATE(latdsw(nlatdsw))
ALLOCATE(IX_dsw(nlondsw))
ALLOCATE(IY_dsw(nlatdsw))
ALLOCATE(work_dsw(nlondsw,nlatdsw))
ALLOCATE(DSW(nlondsw,nlatdsw))
ALLOCATE(loncloud(nloncloud))
ALLOCATE(latcloud(nlatcloud))
ALLOCATE(IX_cloud(nloncloud))
ALLOCATE(IY_cloud(nlatcloud))
ALLOCATE(work_cloud(nloncloud,nlatcloud))
ALLOCATE(CLOUD(nloncloud,nlatcloud))


! Get lons and lats
status_wind = nf90_inq_varid(ncid_wind, "lon",LonVarId)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_wind = nf90_inq_varid(ncid_wind, "lat",LatVarId)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_wind = nf90_get_var(ncid_wind,LonVarId,lonw)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_wind = nf90_get_var(ncid_wind,LatVarId,latw)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_msl = nf90_inq_varid(ncid_msl, "lon",LonVarId)
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
status_msl = nf90_inq_varid(ncid_msl, "lat",LatVarId)
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
status_msl = nf90_get_var(ncid_msl,LonVarId,lonmsl)
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
status_msl = nf90_get_var(ncid_msl,LatVarId,latmsl)
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
status_t2m = nf90_inq_varid(ncid_t2m, "lon",LonVarId)
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
status_t2m = nf90_inq_varid(ncid_t2m, "lat",LatVarId)
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
status_t2m = nf90_get_var(ncid_t2m,LonVarId,lont2m)
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
status_t2m = nf90_get_var(ncid_t2m,LatVarId,latt2m)
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
status_q2m = nf90_inq_varid(ncid_q2m, "lon",LonVarId)
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
status_q2m = nf90_inq_varid(ncid_q2m, "lat",LatVarId)
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
status_q2m = nf90_get_var(ncid_q2m,LonVarId,lonq2m)
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
status_q2m = nf90_get_var(ncid_q2m,LatVarId,latq2m)
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
status_rain = nf90_inq_varid(ncid_rain, "lon",LonVarId)
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
status_rain = nf90_inq_varid(ncid_rain, "lat",LatVarId)
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
status_rain = nf90_get_var(ncid_rain,LonVarId,lonrain)
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
status_rain = nf90_get_var(ncid_rain,LatVarId,latrain)
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
status_dlw = nf90_inq_varid(ncid_dlw, "lon",LonVarId)
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
status_dlw = nf90_inq_varid(ncid_dlw, "lat",LatVarId)
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
status_dlw = nf90_get_var(ncid_dlw,LonVarId,londlw)
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
status_dlw = nf90_get_var(ncid_dlw,LatVarId,latdlw)
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
status_dsw = nf90_inq_varid(ncid_dsw, "lon",LonVarId)
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)
status_dsw = nf90_inq_varid(ncid_dsw, "lat",LatVarId)
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)
status_dsw = nf90_get_var(ncid_dsw,LonVarId,londsw)
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)
status_dsw = nf90_get_var(ncid_dsw,LatVarId,latdsw)
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)

status_cloud = nf90_inq_varid(ncid_cloud, "lon",LonVarId)
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)
status_cloud = nf90_inq_varid(ncid_cloud, "lat",LatVarId)
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)
status_cloud = nf90_get_var(ncid_cloud,LonVarId,loncloud)
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)
status_cloud = nf90_get_var(ncid_cloud,LatVarId,latcloud)
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)

! Get i,j positions of ROMS rho-points on CFSR grid
WHERE(lonw > 180.) lonw = lonw-360.
DO I=1,nlonw
IX_wind(I)= I
ENDDO
DO I=1,nlatw
IY_wind(I)= I
ENDDO
WHERE(lonmsl > 180.) lonmsl = lonmsl-360.
DO I=1,nlonmsl
IX_msl(I)= I
ENDDO
DO I=1,nlatmsl
IY_msl(I)= I
ENDDO
WHERE(lont2m > 180.) lont2m = lont2m-360.
DO I=1,nlont2m
IX_t2m(I)= I
ENDDO
DO I=1,nlatt2m
IY_t2m(I)= I
ENDDO
WHERE(lonq2m > 180.) lonq2m = lonq2m-360.
DO I=1,nlonq2m
IX_q2m(I)= I
ENDDO
DO I=1,nlatq2m
IY_q2m(I)= I
ENDDO
WHERE(lonrain > 180.) lonrain = lonrain-360.
DO I=1,nlonrain
IX_rain(I)= I
ENDDO
DO I=1,nlatrain
IY_rain(I)= I
ENDDO
WHERE(londlw > 180.) londlw = londlw-360.
DO I=1,nlondlw
IX_dlw(I)= I
ENDDO
DO I=1,nlatdlw
IY_dlw(I)= I
ENDDO
WHERE(londsw > 180.) londsw = londsw-360.
DO I=1,nlondsw
IX_dsw(I)= I
ENDDO
DO I=1,nlatdsw
IY_dsw(I)= I
ENDDO
WHERE(loncloud > 180.) loncloud = loncloud-360.
DO I=1,nloncloud
IX_cloud(I)= I
ENDDO
DO I=1,nlatcloud
IY_cloud(I)= I
ENDDO

CALL SSORT (lonw, IX_wind, nlonw, 2)
CALL SSORT (latw, IY_wind, nlatw, 2)
CALL SSORT (lonmsl, IX_msl, nlonmsl, 2)
CALL SSORT (latmsl, IY_msl, nlatmsl, 2)
CALL SSORT (lont2m, IX_t2m, nlont2m, 2)
CALL SSORT (latt2m, IY_t2m, nlatt2m, 2)
CALL SSORT (lonq2m, IX_q2m, nlonq2m, 2)
CALL SSORT (latq2m, IY_q2m, nlatq2m, 2)
CALL SSORT (lonrain, IX_rain, nlonrain, 2)
CALL SSORT (latrain, IY_rain, nlatrain, 2)
CALL SSORT (londlw, IX_dlw, nlondlw, 2)
CALL SSORT (latdlw, IY_dlw, nlatdlw, 2)
CALL SSORT (londsw, IX_dsw, nlondsw, 2)
CALL SSORT (latdsw, IY_dsw, nlatdsw, 2)
CALL SSORT (loncloud, IX_cloud, nloncloud, 2)
CALL SSORT (latcloud, IY_cloud, nlatcloud, 2)

DO j=1,Mp
DO i=1,Lp
DO k=1,nlonw-1
   if((lon_rho(i,j).ge.lonw(k)).AND.(lon_rho(i,j).lt.lonw(k+1)))then
   ipos_wind(i,j) = k
   endif
ENDDO
DO l=1,nlatw-1
   if((lat_rho(i,j).ge.latw(l)).AND.(lat_rho(i,j).lt.latw(l+1)))then
   jpos_wind(i,j) = l
   endif
ENDDO
DO k=1,nlonmsl-1
   if((lon_rho(i,j).ge.lonmsl(k)).AND.(lon_rho(i,j).lt.lonmsl(k+1)))then
   ipos_msl(i,j) = k
   endif
ENDDO
DO l=1,nlatmsl-1
   if((lat_rho(i,j).ge.latmsl(l)).AND.(lat_rho(i,j).lt.latmsl(l+1)))then
   jpos_msl(i,j) = l
   endif
ENDDO
DO k=1,nlont2m-1
   if((lon_rho(i,j).ge.lont2m(k)).AND.(lon_rho(i,j).lt.lont2m(k+1)))then
   ipos_t2m(i,j) = k
   endif
ENDDO
DO l=1,nlatt2m-1
   if((lat_rho(i,j).ge.latt2m(l)).AND.(lat_rho(i,j).lt.latt2m(l+1)))then
   jpos_t2m(i,j) = l
   endif
ENDDO
DO k=1,nlonq2m-1
   if((lon_rho(i,j).ge.lonq2m(k)).AND.(lon_rho(i,j).lt.lonq2m(k+1)))then
   ipos_q2m(i,j) = k
   endif
ENDDO
DO l=1,nlatq2m-1
   if((lat_rho(i,j).ge.latq2m(l)).AND.(lat_rho(i,j).lt.latq2m(l+1)))then
   jpos_q2m(i,j) = l
   endif
ENDDO
DO k=1,nlonrain-1
   if((lon_rho(i,j).ge.lonrain(k)).AND.(lon_rho(i,j).lt.lonrain(k+1)))then
   ipos_rain(i,j) = k
   endif
ENDDO
DO l=1,nlatrain-1
   if((lat_rho(i,j).ge.latrain(l)).AND.(lat_rho(i,j).lt.latrain(l+1)))then
   jpos_rain(i,j) = l
   endif
ENDDO
DO k=1,nlondlw-1
   if((lon_rho(i,j).ge.londlw(k)).AND.(lon_rho(i,j).lt.londlw(k+1)))then
   ipos_dlw(i,j) = k
   endif
ENDDO
DO l=1,nlatdlw-1
   if((lat_rho(i,j).ge.latdlw(l)).AND.(lat_rho(i,j).lt.latdlw(l+1)))then
   jpos_dlw(i,j) = l
   endif
ENDDO
DO k=1,nlondsw-1
   if((lon_rho(i,j).ge.londsw(k)).AND.(lon_rho(i,j).lt.londsw(k+1)))then
   ipos_dsw(i,j) = k
   endif
ENDDO
DO l=1,nlatdsw-1
   if((lat_rho(i,j).ge.latdsw(l)).AND.(lat_rho(i,j).lt.latdsw(l+1)))then
   jpos_dsw(i,j) = l
   endif
ENDDO
DO k=1,nloncloud-1
   if((lon_rho(i,j).ge.loncloud(k)).AND.(lon_rho(i,j).lt.loncloud(k+1)))then
   ipos_cloud(i,j) = k
   endif
ENDDO
DO l=1,nlatcloud-1
   if((lat_rho(i,j).ge.latcloud(l)).AND.(lat_rho(i,j).lt.latcloud(l+1)))then
   jpos_cloud(i,j) = l
   endif
ENDDO
ENDDO
ENDDO
! Find subarea in CFSR grid containing model grid
riminw = 1.e+23
rimaxw =-1.e+23
rjminw = 1.e+23
rjmaxw =-1.e+23
riminmsl = 1.e+23
rimaxmsl =-1.e+23
rjminmsl = 1.e+23
rjmaxmsl =-1.e+23
rimint2m = 1.e+23
rimaxt2m =-1.e+23
rjmint2m = 1.e+23
rjmaxt2m =-1.e+23
riminq2m = 1.e+23
rimaxq2m =-1.e+23
rjminq2m = 1.e+23
rjmaxq2m =-1.e+23
riminrain = 1.e+23
rimaxrain =-1.e+23
rjminrain = 1.e+23
rjmaxrain =-1.e+23
rimindlw = 1.e+23
rimaxdlw =-1.e+23
rjmindlw = 1.e+23
rjmaxdlw =-1.e+23
rimindsw = 1.e+23
rimaxdsw =-1.e+23
rjmindsw = 1.e+23
rjmaxdsw =-1.e+23
rimincloud = 1.e+23
rimaxcloud =-1.e+23
rjmincloud = 1.e+23
rjmaxcloud =-1.e+23
DO j=1,Mp
DO i=1,Lp
   riminw = MIN(riminw,ipos_wind(i,j))
   rimaxw = MAX(rimaxw,ipos_wind(i,j))
   rjminw = MIN(rjminw,jpos_wind(i,j))
   rjmaxw = MAX(rjmaxw,jpos_wind(i,j))
   riminmsl = MIN(riminmsl,ipos_msl(i,j))
   rimaxmsl = MAX(rimaxmsl,ipos_msl(i,j))
   rjminmsl = MIN(rjminmsl,jpos_msl(i,j))
   rjmaxmsl = MAX(rjmaxmsl,jpos_msl(i,j))
   rimint2m = MIN(rimint2m,ipos_t2m(i,j))
   rimaxt2m = MAX(rimaxt2m,ipos_t2m(i,j))
   rjmint2m = MIN(rjmint2m,jpos_t2m(i,j))
   rjmaxt2m = MAX(rjmaxt2m,jpos_t2m(i,j))
   riminq2m = MIN(riminq2m,ipos_q2m(i,j))
   rimaxq2m = MAX(rimaxq2m,ipos_q2m(i,j))
   rjminq2m = MIN(rjminq2m,jpos_q2m(i,j))
   rjmaxq2m = MAX(rjmaxq2m,jpos_q2m(i,j))
   riminrain = MIN(riminrain,ipos_rain(i,j))
   rimaxrain = MAX(rimaxrain,ipos_rain(i,j))
   rjminrain = MIN(rjminrain,jpos_rain(i,j))
   rjmaxrain = MAX(rjmaxrain,jpos_rain(i,j))
   rimindlw = MIN(rimindlw,ipos_dlw(i,j))
   rimaxdlw = MAX(rimaxdlw,ipos_dlw(i,j))
   rjmindlw = MIN(rjmindlw,jpos_dlw(i,j))
   rjmaxdlw = MAX(rjmaxdlw,jpos_dlw(i,j))
   rimindsw = MIN(rimindsw,ipos_dsw(i,j))
   rimaxdsw = MAX(rimaxdsw,ipos_dsw(i,j))
   rjmindsw = MIN(rjmindsw,jpos_dsw(i,j))
   rjmaxdsw = MAX(rjmaxdsw,jpos_dsw(i,j))
   rimincloud = MIN(rimincloud,ipos_cloud(i,j))
   rimaxcloud = MAX(rimaxcloud,ipos_cloud(i,j))
   rjmincloud = MIN(rjmincloud,jpos_cloud(i,j))
   rjmaxcloud = MAX(rjmaxcloud,jpos_cloud(i,j))
ENDDO
ENDDO
i1_wind = FLOOR(riminw)
i2_wind = CEILING(rimaxw)
j1_wind = FLOOR(rjminw)
j2_wind = CEILING(rjmaxw)
i1_msl = FLOOR(riminmsl)
i2_msl = CEILING(rimaxmsl)
j1_msl = FLOOR(rjminmsl)
j2_msl = CEILING(rjmaxmsl)
i1_t2m = FLOOR(rimint2m)
i2_t2m = CEILING(rimaxt2m)
j1_t2m = FLOOR(rjmint2m)
j2_t2m = CEILING(rjmaxt2m)
i1_q2m = FLOOR(riminq2m)
i2_q2m = CEILING(rimaxq2m)
j1_q2m = FLOOR(rjminq2m)
j2_q2m = CEILING(rjmaxq2m)
i1_rain = FLOOR(riminrain)
i2_rain = CEILING(rimaxrain)
j1_rain = FLOOR(rjminrain)
j2_rain = CEILING(rjmaxrain)
i1_dlw = FLOOR(rimindlw)
i2_dlw = CEILING(rimaxdlw)
j1_dlw = FLOOR(rjmindlw)
j2_dlw = CEILING(rjmaxdlw)
i1_dsw = FLOOR(rimindsw)
i2_dsw = CEILING(rimaxdsw)
j1_dsw = FLOOR(rjmindsw)
j2_dsw = CEILING(rjmaxdsw)
i1_cloud = FLOOR(rimincloud)
i2_cloud = CEILING(rimaxcloud)
j1_cloud = FLOOR(rjmincloud)
j2_cloud = CEILING(rjmaxcloud)

IF(i1_wind<1.OR.i2_wind>nlonw.OR.j1_wind<1.OR.j2_wind>nlatw) THEN
   write(*,*) 'min indice to extract in longitude: ',i1_wind
   write(*,*) 'max indice to extract in longitude: ',i2_wind
   write(*,*) 'min indice to extract in latitude: ',j1_wind
   write(*,*) 'max indice to extract in latitude: ',j2_wind
   write(*,*) 'Maximum dimension in longitude :',nlonw
   write(*,*) 'Maximum dimension in latitude :',nlatw
   STOP 'Wind subdomain outside CFSR domain'
ENDIF
IF(i1_msl<1.OR.i2_msl>nlonmsl.OR.j1_msl<1.OR.j2_msl>nlatmsl) THEN
   write(*,*) 'min indice to extract in longitude: ',i1_msl
   write(*,*) 'max indice to extract in longitude: ',i2_msl
   write(*,*) 'min indice to extract in latitude: ',j1_msl
   write(*,*) 'max indice to extract in latitude: ',j2_msl
   write(*,*) 'Maximum dimension in longitude :',nlonmsl
   write(*,*) 'Maximum dimension in latitude :',nlatmsl
   STOP 'MSL subdomain outside CFSR domain'
ENDIF
IF(i1_t2m<1.OR.i2_t2m>nlont2m.OR.j1_t2m<1.OR.j2_t2m>nlatt2m) THEN
   write(*,*) 'min indice to extract in longitude: ',i1_t2m
   write(*,*) 'max indice to extract in longitude: ',i2_t2m
   write(*,*) 'min indice to extract in latitude: ',j1_t2m
   write(*,*) 'max indice to extract in latitude: ',j2_t2m
   write(*,*) 'Maximum dimension in longitude :',nlont2m
   write(*,*) 'Maximum dimension in latitude :',nlatt2m
   STOP 'T2M subdomain outside CFSR domain'
ENDIF
IF(i1_q2m<1.OR.i2_q2m>nlonq2m.OR.j1_q2m<1.OR.j2_q2m>nlatq2m) THEN
   write(*,*) 'min indice to extract in longitude: ',i1_q2m
   write(*,*) 'max indice to extract in longitude: ',i2_q2m
   write(*,*) 'min indice to extract in latitude: ',j1_q2m
   write(*,*) 'max indice to extract in latitude: ',j2_q2m
   write(*,*) 'Maximum dimension in longitude :',nlonq2m
   write(*,*) 'Maximum dimension in latitude :',nlatq2m
   STOP 'Q2M subdomain outside CFSR domain'
ENDIF
IF(i1_rain<1.OR.i2_rain>nlonrain.OR.j1_rain<1.OR.j2_rain>nlatrain) THEN
   write(*,*) 'min indice to extract in longitude: ',i1_rain
   write(*,*) 'max indice to extract in longitude: ',i2_rain
   write(*,*) 'min indice to extract in latitude: ',j1_rain
   write(*,*) 'max indice to extract in latitude: ',j2_rain
   write(*,*) 'Maximum dimension in longitude :',nlonrain
   write(*,*) 'Maximum dimension in latitude :',nlatrain
   STOP 'RAIN subdomain outside CFSR domain'
ENDIF
IF(i1_dlw<1.OR.i2_dlw>nlondlw.OR.j1_dlw<1.OR.j2_dlw>nlatdlw) THEN
   write(*,*) 'min indice to extract in longitude: ',i1_dlw
   write(*,*) 'max indice to extract in longitude: ',i2_dlw
   write(*,*) 'min indice to extract in latitude: ',j1_dlw
   write(*,*) 'max indice to extract in latitude: ',j2_dlw
   write(*,*) 'Maximum dimension in longitude :',nlondlw
   write(*,*) 'Maximum dimension in latitude :',nlatdlw
   STOP 'DLW subdomain outside CFSR domain'
ENDIF
IF(i1_dsw<1.OR.i2_dsw>nlondsw.OR.j1_dsw<1.OR.j2_dsw>nlatdsw) THEN
   write(*,*) 'min indice to extract in longitude: ',i1_dsw
   write(*,*) 'max indice to extract in longitude: ',i2_dsw
   write(*,*) 'min indice to extract in latitude: ',j1_dsw
   write(*,*) 'max indice to extract in latitude: ',j2_dsw
   write(*,*) 'Maximum dimension in longitude :',nlondsw
   write(*,*) 'Maximum dimension in latitude :',nlatdsw
   STOP 'DSW subdomain outside CFSR domain'
ENDIF
IF(i1_cloud<1.OR.i2_cloud>nloncloud.OR.j1_cloud<1.OR.j2_cloud>nlatcloud) THEN
   write(*,*) 'min indice to extract in longitude: ',i1_cloud
   write(*,*) 'max indice to extract in longitude: ',i2_cloud
   write(*,*) 'min indice to extract in latitude: ',j1_cloud
   write(*,*) 'max indice to extract in latitude: ',j2_cloud
   write(*,*) 'Maximum dimension in longitude :',nloncloud
   write(*,*) 'Maximum dimension in latitude :',nlatcloud
   STOP 'CLOUD subdomain outside CFSR domain'
ENDIF

if(yes_no.eq.0)then
i1_wind = i1_wind -5
if(i1_wind.lt.1) i1_wind=1
i2_wind =i2_wind +5
if(i2_wind.gt.nlonw) i2_wind= nlonw
j1_wind= j1_wind-5
if(j1_wind.lt.1) j1_wind=1
j2_wind =j2_wind +5
if(j2_wind.gt.nlatw) j2_wind= nlatw
LP_wind= i2_wind-i1_wind+1
MP_wind= j2_wind-j1_wind+1
i1_t2m = i1_t2m -5
if(i1_t2m.lt.1) i1_t2m=1
i2_t2m =i2_t2m +5
if(i2_t2m.gt.nlont2m) i2_t2m= nlont2m
j1_t2m= j1_t2m-5
if(j1_t2m.lt.1) j1_t2m=1
j2_t2m =j2_t2m +5
if(j2_t2m.gt.nlatt2m) j2_t2m= nlatt2m
LP_t2m= i2_t2m-i1_t2m+1
MP_t2m= j2_t2m-j1_t2m+1
i1_q2m = i1_q2m -5
if(i1_q2m.lt.1) i1_q2m=1
i2_q2m =i2_q2m +5
if(i2_q2m.gt.nlonq2m) i2_q2m= nlonq2m
j1_q2m= j1_q2m-5
if(j1_q2m.lt.1) j1_q2m=1
j2_q2m =j2_q2m +5
if(j2_q2m.gt.nlatq2m) j2_q2m= nlatq2m
LP_q2m= i2_q2m-i1_q2m+1
MP_q2m= j2_q2m-j1_q2m+1
i1_rain = i1_rain -5
if(i1_rain.lt.1) i1_rain=1
i2_rain =i2_rain +5
if(i2_rain.gt.nlonrain) i2_rain= nlonrain
j1_rain= j1_rain-5
if(j1_rain.lt.1) j1_rain=1
j2_rain =j2_rain +5
if(j2_rain.gt.nlatrain) j2_rain= nlatrain
LP_rain= i2_rain-i1_rain+1
MP_rain= j2_rain-j1_rain+1
i1_dlw = i1_dlw -5
if(i1_dlw.lt.1) i1_dlw=1
i2_dlw =i2_dlw +5
if(i2_dlw.gt.nlondlw) i2_dlw= nlondlw
j1_dlw= j1_dlw-5
if(j1_dlw.lt.1) j1_dlw=1
j2_dlw =j2_dlw +5
if(j2_dlw.gt.nlatdlw) j2_dlw= nlatdlw
LP_dlw= i2_dlw-i1_dlw+1
MP_dlw= j2_dlw-j1_dlw+1
i1_dsw = i1_dsw -5
if(i1_dsw.lt.1) i1_dsw=1
i2_dsw =i2_dsw +5
if(i2_dsw.gt.nlondsw) i2_dsw= nlondsw
j1_dsw= j1_dsw-5
if(j1_dsw.lt.1) j1_dsw=1
j2_dsw =j2_dsw +5
if(j2_dsw.gt.nlatdsw) j2_dsw= nlatdsw
LP_dsw= i2_dsw-i1_dsw+1
MP_dsw= j2_dsw-j1_dsw+1
i1_cloud = i1_cloud -5
if(i1_cloud.lt.1) i1_cloud=1
i2_cloud =i2_cloud +5
if(i2_cloud.gt.nloncloud) i2_cloud= nloncloud
j1_cloud= j1_cloud-5
if(j1_cloud.lt.1) j1_cloud=1
j2_cloud =j2_cloud +5
if(j2_cloud.gt.nlatcloud) j2_cloud= nlatcloud
LP_cloud= i2_cloud-i1_cloud+1
MP_cloud= j2_cloud-j1_cloud+1
i1_msl = i1_msl-6
if(i1_msl.lt.1) i1_msl=1
i2_msl =i2_msl+6
if(i2_msl.gt.nlonmsl) i2_msl= nlonmsl
j1_msl= j1_msl-6
if(j1_msl.lt.1) j1_msl=1
j2_msl =j2_msl+6
if(j2_msl.gt.nlatmsl) j2_msl= nlatmsl
LP_msl= i2_msl-i1_msl+1
MP_msl= j2_msl-j1_msl+1
else
LP_wind=LP
MP_wind=MP
LP_t2m=LP
MP_t2m=MP
LP_q2m=LP
MP_q2m=MP
LP_rain=LP
MP_rain=MP
LP_dlw=LP
MP_dlw=MP
LP_dsw=LP
MP_dsw=MP
LP_cloud=LP
MP_cloud=MP
LP_msl=LP
MP_msl=MP
end if

if((LP_wind.eq.LP_t2m).and.(LP_t2m.eq.LP_q2m).and.(LP_q2m.eq.LP_rain).and.(LP_rain.eq.LP_cloud).and.(LP_cloud.eq.LP_dlw).and.(LP_dlw.eq.LP_dsw).and.(LP_dsw.eq.LP_wind))then
LP_dim= LP_wind
else
write(*,*) 'Problem of dimension following the longitude: '
write(*,*) 'Wind dimension: ',LP_wind
write(*,*) 'T2M dimension: ',LP_t2m
write(*,*) 'Q2M dimension: ',LP_q2m
write(*,*) 'RAIN dimension: ',LP_rain
write(*,*) 'DLW dimension: ',LP_dlw
write(*,*) 'DSW dimension: ',LP_dsw
write(*,*) 'CLOUD dimension: ',LP_cloud
STOP
end if

if((MP_wind.eq.MP_t2m).and.(MP_t2m.eq.MP_q2m).and.(MP_q2m.eq.MP_rain).and.(MP_rain.eq.MP_cloud).and.(MP_cloud.eq.MP_dlw).and.(MP_dlw.eq.MP_dsw).and.(MP_dsw.eq.MP_wind))then
MP_dim= MP_wind
else
write(*,*) 'Problem of dimension following the latitude: '
write(*,*) 'Wind dimension: ',MP_wind
write(*,*) 'T2M dimension: ',MP_t2m
write(*,*) 'Q2M dimension: ',MP_q2m
write(*,*) 'RAIN dimension: ',MP_rain
write(*,*) 'DLW dimension: ',MP_dlw
write(*,*) 'DSW dimension: ',MP_dsw
write(*,*) 'CLOUD dimension: ',MP_cloud
STOP
end if

LP_dim_msl = LP_msl
MP_dim_msl = MP_msl


!if(yes_no.eq.0)then
! MSL variable only available on 0.5deg x 0.5deg grid -> interpolation on finer grid needed
! Now letting ROMS do this as using a different file
!ALLOCATE(lonwtemp(Lp_dim))
!ALLOCATE(latwtemp(Mp_dim))
!ALLOCATE(lonmsltemp(Lp_msl))
!ALLOCATE(latmsltemp(Mp_msl))
!ALLOCATE(ipos_msltemp(Lp_dim))
!ALLOCATE(jpos_msltemp(Mp_dim))
!lonmsltemp = lonmsl(i1_msl:i2_msl)
!latmsltemp = latmsl(j1_msl:j2_msl)
!lonwtemp = lonw(i1_wind:i2_wind)
!latwtemp = latw(j1_wind:j2_wind)
!DO i=1,Lp_dim
!DO k=1,Lp_msl-1
!   if((lonwtemp(i).ge.lonmsltemp(k)).AND.(lonwtemp(i).lt.lonmsltemp(k+1)))then
!   ipos_msltemp(i) = k
!   endif
!ENDDO
!ENDDO
!DO j=1,Mp_dim
!DO l=1,Mp_msl-1
!   if((latwtemp(j).ge.latmsltemp(l)).AND.(latwtemp(j).lt.latmsltemp(l+1)))then
!   jpos_msltemp(j) = l
!   endif
!ENDDO
!ENDDO
!riminmsl = 1.e+23
!rimaxmsl =-1.e+23
!rjminmsl = 1.e+23
!rjmaxmsl =-1.e+23
!DO j=1,Mp_dim
!DO i=1,Lp_dim
!   riminmsl = MIN(riminmsl,ipos_msltemp(i))
!   rimaxmsl = MAX(rimaxmsl,ipos_msltemp(i))
!   rjminmsl = MIN(rjminmsl,jpos_msltemp(j))
!   rjmaxmsl = MAX(rjmaxmsl,jpos_msltemp(j))
!ENDDO
!ENDDO
!i1_msltemp = FLOOR(riminmsl)
!i2_msltemp = CEILING(rimaxmsl)
!j1_msltemp = FLOOR(rjminmsl)
!j2_msltemp = CEILING(rjmaxmsl)
!IF(i1_msltemp<1.OR.i2_msltemp>Lp_msl.OR.j1_msltemp<1.OR.j2_msltemp>Mp_msl) THEN
!   write(*,*) 'min indice to extract in longitude: ',i1_msltemp
!   write(*,*) 'max indice to extract in longitude: ',i2_msltemp
!   write(*,*) 'min indice to extract in latitude: ',j1_msltemp
!   write(*,*) 'max indice to extract in latitude: ',j2_msltemp
!   write(*,*) 'Maximum dimension in longitude :',Lp_msl
!   write(*,*) 'Maximum dimension in latitude :',Mp_msl
!   STOP 'MSL subdomain outside of WIND domain'
!ENDIF


!end if

! -------------------------------------------------------------------
! Open output files


!Wind file
statuso = nf90_create(TRIM(fileout)//'wind10m.nc',nf90_clobber,ncidowind)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.1)then
   statuso = nf90_def_dim(ncidowind,"xi_rho",Lp_dim,XiRhoDimIDwind)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidowind,"eta_rho",Mp_dim,EtaRhoDimIDwind)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
else
   statuso = nf90_def_dim(ncidowind,"longitude",Lp_dim,XiRhoDimIDwind)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidowind,"latitude",Mp_dim,EtaRhoDimIDwind)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_def_dim(ncidowind,"time",nf90_unlimited,TimeOutDimIDwind)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   statuso = nf90_put_att(ncidowind,nf90_global,"title", &
        "Atmospheric forcing file - CFSR grid")
else
   statuso = nf90_put_att(ncidowind,nf90_global,"title", &
        "Atmospheric forcing file from CSFR - ROMS grid")
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidowind,nf90_global,"type", &
     "ROMS forcing file")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidowind,nf90_global,"history", &
     "Produced with cfsr2roms.f90")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   ! Latitude variable
   statuso = nf90_def_var(ncidowind,"latitude",nf90_float, &
        (/EtaRhoDimIDwind/), LatVarIdwind)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidowind, LatVarIdwind,"long_Name", &
        "latitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidowind, LatVarIdwind,"units", &
        "degree_east")
   ! Longitude variable
   statuso = nf90_def_var(ncidowind,"longitude",nf90_float, &
        (/ XiRhoDimIDwind /), LonVarIdwind)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidowind, LonVarIdwind,"long_Name", &
        "longitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidowind, LonVarIdwind,"units", &
        "degree_north")
end if
! Time variable
statuso = nf90_def_var(ncidowind,"time",nf90_float, &
     TimeOutDimIDwind, TimeOutVarIdwind)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidowind, TimeOutVarIdwind,"long_Name", &
     "Time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidowind, TimeOutVarIdwind,"units", &
     TRIM(time_units_att))
if(statuso /= nf90_NoErr) call handle_err(statuso)
! UWIND variable
statuso = nf90_def_var(ncidowind,"Uwind",nf90_float, &
     (/ XiRhoDimIDwind, EtaRhoDimIDwind, TimeOutDimIDwind /), UwindVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidowind, UwindVarId,"long_Name", &
     "Xi-component of wind")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidowind, UwindVarId,"units", &
     "meter second-1")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidowind, UwindVarId,"time", &
     "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   statuso = nf90_put_att(ncidowind, UwindVarId,"coordinates", &
        "longitude latitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
! VWIND variable
statuso = nf90_def_var(ncidowind,"Vwind",nf90_float, &
     (/ XiRhoDimIDwind, EtaRhoDimIDwind, TimeOutDimIDwind /), VwindVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidowind, VwindVarId,"long_Name", &
     "Eta-component of wind")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidowind, VwindVarId,"units", &
     "meter second-1")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidowind, VwindVarId,"time", &
     "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   statuso = nf90_put_att(ncidowind, VwindVarId,"coordinates", &
        "longitude latitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_enddef(ncidowind)
if(statuso /= nf90_NoErr) call handle_err(statuso)
WRITE(*,*) 'Output file ', TRIM(fileout)//'wind10m.nc', ' is created.'

!T2M file
statuso = nf90_create(TRIM(fileout)//'t2m.nc',nf90_clobber,ncidot2m)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.1)then
   statuso = nf90_def_dim(ncidot2m,"xi_rho",Lp_dim,XiRhoDimIDt2m)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidot2m,"eta_rho",Mp_dim,EtaRhoDimIDt2m)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
else
   statuso = nf90_def_dim(ncidot2m,"longitude",Lp_dim,XiRhoDimIDt2m)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidot2m,"latitude",Mp_dim,EtaRhoDimIDt2m)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_def_dim(ncidot2m,"time",nf90_unlimited,TimeOutDimIDt2m)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   statuso = nf90_put_att(ncidot2m,nf90_global,"title", &
        "Atmospheric forcing file - CFSR grid")
else
   statuso = nf90_put_att(ncidot2m,nf90_global,"title", &
        "Atmospheric forcing file from CSFR - ROMS grid")
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidot2m,nf90_global,"type", &
     "ROMS forcing file")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidot2m,nf90_global,"history", &
     "Produced with cfsr2roms.f90")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   ! Latitude variable
   statuso = nf90_def_var(ncidot2m,"latitude",nf90_float, &
        (/EtaRhoDimIDt2m/), LatVarIdt2m)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidot2m, LatVarIdt2m,"long_Name", &
        "latitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidot2m, LatVarIdt2m,"units", &
        "degree_east")
   ! Longitude variable
   statuso = nf90_def_var(ncidot2m,"longitude",nf90_float, &
        (/ XiRhoDimIDt2m /), LonVarIdt2m)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidot2m, LonVarIdt2m,"long_Name", &
        "longitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidot2m, LonVarIdt2m,"units", &
        "degree_north")
end if
! Time variable
statuso = nf90_def_var(ncidot2m,"time",nf90_float, &
     TimeOutDimIDt2m, TimeOutVarIdt2m)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidot2m, TimeOutVarIdt2m,"long_Name", &
     "Time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidot2m, TimeOutVarIdt2m,"units", &
     TRIM(time_units_att))
if(statuso /= nf90_NoErr) call handle_err(statuso)

! TAIR variable
statuso = nf90_def_var(ncidot2m,"Tair",nf90_float, &
                 (/ XiRhoDimIDt2m, EtaRhoDimIDt2m, TimeOutDimIDt2m /), TairVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidot2m, TairVarId,"long_Name", &
                 "Air tempcfsrture at 2m")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidot2m, TairVarId,"units", &
                 "degree C")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidot2m, TairVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncidot2m, TairVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if

statuso = nf90_enddef(ncidot2m)
if(statuso /= nf90_NoErr) call handle_err(statuso)
WRITE(*,*) 'Output file ', TRIM(fileout)//'t2m.nc', ' is created.'

!Q2M file
statuso = nf90_create(TRIM(fileout)//'q2m.nc',nf90_clobber,ncidoq2m)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.1)then
   statuso = nf90_def_dim(ncidoq2m,"xi_rho",Lp_dim,XiRhoDimIDq2m)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidoq2m,"eta_rho",Mp_dim,EtaRhoDimIDq2m)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
else
   statuso = nf90_def_dim(ncidoq2m,"longitude",Lp_dim,XiRhoDimIDq2m)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidoq2m,"latitude",Mp_dim,EtaRhoDimIDq2m)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_def_dim(ncidoq2m,"time",nf90_unlimited,TimeOutDimIDq2m)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   statuso = nf90_put_att(ncidoq2m,nf90_global,"title", &
        "Atmospheric forcing file - CFSR grid")
else
   statuso = nf90_put_att(ncidoq2m,nf90_global,"title", &
        "Atmospheric forcing file from CSFR - ROMS grid")
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidoq2m,nf90_global,"type", &
     "ROMS forcing file")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidoq2m,nf90_global,"history", &
     "Produced with cfsr2roms.f90")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   ! Latitude variable
   statuso = nf90_def_var(ncidoq2m,"latitude",nf90_float, &
        (/EtaRhoDimIDq2m/), LatVarIdq2m)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidoq2m, LatVarIdq2m,"long_Name", &
        "latitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidoq2m, LatVarIdq2m,"units", &
        "degree_east")
   ! Longitude variable
   statuso = nf90_def_var(ncidoq2m,"longitude",nf90_float, &
        (/ XiRhoDimIDq2m /), LonVarIdq2m)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidoq2m, LonVarIdq2m,"long_Name", &
        "longitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidoq2m, LonVarIdq2m,"units", &
        "degree_north")
end if
! Time variable
statuso = nf90_def_var(ncidoq2m,"time",nf90_float, &
     TimeOutDimIDq2m, TimeOutVarIdq2m)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidoq2m, TimeOutVarIdq2m,"long_Name", &
     "Time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidoq2m, TimeOutVarIdq2m,"units", &
     TRIM(time_units_att))
if(statuso /= nf90_NoErr) call handle_err(statuso)

! QAIR variable
statuso = nf90_def_var(ncidoq2m,"Qair",nf90_float, &
                 (/ XiRhoDimIDq2m, EtaRhoDimIDq2m, TimeOutDimIDq2m /), QairVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidoq2m, QairVarId,"long_Name", &
                 "Specific humidity at 2m")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidoq2m, QairVarId,"units", &
                 "kg kg-1")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidoq2m, QairVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncidoq2m, QairVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if

statuso = nf90_enddef(ncidoq2m)
if(statuso /= nf90_NoErr) call handle_err(statuso)
WRITE(*,*) 'Output file ', TRIM(fileout)//'q2m.nc', ' is created.'

!Cloud file
statuso = nf90_create(TRIM(fileout)//'cloud.nc',nf90_clobber,ncidocloud)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.1)then
   statuso = nf90_def_dim(ncidocloud,"xi_rho",Lp_dim,XiRhoDimIDcloud)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidocloud,"eta_rho",Mp_dim,EtaRhoDimIDcloud)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
else
   statuso = nf90_def_dim(ncidocloud,"longitude",Lp_dim,XiRhoDimIDcloud)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidocloud,"latitude",Mp_dim,EtaRhoDimIDcloud)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_def_dim(ncidocloud,"time",nf90_unlimited,TimeOutDimIDcloud)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   statuso = nf90_put_att(ncidocloud,nf90_global,"title", &
        "Atmospheric forcing file - CFSR grid")
else
   statuso = nf90_put_att(ncidocloud,nf90_global,"title", &
        "Atmospheric forcing file from CSFR - ROMS grid")
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidocloud,nf90_global,"type", &
     "ROMS forcing file")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidocloud,nf90_global,"history", &
     "Produced with cfsr2roms.f90")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   ! Latitude variable
   statuso = nf90_def_var(ncidocloud,"latitude",nf90_float, &
        (/EtaRhoDimIDcloud/), LatVarIdcloud)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidocloud, LatVarIdcloud,"long_Name", &
        "latitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidocloud, LatVarIdcloud,"units", &
        "degree_east")
   ! Longitude variable
   statuso = nf90_def_var(ncidocloud,"longitude",nf90_float, &
        (/ XiRhoDimIDcloud /), LonVarIdcloud)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidocloud, LonVarIdcloud,"long_Name", &
        "longitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidocloud, LonVarIdcloud,"units", &
        "degree_north")
end if
! Time variable
statuso = nf90_def_var(ncidocloud,"time",nf90_float, &
     TimeOutDimIDcloud, TimeOutVarIdcloud)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidocloud, TimeOutVarIdcloud,"long_Name", &
     "Time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidocloud, TimeOutVarIdcloud,"units", &
     TRIM(time_units_att))
if(statuso /= nf90_NoErr) call handle_err(statuso)

! CLOUD variable
statuso = nf90_def_var(ncidocloud,"cloud",nf90_float, &
                 (/ XiRhoDimIDcloud, EtaRhoDimIDcloud, TimeOutDimIDcloud /), cloudairVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidocloud, cloudairVarId,"long_Name", &
                 "Fraction cloud cover")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidocloud, cloudairVarId,"units", &
                 " ")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, cloudairVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncidocloud, cloudairVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_enddef(ncidocloud)
if(statuso /= nf90_NoErr) call handle_err(statuso)

WRITE(*,*) 'Output file ', TRIM(fileout)//'cloud.nc', ' is created.'

!Rain file
statuso = nf90_create(TRIM(fileout)//'rain.nc',nf90_clobber,ncidorain)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.1)then
   statuso = nf90_def_dim(ncidorain,"xi_rho",Lp_dim,XiRhoDimIDrain)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidorain,"eta_rho",Mp_dim,EtaRhoDimIDrain)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
else
   statuso = nf90_def_dim(ncidorain,"longitude",Lp_dim,XiRhoDimIDrain)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidorain,"latitude",Mp_dim,EtaRhoDimIDrain)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_def_dim(ncidorain,"time",nf90_unlimited,TimeOutDimIDrain)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   statuso = nf90_put_att(ncidorain,nf90_global,"title", &
        "Atmospheric forcing file - CFSR grid")
else
   statuso = nf90_put_att(ncidorain,nf90_global,"title", &
        "Atmospheric forcing file from CSFR - ROMS grid")
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidorain,nf90_global,"type", &
     "ROMS forcing file")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidorain,nf90_global,"history", &
     "Produced with cfsr2roms.f90")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   ! Latitude variable
   statuso = nf90_def_var(ncidorain,"latitude",nf90_float, &
        (/EtaRhoDimIDrain/), LatVarIdrain)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidorain, LatVarIdrain,"long_Name", &
        "latitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidorain, LatVarIdrain,"units", &
        "degree_east")
   ! Longitude variable
   statuso = nf90_def_var(ncidorain,"longitude",nf90_float, &
        (/ XiRhoDimIDrain /), LonVarIdrain)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidorain, LonVarIdrain,"long_Name", &
        "longitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidorain, LonVarIdrain,"units", &
        "degree_north")
end if
! Time variable
statuso = nf90_def_var(ncidorain,"time",nf90_float, &
     TimeOutDimIDrain, TimeOutVarIdrain)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidorain, TimeOutVarIdrain,"long_Name", &
     "Time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidorain, TimeOutVarIdrain,"units", &
     TRIM(time_units_att))
if(statuso /= nf90_NoErr) call handle_err(statuso)

! RAIN variable
statuso = nf90_def_var(ncidorain,"rain",nf90_float, &
                 (/ XiRhoDimIDrain, EtaRhoDimIDrain, TimeOutDimIDrain /), rainairVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidorain, rainairVarId,"long_Name", &
                 "Total precipiation")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidorain, rainairVarId,"units", &
                 "kg meter-2 second-1")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidorain, rainairVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncidorain, rainairVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_enddef(ncidorain)
if(statuso /= nf90_NoErr) call handle_err(statuso)

WRITE(*,*) 'Output file ', TRIM(fileout)//'rain.nc', ' is created.'

!Sward file
statuso = nf90_create(TRIM(fileout)//'dsw.nc',nf90_clobber,ncidodsw)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.1)then
   statuso = nf90_def_dim(ncidodsw,"xi_rho",Lp_dim,XiRhoDimIDdsw)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidodsw,"eta_rho",Mp_dim,EtaRhoDimIDdsw)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
else
   statuso = nf90_def_dim(ncidodsw,"longitude",Lp_dim,XiRhoDimIDdsw)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidodsw,"latitude",Mp_dim,EtaRhoDimIDdsw)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_def_dim(ncidodsw,"time",nf90_unlimited,TimeOutDimIDdsw)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   statuso = nf90_put_att(ncidodsw,nf90_global,"title", &
        "Atmospheric forcing file - CFSR grid")
else
   statuso = nf90_put_att(ncidodsw,nf90_global,"title", &
        "Atmospheric forcing file from CSFR - ROMS grid")
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodsw,nf90_global,"type", &
     "ROMS forcing file")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodsw,nf90_global,"history", &
     "Produced with cfsr2roms.f90")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   ! Latitude variable
   statuso = nf90_def_var(ncidodsw,"latitude",nf90_float, &
        (/EtaRhoDimIDdsw/), LatVarIddsw)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidodsw, LatVarIddsw,"long_Name", &
        "latitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidodsw, LatVarIddsw,"units", &
        "degree_east")
   ! Longitude variable
   statuso = nf90_def_var(ncidodsw,"longitude",nf90_float, &
        (/ XiRhoDimIDdsw /), LonVarIddsw)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidodsw, LonVarIddsw,"long_Name", &
        "longitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidodsw, LonVarIddsw,"units", &
        "degree_north")
end if
! Time variable
statuso = nf90_def_var(ncidodsw,"time",nf90_float, &
     TimeOutDimIDdsw, TimeOutVarIddsw)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodsw, TimeOutVarIddsw,"long_Name", &
     "Time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodsw, TimeOutVarIddsw,"units", &
     TRIM(time_units_att))
if(statuso /= nf90_NoErr) call handle_err(statuso)


! SWRAD variable
statuso = nf90_def_var(ncidodsw,"swrad",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), swradVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, swradVarId,"long_Name", &
                 "Downward shortwave radiation")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodsw, swradVarId,"units", &
                 "Watt meter-2")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodsw, swradVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncidodsw, swradVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_enddef(ncidodsw)
if(statuso /= nf90_NoErr) call handle_err(statuso)

WRITE(*,*) 'Output file ', TRIM(fileout)//'dsw.nc', ' is created.'

!Sward file
statuso = nf90_create(TRIM(fileout)//'dlw.nc',nf90_clobber,ncidodlw)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.1)then
   statuso = nf90_def_dim(ncidodlw,"xi_rho",Lp_dim,XiRhoDimIDdlw)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidodlw,"eta_rho",Mp_dim,EtaRhoDimIDdlw)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
else
   statuso = nf90_def_dim(ncidodlw,"longitude",Lp_dim,XiRhoDimIDdlw)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidodlw,"latitude",Mp_dim,EtaRhoDimIDdlw)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_def_dim(ncidodlw,"time",nf90_unlimited,TimeOutDimIDdlw)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   statuso = nf90_put_att(ncidodlw,nf90_global,"title", &
        "Atmospheric forcing file - CFSR grid")
else
   statuso = nf90_put_att(ncidodlw,nf90_global,"title", &
        "Atmospheric forcing file from CSFR - ROMS grid")
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodlw,nf90_global,"type", &
     "ROMS forcing file")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodlw,nf90_global,"history", &
     "Produced with cfsr2roms.f90")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   ! Latitude variable
   statuso = nf90_def_var(ncidodlw,"latitude",nf90_float, &
        (/EtaRhoDimIDdlw/), LatVarIddlw)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidodlw, LatVarIddlw,"long_Name", &
        "latitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidodlw, LatVarIddlw,"units", &
        "degree_east")
   ! Longitude variable
   statuso = nf90_def_var(ncidodlw,"longitude",nf90_float, &
        (/ XiRhoDimIDdlw /), LonVarIddlw)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidodlw, LonVarIddlw,"long_Name", &
        "longitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidodlw, LonVarIddlw,"units", &
        "degree_north")
end if
! Time variable
statuso = nf90_def_var(ncidodlw,"time",nf90_float, &
     TimeOutDimIDdlw, TimeOutVarIddlw)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodlw, TimeOutVarIddlw,"long_Name", &
     "Time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodlw, TimeOutVarIddlw,"units", &
     TRIM(time_units_att))
if(statuso /= nf90_NoErr) call handle_err(statuso)

! LWRAD variable
statuso = nf90_def_var(ncidodlw,"lwrad_down",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), lwradVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodlw, lwradVarId,"long_Name", &
                 "Downward longwave radiation")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodlw, lwradVarId,"units", &
                 "Watt meter-2")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidodlw, lwradVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncidodlw, lwradVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_enddef(ncidodlw)
if(statuso /= nf90_NoErr) call handle_err(statuso)
WRITE(*,*) 'Output file ', TRIM(fileout)//'dlw.nc', ' is created.'

statuso = nf90_enddef(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

!MSL file
!---- Has different lons and lats dimensions.
statuso = nf90_create(TRIM(fileout)//'msl.nc',nf90_clobber,ncidomsl)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.1)then
   statuso = nf90_def_dim(ncidomsl,"xi_rho",Lp_dim,XiRhoDimIDmsl)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidomsl,"eta_rho",Mp_dim,EtaRhoDimIDmsl)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
else
   statuso = nf90_def_dim(ncidomsl,"longitude",Lp_dim_msl,XiRhoDimIDmsl)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_def_dim(ncidomsl,"latitude",Mp_dim_msl,EtaRhoDimIDmsl)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_def_dim(ncidomsl,"time",nf90_unlimited,TimeOutDimIDmsl)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   statuso = nf90_put_att(ncidomsl,nf90_global,"title", &
        "Atmospheric forcing file - CFSR grid")
else
   statuso = nf90_put_att(ncidomsl,nf90_global,"title", &
        "Atmospheric forcing file from CSFR - ROMS grid")
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidomsl,nf90_global,"type", &
     "ROMS forcing file")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidomsl,nf90_global,"history", &
     "Produced with cfsr2roms.f90")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
   ! Latitude variable
   statuso = nf90_def_var(ncidomsl,"latitude",nf90_float, &
        (/EtaRhoDimID/), LatVarId)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidomsl, LatVarIdmsl,"long_Name", &
        "latitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidomsl, LatVarIdmsl,"units", &
        "degree_east")
   ! Longitude variable
   statuso = nf90_def_var(ncidomsl,"longitude",nf90_float, &
        (/ XiRhoDimID /), LonVarId)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidomsl, LonVarIdmsl,"long_Name", &
        "longitude")
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_att(ncidomsl, LonVarIdmsl,"units", &
        "degree_north")
end if
! Time variable
statuso = nf90_def_var(ncidomsl,"time",nf90_float, &
     TimeOutDimIDmsl, TimeOutVarIdmsl)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidomsl, TimeOutVarIdmsl,"long_Name", &
     "Time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidomsl, TimeOutVarIdmsl,"units", &
     TRIM(time_units_att))
if(statuso /= nf90_NoErr) call handle_err(statuso)

! PAIR variable
statuso = nf90_def_var(ncidomsl,"Pair",nf90_float, &
                 (/ XiRhoDimIDmsl, EtaRhoDimIDmsl, TimeOutDimIDmsl /), PairVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidomsl, PairVarId,"long_Name", &
                 "Atmospheric pressure")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidomsl, PairVarId,"units", &
                 "Pa")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncidomsl, PairVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncidomsl, PairVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if

statuso = nf90_enddef(ncidomsl)
if(statuso /= nf90_NoErr) call handle_err(statuso)
WRITE(*,*) 'Output file ', TRIM(fileout)//'msl.nc', ' is created.'



ALLOCATE(Uwind(Lp_dim,Mp_dim))
ALLOCATE(Vwind(Lp_dim,Mp_dim))
ALLOCATE(Pairtemp(Lp_msl,Mp_msl))
ALLOCATE(Pair(Lp_dim,Mp_dim)) 
ALLOCATE(Tair(Lp_dim,Mp_dim))
ALLOCATE(Qair(Lp_dim,Mp_dim))
ALLOCATE(cloudair(Lp_dim,Mp_dim))
ALLOCATE(rainair(Lp_dim,Mp_dim))
ALLOCATE(swrad(Lp_dim,Mp_dim))
ALLOCATE(lwrad(Lp_dim,Mp_dim))

!Open output file
statuso = nf90_open(trim(fileout)//'wind.nc',nf90_write,ncidowind)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_open(trim(fileout)//'msl.nc',nf90_write,ncidomsl)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_open(trim(fileout)//'q2m.nc',nf90_write,ncidoq2m)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_open(trim(fileout)//'t2m.nc',nf90_write,ncidot2m)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_open(trim(fileout)//'cloud.nc',nf90_write,ncidocloud)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_open(trim(fileout)//'rain.nc',nf90_write,ncidorain)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_open(trim(fileout)//'dsw.nc',nf90_write,ncidodsw)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_open(trim(fileout)//'dlw.nc',nf90_write,ncidodlw)
if(statuso /= nf90_NoErr) call handle_err(statuso)

if(yes_no.eq.0)then
   ! Write latitude and longitude
   statuso = nf90_put_var(ncidowind,LonVarIdwind,lonw(i1_wind:i2_wind))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncidowind,LatVarIdwind,latw(j1_wind:j2_wind))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

   statuso = nf90_put_var(ncidoq2m,LonVarIdq2m,lonq2m(i1_q2m:i2_q2m))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncidoq2m,LatVarIdq2m,latq2m(j1_q2m:j2_q2m))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

   statuso = nf90_put_var(ncidot2m,LonVarIdt2m,lont2m(i1_t2m:i2_t2m))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncidot2m,LatVarIdt2m,latt2m(j1_t2m:j2_t2m))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

   statuso = nf90_put_var(ncidocloud,LonVarIdcloud,loncloud(i1_cloud:i2_cloud))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncidocloud,LatVarIdcloud,latcloud(j1_cloud:j2_cloud))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

   statuso = nf90_put_var(ncidorain,LonVarIdrain,lonrain(i1_rain:i2_rain))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncidorain,LatVarIdrain,latrain(j1_rain:j2_rain))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

   statuso = nf90_put_var(ncidodsw,LonVarIddsw,londsw(i1_dsw:i2_dsw))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncidodsw,LatVarIddsw,latdsw(j1_dsw:j2_dsw))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

   statuso = nf90_put_var(ncidodlw,LonVarIddlw,londlw(i1_dlw:i2_dlw))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncidodlw,LatVarIddlw,latdlw(j1_dlw:j2_dlw))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

! Different lons and lats for MSL
   statuso = nf90_put_var(ncidomsl,LonVarIdmsl,lonmsl(i1_msl:i2_msl))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncidomsl,LatVarIdmsl,latmsl(j1_msl:j2_msl))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

end if

endif

! Get wind time var id
status_wind = nf90_inq_varid(ncid_wind, "time",TimeVarId_wind)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
! Get U10 var id
status_wind = nf90_inq_varid(ncid_wind, "U_GRD_L103",U10VarId)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
! Get V10 var id
status_wind = nf90_inq_varid(ncid_wind, "V_GRD_L103",V10VarId)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
! Get MSL time var id
status_msl = nf90_inq_varid(ncid_msl, "time",TimeVarId_msl)
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
! Get MSL var id
status_msl = nf90_inq_varid(ncid_msl, "PRMSL_L101",MSLVarId)
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
! Get T2M time var id
status_t2m = nf90_inq_varid(ncid_t2m, "time",TimeVarId_t2m)
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
! Get T2M var id
status_t2m = nf90_inq_varid(ncid_t2m, "TMP_L103",T2MVarId)
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
! Get Q2M time var id
status_q2m = nf90_inq_varid(ncid_q2m, "time",TimeVarId_q2m)
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
! Get Q2M var id
status_q2m = nf90_inq_varid(ncid_q2m, "SPF_H_L103",Q2MVarId)
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
! Get RAIN time var id
status_rain = nf90_inq_varid(ncid_rain, "time",TimeVarId_rain)
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
! Get RAIN var id
status_rain = nf90_inq_varid(ncid_rain, "PRATE_L1",RAINVarId)
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
! Get DLW time var id
status_dlw = nf90_inq_varid(ncid_dlw, "time",TimeVarId_dlw)
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
! Get DLW var id
status_dlw = nf90_inq_varid(ncid_dlw, "DLWRF_L1",DLWVarId)
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
! Get DSW time var id
status_dsw = nf90_inq_varid(ncid_dsw, "time",TimeVarId_dsw)
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)
! Get DSW var id
status_dsw = nf90_inq_varid(ncid_dsw, "DSWRF_L1",DSWVarId)
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)
! Get CLOUD time var id
status_cloud = nf90_inq_varid(ncid_cloud, "time",TimeVarId_cloud)
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)
! Get CLOUD var id
status_cloud = nf90_inq_varid(ncid_cloud, "T_CDC_L200",CLOUDVarId)
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)

! Read in Wind Time   
status_wind = nf90_get_var(ncid_wind,TimeVarId_wind,datetime,start = (/itime/))
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
tday_wind= REAL(jdref_cfsr) + datetime/24.d0 - REAL(jdref)
tday_wind= nint(tday_wind*24.d0)*3600.d0
! Read in U10
status_wind = nf90_get_var(ncid_wind,U10VarId,U10,start=(/1, 1,itime /))
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
! Blank out Antarctica
WHERE(ABS(U10)>1.e+7) U10 = 0.
! Read in V10
status_wind = nf90_get_var(ncid_wind,V10VarId,V10,start=(/1, 1,itime /))
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
! Blank out Antarctica
WHERE(ABS(V10)>1.e+7) V10 = 0.

! Read in Time   
status_msl = nf90_get_var(ncid_msl,TimeVarId_msl,datetime,start = (/itime/))
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
tday_msl= REAL(jdref_cfsr) + datetime/24.d0 - REAL(jdref)
tday_msl= nint(tday_msl*24.d0)*3600.d0
! Read in MSL values
status_msl = nf90_get_var(ncid_msl,MSLVarId,MSL,start=(/1, 1,itime /))
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
! Blank out Antarctica
WHERE(ABS(MSL)>1.e+7) MSL = 1.e+5

! Read in Time   
status_t2m = nf90_get_var(ncid_t2m,TimeVarId_t2m,datetime,start = (/itime/))
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
tday_t2m= REAL(jdref_cfsr) + datetime/24.d0 - REAL(jdref)
tday_t2m= nint(tday_t2m*24.d0)*3600.d0
! Read in T2M values
status_t2m = nf90_get_var(ncid_t2m,T2MVarId,T2M,start=(/1, 1,itime /))
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
! Convert from Kelvin to Celsius
T2M = T2M - 273.16
! Blank out Antarctica
WHERE(ABS(T2M)>1.e+7) T2M = 0.

! Read in Time   
status_q2m = nf90_get_var(ncid_q2m,TimeVarId_q2m,datetime,start = (/itime/))
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
tday_q2m= REAL(jdref_cfsr) + datetime/24.d0 - REAL(jdref)
tday_q2m= nint(tday_q2m*24.d0)*3600.d0
! Read in Q2M values
status_q2m = nf90_get_var(ncid_q2m,Q2MVarId,Q2M,start=(/1, 1,itime /))
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
! Blank out Antarctica
WHERE(ABS(Q2M)>1.e+7) Q2M = 1.e-5

! Read in Time   
status_cloud = nf90_get_var(ncid_cloud,TimeVarId_cloud,datetime,start = (/itime/))
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)
tday_cloud= REAL(jdref_cfsr) + datetime/24.d0 - REAL(jdref)
tday_cloud= nint(tday_cloud*24.d0)*3600.d0
! Read in CLOUD values
status_cloud = nf90_get_var(ncid_cloud,CLOUDVarId,CLOUD,start=(/1, 1,itime /))
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)
! Blank out Antarctica
WHERE(ABS(CLOUD)>1.e+7) CLOUD = 0.

! Read in Time   
status_rain = nf90_get_var(ncid_rain,TimeVarId_rain,datetime,start = (/itime/))
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
tday_rain= REAL(jdref_cfsr) + datetime/24.d0 - REAL(jdref)
tday_rain= nint(tday_rain*24.d0)*3600.d0
! Read in RAIN values
status_rain = nf90_get_var(ncid_rain,RAINVarId,RAIN,start=(/1, 1,itime /))
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
! Blank out Antarctica
WHERE(ABS(RAIN)>1.e+7) RAIN = 0.

! Read in Time   
status_dsw = nf90_get_var(ncid_dsw,TimeVarId_dsw,datetime,start = (/itime/))
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)
tday_dsw= REAL(jdref_cfsr) + datetime/24.d0 - REAL(jdref)
tday_dsw= nint(tday_dsw*24.d0)*3600.d0
! Read in DSW values
status_dsw = nf90_get_var(ncid_dsw,DSWVarId,DSW,start=(/1, 1,itime /))
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)
! Blank out Antarctica
WHERE(ABS(DSW)>1.e+7) DSW = 0.

! Read in Time   
status_dlw = nf90_get_var(ncid_dlw,TimeVarId_dlw,datetime,start = (/itime/))
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
tday_dlw= REAL(jdref_cfsr) + datetime/24.d0 - REAL(jdref)
tday_dlw= nint(tday_dlw*24.d0)*3600.d0
! Read in DLW values
status_dlw = nf90_get_var(ncid_dlw,DLWVarId,DLW,start=(/1, 1,itime /))
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
! Blank out Antarctica
WHERE(ABS(DLW)>1.e+7) DLW = 0.


if((tday_wind.eq.tday_msl).and.(tday_msl.eq.tday_t2m).and.(tday_t2m.eq.tday_q2m).and.(tday_q2m.eq.tday_rain).and.(tday_rain.eq.tday_cloud).and.(tday_cloud.eq.tday_dlw).and.(tday_dlw.eq.tday_dsw).and.(tday_dsw.eq.tday_wind))then
else
write(*,*) "Times are not the same for all files"
!write(*,*) 'Problem with the time: '
!write(*,*) 'Wind time: ',tday_wind
!write(*,*) 'MSL time: ',tday_msl
!write(*,*) 'T2M time: ',tday_t2m
!write(*,*) 'Q2M time: ',tday_q2m
!write(*,*) 'RAIN time: ',tday_rain
!write(*,*) 'DLW time: ',tday_dlw
!write(*,*) 'DSW time: ',tday_dsw
!write(*,*) 'CLOUD time: ',tday_cloud
!STOP
end if

work_wind=0.
DO I=1,nlonw
   DO J=1,nlatw
      work_wind(I,J)=U10(IX_wind(I),IY_wind(J))
   ENDDO
ENDDO
U10= work_wind
work_wind=0.
DO I=1,nlonw
   DO J=1,nlatw
      work_wind(I,J)=V10(IX_wind(I),IY_wind(J))
   ENDDO
ENDDO
V10= work_wind
work_msl=0.
DO I=1,nlonmsl
   DO J=1,nlatmsl
      work_msl(I,J)=MSL(IX_msl(I),IY_msl(J))
   ENDDO
ENDDO
MSL= work_msl
work_t2m=0.
DO I=1,nlont2m
   DO J=1,nlatt2m
      work_t2m(I,J)=T2M(IX_t2m(I),IY_t2m(J))
   ENDDO
ENDDO
T2M= work_t2m
work_q2m=0.
DO I=1,nlonq2m
   DO J=1,nlatq2m
      work_q2m(I,J)=Q2M(IX_q2m(I),IY_q2m(J))
   ENDDO
ENDDO
Q2M= work_q2m
work_cloud=0.
DO I=1,nloncloud
   DO J=1,nlatcloud
      work_cloud(I,J)=CLOUD(IX_cloud(I),IY_cloud(J))
   ENDDO
ENDDO
CLOUD= work_cloud
work_rain=0.
DO I=1,nlonrain
   DO J=1,nlatrain
      work_rain(I,J)=RAIN(IX_rain(I),IY_rain(J))
   ENDDO
ENDDO
RAIN= work_rain
work_dlw=0.
DO I=1,nlondlw
   DO J=1,nlatdlw
      work_dlw(I,J)=DLW(IX_dlw(I),IY_dlw(J))
   ENDDO
ENDDO
DLW= work_dlw
work_dsw=0.
DO I=1,nlondsw
   DO J=1,nlatdsw
      work_dsw(I,J)=DSW(IX_dsw(I),IY_dsw(J))
   ENDDO
ENDDO
DSW= work_dsw

if(yes_no.eq.1)then
! Horizontal interpolation
scr1_out = 0.
scr2_out = 0.
Pair=0.
Pairtemp=0.
Tair=0.
Qair=0.
Cloudair=0.
rainair= 0.
DO j=1,Mp
   DO i=1,Lp
      ia = ipos_wind(i,j)
      ja = jpos_wind(i,j)
      rx = (lon_rho(i,j)-lonw(ia))/(lonw(ia+1)-lonw(ia))
      ry = (lat_rho(i,j)-latw(ja))/(latw(ja+1)-latw(ja))
      rxm = 1.0-rx
      rym = 1.0-ry     
      scr1_out(i,j) = rxm*rym*U10(ia,ja)           +   &
           rx*rym*U10(ia+1,ja)                     +   &
           rx*ry*U10(ia+1,ja+1)                    +   &
           rxm*ry*U10(ia,ja+1)      
      scr2_out(i,j) = rxm*rym*V10(ia,ja)           +   &
           rx*rym*V10(ia+1,ja)                     +   &
           rx*ry*V10(ia+1,ja+1)                    +   &
           rxm*ry*V10(ia,ja+1)      
      ia = ipos_msl(i,j)
      ja = jpos_msl(i,j)
      rx = (lon_rho(i,j)-lonmsl(ia))/(lonmsl(ia+1)-lonmsl(ia))
      ry = (lat_rho(i,j)-latmsl(ja))/(latmsl(ja+1)-latmsl(ja))
      rxm = 1.0-rx
      rym = 1.0-ry     
      Pair(i,j)     = rxm*rym*MSL(ia,ja)           +   &
           rx*rym*MSL(ia+1,ja)                     +   &
           rx*ry*MSL(ia+1,ja+1)                    +   &
           rxm*ry*MSL(ia,ja+1)      
      ia = ipos_t2m(i,j)
      ja = jpos_t2m(i,j)
      rx = (lon_rho(i,j)-lont2m(ia))/(lont2m(ia+1)-lont2m(ia))
      ry = (lat_rho(i,j)-latt2m(ja))/(latt2m(ja+1)-latt2m(ja))
      rxm = 1.0-rx
      rym = 1.0-ry     
      Tair(i,j)     = rxm*rym*T2M(ia,ja)           +   &
           rx*rym*T2M(ia+1,ja)                     +   &
           rx*ry*T2M(ia+1,ja+1)                    +   &
           rxm*ry*T2M(ia,ja+1)      
      ia = ipos_q2m(i,j)
      ja = jpos_q2m(i,j)
      rx = (lon_rho(i,j)-lonq2m(ia))/(lonq2m(ia+1)-lonq2m(ia))
      ry = (lat_rho(i,j)-latq2m(ja))/(latq2m(ja+1)-latq2m(ja))
      rxm = 1.0-rx
      rym = 1.0-ry     
      Qair(i,j)     = rxm*rym*Q2M(ia,ja)           +   &
           rx*rym*Q2M(ia+1,ja)                     +   &
           rx*ry*Q2M(ia+1,ja+1)                    +   &
           rxm*ry*Q2M(ia,ja+1)
      ia = ipos_cloud(i,j)
      ja = jpos_cloud(i,j)
      rx = (lon_rho(i,j)-loncloud(ia))/(loncloud(ia+1)-loncloud(ia))
      ry = (lat_rho(i,j)-latcloud(ja))/(latcloud(ja+1)-latcloud(ja))
      rxm = 1.0-rx
      rym = 1.0-ry     
      cloudair(i,j)    = (rxm*rym*CLOUD(ia,ja)           +   &
           rx*rym*CLOUD(ia+1,ja)                     +   &
           rx*ry*CLOUD(ia+1,ja+1)                    +   &
           rxm*ry*CLOUD(ia,ja+1))/100.
      ia = ipos_rain(i,j)
      ja = jpos_rain(i,j)
      rx = (lon_rho(i,j)-lonrain(ia))/(lonrain(ia+1)-lonrain(ia))
      ry = (lat_rho(i,j)-latrain(ja))/(latrain(ja+1)-latrain(ja))
      rxm = 1.0-rx
      rym = 1.0-ry     
      rainair(i,j)     = rxm*rym*RAIN(ia,ja)            +   &
           rx*rym*RAIN(ia+1,ja)                      +   &
           rx*ry*RAIN(ia+1,ja+1)                     +   &
           rxm*ry*RAIN(ia,ja+1)
      ia = ipos_dsw(i,j)
      ja = jpos_dsw(i,j)
      rx = (lon_rho(i,j)-londsw(ia))/(londsw(ia+1)-londsw(ia))
      ry = (lat_rho(i,j)-latdsw(ja))/(latdsw(ja+1)-latdsw(ja))
      rxm = 1.0-rx
      rym = 1.0-ry     
      swrad(i,j)    = rxm*rym*DSW(ia,ja)          +   &
           rx*rym*DSW(ia+1,ja)                    +   &
           rx*ry*DSW(ia+1,ja+1)                   +   &
           rxm*ry*DSW(ia,ja+1)
      ia = ipos_dlw(i,j)
      ja = jpos_dlw(i,j)
      rx = (lon_rho(i,j)-londlw(ia))/(londlw(ia+1)-londlw(ia))
      ry = (lat_rho(i,j)-latdlw(ja))/(latdlw(ja+1)-latdlw(ja))
      rxm = 1.0-rx
      rym = 1.0-ry     
      lwrad(i,j)    = rxm*rym*DLW(ia,ja)          +   &
           rx*rym*DLW(ia+1,ja)                    +   &
           rx*ry*DLW(ia+1,ja+1)                   +   &
           rxm*ry*DLW(ia,ja+1)
   ENDDO
ENDDO

! Rotate wind components from E/W and N/S to Xi and Eta directions

Uwind = scr1_out*cos(angle) + scr2_out*sin(angle)
Vwind = scr2_out*cos(angle) - scr1_out*sin(angle)

else
   Uwind = U10(i1_wind:i2_wind,j1_wind:j2_wind)
   Vwind = V10(i1_wind:i2_wind,j1_wind:j2_wind)               
   Tair  = T2M(i1_t2m:i2_t2m,j1_t2m:j2_t2m)
   Qair  = Q2M(i1_q2m:i2_q2m,j1_q2m:j2_q2m)
   cloudair = CLOUD(i1_cloud:i2_cloud,j1_cloud:j2_cloud)/100.
   rainair  = RAIN(i1_rain:i2_rain,j1_rain:j2_rain)
   swrad = DSW(i1_dsw:i2_dsw,j1_dsw:j2_dsw)
   lwrad = DLW(i1_dlw:i2_dlw,j1_dlw:j2_dlw)     

   Pair = MSL(i1_msl:i2_msl,j1_msl:j2_msl)
! Letting ROMS do the interpolation now, and no need to have MSL on same grid as rest as
! in different files.
!   Pair  = -99999999999.99
! MSL data is only available on 0.5degx0.5deg -> interpolation on the other grid needed
!   Pairtemp  = MSL(i1_msl:i2_msl,j1_msl:j2_msl) 
!   DO I= 1,Lp_dim
!      DO J= 1,Mp_dim
!         ia = ipos_msltemp(i)
!         ja = jpos_msltemp(j)
!         rx = (lonwtemp(i)-lonmsltemp(ia))/(lonmsltemp(ia+1)-lonmsltemp(ia))
!         ry = (latwtemp(j)-latmsltemp(ja))/(latmsltemp(ja+1)-latmsltemp(ja))
!         rxm = 1.0-rx
!         rym = 1.0-ry     
!         Pair(i,j)= rxm*rym*Pairtemp(ia,ja)          +   &
!              rx*rym*Pairtemp(ia+1,ja)               +   &
!              rx*ry*Pairtemp(ia+1,ja+1)              +   &
!              rxm*ry*Pairtemp(ia,ja+1)      
!      ENDDO
!   ENDDO
end if

write(*,*)'Output file record is: ', irec

statuso = nf90_put_var(ncidowind,TimeOutVarId,tday_wind,start=(/ irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_var(ncidoq2m,TimeOutVarId,tday_q2m,start=(/ irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_var(ncidot2m,TimeOutVarId,tday_t2m,start=(/ irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_var(ncidocloud,TimeOutVarId,tday_cloud,start=(/ irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_var(ncidorain,TimeOutVarId,tday_rain,start=(/ irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_var(ncidodsw,TimeOutVarId,tday_dsw,start=(/ irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_var(ncidodlw,TimeOutVarId,tday_dlw,start=(/ irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_var(ncidomsl,TimeOutVarId,tday_msl,start=(/ irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)



statuso = nf90_put_var(ncidowind,UwindVarId,Uwind,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncidowind,VwindVarId,Vwind,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncidomsl,PairVarId,Pair,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncidot2m,TairVarId,Tair,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncidoq2m,QairVarId,Qair,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncidocloud,cloudairVarId,cloudair,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

status_roms = nf90_put_var(ncidorain,rainairVarId,rainair,start=(/1, 1, irec/))
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

status_roms = nf90_put_var(ncidodsw,swradVarId,swrad,start=(/1, 1,  irec/))
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

status_roms = nf90_put_var(ncidodlw,lwradVarId,lwrad,start=(/1, 1, irec/))
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

status_wind = nf90_close(ncid_wind)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_msl = nf90_close(ncid_msl)
if(status_msl /= nf90_NoErr) call handle_err(status_msl)
status_t2m = nf90_close(ncid_t2m)
if(status_t2m /= nf90_NoErr) call handle_err(status_t2m)
status_q2m = nf90_close(ncid_q2m)
if(status_q2m /= nf90_NoErr) call handle_err(status_q2m)
status_rain = nf90_close(ncid_rain)
if(status_rain /= nf90_NoErr) call handle_err(status_rain)
status_cloud = nf90_close(ncid_cloud)
if(status_cloud /= nf90_NoErr) call handle_err(status_cloud)
status_dlw = nf90_close(ncid_dlw)
if(status_dlw /= nf90_NoErr) call handle_err(status_dlw)
status_dsw = nf90_close(ncid_dsw)
if(status_dsw /= nf90_NoErr) call handle_err(status_dsw)

ENDDO

statuso = nf90_close(ncidowind)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_close(ncidomsl)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_close(ncidot2m)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_close(ncidoq2m)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_close(ncidocloud)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_close(ncidorain)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_close(ncidodsw)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_close(ncidodlw)
if(statuso /= nf90_NoErr) call handle_err(statuso)


END PROGRAM cfsr2roms














