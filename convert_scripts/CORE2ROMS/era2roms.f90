! Program to generate ROMS forcing files from ERAI netCDF files
! produced using cvt_gribfiles.bsh

PROGRAM era2roms
use netcdf
implicit none

INTEGER i, j, k, l, imon, itime, ntimes, ia, ja, IFL
INTEGER start_AN, STOP_AN, irec
INTEGER iyear, iday, imonth,ihour, nb_records, indice
INTEGER greg(6)
INTEGER nlons, nlats ,LonDimID, LatDimID, TimeDimId
INTEGER XiDimID, EtaDimID, Lp, Mp, Lp_dim, Mp_dim
INTEGER MaskVarId, LonVarId, LatVarId, DiVarId, DjVarId
INTEGER LonRhoVarId, LatRhoVarId, AngleVarId
INTEGER U10VarId, V10VarId, TimeANVarId, TimeFCVarId
INTEGER MSLVarId, T2MVarId, D2MVarId, TCCVarId, CPVarId
INTEGER LSPVarId, SSRDVarId, STRDVarId
INTEGER year_st, month_st, day_st
INTEGER year_end, month_end, day_end
INTEGER year_ref, month_ref, day_ref

INTEGER XiRhoDimId, EtaRhoDimId, TimeOutDimID
INTEGER TimeOutVarId, UwindVarId,VwindVarId,PairVarId
INTEGER TairVarId, QairVarId,cloudVarId,rainVarId, swradVarId, lwradVarId

INTEGER statuso, ncido
INTEGER status_toto, nctoto
INTEGER status_era, ncid_era
INTEGER status_roms, ncid_roms
INTEGER status_AN, ncid_AN
INTEGER status_FC, ncid_FC

INTEGER ndays, n6hr, n12hr

INTEGER i1, i2, j1, j2
INTEGER jdref, jdref_era
INTEGER yes_no

REAL  Di, Dj, rx, rxm, ry, rym
REAL  rimin, rimax, rjmin, rjmax
REAL, DIMENSION(:,:), allocatable :: ipos
REAL, DIMENSION(:,:), allocatable :: jpos
REAL  delt

REAL U10_scale_factor, U10_add_offset
REAL V10_scale_factor, V10_add_offset
REAL MSL_scale_factor, MSL_add_offset
REAL T2M_scale_factor, T2M_add_offset
REAL D2M_scale_factor, D2M_add_offset
REAL TCC_scale_factor, TCC_add_offset

REAL CP_scale_factor, CP_add_offset
REAL LSP_scale_factor, LSP_add_offset
REAL SSRD_scale_factor, SSRD_add_offset
REAL STRD_scale_factor, STRD_add_offset

REAL, ALLOCATABLE, DIMENSION(:) :: lons
REAL, ALLOCATABLE, DIMENSION(:) :: IX
REAL, ALLOCATABLE, DIMENSION(:) :: lats
REAL, ALLOCATABLE, DIMENSION(:) :: IY
REAL, ALLOCATABLE, DIMENSION(:,:) :: emask
REAL, ALLOCATABLE, DIMENSION(:,:) :: work

REAL, ALLOCATABLE, DIMENSION(:,:) :: scr1_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr2_out

REAL, ALLOCATABLE, DIMENSION(:,:) :: lon_rho
REAL, ALLOCATABLE, DIMENSION(:,:) :: lat_rho
REAL, ALLOCATABLE, DIMENSION(:,:) :: angle

REAL, ALLOCATABLE, DIMENSION(:,:) :: U10
REAL, ALLOCATABLE, DIMENSION(:,:) :: V10
REAL, ALLOCATABLE, DIMENSION(:,:) :: MSL
REAL, ALLOCATABLE, DIMENSION(:,:) :: T2M
REAL, ALLOCATABLE, DIMENSION(:,:) :: D2M
REAL, ALLOCATABLE, DIMENSION(:,:) :: Q2M
REAL, ALLOCATABLE, DIMENSION(:,:) :: evapres
REAL, ALLOCATABLE, DIMENSION(:,:) :: TCC
REAL, ALLOCATABLE, DIMENSION(:,:) :: TP
REAL, ALLOCATABLE, DIMENSION(:,:) :: CP
REAL, ALLOCATABLE, DIMENSION(:,:) :: LSP
REAL, ALLOCATABLE, DIMENSION(:,:) :: SSRD
REAL, ALLOCATABLE, DIMENSION(:,:) :: STRD

REAL, ALLOCATABLE, DIMENSION(:,:) :: Uwind
REAL, ALLOCATABLE, DIMENSION(:,:) :: Vwind
REAL, ALLOCATABLE, DIMENSION(:,:) :: Pair
REAL, ALLOCATABLE, DIMENSION(:,:) :: Tair
REAL, ALLOCATABLE, DIMENSION(:,:) :: Qair
REAL, ALLOCATABLE, DIMENSION(:,:) :: cloud
REAL, ALLOCATABLE, DIMENSION(:,:) :: rain
REAL, ALLOCATABLE, DIMENSION(:,:) :: swrad
REAL, ALLOCATABLE, DIMENSION(:,:) :: lwrad
REAL tday_AN, tday_FC, datetime

REAL, PARAMETER    :: undef = 2.E+35            ! Undefined land value
REAL pi, DTOR

DOUBLE PRECISION date

CHARACTER(len=4) year_char
CHARACTER(len=2) month_char, day_char
CHARACTER(len=50)  ryear, num_temp
CHARACTER(len=160) path_era, maskfile, gridfile, fileout, time_units_att
CHARACTER(len=80)  AN_file
CHARACTER(len=80)  FC_file

interface
     FUNCTION jd(yy,mm,dd)
     implicit none
     INTEGER jd,yy,mm,dd
     end function jd
end interface

pi = ATAN(1.)*4.
DTOR = pi/180.

OPEN(unit=10,file='era2roms.in')
READ(10,*)
READ(10,'(a)') path_era
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
! Get reference Julian day (1900-01-01 00:00:00Z) for input
jdref_era = jd(1900,1,1)


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
ALLOCATE(ipos(Lp,Mp))
ALLOCATE(jpos(Lp,Mp))
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

! Get dimensions of ERAI files
maskfile = TRIM(path_era) // 'landmask.nc'
status_era = nf90_open(TRIM(maskfile),nf90_nowrite,ncid_era)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_inq_dimid(ncid_era, "longitude", LonDimID)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_inq_dimid(ncid_era, "latitude", LatDimID)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_Inquire_Dimension(ncid_era,LonDimID,len = nlons)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_Inquire_Dimension(ncid_era,LatDimID,len = nlats)
if(status_era /= nf90_NoErr) call handle_err(status_era)
ALLOCATE(lons(nlons))
ALLOCATE(lats(nlats))
ALLOCATE(IX(nlons))
ALLOCATE(IY(nlats))
ALLOCATE(emask(nlons,nlats))
ALLOCATE(work(nlons,nlats))
ALLOCATE(U10(nlons,nlats))
ALLOCATE(V10(nlons,nlats))
ALLOCATE(MSL(nlons,nlats))
ALLOCATE(T2M(nlons,nlats))
ALLOCATE(D2M(nlons,nlats))
ALLOCATE(Q2M(nlons,nlats))
ALLOCATE(evapres(nlons,nlats))
ALLOCATE(TCC(nlons,nlats))
ALLOCATE(TP(nlons,nlats))
ALLOCATE(CP(nlons,nlats))
ALLOCATE(LSP(nlons,nlats))
ALLOCATE(SSRD(nlons,nlats))
ALLOCATE(STRD(nlons,nlats))
! Get lons and lats
status_era = nf90_inq_varid(ncid_era, "longitude",LonVarId)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_inq_varid(ncid_era, "latitude",LatVarId)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_get_var(ncid_era,LonVarId,lons)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_get_var(ncid_era,LatVarId,lats)
if(status_era /= nf90_NoErr) call handle_err(status_era)
Di = 0.7
Dj = 0.7
! Get land-sea mask
status_era = nf90_inq_varid(ncid_era, "lsm",MaskVarId)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_get_var(ncid_era,MaskVarId,emask)
if(status_era /= nf90_NoErr) call handle_err(status_era)
emask= emask*1.52597204419215e-05+0.499984740279558
write(*,*) 'nlons,nlats = ',nlons,nlats
status_era = nf90_close(ncid_era)
if(status_era /= nf90_NoErr) call handle_err(status_era)

! Get i,j positions of ROMS rho-points on ERAI grid
WHERE(lons > 180.) lons = lons-360.
DO I=1,nlons
IX(I)= I
ENDDO
DO I=1,nlats
IY(I)= I
ENDDO

CALL SSORT (lons, IX, nlons, 2)
CALL SSORT (lats, IY, nlats, 2)

DO j=1,Mp
DO i=1,Lp
DO k=1,nlons-1
   if((lon_rho(i,j).ge.lons(k)).AND.(lon_rho(i,j).lt.lons(k+1)))then
   ipos(i,j) = k
   endif
ENDDO
DO l=1,nlats-1
   if((lat_rho(i,j).ge.lats(l)).AND.(lat_rho(i,j).lt.lats(l+1)))then
   jpos(i,j) = l
   endif
ENDDO
ENDDO
ENDDO
! Find subarea in ERAI grid containing model grid
rimin = 1.e+23
rimax =-1.e+23
rjmin = 1.e+23
rjmax =-1.e+23
DO j=1,Mp
DO i=1,Lp
   rimin = MIN(rimin,ipos(i,j))
   rimax = MAX(rimax,ipos(i,j))
   rjmin = MIN(rjmin,jpos(i,j))
   rjmax = MAX(rjmax,jpos(i,j))
ENDDO
ENDDO
i1 = FLOOR(rimin)
i2 = CEILING(rimax)
j1 = FLOOR(rjmin)
j2 = CEILING(rjmax)

write(*,*) 'min indice to extract in longitude: ',i1
write(*,*) 'max indice to extract in longitude: ',i2
write(*,*) 'min indice to extract in latitude: ',j1
write(*,*) 'max indice to extract in latitude: ',j2
write(*,*) 'Maximum dimension in longitude :',nlons
write(*,*) 'Maximum dimension in latitude :',nlats

IF(i1<1.OR.i2>nlons.OR.j1<1.OR.j2>nlats) THEN
   STOP 'subdomain outside ERAI domain'
ENDIF

if(yes_no.eq.0)then
i1 = i1 -5
if(i1.lt.1) i1=1
i2 =i2 +5
if(i2.gt.nlons) i2= nlons
j1= j1-5
if(j1.lt.1) j1=1
j2 =j2 +5
if(j2.gt.nlats) j2= nlats
LP_dim= i2-i1+1
MP_dim= j2-j1+1
else
LP_dim=LP
MP_dim=MP
end if


! -------------------------------------------------------------------
! Open output file

statuso = nf90_create(TRIM(fileout),nf90_clobber,ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.1)then
statuso = nf90_def_dim(ncido,"xi_rho",Lp_dim,XiRhoDimID)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_dim(ncido,"eta_rho",Mp_dim,EtaRhoDimID)
if(statuso /= nf90_NoErr) call handle_err(statuso)
else
statuso = nf90_def_dim(ncido,"longitude",Lp_dim,XiRhoDimID)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_dim(ncido,"latitude",Mp_dim,EtaRhoDimID)
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
statuso = nf90_def_dim(ncido,"time",nf90_unlimited,TimeOutDimID)
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncido,nf90_global,"title", &
                  "Atmospheric forcing file - ERAI grid")
else
statuso = nf90_put_att(ncido,nf90_global,"title", &
                  "Atmospheric forcing file - ROMS grid")
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,nf90_global,"type", &
                  "ROMS forcing file")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,nf90_global,"history", &
                  "Produced with era2roms.f90")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
! Latitude variable
statuso = nf90_def_var(ncido,"latitude",nf90_float, &
                 (/EtaRhoDimID/), LatVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, LatVarId,"long_Name", &
                 "latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, LatVarId,"units", &
                 "degree_east")
! Longitude variable
statuso = nf90_def_var(ncido,"longitude",nf90_float, &
                 (/ XiRhoDimID /), LonVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, LatVarId,"long_Name", &
                 "longitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, LonVarId,"units", &
                 "degree_north")
! mask
statuso = nf90_def_var(ncido,"landmask",nf90_int, &
                 (/ XiRhoDimID, EtaRhoDimID/), MaskVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, MaskVarId,"long_Name", &
                 "Land mask")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, MaskVarId,"units", &
                 "land=0/sea=1")
statuso = nf90_put_att(ncido, MaskVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if 
! Time variable
statuso = nf90_def_var(ncido,"time",nf90_float, &
                 TimeOutDimID, TimeOutVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, TimeOutVarId,"long_Name", &
                 "Time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, TimeOutVarId,"units", &
                 TRIM(time_units_att))
if(statuso /= nf90_NoErr) call handle_err(statuso)
! UWIND variable
statuso = nf90_def_var(ncido,"Uwind",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), UwindVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, UwindVarId,"long_Name", &
                 "Xi-component of wind")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, UwindVarId,"units", &
                 "meter second-1")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, UwindVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncido, UwindVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
! VWIND variable
statuso = nf90_def_var(ncido,"Vwind",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), VwindVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, VwindVarId,"long_Name", &
                 "Eta-component of wind")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, VwindVarId,"units", &
                 "meter second-1")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, VwindVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncido, VwindVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
! PAIR variable
statuso = nf90_def_var(ncido,"Pair",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), PairVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, PairVarId,"long_Name", &
                 "Atmospheric pressure")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, PairVarId,"units", &
                 "Pa")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, PairVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncido, PairVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
! TAIR variable
statuso = nf90_def_var(ncido,"Tair",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), TairVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, TairVarId,"long_Name", &
                 "Air temperature at 2m")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, TairVarId,"units", &
                 "degree C")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, TairVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncido, TairVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
! QAIR variable
statuso = nf90_def_var(ncido,"Qair",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), QairVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, QairVarId,"long_Name", &
                 "Specific humidity at 2m")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, QairVarId,"units", &
                 "kg kg-1")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, QairVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncido, QairVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
! CLOUD variable
statuso = nf90_def_var(ncido,"cloud",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), cloudVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, cloudVarId,"long_Name", &
                 "Fraction cloud cover")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, cloudVarId,"units", &
                 " ")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, cloudVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncido, cloudVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
! RAIN variable
statuso = nf90_def_var(ncido,"rain",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), rainVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, rainVarId,"long_Name", &
                 "Total precipiation")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, rainVarId,"units", &
                 "kg meter-2 second-1")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, rainVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncido, rainVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
! SWRAD variable
statuso = nf90_def_var(ncido,"swrad",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), swradVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, swradVarId,"long_Name", &
                 "Downward shortwave radiation")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, swradVarId,"units", &
                 "Watt meter-2")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, swradVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncido, swradVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if
! LWRAD variable
statuso = nf90_def_var(ncido,"lwrad_down",nf90_float, &
                 (/ XiRhoDimID, EtaRhoDimID, TimeOutDimID /), lwradVarId)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, lwradVarId,"long_Name", &
                 "Downward longwave radiation")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, lwradVarId,"units", &
                 "Watt meter-2")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, lwradVarId,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
if(yes_no.eq.0)then
statuso = nf90_put_att(ncido, lwradVarId,"coordinates", &
                 "longitude latitude")
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if

statuso = nf90_enddef(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

WRITE(*,*) 'Output file ', TRIM(fileout), ' is created.'

ALLOCATE(Uwind(Lp_dim,Mp_dim))
ALLOCATE(Vwind(Lp_dim,Mp_dim))
ALLOCATE(Pair(Lp_dim,Mp_dim)) 
ALLOCATE(Tair(Lp_dim,Mp_dim))
ALLOCATE(Qair(Lp_dim,Mp_dim))
ALLOCATE(cloud(Lp_dim,Mp_dim))
ALLOCATE(rain(Lp_dim,Mp_dim))
ALLOCATE(swrad(Lp_dim,Mp_dim))
ALLOCATE(lwrad(Lp_dim,Mp_dim))

!Open output file
statuso = nf90_open(trim(fileout),nf90_write,ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

if(yes_no.eq.0)then
! Write latitude and longitude
statuso = nf90_put_var(ncido,LonVarId,lons(i1:i2))
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_var(ncido,LatVarId,lats(j1:j2))
if(statuso /= nf90_NoErr) call handle_err(statuso)
work=0.
DO I=1,nlons
   DO J=1,nlats
      work(I,J)=emask(IX(I),IY(J))
   ENDDO
ENDDO
emask= INT(work)
statuso = nf90_put_var(ncido,MaskVarId,emask(i1:i2,j1:j2))
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if

nb_records = (jd(year_end,month_end,day_end)-jd(year_st,month_st,day_st)+1)*4

DO irec=1,nb_records

date= REAL(jd(year_st,month_st,day_st))+REAL(irec-1)*6./24.
greg=0
CALL GR(greg,date)
iyear= greg(1)
imonth= greg(2)
iday= greg(3)
ihour= greg(4)
write(*,'(a24,i2,a1,i2,a1,i4,a4,i2,a1)')'Writting data for date: ',greg(3),'/',greg(2),'/',greg(1), ' at ',greg(4),'h'
! First: treatment of the 6hourly files ...
if(ihour.eq.0)then
indice = 1
elseif(ihour.eq.6)then
indice = 2
elseif(ihour.eq.12)then
indice = 3
elseif(ihour.eq.18)then
indice = 4
end if
write(ryear,*) iyear
ryear = adjustl(ryear)
itime= (jd(iyear,imonth,iday)-jd(iyear,1,1))*4 + indice
AN_file = TRIM(path_era)//TRIM(ryear)// '0101AN.nc'
status_AN = nf90_open(TRIM(AN_file),nf90_nowrite,ncid_AN)
if(status_AN /= nf90_NoErr) call handle_err(status_AN)
! Get AN time var id
status_AN = nf90_inq_varid(ncid_AN, "time",TimeANVarId)
if(status_AN /= nf90_NoErr) call handle_err(status_AN)
! Get U10 var id
status_AN = nf90_inq_varid(ncid_AN, "10u",U10VarId)
if(status_AN /= nf90_NoErr) call handle_err(status_AN)
status_AN = nf90_get_att(ncid_AN, U10VarID, "scale_factor", U10_scale_factor)
if (status_AN /= nf90_noerr) call handle_err(status_AN)
status_AN = nf90_get_att(ncid_AN, U10VarID, "add_offset", U10_add_offset)
if (status_AN /= nf90_noerr) call handle_err(status_AN)
! Get V10 var id
status_era = nf90_inq_varid(ncid_AN, "10v",V10VarId)
if(status_AN /= nf90_NoErr) call handle_err(status_AN)
status_AN = nf90_get_att(ncid_AN, V10VarID, "scale_factor", V10_scale_factor)
if (status_AN /= nf90_noerr) call handle_err(status_AN)
status_AN = nf90_get_att(ncid_AN, V10VarID, "add_offset", V10_add_offset)
if (status_AN /= nf90_noerr) call handle_err(status_AN)
! Get MSL var id
status_AN = nf90_inq_varid(ncid_AN, "msl",MSLVarId)
if(status_AN /= nf90_NoErr) call handle_err(status_AN)
status_AN = nf90_get_att(ncid_AN, MSLVarID, "scale_factor", MSL_scale_factor)
if (status_AN /= nf90_noerr) call handle_err(status_AN)
status_AN = nf90_get_att(ncid_AN, MSLVarID, "add_offset", MSL_add_offset)
if (status_AN /= nf90_noerr) call handle_err(status_AN)
! Get T2M var id
status_AN = nf90_inq_varid(ncid_AN, "2t",T2MVarId)
if(status_AN /= nf90_NoErr) call handle_err(status_AN)
status_AN = nf90_get_att(ncid_AN, T2MVarID, "scale_factor", T2M_scale_factor)
if (status_AN /= nf90_noerr) call handle_err(status_AN)
status_AN = nf90_get_att(ncid_AN, T2MVarID, "add_offset", T2M_add_offset)
if (status_AN /= nf90_noerr) call handle_err(status_AN)
! Get D2M var id
status_AN = nf90_inq_varid(ncid_AN, "2d",D2MVarId)
if(status_AN /= nf90_NoErr) call handle_err(status_AN)
status_AN = nf90_get_att(ncid_AN, D2MVarID, "scale_factor", D2M_scale_factor)
if (status_AN /= nf90_noerr) call handle_err(status_AN)
status_AN = nf90_get_att(ncid_AN, D2MVarID, "add_offset", D2M_add_offset)
if (status_AN /= nf90_noerr) call handle_err(status_AN)
! Get TCC var id
status_AN = nf90_inq_varid(ncid_AN, "tcc",TCCVarId)
if(status_AN /= nf90_NoErr) call handle_err(status_AN)
status_AN = nf90_get_att(ncid_AN, TCCVarID, "scale_factor", TCC_scale_factor)
if (status_AN /= nf90_noerr) call handle_err(status_AN)
status_AN = nf90_get_att(ncid_AN, TCCVarID, "add_offset", TCC_add_offset)
if (status_AN /= nf90_noerr) call handle_err(status_AN)
! Read in Time   
status_era = nf90_get_var(ncid_AN,TimeANVarId,datetime,start = (/itime/))
if(status_era /= nf90_NoErr) call handle_err(status_era)
tday_AN= REAL(jdref_era) + REAL(datetime)/24. - REAL(jdref)
tday_AN= tday_AN*24.*3600.
! Read in U10
status_era = nf90_get_var(ncid_AN,U10VarId,U10,start=(/1, 1,itime /))
if(status_era /= nf90_NoErr) call handle_err(status_era)
U10 = U10*U10_scale_factor + U10_add_offset
! Blank out Antarctica
WHERE(ABS(U10)>1.e+7) U10 = 0.
! Read in V10
status_era = nf90_get_var(ncid_AN,V10VarId,V10,start=(/1, 1,itime /))
if(status_era /= nf90_NoErr) call handle_err(status_era)
V10 = V10*V10_scale_factor + V10_add_offset
! Blank out Antarctica
WHERE(ABS(V10)>1.e+7) V10 = 0.
! Read in MSL values
status_era = nf90_get_var(ncid_AN,MSLVarId,MSL,start=(/1, 1,itime /))
if(status_era /= nf90_NoErr) call handle_err(status_era)
MSL = MSL*MSL_scale_factor + MSL_add_offset
! Blank out Antarctica
WHERE(ABS(MSL)>1.e+7) MSL = 1.e+5
! Read in T2M values
status_era = nf90_get_var(ncid_AN,T2MVarId,T2M,start=(/1, 1,itime /))
if(status_era /= nf90_NoErr) call handle_err(status_era)
T2M = T2M*T2M_scale_factor + T2M_add_offset
! Convert from Kelvin to Celsius
T2M = T2M - 273.16
! Blank out Antarctica
WHERE(ABS(T2M)>1.e+7) T2M = 0.
! Read in D2M values
status_era = nf90_get_var(ncid_AN,D2MVarId,D2M,start=(/1, 1,itime /))
if(status_era /= nf90_NoErr) call handle_err(status_era)
D2M = D2M*D2M_scale_factor + D2M_add_offset
! Convert from Kelvin to Celsius
D2M = D2M - 273.16
! Compute specific humidity on ECMWF points
IFL = 1 !Input is dew point temperature
CALL spec_hum(D2M, T2M, MSL, evapres, Q2M,nlons,nlats,IFL)
! Blank out Antarctica
WHERE(ABS(Q2M)>1.e+7) Q2M = 1.e-5
! Read in TCC values
status_era = nf90_get_var(ncid_AN,TCCVarId,TCC,start=(/1, 1,itime /))
if(status_era /= nf90_NoErr) call handle_err(status_era)
TCC = TCC*TCC_scale_factor + TCC_add_offset
! Blank out Antarctica
WHERE(ABS(TCC)>1.e+7) TCC = 0.
status_era = nf90_close(ncid_AN)
if(status_era /= nf90_NoErr) call handle_err(status_era)

! Second: treatment of the 12hourly files ...
delt = 12.*3600.
if(ihour.eq.0)then
indice = 1
elseif(ihour.eq.6)then
indice = 1
elseif(ihour.eq.12)then
indice = 2
elseif(ihour.eq.18)then
indice = 2
end if
write(ryear,*) iyear
ryear = adjustl(ryear)
itime= (jd(iyear,imonth,iday)-jd(iyear,1,1))*2 + indice
FC_file = TRIM(path_era)//TRIM(ryear)// '0101FC.nc'
status_FC = nf90_open(TRIM(FC_file),nf90_nowrite,ncid_FC)
if(status_FC /= nf90_NoErr) call handle_err(status_FC)
! Get FC time var id
status_FC = nf90_inq_varid(ncid_FC, "time",TimeFCVarId)
if(status_FC /= nf90_NoErr) call handle_err(status_FC)
! Get CP var id
status_FC = nf90_inq_varid(ncid_FC, "cp",CPVarId)
if(status_FC /= nf90_NoErr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, CPVarID, "scale_factor", CP_scale_factor)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, CPVarID, "add_offset", CP_add_offset)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
! Get LSP var id
status_FC = nf90_inq_varid(ncid_FC, "lsp",LSPVarId)
if(status_FC /= nf90_NoErr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, LSPVarID, "scale_factor", LSP_scale_factor)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, LSPVarID, "add_offset", LSP_add_offset)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
! Get SSRD var id
status_FC = nf90_inq_varid(ncid_FC, "ssrd",SSRDVarId)
if(status_FC /= nf90_NoErr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, SSRDVarID, "scale_factor", SSRD_scale_factor)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, SSRDVarID, "add_offset", SSRD_add_offset)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
! Get STRD var id
status_FC = nf90_inq_varid(ncid_FC, "strd",STRDVarId)
if(status_FC /= nf90_NoErr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, STRDVarID, "scale_factor", STRD_scale_factor)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, STRDVarID, "add_offset", STRD_add_offset)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
! Read in time value
status_era = nf90_get_var(ncid_FC,TimeFCVarId,datetime,start = (/itime/))
if(status_era /= nf90_NoErr) call handle_err(status_era)
tday_FC = REAL(jdref_era) + REAL(datetime)/24. - REAL(jdref) - 0.25
tday_FC= tday_FC*24.*3600.
! Read in CP values
status_ERA = nf90_get_var(ncid_FC,CPVarId,work,start=(/1, 1,itime /))
if(status_ERA /= nf90_NoErr) call handle_err(status_ERA)
CP = work*CP_scale_factor + CP_add_offset
! Read in LSP values
status_ERA = nf90_get_var(ncid_FC,LSPVarId,work,start=(/1, 1,itime /))
if(status_ERA /= nf90_NoErr) call handle_err(status_ERA)
LSP = work*LSP_scale_factor + LSP_add_offset
TP = CP + LSP
! Read in SSRD values
status_ERA = nf90_get_var(ncid_FC,SSRDVarId,work,start=(/1, 1,itime /))
if(status_ERA /= nf90_NoErr) call handle_err(status_ERA)
SSRD = work*SSRD_scale_factor + SSRD_add_offset
! Read in STRD values
status_ERA = nf90_get_var(ncid_FC,STRDVarId,work,start=(/1, 1,itime /))
if(status_ERA /= nf90_NoErr) call handle_err(status_ERA)
STRD = work*STRD_scale_factor + STRD_add_offset
status_era = nf90_close(ncid_FC)
if(status_era /= nf90_NoErr) call handle_err(status_era)

if((ihour.eq.0).or.(ihour.eq.12))then
date= REAL(jd(year_st,month_st,day_st))+REAL(irec-2)*6./24
greg=0
CALL GR(greg,date)
iyear= greg(1)
imonth= greg(2)
iday= greg(3)
ihour= greg(4)
if(ihour.eq.18)then
indice = 2
elseif(ihour.eq.6)then
indice = 1
end if
write(ryear,*) iyear
ryear = adjustl(ryear)
itime= (jd(iyear,imonth,iday)-jd(iyear,1,1))*2 + indice
FC_file = TRIM(path_era)//TRIM(ryear)// '0101FC.nc'
status_FC = nf90_open(TRIM(FC_file),nf90_nowrite,ncid_FC)
if(status_FC /= nf90_NoErr) call handle_err(status_FC)
! Get FC time var id
status_FC = nf90_inq_varid(ncid_FC, "time",TimeFCVarId)
if(status_FC /= nf90_NoErr) call handle_err(status_FC)
! Get CP var id
status_FC = nf90_inq_varid(ncid_FC, "cp",CPVarId)
if(status_FC /= nf90_NoErr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, CPVarID, "scale_factor", CP_scale_factor)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, CPVarID, "add_offset", CP_add_offset)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
! Get LSP var id
status_FC = nf90_inq_varid(ncid_FC, "lsp",LSPVarId)
if(status_FC /= nf90_NoErr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, LSPVarID, "scale_factor", LSP_scale_factor)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, LSPVarID, "add_offset", LSP_add_offset)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
! Get SSRD var id
status_FC = nf90_inq_varid(ncid_FC, "ssrd",SSRDVarId)
if(status_FC /= nf90_NoErr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, SSRDVarID, "scale_factor", SSRD_scale_factor)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, SSRDVarID, "add_offset", SSRD_add_offset)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
! Get STRD var id
status_FC = nf90_inq_varid(ncid_FC, "strd",STRDVarId)
if(status_FC /= nf90_NoErr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, STRDVarID, "scale_factor", STRD_scale_factor)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
status_FC = nf90_get_att(ncid_FC, STRDVarID, "add_offset", STRD_add_offset)
if (status_FC /= nf90_noerr) call handle_err(status_FC)
! Read in time value
status_era = nf90_get_var(ncid_FC,TimeFCVarId,datetime,start = (/itime/))
if(status_era /= nf90_NoErr) call handle_err(status_era)
datetime = REAL(jdref_era) + REAL(datetime)/24. - REAL(jdref) - 0.25
tday_FC= 0.5*(tday_FC+datetime*24.*3600.)
! Read in CP values
work= 0.
status_ERA = nf90_get_var(ncid_FC,CPVarId,work,start=(/1, 1,itime /))
if(status_ERA /= nf90_NoErr) call handle_err(status_ERA)
work = work*CP_scale_factor + CP_add_offset
CP=0.5*(CP+work)
! Read in LSP values
work = 0.
status_ERA = nf90_get_var(ncid_FC,LSPVarId,work,start=(/1, 1,itime /))
if(status_ERA /= nf90_NoErr) call handle_err(status_ERA)
work = work*LSP_scale_factor + LSP_add_offset
LSP= 0.5*(LSP+work)
TP = CP + LSP
! Read in SSRD values
work= 0.
status_ERA = nf90_get_var(ncid_FC,SSRDVarId,work,start=(/1, 1,itime /))
if(status_ERA /= nf90_NoErr) call handle_err(status_ERA)
work = work*SSRD_scale_factor + SSRD_add_offset
SSRD=0.5*(work+SSRD)
! Read in STRD values
work=0.
status_ERA = nf90_get_var(ncid_FC,STRDVarId,work,start=(/1, 1,itime /))
if(status_ERA /= nf90_NoErr) call handle_err(status_ERA)
work = work*STRD_scale_factor + STRD_add_offset
STRD=0.5*(work+STRD)
status_era = nf90_close(ncid_FC)
if(status_era /= nf90_NoErr) call handle_err(status_era)
end if

work=0.
DO I=1,nlons
   DO J=1,nlats
      work(I,J)=U10(IX(I),IY(J))
   ENDDO
ENDDO
U10= work
work=0.
DO I=1,nlons
   DO J=1,nlats
      work(I,J)=V10(IX(I),IY(J))
   ENDDO
ENDDO
V10= work
work=0.
DO I=1,nlons
   DO J=1,nlats
      work(I,J)=MSL(IX(I),IY(J))
   ENDDO
ENDDO
MSL= work
work=0.
DO I=1,nlons
   DO J=1,nlats
      work(I,J)=T2M(IX(I),IY(J))
   ENDDO
ENDDO
T2M= work
work=0.
DO I=1,nlons
   DO J=1,nlats
      work(I,J)=Q2M(IX(I),IY(J))
   ENDDO
ENDDO
Q2M= work
work=0.
DO I=1,nlons
   DO J=1,nlats
      work(I,J)=TCC(IX(I),IY(J))
   ENDDO
ENDDO
TCC= work
work=0.
DO I=1,nlons
   DO J=1,nlats
      work(I,J)=TP(IX(I),IY(J))
   ENDDO
ENDDO
TP= work
work=0.
DO I=1,nlons
   DO J=1,nlats
      work(I,J)=SSRD(IX(I),IY(J))
   ENDDO
ENDDO
SSRD= work
work=0.
DO I=1,nlons
   DO J=1,nlats
      work(I,J)=STRD(IX(I),IY(J))
   ENDDO
ENDDO
STRD= work

if(yes_no.eq.1)then
! Horizontal interpolation
scr1_out = 0.
scr2_out = 0.
Pair=0.
Tair=0.
Qair=0.
Cloud=0.
rain= 0.
swrad= 0.
lwrad= 0.
DO j=1,Mp
   DO i=1,Lp
      
      ia = ipos(i,j)
      ja = jpos(i,j)
      rx = (lon_rho(i,j)-lons(ia))/(lons(ia+1)-lons(ia))
      ry = (lat_rho(i,j)-lats(ja))/(lats(ja+1)-lats(ja))
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
      Pair(i,j)     = rxm*rym*MSL(ia,ja)           +   &
           rx*rym*MSL(ia+1,ja)                     +   &
           rx*ry*MSL(ia+1,ja+1)                    +   &
           rxm*ry*MSL(ia,ja+1)      
      Tair(i,j)     = rxm*rym*T2M(ia,ja)           +   &
           rx*rym*T2M(ia+1,ja)                     +   &
           rx*ry*T2M(ia+1,ja+1)                    +   &
           rxm*ry*T2M(ia,ja+1)      
      Qair(i,j)     = rxm*rym*Q2M(ia,ja)           +   &
           rx*rym*Q2M(ia+1,ja)                     +   &
           rx*ry*Q2M(ia+1,ja+1)                    +   &
           rxm*ry*Q2M(ia,ja+1)
      cloud(i,j)    = rxm*rym*TCC(ia,ja)           +   &
           rx*rym*TCC(ia+1,ja)                     +   &
           rx*ry*TCC(ia+1,ja+1)                    +   &
           rxm*ry*TCC(ia,ja+1)
      rain(i,j)     = rxm*rym*TP(ia,ja)            +   &
           rx*rym*TP(ia+1,ja)                      +   &
           rx*ry*TP(ia+1,ja+1)                     +   &
           rxm*ry*TP(ia,ja+1)
      swrad(i,j)    = rxm*rym*SSRD(ia,ja)          +   &
           rx*rym*SSRD(ia+1,ja)                    +   &
           rx*ry*SSRD(ia+1,ja+1)                   +   &
           rxm*ry*SSRD(ia,ja+1)
      lwrad(i,j)    = rxm*rym*STRD(ia,ja)          +   &
           rx*rym*STRD(ia+1,ja)                    +   &
           rx*ry*STRD(ia+1,ja+1)                   +   &
           rxm*ry*STRD(ia,ja+1)
   ENDDO
ENDDO

! Rotate wind components from E/W and N/S to Xi and Eta directions
Uwind = scr1_out*cos(angle) + scr2_out*sin(angle)
Vwind = scr2_out*cos(angle) - scr1_out*sin(angle)
rain = rain*1000./delt
swrad = swrad/delt
lwrad = lwrad/delt
else
Uwind = 0.
Vwind = 0.
Pair=0.
Tair=0.
Qair=0.
Cloud=0.
rain= 0.
swrad= 0.
lwrad= 0.
   Uwind = U10(i1:i2,j1:j2)
   Vwind = V10(i1:i2,j1:j2)          
   Pair  = MSL(i1:i2,j1:j2)       
   Tair  = T2M(i1:i2,j1:j2)
   Qair  = Q2M(i1:i2,j1:j2)
   cloud = TCC(i1:i2,j1:j2)
   rain  = TP(i1:i2,j1:j2)*1000./delt 
   swrad = SSRD(i1:i2,j1:j2)/delt
   lwrad = STRD(i1:i2,j1:j2)/delt     
end if

! Write out variables in netCDF file
if(tday_AN.eq.tday_FC)then

write(*,*)'Output file record is: ', irec

statuso = nf90_put_var(ncido,TimeOutVarId,tday_AN,start=(/ irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncido,UwindVarId,Uwind,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncido,VwindVarId,Vwind,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncido,PairVarId,Pair,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncido,TairVarId,Tair,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncido,QairVarId,Qair,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncido,cloudVarId,cloud,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

status_roms = nf90_put_var(ncido,rainVarId,rain,start=(/1, 1, irec/))
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

status_roms = nf90_put_var(ncido,swradVarId,swrad,start=(/1, 1,  irec/))
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

status_roms = nf90_put_var(ncido,lwradVarId,lwrad,start=(/1, 1, irec/))
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

else
write(*,*)'****************************************'
write(*,*)'Problem in time:'
write(*,*)'Time of the AN variables is: ', tday_AN
write(*,*)'Time of the FC variables is: ', tday_FC
write(*,*)'****************************************'
STOP
endif

ENDDO

statuso = nf90_close(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

END PROGRAM era2roms
