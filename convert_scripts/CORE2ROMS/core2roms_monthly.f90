! Program to generate ROMS forcing files from CORE netCDF files
! produced using cvt_gribfiles.bsh

PROGRAM core2roms
use netcdf
implicit none

integer, parameter :: DefInt8 = selected_int_kind(8)


INTEGER i, j, k, l, imon, itime, ntimes, ia, ja, IFL
INTEGER start_AN, STOP_AN, irec
INTEGER iyear, iday, imonth,ihour, nb_records, indice
INTEGER greg(6)
INTEGER nlons, nlats ,LonDimID, LatDimID, TimeDimId
INTEGER XiDimID, EtaDimID, Lp, Mp, Lp_dim, Mp_dim
INTEGER MaskVarId, LonVarId, LatVarId, DiVarId, DjVarId
INTEGER LonRhoVarId, LatRhoVarId, AngleVarId
INTEGER U10VarId, V10VarId, TimeVarId, Time2VarId, Time3VarId
INTEGER MSLVarId, T10VarId, Q10VarId, RAVarId, SNVarID
INTEGER SSRDVarId, STRDVarId
INTEGER year_st, month_st, day_st
INTEGER year_end, month_end, day_end
INTEGER year_ref, month_ref, day_ref
INTEGER year_inref, month_inref, day_inref


INTEGER XiRhoDimId, EtaRhoDimId, TimeOutDimID
INTEGER TimeOutVarId, UwindVarId,VwindVarId,PairVarId
INTEGER TairVarId, QairVarId,rainVarId, raintVarId, snowtVarId

INTEGER statuso, ncido
INTEGER status_toto, nctoto
INTEGER status_core, ncid_core
INTEGER status_roms, ncid_roms
INTEGER status_1, ncid_1, ncid_2,ncid_3,ncid_4,ncid_5,ncid_6,ncid_7

INTEGER ndays, n6hr, n12hr

INTEGER i1, i2, j1, j2
INTEGER(kind = DefInt8) :: jdref, jdref_core
INTEGER yes_no

REAL  Di, Dj, rx, rxm, ry, rym
REAL  rimin, rimax, rjmin, rjmax
REAL, DIMENSION(:,:), allocatable :: ipos
REAL, DIMENSION(:,:), allocatable :: jpos
REAL  delt

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
REAL, ALLOCATABLE, DIMENSION(:,:) :: T10
REAL, ALLOCATABLE, DIMENSION(:,:) :: Q10
REAL, ALLOCATABLE, DIMENSION(:,:) :: evapres
REAL, ALLOCATABLE, DIMENSION(:,:) :: TP
REAL, ALLOCATABLE, DIMENSION(:,:) :: SSRD
REAL, ALLOCATABLE, DIMENSION(:,:) :: STRD

REAL, ALLOCATABLE, DIMENSION(:,:) :: Uwind
REAL, ALLOCATABLE, DIMENSION(:,:) :: Vwind
REAL, ALLOCATABLE, DIMENSION(:,:) :: Pair
REAL, ALLOCATABLE, DIMENSION(:,:) :: Tair
REAL, ALLOCATABLE, DIMENSION(:,:) :: Qair
REAL, ALLOCATABLE, DIMENSION(:,:) :: rain
REAL, ALLOCATABLE, DIMENSION(:,:) :: swrad
REAL, ALLOCATABLE, DIMENSION(:,:) :: lwrad
DOUBLE PRECISION tday_1, tday_2, tday_3
DOUBLE PRECISION datetime

INTEGER year1, year2, year3, month1, month2, month3, day1, day2, day3,hour1

REAL, PARAMETER    :: undef = 2.E+35            ! Undefined land value
REAL pi, DTOR

DOUBLE PRECISION date, date2, date3, date4

CHARACTER(len=4) year_char
CHARACTER(len=2) month_char, day_char
CHARACTER(len=50)  ryear, num_temp
CHARACTER(len=160) path_core, maskfile, gridfile, fileout, time_units_att
CHARACTER(len=80)  file1,file2,file3,file4,file5,file6,file7

interface
     FUNCTION jd(yy,mm,dd)
     implicit none
     INTEGER jd,yy,mm,dd
     end function jd
end interface

pi = ATAN(1.)*4.
DTOR = pi/180.

OPEN(unit=10,file='core2roms_monthly.in')
READ(10,*)
READ(10,'(a)') path_core
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
READ(10,*)! reference time for input files
READ(10,*) year_inref, month_inref, day_inref
CLOSE(10)

! Get reference Julian day for output
jdref = jd(year_ref,month_ref,day_ref)
write(*,*) 'jdref is ', jdref
! Get reference Julian day (1948-01-01 00:00:00Z) for input
jdref_core = jd(year_inref,month_inref,day_inref)
write(*,*) 'jdref_core is ', jdref_core

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

! Get dimensions of CORE files
maskfile = TRIM(path_core) // 'land.sfc.gauss.nc'
status_core = nf90_open(TRIM(maskfile),nf90_nowrite,ncid_core)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_inq_dimid(ncid_core, "lon", LonDimID)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_inq_dimid(ncid_core, "lat", LatDimID)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_Inquire_Dimension(ncid_core,LonDimID,len = nlons)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_Inquire_Dimension(ncid_core,LatDimID,len = nlats)
if(status_core /= nf90_NoErr) call handle_err(status_core)
ALLOCATE(lons(nlons))
ALLOCATE(lats(nlats))
ALLOCATE(IX(nlons))
ALLOCATE(IY(nlats))
ALLOCATE(emask(nlons,nlats))
ALLOCATE(work(nlons,nlats))
ALLOCATE(U10(nlons,nlats))
ALLOCATE(V10(nlons,nlats))
ALLOCATE(MSL(nlons,nlats))
ALLOCATE(T10(nlons,nlats))
ALLOCATE(Q10(nlons,nlats))
ALLOCATE(evapres(nlons,nlats))
ALLOCATE(TP(nlons,nlats))
ALLOCATE(SSRD(nlons,nlats))
ALLOCATE(STRD(nlons,nlats))
! Get lons and lats
status_core = nf90_inq_varid(ncid_core, "lon",LonVarId)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_inq_varid(ncid_core, "lat",LatVarId)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_get_var(ncid_core,LonVarId,lons)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_get_var(ncid_core,LatVarId,lats)
if(status_core /= nf90_NoErr) call handle_err(status_core)
Di = 0.7
Dj = 0.7

status_core = nf90_inq_varid(ncid_core, "land",MaskVarId)
if(status_core /= nf90_NoErr) call handle_err(status_core)
status_core = nf90_get_var(ncid_core,MaskVarId,emask)

if(status_core /= nf90_NoErr) call handle_err(status_core)
! Don't need following line for CORE data - already 0 or 1.
!emask= emask*1.52597204419215e-05+0.499984740279558
write(*,*) 'nlons,nlats = ',nlons,nlats
status_core = nf90_close(ncid_core)
if(status_core /= nf90_NoErr) call handle_err(status_core)

! Get i,j positions of ROMS rho-points on COREI grid
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
! Find subarea in COREI grid containing model grid
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
   STOP 'subdomain outside CORE domain'
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
                  "Atmospheric forcing file - CORE grid")
else
statuso = nf90_put_att(ncido,nf90_global,"title", &
                  "Atmospheric forcing file - ROMS grid")
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,nf90_global,"type", &
                  "ROMS forcing file")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,nf90_global,"history", &
                  "Produced with core2roms.f90")
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

statuso = nf90_enddef(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

WRITE(*,*) 'Output file ', TRIM(fileout), ' is created.'

ALLOCATE(rain(Lp_dim,Mp_dim))

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

!Number of monthly records
nb_records = (year_end - year_st) * 12 + (month_end - month_st) + 1
write(*,*) 'monthly', nb_records

DO irec=1,nb_records

write(*,*) year_st+FLOOR((irec-1)/12.0)
write(*,*) month_st+(irec-1)-12*FLOOR((irec-1)/12.0)
date= REAL(jd(year_st+FLOOR((irec-1)/12.0),month_st+(irec-1)-12*FLOOR((irec-1)/12.0),15))

greg=0
CALL GR(greg,date)
iyear= greg(1)
imonth= greg(2)
iday= greg(3)
ihour= greg(4)

write(*,'(a24,i2,a1,i2,a1,i4,a4,i2,a1)')'Writing data for date: ',greg(3),'/',greg(2),'/',greg(1), ' at ',greg(4),'h'

!treatment of the monthly precip files
delt = 12.*3600.
write(ryear,*) iyear
ryear = adjustl(ryear)
itime = (iyear - year_inref) * 12 + (imonth - month_inref) + 1
write(*,*) 'monthly itime',itime
file6 = TRIM(path_core)//'precip.1948-2009.nc'
status_1 = nf90_open(TRIM(file6),nf90_nowrite,ncid_6)
if(status_1 /= nf90_NoErr) call handle_err(status_1)
! Get FC time var id
status_1 = nf90_inq_varid(ncid_6, "TIME",Time2VarId)
if(status_1 /= nf90_NoErr) call handle_err(status_1)
! Get SSRD var id
status_1 = nf90_inq_varid(ncid_6, "RAIN",RaintVarId)
if(status_1 /= nf90_NoErr) call handle_err(status_1)
! Get STRD var id
status_1 = nf90_inq_varid(ncid_6, "SNOW",SnowtVarId)
if(status_1 /= nf90_NoErr) call handle_err(status_1)
! Read in time value
status_core = nf90_get_var(ncid_6,Time2VarId,datetime,start = (/itime/))
if(status_core /= nf90_NoErr) call handle_err(status_core)
tday_2 = REAL(jdref_core) + datetime - REAL(jdref)
date3 =  REAL(jdref_core) + datetime
CALL GR(greg,date3)
year2= greg(1)
month2= greg(2)
day2= greg(3)
! convert tday to seconds
tday_2= tday_2*24.*3600.
write(*,*) year2, month2, day2


! Read in RAIN values
status_CORE = nf90_get_var(ncid_6,RaintVarId,work,start=(/1, 1,itime /))
if(status_CORE /= nf90_NoErr) call handle_err(status_CORE)
!SSRD = work*SSRD_scale_factor + SSRD_add_offset
TP = work
! Read in STRD values
status_CORE = nf90_get_var(ncid_6,SnowtVarId,work,start=(/1, 1,itime /))
if(status_CORE /= nf90_NoErr) call handle_err(status_CORE)
!STRD = work*STRD_scale_factor + STRD_add_offset
TP = TP + work

status_core = nf90_close(ncid_6)
if(status_core /= nf90_NoErr) call handle_err(status_core)

work=0.
DO I=1,nlons
   DO J=1,nlats
      work(I,J)=TP(IX(I),IY(J))
   ENDDO
ENDDO
TP= work
work=0.

if(yes_no.eq.1)then
! Horizontal interpolation
rain= 0.
DO j=1,Mp
   DO i=1,Lp
      
      ia = ipos(i,j)
      ja = jpos(i,j)
      rx = (lon_rho(i,j)-lons(ia))/(lons(ia+1)-lons(ia))
      ry = (lat_rho(i,j)-lats(ja))/(lats(ja+1)-lats(ja))
      rxm = 1.0-rx
      rym = 1.0-ry
      
      rain(i,j)     = rxm*rym*TP(ia,ja)            +   &
           rx*rym*TP(ia+1,ja)                      +   &
           rx*ry*TP(ia+1,ja+1)                     +   &
           rxm*ry*TP(ia,ja+1)
   ENDDO
ENDDO

! Convert rain from m to kg/m2s - no need for CORE
!rain = rain*1000./delt
! Convert radiation from W/m2s to W/m2 - no need for CORE
!swrad = swrad/delt
!lwrad = lwrad/delt
else
!   rain  = TP(i1:i2,j1:j2)*1000./delt 
!   swrad = SSRD(i1:i2,j1:j2)/delt
!   lwrad = STRD(i1:i2,j1:j2)/delt     

end if

! Write out variables in netCDF file
if(iyear.eq.year2 .AND. imonth .eq. month2) then

write(*,*)'Output file record is: ', irec

statuso = nf90_put_var(ncido,TimeOutVarId,tday_2,start=(/ irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

status_roms = nf90_put_var(ncido,rainVarId,rain,start=(/1, 1,  irec/))
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

else
write(*,*)'****************************************'
write(*,*)'Problem in time:'
write(*,*)'Time of the 6hrly variables is: ', iyear, imonth 
write(*,*)'Time of the monthly variables is: ', year2, month2
write(*,*)'****************************************'
STOP
endif

ENDDO

statuso = nf90_close(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

END PROGRAM core2roms
