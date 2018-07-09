! Program to generate ROMS forcing files from CFSR netCDF files
! produced using cvt_gribfiles.bsh

PROGRAM cfsr_winds2roms
use netcdf
implicit none

INTEGER i, j, k, l, imon, itime, ntimes, ia, ja, IFL
INTEGER irec
INTEGER iyear, iday, imonth,ihour, nb_records, indice
INTEGER greg(6)
INTEGER LonDimID, LatDimID, TimeDimId
INTEGER nlonw,nlatw
INTEGER XiDimID, EtaDimID, Lp, Mp, Lp_dim, Mp_dim
INTEGER Lp_wind,Mp_wind 
INTEGER LonVarId, LatVarId
INTEGER LonRhoVarId, LatRhoVarId, AngleVarId
INTEGER U10VarId, V10VarId
INTEGER TimeVarId_wind
INTEGER year_st, month_st, day_st
INTEGER year_end, month_end, day_end
INTEGER year_ref, month_ref, day_ref

INTEGER XiRhoDimId, EtaRhoDimId, TimeOutDimID
INTEGER TimeOutVarId, UwindVarId,VwindVarId

INTEGER statuso, ncido
INTEGER status_cfsr, ncid_cfsr
INTEGER status_roms, ncid_roms
INTEGER status_wind, ncid_wind

INTEGER i1_wind, i2_wind, j1_wind, j2_wind

INTEGER jdref, jdref_cfsr
INTEGER yes_no

REAL  rx, rxm, ry, rym
REAL  riminw, rimaxw, rjminw, rjmaxw

REAL, DIMENSION(:,:), allocatable :: ipos_wind
REAL, DIMENSION(:,:), allocatable :: jpos_wind

REAL, ALLOCATABLE, DIMENSION(:) :: lonw,latw
REAL, ALLOCATABLE, DIMENSION(:) :: lonwtemp,latwtemp
REAL, ALLOCATABLE, DIMENSION(:) :: IX_wind,IY_wind
REAL, ALLOCATABLE, DIMENSION(:,:) :: work_wind

REAL, ALLOCATABLE, DIMENSION(:,:) :: scr1_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr2_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: lon_rho
REAL, ALLOCATABLE, DIMENSION(:,:) :: lat_rho
REAL, ALLOCATABLE, DIMENSION(:,:) :: angle

REAL, ALLOCATABLE, DIMENSION(:,:) :: U10
REAL, ALLOCATABLE, DIMENSION(:,:) :: V10

REAL, ALLOCATABLE, DIMENSION(:,:) :: Uwind
REAL, ALLOCATABLE, DIMENSION(:,:) :: Vwind

DOUBLE PRECISION tday_wind, tday_msl, tday_t2m, tday_q2m, tday_rain, tday_cloud, tday_dlw, tday_dsw, tday
REAL datetime

REAL, PARAMETER    :: undef = 2.E+35            ! Undefined land value
REAL pi, DTOR

DOUBLE PRECISION date

CHARACTER(len=4) year_char
CHARACTER(len=2) month_char, day_char
CHARACTER(len=50)  ryear, num_temp, rmonth
CHARACTER(len=160) path_cfsr, maskfile, gridfile, fileout, time_units_att
CHARACTER(len=80) wind10m_file

interface
     FUNCTION jd(yy,mm,dd)
     implicit none
     INTEGER jd,yy,mm,dd
     end function jd
end interface

pi = ATAN(1.)*4.
DTOR = pi/180.

OPEN(unit=10,file='cfsr_winds2roms.in')
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

date= REAL(jd(year_st,month_st,day_st))+REAL(irec-1)*(6.d0/24.d0)+(1.d0/24.d0)+ 0.00001
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

if(irec.eq.1)then
status_wind = nf90_inq_dimid(ncid_wind, "lon", LonDimID)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_wind = nf90_inq_dimid(ncid_wind, "lat", LatDimID)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_wind = nf90_Inquire_Dimension(ncid_wind,LonDimID,len = nlonw)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_wind = nf90_Inquire_Dimension(ncid_wind,LatDimID,len = nlatw)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)

ALLOCATE(lonw(nlonw))
ALLOCATE(latw(nlatw))
ALLOCATE(IX_wind(nlonw))
ALLOCATE(IY_wind(nlatw))
ALLOCATE(work_wind(nlonw,nlatw))
ALLOCATE(U10(nlonw,nlatw))
ALLOCATE(V10(nlonw,nlatw))

! Get lons and lats
status_wind = nf90_inq_varid(ncid_wind, "lon",LonVarId)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_wind = nf90_inq_varid(ncid_wind, "lat",LatVarId)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_wind = nf90_get_var(ncid_wind,LonVarId,lonw)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
status_wind = nf90_get_var(ncid_wind,LatVarId,latw)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)

! Get i,j positions of ROMS rho-points on CFSR grid
WHERE(lonw > 180.) lonw = lonw-360.
DO I=1,nlonw
IX_wind(I)= I
ENDDO
DO I=1,nlatw
IY_wind(I)= I
ENDDO

CALL SSORT (lonw, IX_wind, nlonw, 2)
CALL SSORT (latw, IY_wind, nlatw, 2)

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
ENDDO
ENDDO
! Find subarea in CFSR grid containing model grid
riminw = 1.e+23
rimaxw =-1.e+23
rjminw = 1.e+23
rjmaxw =-1.e+23
DO j=1,Mp
DO i=1,Lp
   riminw = MIN(riminw,ipos_wind(i,j))
   rimaxw = MAX(rimaxw,ipos_wind(i,j))
   rjminw = MIN(rjminw,jpos_wind(i,j))
   rjmaxw = MAX(rjmaxw,jpos_wind(i,j))
ENDDO
ENDDO
i1_wind = FLOOR(riminw)
i2_wind = CEILING(rimaxw)
j1_wind = FLOOR(rjminw)
j2_wind = CEILING(rjmaxw)

IF(i1_wind<1.OR.i2_wind>nlonw.OR.j1_wind<1.OR.j2_wind>nlatw) THEN
   write(*,*) 'min indice to extract in longitude: ',i1_wind
   write(*,*) 'max indice to extract in longitude: ',i2_wind
   write(*,*) 'min indice to extract in latitude: ',j1_wind
   write(*,*) 'max indice to extract in latitude: ',j2_wind
   write(*,*) 'Maximum dimension in longitude :',nlonw
   write(*,*) 'Maximum dimension in latitude :',nlatw
   STOP 'Wind subdomain outside CFSR domain'
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
else
LP_wind=LP
MP_wind=MP
end if

LP_dim= LP_wind
MP_dim= MP_wind

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
                  "Atmospheric forcing file - CFSR grid")
else
statuso = nf90_put_att(ncido,nf90_global,"title", &
                  "Atmospheric forcing file from CSFR - ROMS grid")
end if
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,nf90_global,"type", &
                  "ROMS forcing file")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,nf90_global,"history", &
                  "Produced with cfsr2roms.f90")
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
statuso = nf90_enddef(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

WRITE(*,*) 'Output file ', TRIM(fileout), ' is created.'

ALLOCATE(Uwind(Lp_dim,Mp_dim))
ALLOCATE(Vwind(Lp_dim,Mp_dim))

!Open output file
statuso = nf90_open(trim(fileout),nf90_write,ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

if(yes_no.eq.0)then
! Write latitude and longitude
statuso = nf90_put_var(ncido,LonVarId,lonw(i1_wind:i2_wind))
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_var(ncido,LatVarId,latw(j1_wind:j2_wind))
if(statuso /= nf90_NoErr) call handle_err(statuso)
end if

endif


! Get wind time var id
status_wind = nf90_inq_varid(ncid_wind, "time",TimeVarId_wind)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
! Get U10 var id
status_wind = nf90_inq_varid(ncid_wind, "U_GRD_103_10",U10VarId)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)
! Get V10 var id
status_wind = nf90_inq_varid(ncid_wind, "V_GRD_103_10",V10VarId)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)

! Read in Time   
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

tday= tday_wind

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

if(yes_no.eq.1)then
! Horizontal interpolation
scr1_out = 0.
scr2_out = 0.
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
   ENDDO
ENDDO

! Rotate wind components from E/W and N/S to Xi and Eta directions

Uwind = scr1_out*cos(angle) + scr2_out*sin(angle)
Vwind = scr2_out*cos(angle) - scr1_out*sin(angle)

else
   Uwind = U10(i1_wind:i2_wind,j1_wind:j2_wind)
   Vwind = V10(i1_wind:i2_wind,j1_wind:j2_wind)               
end if

write(*,*)'Output file record is: ', irec

statuso = nf90_put_var(ncido,TimeOutVarId,tday,start=(/ irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncido,UwindVarId,Uwind,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncido,VwindVarId,Vwind,start=(/1, 1, irec /))
if(statuso /= nf90_NoErr) call handle_err(statuso)

status_wind = nf90_close(ncid_wind)
if(status_wind /= nf90_NoErr) call handle_err(status_wind)

ENDDO

statuso = nf90_close(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

END PROGRAM cfsr_winds2roms














