PROGRAM hycom2stations
! Program to obtain vertical interpolation of T,S,u,v from z to s levels on rho grid

use netcdf
USE, INTRINSIC :: IEEE_ARITHMETIC
implicit none

interface
     FUNCTION jd(yy,mm,dd)
     implicit none
     INTEGER jd,yy,mm,dd
     end function jd
end interface


INTEGER i, j, k, l,t,itime
INTEGER imin,imax,jmin,jmax,yes_no
INTEGER yy,mm,dd,nb_day,start_record,start_check,record
INTEGER yy_ref,mm_ref,dd_ref
INTEGER ia,ja,nb_stat
REAL rx,ry,rxm,rym

CHARACTER*160 station_list, fileout,num_temp
CHARACTER*160 pathin, pathout, filein_salt, filein_temp, filein_u, filein_v
CHARACTER*160 PHCfile_salt, PHCfile_temp, PHCfile_u, PHCfile_v, PHCfile_ssh
CHARACTER*80 unit_name
CHARACTER*4 year_char
CHARACTER*2 month_char,day_char

INTEGER XiDimID, EtaDimID
INTEGER LonDimId,LatDimId,LevelDimId,nlon,nlat,nlev
INTEGER LonRhoVarId, LatRhoVarId, AngleVarId
INTEGER status_roms,ncid_roms
INTEGER statuso,ncido

INTEGER, ALLOCATABLE, DIMENSION(:)     :: ipos
INTEGER, ALLOCATABLE, DIMENSION(:)     :: jpos
INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: IX_pol
INTEGER, ALLOCATABLE, DIMENSION(:)     :: IX

REAL, ALLOCATABLE, DIMENSION(:) :: lon_rho
REAL, ALLOCATABLE, DIMENSION(:) :: lat_rho

REAL, ALLOCATABLE, DIMENSION(:)     :: lon_hycom
REAL, ALLOCATABLE, DIMENSION(:)     :: lat_hycom
REAL, ALLOCATABLE, DIMENSION(:,:)   :: lon_hycom_pol
REAL, ALLOCATABLE, DIMENSION(:,:)   :: lat_hycom_pol
REAL, ALLOCATABLE, DIMENSION(:)     :: depth_hycom
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SaltPHC
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: TempPHC
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: UPHC
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: VPHC
REAL, ALLOCATABLE, DIMENSION(:,:)   :: SSHPHC
REAL, ALLOCATABLE, DIMENSION(:,:)   :: work2d
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: work

DOUBLE PRECISION TimePHC, time_offset, time

INTEGER ncid_PHC_ssh, ncid_PHC_salt, ncid_PHC_temp, ncid_PHC_u, ncid_PHC_v, status_PHC
INTEGER SSHPHCVarId, TimePHCVarId, LevelPHCVarId, LonPHCVarId, LatPHCVarId
INTEGER UPHCVarId, VPHCVarId, ZetaPHCVarId, SaltPHCVarId, TempPHCVarId

INTEGER XiRhoDimIdo,EtaRhoDimIdo,DepthDimIdo,TimeDimIdo
INTEGER XiUDimIdo,EtaUDimIdo,XiVDimIdo,EtaVDimIdo
INTEGER ZetaVarIdo, TimeVarIdo, DepthVarIdo,UVarIdo, VVarIdo, SaltVarIdo, TempVarIdo

REAL, ALLOCATABLE, DIMENSION(:,:) :: salt_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: u_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: v_out
REAL, ALLOCATABLE, DIMENSION(:)   :: zeta_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr3d_u
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr3d_v
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr3d_uout
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr3d_vout


CHARACTER*160, ALLOCATABLE, DIMENSION(:) :: path_final
CHARACTER*4,   ALLOCATABLE, DIMENSION(:) :: num_file

open(unit=10,file='hycom2stations.in')
read(10,*) ! File containing the number of stations and the coordinates (lon lat)
read(10,'(a)') station_list
read(10,*) ! Path of the folder containing the HYCOM data
read(10,'(a)') pathin
read(10,*) ! Path and name of the Output NetCDF file
read(10,'(a)') pathout
read(10,*) ! Year Month and Day of the first day to extract and number of days to extract
read(10,*) yy,mm,dd,nb_day
read(10,*) ! Reference date of the netcdf file
read(10,*) yy_ref,mm_ref,dd_ref
read(10,*) ! Polar coordinates YES=1/NO=0
read(10,*) yes_no   ! polar coordinates in HYCOM
close(10)

!----------------------------------------------------------------------
write(*,*) 'station file is ',station_list
open(unit=10,file= station_list)
read(10,*) nb_stat
ALLOCATE(lon_rho(nb_stat))
ALLOCATE(lat_rho(nb_stat))
ALLOCATE(ipos(nb_stat))
ALLOCATE(jpos(nb_stat))
do i=1,nb_stat
read(10,*) lon_rho(I),lat_rho(I)
end do
close(10)

write(*,*) 'The number of stations is ',nb_stat

!-------------------------------------------------------------------------------
! Get time offset
!   HYCOM times are in seconds referenced to 1970-01-01 00:00:00
time_offset = DBLE(jd(1970,1,1)-jd(yy_ref,mm_ref,dd_ref))*24.d0*3600.d0
! ------------------------------------------------------------------------------
write(*,*) 'Define HYCOM files to horizontally interpolate into the ROMS gridfile'

ALLOCATE(path_final(nb_day))
ALLOCATE(num_file(nb_day))

if (((yy.eq.2003).and.(mm.lt.11)).or.((yy.eq. 2003).and.(mm.eq.11).and.(dd.lt.3)))then
start_record=INT(jd(yy,mm,dd)-jd(2003,1,2))+1
start_check=1
else
start_record=INT(jd(yy,mm,dd)-jd(2003,11,3))+1
start_check=0
endif
if (start_check.eq.1)then
    do i=1,nb_day
        record =  INT(jd(yy,mm,dd)-jd(2003,1,2))+ i
        if (record.lt.298)then
            path_final(i)= TRIM(pathin)//'2003-01-02_to_2003-11-02/'
        else
            path_final(i)= TRIM(pathin)//'2003-11-03_to_nowdays/'
            record= record -298
        endif
        write(num_temp,*) record
        num_temp=ADJUSTL(num_temp)
        if (len_trim(num_temp).EQ.1) num_file(i)= '000'//TRIM(num_temp)
        if (len_trim(num_temp).EQ.2) num_file(i)= '00'//TRIM(num_temp)
        if (len_trim(num_temp).EQ.3) num_file(i)= '0'//TRIM(num_temp)
        if (len_trim(num_temp).EQ.4) num_file(i)= TRIM(num_temp)
    enddo
else
    do i=1,nb_day
        record = INT(jd(yy,mm,dd)-jd(2003,11,3))+ i
        path_final(i)= TRIM(pathin)//'2003-11-03_to_nowdays/'
        write(num_temp,*) record
        num_temp=ADJUSTL(num_temp)
        if (len_trim(num_temp).EQ.1) num_file(i)= '000'//TRIM(num_temp)
        if (len_trim(num_temp).EQ.2) num_file(i)= '00'//TRIM(num_temp)
        if (len_trim(num_temp).EQ.3) num_file(i)= '0'//TRIM(num_temp)
        if (len_trim(num_temp).EQ.4) num_file(i)= TRIM(num_temp)
    enddo    
endif

! ------------------------------------------------------------------------------
! Get dimensions of HYCOM files

PHCfile_salt = TRIM(path_final(1))//'salt_opendap_'//TRIM(num_file(1))//'.nc'
status_PHC = nf90_open(TRIM(PHCfile_salt),nf90_nowrite,ncid_PHC_salt)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_dimid(ncid_PHC_salt, "X", LonDimId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_dimid(ncid_PHC_salt, "Y", LatDimId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_dimid(ncid_PHC_salt, "depth",LevelDimId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_Inquire_Dimension(ncid_PHC_salt,LonDimID,len = nlon)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_Inquire_Dimension(ncid_PHC_salt,LatDimID,len = nlat)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_Inquire_Dimension(ncid_PHC_salt,LevelDimID,len = nlev)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

if(yes_no.eq.0)then
ALLOCATE(lon_hycom(nlon))
ALLOCATE(lat_hycom(nlat))
ALLOCATE(IX(nlon))
elseif(yes_no.eq.1)then
ALLOCATE(lon_hycom_pol(nlon,nlat))
ALLOCATE(lat_hycom_pol(nlon,nlat))
ALLOCATE(IX_pol(nlon,nlat))
endif
ALLOCATE(depth_hycom(nlev))

status_PHC = nf90_inq_varid(ncid_PHC_salt, "longitude",LonPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_PHC_salt, "latitude",LatPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_PHC_salt, "depth",LevelPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
if(yes_no.eq.0)then
status_PHC = nf90_get_var(ncid_PHC_salt,LonPHCVarId,lon_hycom)
if(status_PHC /= nf90_NoErr) call handle_err(status_roms)
status_PHC = nf90_get_var(ncid_PHC_salt,LatPHCVarId,lat_hycom)
if(status_PHC /= nf90_NoErr) call handle_err(status_roms)
DO i=1,nlon
IX(I)= I
if(lon_hycom(i).gt.180.)then
  lon_hycom(i)=lon_hycom(i)-360.d0
endif
ENDDO
CALL SSORT (lon_hycom, IX, nlon, 2)
elseif(yes_no.eq.1)then
status_PHC = nf90_get_var(ncid_PHC_salt,LonPHCVarId,lon_hycom_pol)
if(status_PHC /= nf90_NoErr) call handle_err(status_roms)
status_PHC = nf90_get_var(ncid_PHC_salt,LatPHCVarId,lat_hycom_pol)
if(status_PHC /= nf90_NoErr) call handle_err(status_roms)
DO j=1,nlat
DO i=1,nlon
IX_pol(I,J)= I
if(lon_hycom_pol(i,j).gt.180.)then
  lon_hycom_pol(i,j)=lon_hycom_pol(i,j)-360.d0
endif
ENDDO
CALL SSORT (lon_hycom_pol(:,J), IX_pol(:,J), nlon, 2)
ENDDO
endif
status_PHC = nf90_get_var(ncid_PHC_salt,LevelPHCVarId,depth_hycom)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

status_PHC = nf90_close(ncid_PHC_salt)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

! ------------------------------------------------------------------------------
! grid specification for input and output

ALLOCATE(SaltPHC(nlev,nlon,nlat))
ALLOCATE(TempPHC(nlev,nlon,nlat))
ALLOCATE(UPHC(nlev,nlon,nlat))
ALLOCATE(VPHC(nlev,nlon,nlat))
ALLOCATE(SSHPHC(nlon,nlat))
ALLOCATE(work(nlev,nlon,nlat))
ALLOCATE(work2d(nlon,nlat))

ALLOCATE(salt_out(nb_stat,nlev))
ALLOCATE(temp_out(nb_stat,nlev))
ALLOCATE(u_out(nb_stat,nlev))
ALLOCATE(v_out(nb_stat,nlev))
ALLOCATE(zeta_out(nb_stat))

! ------------------------------------------------------------------------------
! Get i,j positions of ROMS rho-points on HYCOM grid

if(yes_no.eq.0)then
DO i=1,nb_stat
DO k=1,nlon-1
   if((lon_rho(i).ge.lon_hycom(k)).AND.(lon_rho(i).lt.lon_hycom(k+1)))then
   ipos(i) = k
   endif
ENDDO
DO l=1,nlat-1
   if((lat_rho(i).ge.lat_hycom(l)).AND.(lat_rho(i).lt.lat_hycom(l+1)))then
   jpos(i) = l
   endif
ENDDO
ENDDO
elseif(yes_no.eq.1)then
DO i=1,nb_stat
DO k=1,nlon-1
DO l=1,nlat-1
   if( &
   (((lon_rho(i).ge.lon_hycom_pol(k,l)).AND.(lon_rho(i).lt.lon_hycom_pol(k+1,l)))   &
   .AND.(((lat_rho(i).ge.lat_hycom_pol(k,l)).AND.(lat_rho(i).lt.lat_hycom_pol(k,l+1))))) &
   .OR. &
   (((lon_rho(i).ge.lon_hycom_pol(k,l+1)).AND.(lon_rho(i).lt.lon_hycom_pol(k+1,l)))   &
   .AND.(((lat_rho(i).ge.lat_hycom_pol(k,l)).AND.(lat_rho(i).lt.lat_hycom_pol(k+1,l+1))))) &
   .OR.  &
   (((lon_rho(i).ge.lon_hycom_pol(k,l)).AND.(lon_rho(i).lt.lon_hycom_pol(k+1,l+1)))   &
   .AND.(((lat_rho(i).ge.lat_hycom_pol(k+1,l)).AND.(lat_rho(i).lt.lat_hycom_pol(k,l+1))))) &
   )then
   ipos(i) = k
   jpos(i) = l
   endif
ENDDO
ENDDO
if(ipos(i).eq.0.or.jpos(i).eq.0) write(*,*) 'Problem with station ',i,' at longitude= ',lon_rho(i),' and latitude= ',lat_rho(i)
ENDDO
endif

! Find subarea in HYCOM grid containing model grid
imin =  1000000
imax = -1000000
jmin =  1000000
jmax = -1000000
DO i=1,nb_stat
   imin = MIN(imin,ipos(i))
   imax = MAX(imax,ipos(i))
   jmin = MIN(jmin,jpos(i))
   jmax = MAX(jmax,jpos(i))
ENDDO

IF(imin<1.OR.imax>nlon.OR.jmin<1.OR.jmax>nlat) THEN
   write(*,*) 'imin,imax,jmin,jmax = ',imin,imax,jmin,jmax
   write(*,*) 'nlon, nlat = ',nlon,nlat
   STOP 'station outside HYCOM domain'
ELSE
   write(*,*) 'imin,imax,jmin,jmax = ',imin,imax,jmin,jmax
ENDIF

! ------------------------------------------------------
! Generate and open output netCDF files

fileout = TRIM(pathout)
statuso = nf90_create(TRIM(fileout),nf90_clobber,ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)
! Dimensions
statuso = nf90_def_dim(ncido,"stations",nb_stat,XiRhoDimIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_dim(ncido,"depth",nlev,DepthDimIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_dim(ncido,"time",nf90_unlimited,TimeDimIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
! Global attributes
statuso = nf90_put_att(ncido,nf90_global,"title", &
                  "Stations from HYCOM z-levels")
if(statuso /= nf90_NoErr) call handle_err(status_roms)
statuso = nf90_put_att(ncido,nf90_global,"type", &
                  "HYCOM station file")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,nf90_global,"history", &
                  "Produced with hycom2stations.f90")
if(statuso /= nf90_NoErr) call handle_err(statuso)
!Salinity
statuso = nf90_def_var(ncido,"salinity",nf90_float, &
                 (/ XiRhoDimIdo, DepthDimIdo,TimeDimIdo /), SaltVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, SaltVarIdo,"long_Name", &
                 "Salinity")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, SaltVarIdo,"units", &
                 "PSU")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, SaltVarIdo,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
!Temperature
statuso = nf90_def_var(ncido,"temperature",nf90_float, &
                 (/ XiRhoDimIdo, DepthDimIdo,TimeDimIdo /), TempVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, TempVarIdo,"long_Name", &
                 "Temperature")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, TempVarIdo,"units", &
                 "Degree C")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, TempVarIdo,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
!U velocity
statuso = nf90_def_var(ncido,"u",nf90_float, &
                 (/ XiRhoDimIdo, DepthDimIdo,TimeDimIdo /), UVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, UVarIdo,"long_Name", &
                 "U velocity")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, UVarIdo,"units", &
                 "meter.second-1")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, UVarIdo,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
!V velocity
statuso = nf90_def_var(ncido,"v",nf90_float, &
                 (/ XiRhoDimIdo, DepthDimIdo,TimeDimIdo /), VVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, VVarIdo,"long_Name", &
                 "V velocity")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, VVarIdo,"units", &
                 "meter.second-1")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, VVarIdo,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
!ssh
statuso = nf90_def_var(ncido,"ssh",nf90_float, &
                 (/ XiRhoDimIdo, TimeDimIdo /), ZetaVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, ZetaVarIdo,"long_Name", &
                 "Sea surface height")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, ZetaVarIdo,"units", &
                 "meter")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, ZetaVarIdo,"time", &
                 "time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
! depth
statuso = nf90_def_var(ncido,"depth",nf90_float, &
                 DepthDimIdo, DepthVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, DepthVarIdo,"long_Name", &
                 "Depth of the z-levels")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, DepthVarIdo,"units", &
                 "meter")
if(statuso /= nf90_NoErr) call handle_err(statuso)
! Time
statuso = nf90_def_var(ncido,"time",nf90_float, &
                 TimeDimIdo, TimeVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
!
write(num_temp,*) yy_ref
year_char=ADJUSTL(num_temp)
write(num_temp,*) mm_ref
num_temp=ADJUSTL(num_temp)
if (len_trim(num_temp).EQ.1) month_char= '0'//TRIM(num_temp)
if (len_trim(num_temp).EQ.2) month_char= TRIM(num_temp)
write(num_temp,*) dd_ref
num_temp=ADJUSTL(num_temp)
if (len_trim(num_temp).EQ.1) day_char= '0'//TRIM(num_temp)
if (len_trim(num_temp).EQ.2) day_char= TRIM(num_temp)
unit_name='seconds since '//TRIM(year_char)//'-'//TRIM(month_char)//'-'//TRIM(day_char)//' 00:00:00'
!
statuso = nf90_put_att(ncido, TimeVarIdo,"long_Name", &
                 "Climatology time")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, TimeVarIdo,"units", &
                 TRIM(unit_name))
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, TimeVarIdo,"calendar", &
                 "gregorian")
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido, TimeVarIdo,"field", &
                 "time, scalar, series")
if(statuso /= nf90_NoErr) call handle_err(statuso)
! End of definition
statuso = nf90_enddef(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

!Open and read output file

statuso = nf90_open(trim(fileout),nf90_write,ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_inq_varid(ncido, "salinity",SaltVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output salt variable'
statuso = nf90_inq_varid(ncido, "temperature",TempVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output temp variable'
statuso = nf90_inq_varid(ncido, "u",UVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output u variable'
statuso = nf90_inq_varid(ncido, "v",VVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output v variable'
statuso = nf90_inq_varid(ncido, "ssh",ZetaVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output zeta variable'
statuso = nf90_inq_varid(ncido, "depth",DepthVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output depth variable'
statuso = nf90_inq_varid(ncido, "time",TimeVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output clim_time variable'

write(*,*) 'Finished opening output file: ', fileout

! Output depth to netCDF file
   statuso = nf90_put_var(ncido,DepthVarIdo,depth_hycom)
   if(statuso /= nf90_NoErr) call handle_err(statuso)

! -----------------------------------------------------------------------
! Generate the climatology file on the z-levels

DO T=1,nb_day

itime=t

! Open PHC file and find variables
! Open input averages file

PHCfile_salt = TRIM(path_final(T))//'salt_opendap_'//TRIM(num_file(T))//'.nc'
PHCfile_temp = TRIM(path_final(T))//'temp_opendap_'//TRIM(num_file(T))//'.nc'
PHCfile_u =    TRIM(path_final(T))//'uvel_opendap_'//TRIM(num_file(T))//'.nc'
PHCfile_v =    TRIM(path_final(T))//'vvel_opendap_'//TRIM(num_file(T))//'.nc'
PHCfile_ssh =  TRIM(path_final(T))//'ssh_opendap_' //TRIM(num_file(T))//'.nc'
status_PHC = nf90_open(trim(PHCfile_salt),nf90_nowrite,ncid_PHC_salt)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Opened input salinity file'
status_PHC = nf90_inq_varid(ncid_PHC_salt, "time",TimePHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable time'
status_PHC = nf90_inq_varid(ncid_PHC_salt, "salinity",SaltPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable salinity'
status_PHC = nf90_open(trim(PHCfile_temp),nf90_nowrite,ncid_PHC_temp)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Opened input temperature file'
status_PHC = nf90_inq_varid(ncid_PHC_temp, "temperature",TempPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable temperature'
status_PHC = nf90_open(trim(PHCfile_u),nf90_nowrite,ncid_PHC_u)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Opened input u file'
status_PHC = nf90_inq_varid(ncid_PHC_u, "u",UPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable u'
status_PHC = nf90_open(trim(PHCfile_v),nf90_nowrite,ncid_PHC_v)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Opened input v file'
status_PHC = nf90_inq_varid(ncid_PHC_v, "v",VPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable v'
status_PHC = nf90_open(trim(PHCfile_ssh),nf90_nowrite,ncid_PHC_ssh)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Opened input ssh file'
status_PHC = nf90_inq_varid(ncid_PHC_ssh, "ssh",SSHPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable ssh'

!read in time
   status_PHC = nf90_get_var(ncid_PHC_salt,TimePHCVarId,TimePHC)
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
   write(*,*) 'Time read in'

   time = DBLE(time_offset + TimePHC)

!-------------------------------------------------------
! Do salinity

   status_PHC = nf90_get_var(ncid_PHC_salt,SaltPHCVarId,SaltPHC)
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

! Horizontal interpolation
   salt_out = 0.
   work=0.
   DO I=1,nlon
      DO J=1,nlat
         DO K=1,nlev
            if(yes_no.eq.0)then
               work(K,I,J)=SaltPHC(K,IX(I),J)
            elseif(yes_no.eq.1)then
               work(K,I,J)=SaltPHC(K,IX_pol(I,J),J)
            endif
         ENDDO
      ENDDO
   ENDDO
   SaltPHC = work

   DO i=1,nb_stat
      ia=ipos(i)
      ja=jpos(i)
      if(yes_no.eq.0)then
      rx =  (lon_rho(i)-lon_hycom(ia))/(lon_hycom(ia+1)-lon_hycom(ia))
      ry =  (lat_rho(i)-lat_hycom(ja))/(lat_hycom(ja+1)-lat_hycom(ja))
      elseif(yes_no.eq.1)then
      rx =  (lon_rho(i)-lon_hycom_pol(ia,ja))/(lon_hycom_pol(ia+1,ja)-lon_hycom_pol(ia,ja))
      ry =  (lat_rho(i)-lat_hycom_pol(ia,ja))/(lat_hycom_pol(ia,ja+1)-lat_hycom_pol(ia,ja))
      endif
      rxm = 1.0-rx
      rym = 1.0-ry
   DO k=1,nlev
      salt_out(i,k)   = rxm*rym*SaltPHC(k,ia,ja)         +   &
                        rx*rym*SaltPHC(k,ia+1,ja)        +   &
                        rx*ry*SaltPHC(k,ia+1,ja+1)       +   &
                        rxm*ry*SaltPHC(k,ia,ja+1)
   ENDDO
   ENDDO
 

 ! Output to netCDF file
   statuso = nf90_put_var(ncido,SaltVarIdo,salt_out,(/1,1,itime/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'Salinity done at time=',time

! ---------------------------------------------------------------
! Do temperature

   status_PHC = nf90_get_var(ncid_PHC_temp,TempPHCVarId,TempPHC)
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

! Horizontal interpolation
   temp_out = 0.
   work=0.
   DO I=1,nlon
      DO J=1,nlat
         DO K=1,nlev
            if(yes_no.eq.0)then
               work(K,I,J)=TempPHC(K,IX(I),J)
            elseif(yes_no.eq.1)then
               work(K,I,J)=TempPHC(K,IX_pol(I,J),J)
            endif
         ENDDO
      ENDDO
   ENDDO
   SaltPHC = work

   DO i=1,nb_stat
      ia=ipos(i)
      ja=jpos(i)
      if(yes_no.eq.0)then
      rx =  (lon_rho(i)-lon_hycom(ia))/(lon_hycom(ia+1)-lon_hycom(ia))
      ry =  (lat_rho(i)-lat_hycom(ja))/(lat_hycom(ja+1)-lat_hycom(ja))
      elseif(yes_no.eq.1)then
      rx =  (lon_rho(i)-lon_hycom_pol(ia,ja))/(lon_hycom_pol(ia+1,ja)-lon_hycom_pol(ia,ja))
      ry =  (lat_rho(i)-lat_hycom_pol(ia,ja))/(lat_hycom_pol(ia,ja+1)-lat_hycom_pol(ia,ja))
      endif
      rxm = 1.0-rx
      rym = 1.0-ry
   DO k=1,nlev
      temp_out(i,k) = rxm*rym*TempPHC(k,ia,ja)         +   &
                        rx*rym*TempPHC(k,ia+1,ja)        +   &
                        rx*ry*TempPHC(k,ia+1,ja+1)       +   &
                        rxm*ry*TempPHC(k,ia,ja+1)
   ENDDO
   ENDDO

! Output to netCDF file
   statuso = nf90_put_var(ncido,TempVarIdo,temp_out,(/1,1,itime/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'Temperature done at time=',time

! ---------------------------------------------------------------
! Do U, V velocity components

! Get U
   status_PHC = nf90_get_var(ncid_PHC_u,UPHCVarId,UPHC)
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
! Get V
   status_PHC = nf90_get_var(ncid_PHC_v,VPHCVarId,VPHC)
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

! Horizontal interpolation
   u_out = 0.
   v_out=  0.
   work=0.
   DO I=1,nlon
      DO J=1,nlat
         DO K=1,nlev
            if(yes_no.eq.0)then
               work(K,I,J)=UPHC(K,IX(I),J)
            elseif(yes_no.eq.1)then
               work(K,I,J)=UPHC(K,IX_pol(I,J),J)
            endif
         ENDDO
      ENDDO
   ENDDO
   UPHC = work
   work=0.
   DO I=1,nlon
      DO J=1,nlat
         DO K=1,nlev
            if(yes_no.eq.0)then
               work(K,I,J)=VPHC(K,IX(I),J)
            elseif(yes_no.eq.1)then
               work(K,I,J)=VPHC(K,IX_pol(I,J),J)
            endif
         ENDDO
      ENDDO
   ENDDO
   VPHC = work
   DO i=1,nb_stat
      ia=ipos(i)
      ja=jpos(i)
      if(yes_no.eq.0)then
      rx =  (lon_rho(i)-lon_hycom(ia))/(lon_hycom(ia+1)-lon_hycom(ia))
      ry =  (lat_rho(i)-lat_hycom(ja))/(lat_hycom(ja+1)-lat_hycom(ja))
      elseif(yes_no.eq.1)then
      rx =  (lon_rho(i)-lon_hycom_pol(ia,ja))/(lon_hycom_pol(ia+1,ja)-lon_hycom_pol(ia,ja))
      ry =  (lat_rho(i)-lat_hycom_pol(ia,ja))/(lat_hycom_pol(ia,ja+1)-lat_hycom_pol(ia,ja))
      endif
      rxm = 1.0-rx
      rym = 1.0-ry
   DO k=1,nlev
      u_out(i,k)     =  rxm*rym*UPHC(k,ia,ja)         +   &
                        rx*rym*UPHC(k,ia+1,ja)        +   &
                        rx*ry*UPHC(k,ia+1,ja+1)       +   &
                        rxm*ry*UPHC(k,ia,ja+1)

      v_out(i,k)    =   rxm*rym*VPHC(k,ia,ja)         +   &
                        rx*rym*VPHC(k,ia+1,ja)        +   &
                        rx*ry*VPHC(k,ia+1,ja+1)       +   &
                        rxm*ry*VPHC(k,ia,ja+1)
   ENDDO
   ENDDO

! Output to netCDF file
   statuso = nf90_put_var(ncido,UVarIdo,u_out,(/1,1,itime/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncido,VVarIdo,v_out,(/1,1,itime/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'u,v done at time=',time

! ----------------------------------------------------------------
! Do sea level height

   status_PHC = nf90_get_var(ncid_PHC_ssh,SSHPHCVarId,SSHPHC)
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

! Horizontal interpolation
   zeta_out = 0.
   work2d=0.
   DO I=1,nlon
      DO J=1,nlat
         if(yes_no.eq.0)then
            work2d(I,J)=SSHPHC(IX(I),J)
         elseif(yes_no.eq.1)then
            work2d(I,J)=SSHPHC(IX_pol(I,J),J)
         endif
      ENDDO
   ENDDO
   SSHPHC = work2d
   DO i=1,nb_stat
      ia=ipos(i)
      ja=jpos(i)
      if(yes_no.eq.0)then
      rx =  (lon_rho(i)-lon_hycom(ia))/(lon_hycom(ia+1)-lon_hycom(ia))
      ry =  (lat_rho(i)-lat_hycom(ja))/(lat_hycom(ja+1)-lat_hycom(ja))
      elseif(yes_no.eq.1)then
      rx =  (lon_rho(i)-lon_hycom_pol(ia,ja))/(lon_hycom_pol(ia+1,ja)-lon_hycom_pol(ia,ja))
      ry =  (lat_rho(i)-lat_hycom_pol(ia,ja))/(lat_hycom_pol(ia,ja+1)-lat_hycom_pol(ia,ja))
      endif
      rxm = 1.0-rx
      rym = 1.0-ry
      zeta_out(i)   =   rxm*rym*SSHPHC(ia,ja)         +   &
                        rx*rym*SSHPHC(ia+1,ja)        +   &
                        rx*ry*SSHPHC(ia+1,ja+1)       +   &
                        rxm*ry*SSHPHC(ia,ja+1)
   ENDDO
   
! Output to netCDF file
   statuso = nf90_put_var(ncido,ZetaVarIdo,zeta_out,(/1,itime/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'zeta done at time=',time

! ------------------------------------------------------------------
! write time
   statuso = nf90_put_var(ncido,TimeVarIdo,time,(/itime/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
! ------------------------------------------------------------------
status_PHC = nf90_close(ncid_PHC_salt)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_close(ncid_PHC_temp)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_close(ncid_PHC_u)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_close(ncid_PHC_v)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_close(ncid_PHC_ssh)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

END DO

statuso = nf90_close(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)


END PROGRAM hycom2stations
