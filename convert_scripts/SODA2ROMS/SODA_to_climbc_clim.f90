PROGRAM SODA_to_climbc
! Program to obtain climatology from SODA monthly T,S,U,V,zeta climatology

use netcdf
use typesizes

implicit none

INTEGER, PARAMETER :: nlev = 35, nlat = 269, nlon = 271
INTEGER :: nlonp, nlatp
REAL, DIMENSION(1) :: TimePHC
REAL, DIMENSION(nlev) :: LevelPHC
REAL, DIMENSION(nlon) :: LonPHC
REAL, DIMENSION(nlat) :: LatPHC
REAL, DIMENSION(nlon+2) :: lonp
REAL, DIMENSION(nlat+2) :: latp
REAL, DIMENSION(nlon,nlat,nlev) :: SaltPHC,TempPHC,UPHC,VPHC
REAL, DIMENSION(nlon,nlat) :: ZetaPHC
REAL, DIMENSION(nlon,nlat) :: scr
REAL, DIMENSION(nlon,nlat) :: scr1
REAL, DIMENSION(nlon,nlat) :: scra
REAL, DIMENSION(nlon,nlat) :: scrb
REAL, DIMENSION(nlon+2,nlat+2) :: scrup
REAL, DIMENSION(nlon+2,nlat+2) :: scrvp
REAL, DIMENSION(nlon+2,nlat+2) :: scrap
REAL, DIMENSION(nlon+2,nlat+2) :: scrbp
REAL, DIMENSION(nlon+2,nlat+2) :: scrp

INTEGER year, month, day, jd, jday1, jday2
REAL tday


INTEGER i, j, k, idate, ict
INTEGER No
REAL theta_so, theta_bo, Tclineo
CHARACTER*80 gridfile, PHCfile, outfile, PHCpath, outpath, filelst
CHARACTER*80 clim_root, avgfile, avgfilei

INTEGER XiDimID, EtaDimID, Lp, Mp, L, M
INTEGER XiRhoDimId, EtaRhoDimId
INTEGER LonRhoVarId, LatRhoVarId, HVarId, AngleVarId
INTEGER status_roms, ncid_roms
REAL, ALLOCATABLE, DIMENSION(:,:) :: lon_rho
REAL, ALLOCATABLE, DIMENSION(:,:) :: lat_rho
REAL, ALLOCATABLE, DIMENSION(:,:) :: h
REAL, ALLOCATABLE, DIMENSION(:,:) :: angle
REAL, ALLOCATABLE, DIMENSION(:,:) :: work
REAL, ALLOCATABLE, DIMENSION(:,:) :: error
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: z_ro
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: z_wo
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Hzo
REAL, ALLOCATABLE, DIMENSION(:) :: Cs_ro, sc_ro
REAL cffo, cff1o, cff2o, hco, hmino, cff_ro, cff1_ro, cff2_ro


INTEGER ncid_PHC, status_PHC
INTEGER TimePHCVarId, LevelPHCVarId, LonPHCVarId, LatPHCVarId
INTEGER UPHCVarId, VPHCVarId, ZetaPHCVarId, SaltPHCVarId, TempPHCVarId
INTEGER statuso, ncido, dim_xi_rhoo, dim_eta_rhoo, dim_timeo, dim_s_rhoo
INTEGER dim_xi_uo, dim_eta_uo, dim_xi_vo, dim_eta_vo
INTEGER UVarIdo, VVarIdo, TimeVarIdo, SaltVarIdo, TempVarIdo
INTEGER UbarVarIdo, VbarVarIdo, ZetaVarIdo
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: salt_out
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp_out
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: u_out
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: v_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: zeta_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: ubar_out
REAL, ALLOCATABLE, DIMENSION(:,:) :: vbar_out
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: u_rhoo
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: v_rhoo
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: u_rhoz
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: v_rhoz
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: scr3d
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr1o
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr2o



REAL dlon, dlat, lon0, lat0, lonr, latr
REAL rx, ry, rxm, rym, rz1, rz2
INTEGER ia, ja, kT
REAL, ALLOCATABLE, DIMENSION(:,:) :: ipos, jpos

REAL sumb, sumt
REAL rimin, rimax, rjmin, rjmax
INTEGER i1, i2, j1, j2
REAL, PARAMETER :: undef = 2.E+35 ! Undefined land value
REAL pi, DTOR
INTEGER nvalue
REAL tx, critx, cor
INTEGER mxs

INTEGER feof, nc
CHARACTER*6 numext

pi = ATAN(1.)*4.
DTOR = pi/180.
tx = 0.9*undef
critx = 0.01
cor = 1.6
mxs = 100

nlonp = nlon+2
nlatp = nlat+2

open(unit=10,file='SODA_clim.in')
read(*,'(a)') gridfile
read(*,'(a)') filelst
read(*,'(a)') PHCpath
read(*,'(a)') outpath
read(*,'(a)') clim_root
read(*,*) No, theta_so, theta_bo, Tclineo


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

L = Lp-1
M = Mp-1

ALLOCATE(work(nlonp,nlatp))
ALLOCATE(error(nlonp,nlatp))


ALLOCATE(lon_rho(Lp,Mp))
ALLOCATE(lat_rho(Lp,Mp))
ALLOCATE(h(Lp,Mp))
ALLOCATE(angle(Lp,Mp))

! Get lons, lats and h
status_roms = nf90_inq_varid(ncid_roms, "lon_rho",LonRhoVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_inq_varid(ncid_roms, "lat_rho",LatRhoVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_inq_varid(ncid_roms, "h",HVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_inq_varid(ncid_roms, "angle",AngleVarId)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

status_roms = nf90_get_var(ncid_roms,LonRhoVarId,lon_rho)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_get_var(ncid_roms,LatRhoVarId,lat_rho)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_get_var(ncid_roms,HVarId,h)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_roms = nf90_get_var(ncid_roms,AngleVarId,angle)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)


status_roms = nf90_close(ncid_roms)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

hmino = 1.0e+23
do j=1,Mp
do i=1,Lp
   hmino = min(hmino,h(i,j))
enddo
enddo

! Vertical grid specification for output
ALLOCATE(sc_ro(No))
ALLOCATE(Cs_ro(No))
ALLOCATE(z_ro(Lp,Mp,No))
ALLOCATE(z_wo(Lp,Mp,No+1))
ALLOCATE(Hzo(Lp,Mp,No))
ALLOCATE(scr1o(Lp,Mp))
ALLOCATE(scr2o(Lp,Mp))
ALLOCATE(scr3d(Lp,Mp,nlev))
ALLOCATE(salt_out(Lp,Mp,No))
ALLOCATE(temp_out(Lp,Mp,No))
ALLOCATE(u_out(L,Mp,No))
ALLOCATE(v_out(Lp,M,No))
ALLOCATE(zeta_out(Lp,Mp))
ALLOCATE(ubar_out(L,Mp))
ALLOCATE(vbar_out(Lp,M))
ALLOCATE(u_rhoo(Lp,Mp,No))
ALLOCATE(v_rhoo(Lp,Mp,No))
ALLOCATE(u_rhoz(Lp,Mp,nlev))
ALLOCATE(v_rhoz(Lp,Mp,nlev))

hco=MIN(hmino,Tclineo)
IF (theta_so.ne.0.0) THEN
  cff1o=1.0/SINH(theta_so)
  cff2o=0.5/TANH(0.5*theta_so)
END IF
cffo=1.0/REAL(No)
DO k=1,No
   sc_ro(k)=cffo*(REAL(k-No)-0.5)
   IF (theta_so.ne.0.0) THEN
       Cs_ro(k)=(1.0-theta_bo)* &
          cff1o*SINH(theta_so*sc_ro(k))+ &
          theta_bo* &
          (cff2o*TANH(theta_so* &
          (sc_ro(k)+0.5))- 0.5)
   ELSE
      Cs_ro(k)=sc_ro(k)
   END IF
END DO
      DO j=1,Mp
      DO i=1,Lp
        DO k=1,No
          cff_ro=hco*(sc_ro(k)-Cs_ro(k))
          cff1_ro=Cs_ro(k)
          cff2_ro=sc_ro(k)+1.0
          z_ro(i,j,k)=cff_ro+cff1_ro*h(i,j)
        END DO
      END DO
      END DO
      DO k=2,No
         z_wo(:,:,k) = 0.5*(z_ro(:,:,k)+z_ro(:,:,k+1))
      ENDDO
      z_wo(:,:,1) = -h(:,:)
      z_wo(:,:,No+1) = 0
      DO k=1,No
         Hzo(:,:,k) = z_wo(:,:,k+1)-z_wo(:,:,k)
      ENDDO

! ..........................
! Open list of input files

write(*,'(a)') filelst
open(unit=55,file=filelst,form='formatted')

ict = 0
feof = 0
DO WHILE (feof==0)

read(55,'(a)',end=999) avgfile

ict = ict+1

avgfilei = TRIM(PHCpath)

nc = LEN_TRIM(avgfilei)
numext = avgfilei(nc-4:nc-3)
outfile = TRIM(outpath)


read(numext,'(i2)') month

tday = 0.5*365.25/12 + (month-1)*365.25/12

write(*,*) 'tday = ',tday

write(*,'(a)') TRIM(avgfilei)
write(*,'(a)') TRIM(outfile)


! ------------------------------------------------------
! Prepare output netCDF file
statuso = nf90_create(trim(outfile),nf90_clobber,ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

!
! Define dimensions
!
statuso = nf90_def_dim(ncido,'xi_rho',Lp,dim_xi_rhoo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_dim(ncido,'eta_rho',Mp,dim_eta_rhoo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_dim(ncido,'xi_u',L,dim_xi_uo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_dim(ncido,'eta_u',Mp,dim_eta_uo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_dim(ncido,'xi_v',Lp,dim_xi_vo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_dim(ncido,'eta_v',M,dim_eta_vo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_dim(ncido,'s_rho',No,dim_s_rhoo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_dim(ncido,'time',0,dim_timeo)
if(statuso /= nf90_NoErr) call handle_err(statuso)

!
! Define variables
!
statuso = nf90_def_var(ncido,'clim_time',nf90_double,dim_timeo,TimeVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_var(ncido,'salt',nf90_float, &
         (/dim_xi_rhoo, dim_eta_rhoo, dim_s_rhoo, dim_timeo/),SaltVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_var(ncido,'temp',nf90_float, &
         (/dim_xi_rhoo, dim_eta_rhoo, dim_s_rhoo, dim_timeo/),TempVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_var(ncido,'u',nf90_float, &
         (/dim_xi_uo, dim_eta_uo, dim_s_rhoo, dim_timeo/),UVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_var(ncido,'v',nf90_float, &
         (/dim_xi_vo, dim_eta_vo, dim_s_rhoo, dim_timeo/),VVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_var(ncido,'zeta',nf90_float, &
         (/dim_xi_rhoo, dim_eta_rhoo, dim_timeo/),ZetaVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_var(ncido,'ubar',nf90_float, &
         (/dim_xi_uo, dim_eta_uo, dim_timeo/),UbarVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_def_var(ncido,'vbar',nf90_float, &
         (/dim_xi_vo, dim_eta_vo, dim_timeo/),VbarVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)

!
! Include variable attributes
!
statuso = nf90_put_att(ncido,SaltVarIdo,'long_name','salinity')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,SaltVarIdo,'units','PSU')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,SaltVarIdo,'time','clim_time')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,SaltVarIdo,'field','salinity, scalar, series')
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_att(ncido,TempVarIdo,'long_name','potential temperature')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,TempVarIdo,'units','Celsius')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,TempVarIdo,'time','clim_time')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,TempVarIdo,'field','temperature, scalar, series')
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_att(ncido,UVarIdo,'long_name','xi velocity component')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,UVarIdo,'units','m sec-1')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,UVarIdo,'time','clim_time')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,UVarIdo,'field','u, scalar, series')
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_att(ncido,VVarIdo,'long_name','eta velocity component')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,VVarIdo,'units','m sec-1')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,VVarIdo,'time','clim_time')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,VVarIdo,'field','v, scalar, series')
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_att(ncido,UbarVarIdo,'long_name','xi barotropic velocity component')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,UbarVarIdo,'units','m sec-1')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,UbarVarIdo,'time','clim_time')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,UbarVarIdo,'field','ubar, scalar, series')
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_att(ncido,VbarVarIdo,'long_name','eta barotropic velocity component')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,VbarVarIdo,'units','m sec-1')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,VbarVarIdo,'time','clim_time')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,VbarVarIdo,'field','vbar, scalar, series')
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_att(ncido,ZetaVarIdo,'long_name','free surface elevation')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,ZetaVarIdo,'units','m')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,ZetaVarIdo,'time','clim_time')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,ZetaVarIdo,'field','zeta, scalar, series')
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_att(ncido,TimeVarIdo,'long_name','time')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,TimeVarIdo,'units','day of year')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,TimeVarIdo,'field','time, scalar, series')
if(statuso /= nf90_NoErr) call handle_err(statuso)
!statuso = nf90_put_att(ncido,TimeVarIdo,'cycle_length',365.25)
!statuso = nf90_put_att(ncido,TimeVarIdo,'cycle_length',365)
if(statuso /= nf90_NoErr) call handle_err(statuso)


statuso = nf90_put_att(ncido,nf90_global,'type','ROMS/TOMS climatology file')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,nf90_global,'title','ROMS 3.0 Global 20km')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,nf90_global,'history','Climatology from SODA')
if(statuso /= nf90_NoErr) call handle_err(statuso)

!
! Exit definition mode
!
statuso = nf90_enddef(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'Finished file definition'

! -----------------------------------------------------------------------

write(*,*) 'writing time = ',tday
! Output to time netCDF file
   statuso = nf90_put_var(ncido,TimeVarIdo,tday)
   if(statuso /= nf90_NoErr) call handle_err(statuso)

! Open PHC file and find variables
! Open input averages file

status_PHC = nf90_open(trim(avgfilei),nf90_nowrite,ncid_PHC)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_PHC, "time",TimePHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_PHC, "depth",LevelPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_PHC, "lon",LonPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_PHC, "lat",LatPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_PHC, "salt",SaltPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_PHC, "temp",TempPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_PHC, "u",UPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_PHC, "v",VPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_PHC, "ssh",ZetaPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_get_var(ncid_PHC,TimePHCVarId,TimePHC)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'got time'
status_PHC = nf90_get_var(ncid_PHC,LonPHCVarId,LonPHC)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'got lon'
status_PHC = nf90_get_var(ncid_PHC,LatPHCVarId,LatPHC)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'got lat'
status_PHC = nf90_get_var(ncid_PHC,LevelPHCVarId,LevelPHC)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'got depth'

!-------------------------------------------------------

IF(ict.eq.1) THEN
write(*,*) 'before get ipos,jpos'

ALLOCATE(ipos(Lp,Mp))
ALLOCATE(jpos(Lp,Mp))

dlon = LonPHC(2)-LonPHC(1)
dlat = LatPHC(2)-LatPHC(1)
lon0 = LonPHC(1)-dlon
lat0 = LatPHC(1)
lonp(2:nlon+1) = LonPHC(:)
lonp(1) =lonp(2)-dlon
lonp(nlonp) = lonp(nlonp-1)+dlon
latp(2:nlat+1) = LatPHC(:)
latp(1) = -90.
latp(nlatp) = 90.
DO j=1,Mp
DO i=1,Lp
   lonr = lon_rho(i,j)
   IF(lonr<0.) lonr = lonr+360.
   ipos(i,j) = (lonr-lon0)/dlon + 1.0
   latr = lat_rho(i,j)
   IF(latr<LatPHC(1)) THEN
      jpos(i,j) = (latr+90.)/(LatPHC(1)+90.) + 1.0
   ELSEIF(latr>LatPHC(nlat)) THEN
      jpos(i,j) = REAL(nlat+1) + (latr-LatPHC(nlat))/(0.75)
   ELSE
      jpos(i,j) = (latr-LatPHC(1))/dlat + 2.0
   ENDIF
ENDDO
ENDDO

! Find subarea in SODA grid containing model grid
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

IF(i1<1.OR.i2>nlon+2.OR.j1<1.OR.j2>nlat+2) THEN
   write(*,*) 'i1,i2,j1,j2 = ',i1,i2,j1,j2
   write(*,*) 'nlon+2, nlat+2 = ',nlon+1,nlat+1
   STOP 'subdomain outside SODA domain'
ENDIF

write(*,*) 'after get ipos,jpos'
ENDIF



! ---------------------------------------------------------------
! Do salinity

   status_PHC = nf90_inq_varid(ncid_PHC,'salt',SaltPHCVarId)
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
   status_PHC = nf90_get_var(ncid_PHC,SaltPHCVarId,SaltPHC, &
 & start=(/ 1, 1, 1, 1/))
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

   where(abs(SaltPHC)>1.e4) SaltPHC=undef

   scr3d = 0.

   DO k=1,nlev

! Fill in masked-out values
      scr = 0.
      scrp = 0.
      scr(1:nlon,1:nlat) = SaltPHC(1:nlon,1:nlat,k)
      sumb = 0.
      sumt = 0.
      DO i=1,nlon
         sumb = sumb + scr(i,1)
         sumt = sumt + scr(i,nlat)
      ENDDO
      sumb = sumb/REAL(nlon)
      sumt = sumt/REAL(nlon)
      scrp(1:nlon+2,1) = sumb
      scrp(1:nlon+2,nlat+2) = sumt
      scrp(2:nlon+1,2:nlat+1) = scr(1:nlon,1:nlat)
      scrp(1,2:nlat+1) = scr(nlon,1:nlat)
      scrp(nlon+2,2:nlat+1) = scr(1,1:nlat)
      WHERE(ABS(scrp)>1.e+4) scrp = undef
      CALL fill(nlonp,nlatp,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work,error,nvalue)
      DO j=1,Mp
      DO i=1,Lp
         ia = ifix(ipos(i,j))
         ja = ifix(jpos(i,j))
         rx = (ipos(i,j)-ia)
         ry = (jpos(i,j)-ja)
         rxm = 1.0-rx
         rym = 1.0-ry
         scr3d(i,j,k) = rxm*rym*scrp(ia,ja) + &
                       rx*rym*scrp(ia+1,ja) + &
                       rx*ry*scrp(ia+1,ja+1) + &
                       rxm*ry*scrp(ia,ja+1)
      ENDDO
      ENDDO
   ENDDO

   DO k=2,nlev
      WHERE(ABS(scr3d(:,:,k))>1.e4) scr3d(:,:,k) = scr3d(:,:,k-1)
   ENDDO

! Vertical interpolation
   DO j=1,Mp
   DO i=1,Lp
      DO k=1,No
         IF(-z_ro(i,j,k).LT.LevelPHC(1)) THEN
            salt_out(i,j,k) = scr3d(i,j,1)
         ELSEIF(-z_ro(i,j,k).GT.LevelPHC(nlev)) THEN
            salt_out(i,j,k) = scr3d(i,j,nlev)
         ELSE
            DO kT=1,nlev
               IF(-z_ro(i,j,k).LT.LevelPHC(kT+1).AND. &
                  -z_ro(i,j,k).GE.LevelPHC(kT)) THEN
                  rz2 = (-z_ro(i,j,k)-LevelPHC(kT))/ &
                          (LevelPHC(kT+1)-LevelPHC(kT))
                  rz1 = 1.0-rz2
                  salt_out(i,j,k) = rz1*scr3d(i,j,kT) + rz2*scr3d(i,j,kT+1)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   ENDDO

! Output to netCDF file
   statuso = nf90_put_var(ncido,SaltVarIdo,salt_out,(/1,1,1,1/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)


! ---------------------------------------------------------------
! Do temperature

   status_PHC = nf90_inq_varid(ncid_PHC,'temp',TempPHCVarId)
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
   status_PHC = nf90_get_var(ncid_PHC,TempPHCVarId,TempPHC, &
 & start=(/ 1, 1, 1, 1/))
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

   where(abs(TempPHC)>1.e4) TempPHC=undef

   scr3d = 0.
   DO k=1,nlev

! Fill in masked-out values
      scr = 0.
      scrp = 0.
      scr(1:nlon,1:nlat) = TempPHC(1:nlon,1:nlat,k)
      sumb = 0.
      sumt = 0.
      DO i=1,nlon
         sumb = sumb + scr(i,1)
         sumt = sumt + scr(i,nlat)
      ENDDO
      sumb = sumb/REAL(nlon)
      sumt = sumt/REAL(nlon)
      scrp(1:nlon+2,1) = sumb
      scrp(1:nlon+2,nlat+2) = sumt
      scrp(2:nlon+1,2:nlat+1) = scr(1:nlon,1:nlat)
      scrp(1,2:nlat+1) = scr(nlon,1:nlat)
      scrp(nlon+2,2:nlat+1) = scr(1,1:nlat)
      WHERE(ABS(scrp)>1.e+4) scrp = undef
      CALL fill(nlonp,nlatp,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work,error,nvalue)

      DO j=1,Mp
      DO i=1,Lp
         ia = ifix(ipos(i,j))
         ja = ifix(jpos(i,j))
         rx = (ipos(i,j)-ia)
         ry = (jpos(i,j)-ja)
         rxm = 1.0-rx
         rym = 1.0-ry
         scr3d(i,j,k) = rxm*rym*scrp(ia,ja) + &
                       rx*rym*scrp(ia+1,ja) + &
                       rx*ry*scrp(ia+1,ja+1) + &
                       rxm*ry*scrp(ia,ja+1)
      ENDDO
      ENDDO
   ENDDO

   DO k=2,nlev
      WHERE(ABS(scr3d(:,:,k))>1.e4) scr3d(:,:,k) = scr3d(:,:,k-1)
   ENDDO

! Vertical interpolation
   DO j=1,Mp
   DO i=1,Lp
      DO k=1,No
         IF(-z_ro(i,j,k).LT.LevelPHC(1)) THEN
            temp_out(i,j,k) = scr3d(i,j,1)
         ELSEIF(-z_ro(i,j,k).GT.LevelPHC(nlev)) THEN
            temp_out(i,j,k) = scr3d(i,j,nlev)
         ELSE
            DO kT=1,nlev
               IF(-z_ro(i,j,k).LT.LevelPHC(kT+1).AND. &
                  -z_ro(i,j,k).GE.LevelPHC(kT)) THEN
                  rz2 = (-z_ro(i,j,k)-LevelPHC(kT))/ &
                          (LevelPHC(kT+1)-LevelPHC(kT))
                  rz1 = 1.0-rz2
                  temp_out(i,j,k) = rz1*scr3d(i,j,kT) + rz2*scr3d(i,j,kT+1)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   ENDDO

! Output to netCDF file
   statuso = nf90_put_var(ncido,TempVarIdo,temp_out,(/1,1,1,1/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

! ---------------------------------------------------------------
! Do U, V velocity components

! Get U
   status_PHC = nf90_inq_varid(ncid_PHC,'u',UPHCVarId)
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
   status_PHC = nf90_get_var(ncid_PHC,UPHCVarId,UPHC, &
 & start=(/ 1, 1, 1, 1/))
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
! Get V
   status_PHC = nf90_inq_varid(ncid_PHC,'v',VPHCVarId)
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
   status_PHC = nf90_get_var(ncid_PHC,VPHCVarId,VPHC, &
 & start=(/ 1, 1, 1, 1/))
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

   where(abs(UPHC)>1.e4) UPHC=undef
   where(abs(VPHC)>1.e4) VPHC=undef

   u_rhoz = 0.
   v_rhoz = 0.

   DO k=1,nlev
! Rotate components to polar stereographic
    DO i=1,nlon
      scra(i,:) = UPHC(i,:,k)*COS(DTOR*LonPHC(i)) - VPHC(i,:,k)*SIN(DTOR*LonPHC(i))
      scrb(i,:) = VPHC(i,:,k)*COS(DTOR*LonPHC(i)) + UPHC(i,:,k)*SIN(DTOR*LonPHC(i))
    ENDDO
! Fill in masked-out values
      scr = 0.
      scrp = 0.
      scr(1:nlon,1:nlat) = scra(1:nlon,1:nlat)
      sumb = 0.
      sumt = 0.
      DO i=1,nlon
         sumb = sumb + scr(i,1)
         sumt = sumt + scr(i,nlat)
      ENDDO
      sumb = sumb/REAL(nlon)
      sumt = sumt/REAL(nlon)
      scrp(1:nlon+2,1) = sumb
      scrp(1:nlon+2,nlat+2) = sumt
      scrp(2:nlon+1,2:nlat+1) = scr(1:nlon,1:nlat)
      scrp(1,2:nlat+1) = scr(nlon,1:nlat)
      scrp(nlon+2,2:nlat+1) = scr(1,1:nlat)
      scrup = scrp

      scr = 0.
      scrp = 0.
      scr(1:nlon,1:nlat) = scrb(1:nlon,1:nlat)
      sumb = 0.
      sumt = 0.
      DO i=1,nlon
         sumb = sumb + scr(i,1)
         sumt = sumt + scr(i,nlat)
      ENDDO
      sumb = sumb/REAL(nlon)
      sumt = sumt/REAL(nlon)
      scrp(1:nlon+2,1) = sumb
      scrp(1:nlon+2,nlat+2) = sumt
      scrp(2:nlon+1,2:nlat+1) = scr(1:nlon,1:nlat)
      scrp(1,2:nlat+1) = scr(nlon,1:nlat)
      scrp(nlon+2,2:nlat+1) = scr(1,1:nlat)
      scrvp = scrp

    scrp = scrup
    WHERE(ABS(scrp)>1.e+4) scrp = undef
    CALL fill(nlonp,nlatp,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work,error,nvalue)
    scrup = scrp

    scrp = 0.
    scrp = scrvp
    WHERE(ABS(scrp)>1.e+4) scrp = undef
    CALL fill(nlonp,nlatp,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work,error,nvalue)
    scrvp = scrp


    DO i=1,nlonp
      scrap(i,:) = scrup(i,:)*COS(DTOR*lonp(i)) + scrvp(i,:)*SIN(DTOR*lonp(i))
      scrbp(i,:) = scrvp(i,:)*COS(DTOR*lonp(i)) - scrup(i,:)*SIN(DTOR*lonp(i))
    ENDDO
! Horizontal interpolation
      DO j=1,Mp
      DO i=1,Lp
         ia = ifix(ipos(i,j))
         ja = ifix(jpos(i,j))
         rx = (ipos(i,j)-ia)
         ry = (jpos(i,j)-ja)
         rxm = 1.0-rx
         rym = 1.0-ry
         u_rhoz(i,j,k) = rxm*rym*scrap(ia,ja) + &
                       rx*rym*scrap(ia+1,ja) + &
                       rx*ry*scrap(ia+1,ja+1) + &
                       rxm*ry*scrap(ia,ja+1)
      ENDDO
      ENDDO
      DO j=1,Mp
      DO i=1,Lp
         ia = ifix(ipos(i,j))
         ja = ifix(jpos(i,j))
         rx = (ipos(i,j)-ia)
         ry = (jpos(i,j)-ja)
         rxm = 1.0-rx
         rym = 1.0-ry
         v_rhoz(i,j,k) = rxm*rym*scrbp(ia,ja) + &
                       rx*rym*scrbp(ia+1,ja) + &
                       rx*ry*scrbp(ia+1,ja+1) + &
                       rxm*ry*scrbp(ia,ja+1)
      ENDDO
      ENDDO
   ENDDO
   DO k=2,nlev
      WHERE(ABS(u_rhoz(:,:,k))>1.e4) u_rhoz(:,:,k) = u_rhoz(:,:,k-1)
      WHERE(ABS(v_rhoz(:,:,k))>1.e4) v_rhoz(:,:,k) = v_rhoz(:,:,k-1)
   ENDDO
! Vertical interpolation
   DO j=1,Mp
   DO i=1,Lp
      DO k=1,No
         IF(-z_ro(i,j,k).LT.LevelPHC(1)) THEN
            scr1o(i,j) = u_rhoz(i,j,1)
            scr2o(i,j) = v_rhoz(i,j,1)
         ELSEIF(-z_ro(i,j,k).GT.LevelPHC(nlev)) THEN
            scr1o(i,j) = u_rhoz(i,j,nlev)
            scr2o(i,j) = v_rhoz(i,j,nlev)
         ELSE
            DO kT=1,nlev
               IF(-z_ro(i,j,k).LT.LevelPHC(kT+1).AND. &
                  -z_ro(i,j,k).GE.LevelPHC(kT)) THEN
                  rz2 = (-z_ro(i,j,k)-LevelPHC(kT))/ &
                          (LevelPHC(kT+1)-LevelPHC(kT))
                  rz1 = 1.0-rz2
                  scr1o(i,j) = rz1*u_rhoz(i,j,kT) + rz2*u_rhoz(i,j,kT+1)
                  scr2o(i,j) = rz1*v_rhoz(i,j,kT) + rz2*v_rhoz(i,j,kT+1)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         u_rhoo(i,j,k) = scr1o(i,j)*cos(angle(i,j)) + scr2o(i,j)*sin(angle(i,j))
         v_rhoo(i,j,k) = scr2o(i,j)*cos(angle(i,j)) - scr1o(i,j)*sin(angle(i,j))
      ENDDO
   ENDDO
   ENDDO
! Put on staggered grid and convert cm/s to m/s
   DO k=1,No
   DO j=1,Mp
   DO i=1,L
      u_out(i,j,k) = 0.005*(u_rhoo(i,j,k) + u_rhoo(i+1,j,k))
   ENDDO
   ENDDO
   ENDDO
   DO k=1,No
   DO j=1,M
   DO i=1,Lp
      v_out(i,j,k) = 0.005*(v_rhoo(i,j,k) + v_rhoo(i,j+1,k))
   ENDDO
   ENDDO
   ENDDO
! Output to netCDF file
   statuso = nf90_put_var(ncido,UVarIdo,u_out,(/1,1,1,1/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncido,VVarIdo,v_out,(/1,1,1,1/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
! ---------------------------------------------------------------
! Do free surface elevation
   status_PHC = nf90_inq_varid(ncid_PHC,'ssh',ZetaPHCVarId)
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
   status_PHC = nf90_get_var(ncid_PHC,ZetaPHCVarId,ZetaPHC, &
 & start=(/ 1, 1, 1/))
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
   where(abs(ZetaPHC)>1.e4) ZetaPHC=undef
! Fill in masked-out values
      scr = 0.
      scrp = 0.
      scr(1:nlon,1:nlat) = ZetaPHC(1:nlon,1:nlat)
      sumb = 0.
      sumt = 0.
      DO i=1,nlon
         sumb = sumb + scr(i,1)
         sumt = sumt + scr(i,nlat)
      ENDDO
      sumb = sumb/REAL(nlon)
      sumt = sumt/REAL(nlon)
      scrp(1:nlon+2,1) = sumb
      scrp(1:nlon+2,nlat+2) = sumt
      scrp(2:nlon+1,2:nlat+1) = scr(1:nlon,1:nlat)
      scrp(1,2:nlat+1) = scr(nlon,1:nlat)
      scrp(nlon+2,2:nlat+1) = scr(1,1:nlat)
      WHERE(abs(scrp)>1.e+4) scrp = undef
      CALL fill(nlonp,nlatp,i1,i2,j1,j2,scrp,tx,critx,cor,mxs,work,error,nvalue)
      DO j=1,Mp
      DO i=1,Lp
         ia = ifix(ipos(i,j))
         ja = ifix(jpos(i,j))
         rx = (ipos(i,j)-ia)
         ry = (jpos(i,j)-ja)
         rxm = 1.0-rx
         rym = 1.0-ry
         zeta_out(i,j) = rxm*rym*scrp(ia,ja) + &
                       rx*rym*scrp(ia+1,ja) + &
                       rx*ry*scrp(ia+1,ja+1) + &
                       rxm*ry*scrp(ia,ja+1)
      ENDDO
      ENDDO
! Convert from cm to m
      zeta_out = 0.01*zeta_out
! Output to netCDF file
   statuso = nf90_put_var(ncido,ZetaVarIdo,zeta_out,(/1,1,1/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
! ---------------------------------------------------------------
! Do barotropic velocities
      ubar_out = 0.
      vbar_out = 0.
      DO k=1,No
         ubar_out(:,:) = ubar_out(:,:) + u_out(:,:,k)*Hzo(:,:,k)
         vbar_out(:,:) = vbar_out(:,:) + v_out(:,:,k)*Hzo(:,:,k)
      ENDDO
      ubar_out = ubar_out/h
      vbar_out = vbar_out/h
! Output to netCDF file
   statuso = nf90_put_var(ncido,UbarVarIdo,ubar_out,(/1,1,1/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncido,VbarVarIdo,vbar_out,(/1,1,1/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
status_PHC = nf90_close(ncid_PHC)
if(status_roms /= nf90_NoErr) call handle_err(status_PHC)
statuso = nf90_close(ncido)
if(status_roms /= nf90_NoErr) call handle_err(statuso)
ENDDO
999 CONTINUE
END PROGRAM SODA_to_climbc
