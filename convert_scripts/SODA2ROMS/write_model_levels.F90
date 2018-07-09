PROGRAM write_model_levels
! Program to obtain climatology from SODA monthly T,S,U,V,zeta climatology

#undef  PST_OUT
#define GEO_OUT

use netcdf

implicit none

INTEGER, PARAMETER :: nlev = 40, nlat = 330, nlon = 720
INTEGER :: nlonp, nlatp,ounit
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

INTEGER refyear_new, refmonth_new, refday_new
REAL refdate1,refdate2,curdate,newdate
CHARACTER*30 datestring
CHARACTER*4 refyear_char
CHARACTER*2 refmonth_char,refday_char

INTEGER i, j, k, idate, ict
INTEGER No, Vstretching, Vtransform
REAL theta_so, theta_bo, Tclineo
REAL Aweight, Bweight, Cweight, ds, ssc_w, ssc_r, Csur, Cbot
REAL cff_wo, cff1_wo, cff2_wo, hinv, hwater
CHARACTER*160 gridfile, PHCfile, outfile,PHCpath, outpath, filelst
CHARACTER*160 clim_root, avgfile, avgfilei

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
REAL, ALLOCATABLE, DIMENSION(:) :: Cs_ro, sc_ro, Cs_wo, sc_wo
REAL cffo, cff1o, cff2o, hco, hmino, cff_ro, cff1_ro, cff2_ro


INTEGER ncid_PHC, status_PHC
INTEGER TimePHCVarId, LevelPHCVarId, LonPHCVarId, LatPHCVarId
INTEGER UPHCVarId, VPHCVarId, ZetaPHCVarId, SaltPHCVarId, TempPHCVarId
INTEGER statuso, ncido, dim_xi_rhoo, dim_eta_rhoo, dim_timeo, dim_s_rhoo
INTEGER dim_xi_uo, dim_eta_uo, dim_xi_vo, dim_eta_vo 
INTEGER UVarIdo, VVarIdo, TimeVarIdo, SaltVarIdo, TempVarIdo, MLVarIdo, HVarIdo
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
REAL, PARAMETER    :: undef = 2.E+35            ! Undefined land value
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

write(*,*) 'reading input file'
read(*,'(a)') gridfile
write(*,*) gridfile
read(*,'(a)') outfile
read(*,*) Vtransform, Vstretching
read(*,*) No, theta_so, theta_bo, Tclineo

write(*,*) Vtransform, Vstretching, theta_so, theta_bo

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

write(*,*) 'read gridfile'
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

write(*,*) 'finished reading gridfile'

hmino = 1.0e+23
do j=1,Mp
do i=1,Lp
   hmino = min(hmino,h(i,j))
enddo
enddo

! Vertical grid specification for output
ALLOCATE(sc_ro(No))
ALLOCATE(Cs_ro(No))
ALLOCATE(sc_wo(No+1))
ALLOCATE(Cs_wo(No+1))
ALLOCATE(z_ro(Lp,Mp,No))
ALLOCATE(z_wo(Lp,Mp,No))
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

IF(Vtransform.EQ.1.AND.Vstretching.EQ.1) THEN
hco=MIN(hmino,Tclineo)
IF (theta_so.ne.0.0) THEN
  cff1o=1.0/SINH(theta_so)
  cff2o=0.5/TANH(0.5*theta_so)
END IF
cffo=1.0/REAL(No)
DO k=1,No
   sc_ro(k)=cffo*(REAL(k-No)-0.5)
   IF (theta_so.ne.0.0) THEN
       Cs_ro(k)=(1.0-theta_bo)*                  &
          cff1o*SINH(theta_so*sc_ro(k))+         &
          theta_bo*                              &
          (cff2o*TANH(theta_so*                  &
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

ELSEIF(Vtransform.EQ.2.AND.Vstretching.EQ.2) THEN
!                                                                                              
        Aweight=1.0
        Bweight=1.0
        ds=1.0/REAL(No)
!                                                                                              
        sc_wo(No+1)=0.0
        Cs_wo(No+1)=0.0
        DO k=No,2,-1
          ssc_w=ds*REAL(k-No-1)
          sc_wo(k)=ssc_w
          IF (theta_so.gt.0.0) THEN
            Csur=(1.0-COSH(theta_so*ssc_w))/                       &
     &           (COSH(theta_so)-1.0)
            IF (theta_bo.gt.0.0) THEN
              Cbot=SINH(theta_bo*(ssc_w+1.0))/                     &
     &             SINH(theta_bo)-1.0
              Cweight=(ssc_w+1.0)**Aweight*                           &
     &                (1.0+(Aweight/Bweight)*                        &
     &                        (1.0-(ssc_w+1.0)**Bweight))
              Cs_wo(k)=Cweight*Csur+(1.0-Cweight)*Cbot
            ELSE
              Cs_wo(k)=Csur
            END IF
          ELSE
            Cs_wo(k)=ssc_w
          END IF
        END DO
        sc_wo(1)=-1.0
        Cs_wo(1)=-1.0
!                                                                                              
        DO k=1,No
          ssc_r=ds*(REAL(k-No)-0.5)
          sc_ro(k)=ssc_r
          IF (theta_so.gt.0.0) THEN
            Csur=(1.0-COSH(theta_so*ssc_r))/                       &
     &           (COSH(theta_so)-1.0)
            IF (theta_bo.gt.0.0) THEN
              Cbot=SINH(theta_bo*(ssc_r+1.0))/                     &
     &             SINH(theta_bo)-1.0
              Cweight=(ssc_r+1.0)**Aweight*                           &
     &                (1.0+(Aweight/Bweight)*                        &
     &                        (1.0-(ssc_r+1.0)**Bweight))
              Cs_ro(k)=Cweight*Csur+(1.0-Cweight)*Cbot
            ELSE
              Cs_ro(k)=Csur
            END IF
          ELSE
            Cs_ro(k)=ssc_r
          END IF
        END DO

      hco = Tclineo
      DO j=1,Mp
      DO i=1,Lp
          z_wo(i,j,1)=-h(i,j)
      ENDDO
      ENDDO
      DO j=1,Mp
         DO i=1,Lp
            DO k=1,No-1
               cff_ro=hco*sc_ro(k)
               cff_wo=hco*sc_wo(k+1)
               cff1_ro=Cs_ro(k)
               cff1_wo=Cs_wo(k+1)
              hwater=h(i,j)
              hinv=1.0/(hco+hwater)
              cff2_ro=(cff_ro+cff1_ro*hwater)*hinv
              cff2_wo=(cff_wo+cff1_wo*hwater)*hinv
              z_wo(i,j,k+1)=hwater*cff2_wo
             END DO
          END DO
        END DO
!
ELSEIF(Vtransform.EQ.2.AND.Vstretching.EQ.3) THEN
!
!  This stretching function is intended for very shallow coastal
!  applications, like gravity sediment flows.
!
!  At the surface, C(s=0)=0
!
!      C(s) = - LOG(COSH(Hscale * ABS(s) ** alpha)) /
!               LOG(COSH(Hscale))
!
!  At the bottom, C(s=-1)=-1
!
!      C(s) = LOG(COSH(Hscale * (s + 1) ** beta)) /
!             LOG(COSH(Hscale)) - 1
!
!  where
!
!       Hscale : scale value for all hyperbolic functions
!                  Hscale = 3.0    set internally here
!        alpha : surface stretching exponent
!                  alpha = 0.65   minimal increase of surface resolution
!                          1.0    significant amplification
!         beta : bottoom stretching exponent
!                  beta  = 0.58   no amplification
!                          1.0    significant amplification
!                          3.0    super-high bottom resolution
!            s : stretched vertical coordinate, -1 <= s <= 0
!                  s(k) = (k-N)/N       k=0:N,    W-points  (s_w)
!                  s(k) = (k-N-0.5)/N   k=1:N,  RHO-points  (s_rho)
!
!  The stretching exponents (alpha, beta) are specify at input:
!
!         alpha = theta_s
!         beta  = theta_b
!
   ds=1.0/REAL(No)
!
   sc_wo(No+1)=0.0
   Cs_wo(No+1)=0.0
   DO k=No,2,-1
      sc_wo(k)=ds*REAL(k-No-1)
      Cbot= LOG(COSH(3.0*(sc_wo(k)+1.0)**theta_bo))/               &
           LOG(COSH(3.0))-1.0
      Csur=-LOG(COSH(3.0*ABS(sc_wo(k))**theta_so))/                   &
           LOG(COSH(3.0))
      Cweight=0.5*(1.0-TANH(3.0*(sc_wo(k)+0.5)))
      Cs_wo(k)=Cweight*Cbot+(1.0-Cweight)*Csur
   END DO
   sc_wo(1)=-1.0
   Cs_wo(1)=-1.0
!
   DO k=1,No
      sc_ro(k)=ds*(REAL(k-No)-0.5)
      Cbot= LOG(COSH(3.0*(sc_ro(k)+1.0)**theta_bo))/               &
           LOG(COSH(3.0))-1.0
      Csur=-LOG(COSH(3.0*ABS(sc_ro(k))**theta_so))/                   &
           LOG(COSH(3.0))
      Cweight=0.5*(1.0-TANH(3.0*(sc_ro(k)+0.5)))
      Cs_ro(k)=Cweight*Cbot+(1.0-Cweight)*Csur
   END DO
! 
  hco = Tclineo
   DO j=1,Mp
      DO i=1,Lp
         z_wo(i,j,1)=-h(i,j)
      ENDDO
   ENDDO
   DO j=1,Mp
      DO i=1,Lp
         DO k=1,No-1
            cff_ro=hco*sc_ro(k)
            cff_wo=hco*sc_wo(k+1)
            cff1_ro=Cs_ro(k)
            cff1_wo=Cs_wo(k+1)
            hwater=h(i,j)
            hinv=1.0/(hco+hwater)
            cff2_ro=(cff_ro+cff1_ro*hwater)*hinv
            cff2_wo=(cff_wo+cff1_wo*hwater)*hinv
            z_wo(i,j,k+1)=hwater*cff2_wo
          END DO
      END DO
   END DO
!
ELSEIF(Vtransform.EQ.2.AND.Vstretching.EQ.4) THEN
!-----------------------------------------------------------------------
!  A. Shchepetkin improved double vertical stretching functions with
!  bottom refiment.
!-----------------------------------------------------------------------
!
!  The range of meaningful values for the control parameters are:
!
!       0 <  theta_s <= 10.0
!       0 <= theta_b <=  3.0
!
!  Users need to pay attention to extreme r-factor (rx1) values near
!  the bottom.
!
!  This vertical stretching function is defined, in the simplest form,
!  as:
!
!      C(s) = [1.0 - COSH(theta_s * s)] / [COSH(theta_s) - 1.0]
!
!  it is similar in meaning to the original vertical stretcing function
!  (Song and Haidvogel, 1994), but note that hyperbolic functions are
!  COSH, and not SINH.
!
!  Note that the above definition results in
!
!         -1 <= C(s) <= 0
!
!  as long as
!
!         -1 <= s <= 0
!
!  and
!
!         d[C(s)]/ds  -->  0      if  s -->  0
!
!  For the purpose of bottom boundary layer C(s) is further modified
!  to allow near-bottom refinement by using a continuous, second
!  stretching function
!
!         C(s) = [EXP(theta_b * C(s)) - 1.0] / [1.0 - EXP(-theta_b)]
!
!  This double transformation is continuous with respect to "theta_s"
!  and "theta_b", as both values approach to zero.
!
   ds=1.0/REAL(No)
!
   sc_wo(No+1)=0.0
   Cs_wo(No+1)=0.0
   DO k=No,2,-1
      sc_wo(k)=ds*REAL(k-No-1)
      IF (theta_so.gt.0.0) THEN
         Csur=(1.0-COSH(theta_so*sc_wo(k)))/ &
              (COSH(theta_so)-1.0)
      ELSE
         Csur=-sc_wo(k)**2
      END IF
      IF (theta_bo.gt.0.0) THEN
         Cbot=(EXP(theta_bo*Csur)-1.0)/      &
              (1.0-EXP(-theta_bo))
         Cs_wo(k)=Cbot
      ELSE
         Cs_wo(k)=Csur
      END IF
   END DO
   sc_wo(1)=-1.0
   Cs_wo(1)=-1.0
!
   DO k=1,No
      sc_ro(k)=ds*(REAL(k-No)-0.5)
      IF (theta_so.gt.0.0) THEN
         Csur=(1.0-COSH(theta_so*sc_ro(k)))/                       &
              (COSH(theta_so)-1.0)
      ELSE
         Csur=-sc_ro(k)**2
      END IF
      IF (theta_bo.gt.0.0) THEN
         Cbot=(EXP(theta_bo*Csur)-1.0)/                        &
              (1.0-EXP(-theta_bo))
         Cs_ro(k)=Cbot
      ELSE
         Cs_ro(k)=Csur
      END IF
   END DO
! 
  hco = Tclineo
   DO j=1,Mp
      DO i=1,Lp
         z_wo(i,j,1)=-h(i,j)
      ENDDO
   ENDDO
   DO j=1,Mp
      DO i=1,Lp
         DO k=1,No-1
            cff_ro=hco*sc_ro(k)
            cff_wo=hco*sc_wo(k+1)
            cff1_ro=Cs_ro(k)
            cff1_wo=Cs_wo(k+1)
            hwater=h(i,j)
            hinv=1.0/(hco+hwater)
            cff2_ro=(cff_ro+cff1_ro*hwater)*hinv
            cff2_wo=(cff_wo+cff1_wo*hwater)*hinv
            z_wo(i,j,k+1)=hwater*cff2_wo
          END DO
      END DO
   END DO 
   DO k = 1,No
      write(*,*) z_wo(:,1,k)
   END DO

ELSE
   write(*,*) 'Not programmed these vertical grid options.'
   STOP
ENDIF

! ------------------------------------------------------
! Prepare output netCDF file
write(*,*) 'creating outfile'
statuso = nf90_create(trim(outfile),nf90_clobber,ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'outfile created'
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

write(*,*) 'dims defined'

!
! Define variables
!
statuso = nf90_def_var(ncido,'level_depths',nf90_float,                      &
         (/dim_xi_rhoo, dim_eta_rhoo, dim_s_rhoo/),MLVarIdo)
statuso = nf90_def_var(ncido,'h',nf90_float,                      &
         (/dim_xi_rhoo, dim_eta_rhoo/),HVarIdo)
write(*,*) 'var defined'

!
! Include variable attributes
!
statuso = nf90_put_att(ncido,MLVarIdo,'long_name','model level depths')
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_att(ncido,MLVarIdo,'units','m')
if(statuso /= nf90_NoErr) call handle_err(statuso)

!
! Exit definition mode
!
statuso = nf90_enddef(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'Finished file definition'

! -----------------------------------------------------------------------

! Output to time netCDF file

   statuso = nf90_put_var(ncido,MLVarIdo,z_wo)
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncido,HVarIdo,h)
   if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_close(ncido)
if(status_roms /= nf90_NoErr) call handle_err(statuso)

999 CONTINUE

END PROGRAM write_model_levels
