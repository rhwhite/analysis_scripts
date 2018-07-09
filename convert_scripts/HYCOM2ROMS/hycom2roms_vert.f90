PROGRAM hycom2roms_vert
! Program to obtain vertical interpolation of T,S,u,v from z to s levels on rho grid

use netcdf
USE, INTRINSIC :: IEEE_ARITHMETIC
implicit none


INTEGER i, j, k,t,Kt
INTEGER No, Vtransform, Vstretching
REAL theta_so, theta_bo, Tclineo
REAL Aweight, Bweight, Cweight, ds, ssc_w, ssc_r, Csur, Cbot
REAL cff_wo, cff1_wo, cff2_wo, hinv, hwater
CHARACTER*160 gridfile, inputfile, outputfile
CHARACTER*160 PHCfile,num_temp,line_temp,unit_time_att

INTEGER XiDimID, EtaDimID, Lp, Mp, L, M,len_num
INTEGER XiRhoDimId, EtaRhoDimId, DepthDimId, TimeDimId
INTEGER LonRhoVarId, LatRhoVarId, HVarId, AngleVarId
INTEGER status_roms, ncid_roms

CHARACTER*160, DIMENSION(58) :: line_clim_file


REAL, ALLOCATABLE, DIMENSION(:)     :: LevelPHC
REAL, ALLOCATABLE, DIMENSION(:,:)   :: lon_rho
REAL, ALLOCATABLE, DIMENSION(:,:)   :: lat_rho
REAL, ALLOCATABLE, DIMENSION(:,:)   :: h
REAL, ALLOCATABLE, DIMENSION(:,:)   :: h_u
REAL, ALLOCATABLE, DIMENSION(:,:)   :: h_v
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: z_ro
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: z_ro_u
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: z_ro_v
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: z_wo
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Hzo
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Hzo_u
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: Hzo_v
REAL, ALLOCATABLE, DIMENSION(:)     :: Cs_ro, sc_ro, Cs_wo, sc_wo
REAL cffo, cff1o, cff2o, hco, hmino, cff_ro, cff1_ro, cff2_ro
REAL rz1, rz2

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SaltPHC
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: TempPHC
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: UPHC
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: VPHC
REAL, ALLOCATABLE, DIMENSION(:,:)   :: SSHPHC
DOUBLE PRECISION TimePHC


INTEGER nlev,NT
INTEGER ncid_PHC, status_PHC
INTEGER SSHPHCVarId,TimePHCVarId, LevelPHCVarId
INTEGER UPHCVarId, VPHCVarId, ZetaPHCVarId, SaltPHCVarId, TempPHCVarId,TimeVarId
INTEGER statuso, ncido
INTEGER UVarIdo, VVarIdo, TimeVarIdo, SaltVarIdo, TempVarIdo
INTEGER UbarVarIdo, VbarVarIdo, ZetaVarIdo
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: salt_out
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp_out
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: u_out
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: v_out
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ubar_out
REAL, ALLOCATABLE, DIMENSION(:,:)   :: vbar_out
REAL, ALLOCATABLE, DIMENSION(:,:)   :: zeta_out

OPEN(unit=10,file='hycom2roms_vert.in')
read(10,'(a)') gridfile
read(10,'(a)') inputfile
read(10,'(a)') outputfile
read(10,*) No, theta_so, theta_bo, Tclineo
read(10,*) Vtransform, Vstretching
CLOSE(10)


write(*,*) 'gridfile is ',gridfile
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
ALLOCATE(h(Lp,Mp))
ALLOCATE(h_u(Lp-1,Mp))
ALLOCATE(h_v(Lp,Mp-1))

write(*,*) 'Grid is Lp=',Lp,' Mp=',Mp,' N=',No
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


status_roms = nf90_close(ncid_roms)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)


! Get info on ROMS climatology on z-levels
status_PHC = nf90_open(TRIM(inputfile),nf90_nowrite,ncid_PHC)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_dimid(ncid_PHC, "depth", DepthDimId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_dimid(ncid_roms, "time", TimeDimId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_roms, "time",TimeVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

status_PHC = nf90_Inquire_Dimension(ncid_PHC,DepthDimID,len = nlev)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_Inquire_Dimension(ncid_PHC,TimeDimID,len = NT)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)
status_PHC = nf90_get_att(ncid_PHC,TimeVarId,'units',unit_time_att)

if(status_PHC /= nf90_NoErr) call handle_err(status_roms)

ALLOCATE(LevelPHC(nlev))

write(*,*) 'Number of z-levels nlev=',nlev,' and number of time NT=',NT

status_PHC = nf90_close(ncid_PHC)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)


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
ALLOCATE(z_ro_u(Lp-1,Mp,No))
ALLOCATE(z_ro_v(Lp,Mp-1,No))

ALLOCATE(z_wo(Lp,Mp,No+1))
ALLOCATE(Hzo(Lp,Mp,No))
ALLOCATE(Hzo_u(Lp-1,Mp,No))
ALLOCATE(Hzo_v(Lp,Mp-1,No))

ALLOCATE(salt_out(Lp,Mp,No))
ALLOCATE(temp_out(Lp,Mp,No))
ALLOCATE(u_out(Lp-1,Mp,No))
ALLOCATE(v_out(Lp,Mp-1,No))
ALLOCATE(ubar_out(Lp-1,Mp))
ALLOCATE(vbar_out(Lp,Mp-1))
ALLOCATE(zeta_out(Lp,Mp))

ALLOCATE(SaltPHC(Lp,Mp,nlev))
ALLOCATE(TempPHC(Lp,Mp,nlev))
ALLOCATE(UPHC(Lp-1,Mp,nlev))
ALLOCATE(VPHC(Lp,Mp-1,nlev))
ALLOCATE(SSHPHC(Lp,Mp))


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
      z_wo(:,:,k) = 0.5*(z_ro(:,:,k-1)+z_ro(:,:,k))
   ENDDO
   z_wo(:,:,1) = -h(:,:)
   z_wo(:,:,No+1) = 0.
   DO k=1,No
      Hzo(:,:,k) = z_wo(:,:,k+1)-z_wo(:,:,k)
   ENDDO
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
              (COSH(theta_so)-1.0)
         IF (theta_bo.gt.0.0) THEN
            Cbot=SINH(theta_bo*(ssc_r+1.0))/                    &
                 SINH(theta_bo)-1.0
            Cweight=(ssc_r+1.0)**Aweight*                       &
                 (1.0+(Aweight/Bweight)*                        &
                 (1.0-(ssc_r+1.0)**Bweight))
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
         DO k=1,No
            cff_ro=hco*sc_ro(k)
            cff_wo=hco*sc_wo(k+1)
            cff1_ro=Cs_ro(k)
            cff1_wo=Cs_wo(k+1)
            hwater=h(i,j)
            hinv=1.0/(hco+hwater)
            cff2_ro=(cff_ro+cff1_ro*hwater)*hinv
            cff2_wo=(cff_wo+cff1_wo*hwater)*hinv
            z_wo(i,j,k+1)=hwater*cff2_wo
            z_ro(i,j,k)=hwater*cff2_ro
            Hzo(i,j,k)=z_wo(i,j,k+1)-z_wo(i,j,k)
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
         DO k=1,No
            cff_ro=hco*sc_ro(k)
            cff_wo=hco*sc_wo(k+1)
            cff1_ro=Cs_ro(k)
            cff1_wo=Cs_wo(k+1)
            hwater=h(i,j)
            hinv=1.0/(hco+hwater)
            cff2_ro=(cff_ro+cff1_ro*hwater)*hinv
            cff2_wo=(cff_wo+cff1_wo*hwater)*hinv
            z_wo(i,j,k+1)=hwater*cff2_wo
            z_ro(i,j,k)=hwater*cff2_ro
            Hzo(i,j,k)=z_wo(i,j,k+1)-z_wo(i,j,k)
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
         DO k=1,No
            cff_ro=hco*sc_ro(k)
            cff_wo=hco*sc_wo(k+1)
            cff1_ro=Cs_ro(k)
            cff1_wo=Cs_wo(k+1)
            hwater=h(i,j)
            hinv=1.0/(hco+hwater)
            cff2_ro=(cff_ro+cff1_ro*hwater)*hinv
            cff2_wo=(cff_wo+cff1_wo*hwater)*hinv
            z_wo(i,j,k+1)=hwater*cff2_wo
            z_ro(i,j,k)=hwater*cff2_ro
            Hzo(i,j,k)=z_wo(i,j,k+1)-z_wo(i,j,k)
         END DO
      END DO
   END DO      
ELSE
   write(*,*) 'Not programmed these vertical grid options.'
   STOP
ENDIF

z_ro_u= 0.5*(z_ro(1:Lp-1,:,:)+z_ro(2:Lp,:,:))
z_ro_v= 0.5*(z_ro(:,1:Mp-1,:)+z_ro(:,2:Mp,:))
Hzo_u=  0.5*(Hzo(1:Lp-1,:,:)+Hzo(2:Lp,:,:))
Hzo_v=  0.5*(Hzo(:,1:Mp-1,:)+Hzo(:,2:Mp,:))
h_u=    0.5*(h(1:Lp-1,:)+h(2:Lp,:))
h_v=    0.5*(h(:,1:Mp-1)+h(:,2:Mp))

! ..........................
write(*,*) 'Finished computing z-values in ROMS grid'

OPEN(unit=10,file= 'clim_blank.cdl')
DO I=1,58
READ(10,'(a)') line_clim_file(I)
ENDDO
CLOSE(10)

! s_rho
write(num_temp,*) No
num_temp=ADJUSTL(num_temp)
len_num= len_trim(num_temp)
line_temp= line_clim_file(4)
line_temp(10:11+len_num) = trim(num_temp)//' ;'
line_clim_file(4)=trim(line_temp)
! eta_rho
write(num_temp,*) MP
num_temp=ADJUSTL(num_temp)
len_num= len_trim(num_temp)
line_temp= line_clim_file(5)
line_temp(12:13+len_num) = trim(num_temp)//' ;'
line_clim_file(5)=trim(line_temp)
! xi_rho
write(num_temp,*) LP
num_temp=ADJUSTL(num_temp)
len_num= len_trim(num_temp)
line_temp= line_clim_file(6)
line_temp(11:12+len_num) = trim(num_temp)//' ;'
line_clim_file(6)=trim(line_temp)
!eta_u
write(num_temp,*) MP
num_temp=ADJUSTL(num_temp)
len_num= len_trim(num_temp)
line_temp= line_clim_file(7)
line_temp(10:11+len_num) = trim(num_temp)//' ;'
line_clim_file(7)=trim(line_temp)
!xi_u
write(num_temp,*) LP-1
num_temp=ADJUSTL(num_temp)
len_num= len_trim(num_temp)
line_temp= line_clim_file(8)
line_temp(9:10+len_num) = trim(num_temp)//' ;'
line_clim_file(8)=trim(line_temp)
!eta_v
write(num_temp,*) MP-1
num_temp=ADJUSTL(num_temp)
len_num= len_trim(num_temp)
line_temp= line_clim_file(9)
line_temp(10:11+len_num) = trim(num_temp)//' ;'
line_clim_file(9)=trim(line_temp)
!xi_v
write(num_temp,*) LP
num_temp=ADJUSTL(num_temp)
len_num= len_trim(num_temp)
line_temp= line_clim_file(10)
line_temp(9:10+len_num) = trim(num_temp)//' ;'
line_clim_file(10)=trim(line_temp)
! Time attribute
unit_time_att= '"'//trim(unit_time_att)//'" ;'
len_num= len_trim(unit_time_att)
line_temp= line_clim_file(24)
line_temp(16:16+len_num) = TRIM(unit_time_att)
line_clim_file(24)=trim(line_temp)
!Title
num_temp='"ROMS Climatology file generated for the domain '//TRIM(gridfile)//'" ;'
len_num= len_trim(num_temp)
line_temp= line_clim_file(55)
line_temp(11:11+len_num) = num_temp
line_clim_file(55)=trim(line_temp)

OPEN(unit=11,file='clim_final.cdl')
DO I=1,54
WRITE(11,'(a60)') line_clim_file(I)
ENDDO
WRITE(11,'(a150)') line_clim_file(55)
WRITE(11,'(a60)') line_clim_file(56)
WRITE(11,'(a100)') line_clim_file(57)
WRITE(11,'(a60)') line_clim_file(58)
CLOSE(11)

call system('ncgen -o '//TRIM(outputfile)// ' ' // 'clim_final.cdl')
! ------------------------------------------------------
! Open output netCDF files
statuso = nf90_open(trim(outputfile),nf90_write,ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_inq_varid(ncido, "salt",SaltVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output salt variable'
statuso = nf90_inq_varid(ncido, "temp",TempVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output temp variable'
statuso = nf90_inq_varid(ncido, "u",UVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output u variable'
statuso = nf90_inq_varid(ncido, "v",VVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output v variable'

statuso = nf90_inq_varid(ncido, "ubar",UbarVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output ubar variable'
statuso = nf90_inq_varid(ncido, "vbar",VbarVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output vbar variable'
statuso = nf90_inq_varid(ncido, "zeta",ZetaVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output zeta variable'

statuso = nf90_inq_varid(ncido, "time",TimeVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output time variable'

write(*,*) 'Finished opening output files'

! Open PHC file and find variables
! Open input averages file
status_PHC = nf90_open(trim(inputfile),nf90_nowrite,ncid_PHC)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
status_PHC = nf90_inq_varid(ncid_PHC, "depth",LevelPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable depth'
status_PHC = nf90_inq_varid(ncid_PHC, "time",TimePHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable time'
status_PHC = nf90_inq_varid(ncid_PHC, "salinity",SaltPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable salinity'
status_PHC = nf90_inq_varid(ncid_PHC, "temperature",TempPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable temperature'
status_PHC = nf90_inq_varid(ncid_PHC, "u",UPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable u'
status_PHC = nf90_inq_varid(ncid_PHC, "v",VPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable v'
status_PHC = nf90_inq_varid(ncid_PHC, "ssh",SSHPHCVarId)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
write(*,*) 'Found input variable ssh'

!Read in depths
   status_PHC = nf90_get_var(ncid_PHC,LevelPHCVarId,LevelPHC)
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
   write(*,*) 'Depth levels read in'

DO T=1,NT

!read in time
   status_PHC = nf90_get_var(ncid_PHC,TimePHCVarId,TimePHC,(/t/))
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)
   write(*,*) 'Time = ', TimePHC

!-------------------------------------------------------
! Do salinity

   status_PHC = nf90_get_var(ncid_PHC,SaltPHCVarId,SaltPHC,(/1,1,1,t/))
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

 
 WHERE(isnan(SaltPHC)) SaltPHC=35.0

! Vertical interpolation
   DO j=1,Mp
   DO i=1,Lp
      DO k=1,No
         IF(-z_ro(i,j,k).LT.LevelPHC(1)) THEN
            salt_out(i,j,k) = SaltPHC(i,j,1)
         ELSEIF(-z_ro(i,j,k).GT.LevelPHC(nlev)) THEN
            salt_out(i,j,k) = SaltPHC(i,j,nlev)
         ELSE
            DO kT=1,nlev
               IF(-z_ro(i,j,k).LT.LevelPHC(kT+1).AND.              &
                  -z_ro(i,j,k).GE.LevelPHC(kT)) THEN
                  rz2 = (-z_ro(i,j,k)-LevelPHC(kT))/               &
                          (LevelPHC(kT+1)-LevelPHC(kT))
                  rz1 = 1.0-rz2 
                  salt_out(i,j,k) = rz1*SaltPHC(i,j,kT) + rz2*SaltPHC(i,j,kT+1)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   ENDDO

! Output to netCDF file
   statuso = nf90_put_var(ncido,SaltVarIdo,salt_out,(/1,1,1,t/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'Vertical interpolation of Salinity done at time = ',TimePHC

! ---------------------------------------------------------------
! Do temperature

   status_PHC = nf90_get_var(ncid_PHC,TempPHCVarId,TempPHC,(/1,1,1,t/))
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

 
 WHERE(isnan(TempPHC)) TempPHC=2.0

! Vertical interpolation
   DO j=1,Mp
   DO i=1,Lp
      DO k=1,No
         IF(-z_ro(i,j,k).LT.LevelPHC(1)) THEN
            temp_out(i,j,k) = TempPHC(i,j,1)
         ELSEIF(-z_ro(i,j,k).GT.LevelPHC(nlev)) THEN
            temp_out(i,j,k) = TempPHC(i,j,nlev)
         ELSE
            DO kT=1,nlev
               IF(-z_ro(i,j,k).LT.LevelPHC(kT+1).AND.              &
                  -z_ro(i,j,k).GE.LevelPHC(kT)) THEN
                  rz2 = (-z_ro(i,j,k)-LevelPHC(kT))/               &
                          (LevelPHC(kT+1)-LevelPHC(kT))
                  rz1 = 1.0-rz2 
                  temp_out(i,j,k) = rz1*TempPHC(i,j,kT) + rz2*TempPHC(i,j,kT+1)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   ENDDO

! Output to netCDF file
   statuso = nf90_put_var(ncido,TempVarIdo,temp_out,(/1,1,1,t/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'Vertical interpolation of Temperature done at time= ',TimePHC

! ---------------------------------------------------------------
! Do U velocity component

   status_PHC = nf90_get_var(ncid_PHC,UPHCVarId,UPHC,(/1,1,1,t/))
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

 WHERE(isnan(UPHC)) UPHC=0.0

! Vertical interpolation
   DO j=1,Mp
   DO i=1,Lp-1
      DO k=1,No
         IF(-z_ro_u(i,j,k).LT.LevelPHC(1)) THEN
            u_out(i,j,k) = UPHC(i,j,1)
         ELSEIF(-z_ro_u(i,j,k).GT.LevelPHC(nlev)) THEN
            u_out(i,j,k) = UPHC(i,j,nlev)
         ELSE
            DO kT=1,nlev
               IF(-z_ro_u(i,j,k).LT.LevelPHC(kT+1).AND.              &
                  -z_ro_u(i,j,k).GE.LevelPHC(kT)) THEN
                  rz2 = (-z_ro_u(i,j,k)-LevelPHC(kT))/               &
                          (LevelPHC(kT+1)-LevelPHC(kT))
                  rz1 = 1.0-rz2 
                  u_out(i,j,k) = rz1*UPHC(i,j,kT) + rz2*UPHC(i,j,kT+1)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   ENDDO

! Output to netCDF file
   statuso = nf90_put_var(ncido,UVarIdo,u_out,(/1,1,1,t/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'Vertical interpolation of U velocity done at time= ',TimePHC

! ---------------------------------------------------------------
! Do V velocity component

   status_PHC = nf90_get_var(ncid_PHC,VPHCVarId,VPHC,(/1,1,1,t/))
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

 WHERE(isnan(VPHC)) VPHC=0.0

! Vertical interpolation
   DO j=1,Mp-1
   DO i=1,Lp
      DO k=1,No
         IF(-z_ro_v(i,j,k).LT.LevelPHC(1)) THEN
            v_out(i,j,k) = VPHC(i,j,1)
         ELSEIF(-z_ro_u(i,j,k).GT.LevelPHC(nlev)) THEN
            v_out(i,j,k) = VPHC(i,j,nlev)
         ELSE
            DO kT=1,nlev
               IF(-z_ro_v(i,j,k).LT.LevelPHC(kT+1).AND.              &
                  -z_ro_v(i,j,k).GE.LevelPHC(kT)) THEN
                  rz2 = (-z_ro_v(i,j,k)-LevelPHC(kT))/               &
                          (LevelPHC(kT+1)-LevelPHC(kT))
                  rz1 = 1.0-rz2 
                  v_out(i,j,k) = rz1*VPHC(i,j,kT) + rz2*VPHC(i,j,kT+1)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   ENDDO

! Output to netCDF file
   statuso = nf90_put_var(ncido,VVarIdo,v_out,(/1,1,1,t/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'Vertical interpolation of V velocity done at time = ',TimePHC


! ---------------------------------------------------------------
! Do barotropic velocities

      ubar_out = 0.
      vbar_out = 0.

      DO k=1,No
         ubar_out(:,:) = ubar_out(:,:) + u_out(:,:,k)*Hzo_u(:,:,k)
         vbar_out(:,:) = vbar_out(:,:) + v_out(:,:,k)*Hzo_v(:,:,k)
      ENDDO
      ubar_out = ubar_out/h_u
      vbar_out = vbar_out/h_v

! Output to netCDF file
   statuso = nf90_put_var(ncido,UbarVarIdo,ubar_out,(/1,1,t/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)
   statuso = nf90_put_var(ncido,VbarVarIdo,vbar_out,(/1,1,t/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

   write(*,*) 'ubar,vbar done at time= ',TimePHC

! ----------------------------------------------------------------
! Do sea level height

   status_PHC = nf90_get_var(ncid_PHC,SSHPHCVarId,SSHPHC,(/1,1,t/))
   if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

   WHERE(isnan(SSHPHC)) SSHPHC=0.0
   zeta_out = SSHPHC
   
! Output to netCDF file
   statuso = nf90_put_var(ncido,ZetaVarIdo,zeta_out,(/1,1,t/))
   if(statuso /= nf90_NoErr) call handle_err(statuso)

   write(*,*) 'zeta done at time = ', TimePHC

! ------------------------------------------------------------------
! write time
statuso = nf90_put_var(ncido,TimeVarIdo,TimePHC,(/t/))
if(statuso /= nf90_NoErr) call handle_err(statuso)
! ------------------------------------------------------------------
END DO

status_PHC = nf90_close(ncid_PHC)
if(status_PHC /= nf90_NoErr) call handle_err(status_PHC)

statuso = nf90_close(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

END PROGRAM hycom2roms_vert
