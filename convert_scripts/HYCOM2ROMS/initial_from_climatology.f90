PROGRAM initial_from_climatology
! Program to obtain the initial condition from the climatology file

use netcdf
USE, INTRINSIC :: IEEE_ARITHMETIC
implicit none

interface
     FUNCTION jd(yy,mm,dd)
     implicit none
     INTEGER jd,yy,mm,dd
     end function jd
end interface

INTEGER i, j, k,t, interp_index
INTEGER year_ini,month_ini,day_ini,hour_ini,min_ini,sec_ini
INTEGER year_ref,month_ref,day_ref,hour_ref,min_ref,sec_ref
CHARACTER*160 inputfile, outputfile
CHARACTER*160 PHCfile,num_temp,line_temp,unit_time_att,gridfile

INTEGER XiDimID, EtaDimID, Lp, Mp,No,len_num
INTEGER XiRhoDimId, EtaRhoDimId, SDimId, TimeDimId
INTEGER LonRhoVarId, LatRhoVarId, HVarId, AngleVarId
INTEGER status_roms, ncid_roms

CHARACTER*160, DIMENSION(64) :: line_clim_file

REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SaltPHC
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: TempPHC
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: UPHC
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: VPHC
REAL, ALLOCATABLE, DIMENSION(:,:)   :: SSHPHC
REAL, ALLOCATABLE, DIMENSION(:,:)   :: UBARPHC
REAL, ALLOCATABLE, DIMENSION(:,:)   :: VBARPHC
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Vect_TimePHC
DOUBLE PRECISION, DIMENSION(2) :: time_temp
DOUBLE PRECISION TimePHC, tinc,tincm, time_ini
 
INTEGER ntime
INTEGER nc_clim, status
INTEGER SSHPHCVarId,TimePHCVarId, LevelPHCVarId
INTEGER UPHCVarId,VPHCVarId,ZetaPHCVarId,SaltPHCVarId,TempPHCVarId,TimeVarId,UBARPHCVarId,VBARPHCVarId
INTEGER statuso, ncido
INTEGER UVarIdo, VVarIdo, TimeVarIdo,ocean_TimeVarIdo, SaltVarIdo, TempVarIdo
INTEGER UbarVarIdo, VbarVarIdo, ZetaVarIdo
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: salt_out
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp_out
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: u_out
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: v_out
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ubar_out
REAL, ALLOCATABLE, DIMENSION(:,:)   :: vbar_out
REAL, ALLOCATABLE, DIMENSION(:,:)   :: zeta_out

REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: salt_temp
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: temp_temp
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: u_temp
REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: v_temp
REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: ubar_temp
REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: vbar_temp
REAL, ALLOCATABLE, DIMENSION(:,:,:)   :: zeta_temp

OPEN(unit=10,file='initial_from_climatology.in')
read(10,'(a)') gridfile
read(10,'(a)') inputfile
read(10,'(a)') outputfile
read(10,*) year_ini,month_ini,day_ini,hour_ini,min_ini,sec_ini
CLOSE(10)

! Get info on Climatology grid
status = nf90_open(TRIM(inputfile),nf90_nowrite,nc_clim)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_dimid(nc_clim, "xi_rho", XiDimID)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_dimid(nc_clim, "eta_rho", EtaDimID)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_Inquire_Dimension(nc_clim,XiDimID,len = Lp)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_Inquire_Dimension(nc_clim,EtaDimID,len = Mp)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_dimid(nc_clim, "s_rho", SDimId)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_dimid(nc_clim, "time", TimeDimId)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_inq_varid(nc_clim, "time",TimeVarId)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_Inquire_Dimension(nc_clim,SDimID,len = No)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_Inquire_Dimension(nc_clim,TimeDimID,len = ntime)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_get_att(nc_clim,TimeVarId,'units',unit_time_att)
if(status /= nf90_NoErr) call handle_err(status)
ALLOCATE(Vect_TimePHC(ntime))
status = nf90_get_var(nc_clim,TimeVarId,Vect_TimePHC)
if(status /= nf90_NoErr) call handle_err(status)
status = nf90_close(nc_clim)
if(status /= nf90_NoErr) call handle_err(status)    


write(*,*) 'Define HYCOM files to horizontally interpolate into the ROMS gridfile'

write(*,*) 'Lp = ',Lp,',  Mp = ',Mp,',  N = ',No


ALLOCATE(salt_out(Lp,Mp,No))
ALLOCATE(temp_out(Lp,Mp,No))
ALLOCATE(u_out(Lp-1,Mp,No))
ALLOCATE(v_out(Lp,Mp-1,No))
ALLOCATE(ubar_out(Lp-1,Mp))
ALLOCATE(vbar_out(Lp,Mp-1))
ALLOCATE(zeta_out(Lp,Mp))

ALLOCATE(salt_temp(Lp,Mp,No,2))
ALLOCATE(temp_temp(Lp,Mp,No,2))
ALLOCATE(u_temp(Lp-1,Mp,No,2))
ALLOCATE(v_temp(Lp,Mp-1,No,2))
ALLOCATE(ubar_temp(Lp-1,Mp,2))
ALLOCATE(vbar_temp(Lp,Mp-1,2))
ALLOCATE(zeta_temp(Lp,Mp,2))

ALLOCATE(SaltPHC(Lp,Mp,No))
ALLOCATE(TempPHC(Lp,Mp,No))
ALLOCATE(UPHC(Lp-1,Mp,No))
ALLOCATE(VPHC(Lp,Mp-1,No))
ALLOCATE(SSHPHC(Lp,Mp))
ALLOCATE(UBARPHC(Lp-1,Mp))
ALLOCATE(VBARPHC(Lp,Mp-1))

read(unit_time_att(15:18),*)year_ref
read(unit_time_att(20:21),*)month_ref
read(unit_time_att(23:24),*)day_ref
read(unit_time_att(26:27),*)min_ref
read(unit_time_att(29:30),*)hour_ref
read(unit_time_att(32:33),*)sec_ref

!-------------------------------------------------------------------------------
! Get time ini
time_ini = DBLE(jd(year_ini,month_ini,day_ini)-jd(year_ref,month_ref,day_ref))*24.d0*3600.d0 &
                   +  DBLE(hour_ini-hour_ref)*3600.d0                                        &
                   +  DBLE(min_ini-min_ref)*60.d0                                            &
                   +  DBLE(sec_ini-sec_ref)                                                  
! ------------------------------------------------------------------------------
write(*,*) 'Initial time is: ', time_ini, trim(unit_time_att)

DO i=1,ntime-1
IF((time_ini.ge.Vect_TimePHC(i)).and.(time_ini.lt.Vect_TimePHC(i+1)) ) then
interp_index= i
ENDIF
END DO

OPEN(unit=10,file= 'initial_blank.cdl')
DO I=1,64
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
line_temp= line_clim_file(14)
line_temp(16:16+len_num) = TRIM(unit_time_att)
line_clim_file(14)=trim(line_temp)
! Ocean Time attribute
line_temp= line_clim_file(19)
line_temp(21:21+len_num) = TRIM(unit_time_att)
line_clim_file(19)=trim(line_temp)
!Title
num_temp='"ROMS ini file for the domain '//TRIM(gridfile)//'" ;'
len_num= len_trim(num_temp)
line_temp= line_clim_file(61)
line_temp(11:11+len_num) = num_temp
line_clim_file(61)=trim(line_temp)

OPEN(unit=11,file= 'initial_final.cdl')
DO I=1,59
WRITE(11,'(a60)') line_clim_file(I)
ENDDO
WRITE(11,'(a70)') line_clim_file(60)
WRITE(11,'(a125)') line_clim_file(61)
WRITE(11,'(a60)') line_clim_file(62)
WRITE(11,'(a100)') line_clim_file(63)
WRITE(11,'(a60)') line_clim_file(64)
CLOSE(11)

call system('ncgen -o '//TRIM(outputfile)// ' ' // 'initial_final.cdl')
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

statuso = nf90_inq_varid(ncido, "ocean_time",ocean_TimeVarIdo)
if(statuso /= nf90_NoErr) call handle_err(statuso)
write(*,*) 'Found output ocean_time variable'

write(*,*) 'Finished opening output files'

! Open PHC file and find variables
! Open input averages file
status = nf90_open(TRIM(inputfile),nf90_nowrite,nc_clim)
if(status /= nf90_NoErr) call handle_err(status)

status = nf90_inq_varid(nc_clim, "time",TimePHCVarId)
if(status /= nf90_NoErr) call handle_err(status)
write(*,*) 'Found input variable time'
status = nf90_inq_varid(nc_clim, "salt",SaltPHCVarId)
if(status /= nf90_NoErr) call handle_err(status)
write(*,*) 'Found input variable salinity'
status = nf90_inq_varid(nc_clim, "temp",TempPHCVarId)
if(status /= nf90_NoErr) call handle_err(status)
write(*,*) 'Found input variable temperature'
status = nf90_inq_varid(nc_clim, "u",UPHCVarId)
if(status /= nf90_NoErr) call handle_err(status)
write(*,*) 'Found input variable u'
status = nf90_inq_varid(nc_clim, "v",VPHCVarId)
if(status /= nf90_NoErr) call handle_err(status)
write(*,*) 'Found input variable v'
status = nf90_inq_varid(nc_clim, "zeta",SSHPHCVarId)
if(status /= nf90_NoErr) call handle_err(status)
write(*,*) 'Found input variable zeta'
status = nf90_inq_varid(nc_clim, "ubar",UBARPHCVarId)
if(status /= nf90_NoErr) call handle_err(status)
write(*,*) 'Found input variable ubar'
status = nf90_inq_varid(nc_clim, "vbar",VBARPHCVarId)
if(status /= nf90_NoErr) call handle_err(status)
write(*,*) 'Found input variable vbar'


i= interp_index

! Time
status = nf90_get_var(nc_clim,TimePHCVarId,TimePHC,(/i/))
if(status /= nf90_NoErr) call handle_err(status)
time_temp(1)=TimePHC
status = nf90_get_var(nc_clim,TimePHCVarId,TimePHC,(/i+1/))
if(status /= nf90_NoErr) call handle_err(status)
time_temp(2)=TimePHC

tinc= (time_ini-time_temp(1))/(time_temp(2)-time_temp(1))
tincm= 1-tinc

! Output to netCDF file

statuso = nf90_put_var(ncido,TimeVarIdo,time_ini,(/1/))
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_put_var(ncido,ocean_TimeVarIdo,time_ini,(/1/))
if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'time done'


! Salinity
status = nf90_get_var(nc_clim,SaltPHCVarId,SaltPHC,(/1,1,1,i/))
if(status /= nf90_NoErr) call handle_err(status)
salt_temp(:,:,:,1)=SaltPHC
status = nf90_get_var(nc_clim,SaltPHCVarId,SaltPHC,(/1,1,1,i+1/))
if(status /= nf90_NoErr) call handle_err(status)
salt_temp(:,:,:,2)=SaltPHC

salt_out= tincm*salt_temp(:,:,:,1)+tinc*salt_temp(:,:,:,2)

! Output to netCDF file

statuso = nf90_put_var(ncido,SaltVarIdo,salt_out,(/1,1,1,1/))
if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'salt done'

! Temperature
status = nf90_get_var(nc_clim,TempPHCVarId,TempPHC,(/1,1,1,i/))
if(status /= nf90_NoErr) call handle_err(status)
temp_temp(:,:,:,1)=TempPHC
status = nf90_get_var(nc_clim,TempPHCVarId,TempPHC,(/1,1,1,i+1/))
if(status /= nf90_NoErr) call handle_err(status)
temp_temp(:,:,:,2)=TempPHC

temp_out= tincm*temp_temp(:,:,:,1)+tinc*temp_temp(:,:,:,2)

! Output to netCDF file

statuso = nf90_put_var(ncido,TempVarIdo,temp_out,(/1,1,1,1/))
if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'temp done'

! U velocity component
status = nf90_get_var(nc_clim,UPHCVarId,UPHC,(/1,1,1,i/))
if(status /= nf90_NoErr) call handle_err(status)
u_temp(:,:,:,1)=UPHC
status = nf90_get_var(nc_clim,UPHCVarId,UPHC,(/1,1,1,i+1/))
if(status /= nf90_NoErr) call handle_err(status)
u_temp(:,:,:,2)=UPHC

u_out= tincm*u_temp(:,:,:,1)+tinc*u_temp(:,:,:,2)
 
! V velocity component
status = nf90_get_var(nc_clim,VPHCVarId,VPHC,(/1,1,1,i/))
if(status /= nf90_NoErr) call handle_err(status)
v_temp(:,:,:,1)=VPHC
status = nf90_get_var(nc_clim,VPHCVarId,VPHC,(/1,1,1,i+1/))
if(status /= nf90_NoErr) call handle_err(status)
v_temp(:,:,:,2)=VPHC

v_out= tincm*v_temp(:,:,:,1)+tinc*v_temp(:,:,:,2)

! Output to netCDF file

statuso = nf90_put_var(ncido,UVarIdo,u_out,(/1,1,1,1/))
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_var(ncido,VVarIdo,v_out,(/1,1,1,1/))
if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'u,v done'


! Ubar velocity component
status = nf90_get_var(nc_clim,UBARPHCVarId,UBARPHC,(/1,1,i/))
if(status /= nf90_NoErr) call handle_err(status)
ubar_temp(:,:,1)=UBARPHC
status = nf90_get_var(nc_clim,UBARPHCVarId,UBARPHC,(/1,1,i+1/))
if(status /= nf90_NoErr) call handle_err(status)
ubar_temp(:,:,2)=UBARPHC
 
ubar_out= tincm*ubar_temp(:,:,1)+tinc*ubar_temp(:,:,2)

! Vbar velocity component
status = nf90_get_var(nc_clim,VBARPHCVarId,VBARPHC,(/1,1,i/))
if(status /= nf90_NoErr) call handle_err(status)
vbar_temp(:,:,1)=VBARPHC
status = nf90_get_var(nc_clim,VBARPHCVarId,VBARPHC,(/1,1,i+1/))
if(status /= nf90_NoErr) call handle_err(status)
vbar_temp(:,:,2)=VBARPHC

vbar_out= tincm*vbar_temp(:,:,1)+tinc*vbar_temp(:,:,2)

! Output to netCDF file

statuso = nf90_put_var(ncido,UbarVarIdo,ubar_out,(/1,1,1/))
if(statuso /= nf90_NoErr) call handle_err(statuso)
statuso = nf90_put_var(ncido,VbarVarIdo,vbar_out,(/1,1,1/))
if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'ubar,vbar done'


! Sea level height
status = nf90_get_var(nc_clim,SSHPHCVarId,SSHPHC,(/1,1,i/))
if(status /= nf90_NoErr) call handle_err(status)
zeta_temp(:,:,1)=SSHPHC
status = nf90_get_var(nc_clim,SSHPHCVarId,SSHPHC,(/1,1,i+1/))
if(status /= nf90_NoErr) call handle_err(status)
zeta_temp(:,:,2)=SSHPHC

zeta_out= tincm*zeta_temp(:,:,1)+tinc*zeta_temp(:,:,2)

! Output to netCDF file

statuso = nf90_put_var(ncido,ZetaVarIdo,zeta_out,(/1,1,1/))
if(statuso /= nf90_NoErr) call handle_err(statuso)

write(*,*) 'zeta done'

status = nf90_close(nc_clim)
if(status /= nf90_NoErr) call handle_err(status)    
statuso = nf90_close(ncido)
if(statuso /= nf90_NoErr) call handle_err(statuso)

END PROGRAM initial_from_climatology
