PROGRAM SODA_change_time
! Program to obtain climatology from SODA monthly T,S,U,V,zeta climatology

#undef  PST_OUT
#define GEO_OUT

use netcdf

implicit none


INTEGER year, month, day, jd, jday1, jday2, refyear_new, refmonth_new, refday_new
REAL tday, refdate1,refdate2,curdate,newdate

CHARACTER*160 outpath1,outpath2,outpath3, filelst
CHARACTER*160 clim_root, avgfile, avgfilei

INTEGER L1,L2,L3
INTEGER time1_dimid,time1_varid,tt
INTEGER status_roms, ncid1

INTEGER statuso


REAL*8, ALLOCATABLE, DIMENSION(:) :: times
REAL*8, ALLOCATABLE, DIMENSION(:) :: times_final


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

read(*,'(a)') outpath1
write(*,*) outpath1
read(*,'(a)') outpath2
write(*,*) outpath2
read(*,'(a)') outpath3
write(*,*) outpath3
read(*,*) refyear_new, refmonth_new, refday_new
write(*,*) refyear_new, refmonth_new, refday_new

!Julian Day for 1900-01-01
refdate1 = JD(1900,1,1)
!Julian Dat for reference time
refdate2 = JD(refyear_new,refmonth_new,refday_new)


! Change times for outpath1

status_roms = nf90_open(TRIM(outpath1),nf90_write,ncid1)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

statuso = nf90_inq_dimid(ncid1,'ocean_time',time1_dimid)
if(statuso /= nf90_NoErr) call handle_err(statuso)
status_roms = nf90_Inquire_Dimension(ncid1,time1_dimid,len = L1)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

write(*,*) 'L1 is', L1
ALLOCATE(times(L1))
ALLOCATE(times_final(L1))


statuso = nf90_inq_varid(ncid1,'ocean_time',time1_varid)
if(statuso /= nf90_NoErr) call handle_err(statuso)

statuso = nf90_get_var(ncid1,time1_varid,times)
if(statuso /= nf90_NoErr) call handle_err(statuso)

!Time currently in days since 1900-01-01 00:00:00
do tt=1,L1
   curdate = refdate1 + times(tt)
! change to seconds
   times_final(tt) = (curdate - refdate2) * 60. * 60. * 24
end do

statuso = nf90_put_var(ncid1,time1_varid,times_final)
if(statuso /= nf90_NoErr) call handle_err(statuso)

status_roms = nf90_close(ncid1)
if(status_roms /= nf90_NoErr) call handle_err(status_roms)

DEALLOCATE(times)
DEALLOCATE(times_final)

END PROGRAM SODA_change_time
