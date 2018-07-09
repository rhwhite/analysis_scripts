PROGRAM smooth_FC_fields
! ______________________________________________
use netcdf
USE, INTRINSIC :: IEEE_ARITHMETIC
implicit none

interface
     FUNCTION jd(yy,mm,dd)
     implicit none
     INTEGER jd,yy,mm,dd
     end function jd
end interface

INTEGER iyear, itime
INTEGER year_st, year_end, iday1, iday2, ndays, n12hr
INTEGER ncid_era, status_era
INTEGER LonDimID, LatDimID, MaskVarId, CPVarId, LSPVarId
INTEGER SSRDVarId, STRDVarId
INTEGER nlon, nlat
INTEGER nvalue
REAL, ALLOCATABLE, DIMENSION(:,:) :: scr, scrp, emask
REAL, ALLOCATABLE, DIMENSION(:,:) :: error, work
REAL :: tx
CHARACTER(len=160) pathin, pathout,filein,fileout,maskfile
CHARACTER(len= 50):: ryear
REAL, PARAMETER :: undef = 2.E+35            ! Undefined land value
REAL, PARAMETER :: critx=0.01
REAL, PARAMETER :: cor=1.6
INTEGER, PARAMETER :: mxs=100

tx=0.9*undef

OPEN(unit=10,file='smooth_FC_fields.in')
READ(10,'(a)') pathin
READ(10,'(a)') pathout
READ(10,*) year_st, year_end
CLOSE(10)

maskfile = TRIM(pathout)//'landmask.nc'

! Get dimensions of ERAI files
status_era = nf90_open(TRIM(maskfile),nf90_nowrite,ncid_era)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_inq_dimid(ncid_era, "longitude", LonDimID)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_inq_dimid(ncid_era, "latitude", LatDimID)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_Inquire_Dimension(ncid_era,LonDimID,len = nlon)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_Inquire_Dimension(ncid_era,LatDimID,len = nlat)
if(status_era /= nf90_NoErr) call handle_err(status_era)

ALLOCATE(scr(nlon,nlat))
ALLOCATE(scrp(nlon,nlat))
ALLOCATE(emask(nlon,nlat))
ALLOCATE(work(nlon,nlat))
ALLOCATE(error(nlon,nlat))

! Get land-sea mask
status_era = nf90_inq_varid(ncid_era, "lsm",MaskVarId)
if(status_era /= nf90_NoErr) call handle_err(status_era)
status_era = nf90_get_var(ncid_era,MaskVarId,scr)
if(status_era /= nf90_NoErr) call handle_err(status_era)

! Close netcdf file
status_era = nf90_close(ncid_era)
if(status_era /= nf90_NoErr) call handle_err(status_era)

DO iyear= year_st,year_end

write(ryear,*) iyear
ryear = adjustl(ryear)

write(*,*) 'Treament of year ',TRIM(ryear), ' in progress ...'

filein=TRIM(pathin)//TRIM(ryear)// '0101FC.nc'
fileout=TRIM(pathout)//TRIM(ryear)// '0101FC.nc'

CALL system('mv ' // TRIM(filein) // ' ' // TRIM(fileout))

status_era = nf90_open(TRIM(fileout),nf90_write,ncid_era)
if(status_era /= nf90_NoErr) call handle_err(status_era)

! Get number of 6-hourly records to process
iday1 = jd(iyear,1,1)
iday2 = jd(iyear+1,1,1)
ndays = iday2-iday1
n12hr = ndays*2

DO itime= 1,n12hr
! _________________________________________________________________
! Read in CP value
! Get CP var id
    status_era = nf90_inq_varid(ncid_era, "cp",CPVarId)
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    status_era = nf90_get_var(ncid_era,CPVarId,scr,start=(/1, 1,itime /))
    if(status_era /= nf90_NoErr) call handle_err(status_era)
! Fill in masked-out values 
    scrp = 0.
    scrp = scr
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
    scr = scrp
! Smooth zonally
    CALL smooth_polar(scrp, scr, nlon, nlat)
! Write out in netCDF file
    status_era = nf90_put_var(ncid_era,CPVarId,scr,(/1,1,itime/))
    if(status_era /= nf90_NoErr) call handle_err(status_era)   
! _________________________________________________________________
! Read in LSP value
! Get LSP var id
    status_era = nf90_inq_varid(ncid_era, "lsp",LSPVarId)
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    status_era = nf90_get_var(ncid_era,LSPVarId,scr,start=(/1, 1,itime /))
    if(status_era /= nf90_NoErr) call handle_err(status_era)
! Fill in masked-out values 
    scrp = 0.
    scrp = scr
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
    scr = scrp
! Smooth zonally
    CALL smooth_polar(scrp, scr, nlon, nlat)
! Write out in netCDF file
    status_era = nf90_put_var(ncid_era,LSPVarId,scr,(/1,1,itime/))
    if(status_era /= nf90_NoErr) call handle_err(status_era)   
! _________________________________________________________________
! Read in SSRD values
! Get SSRD var id
    status_era = nf90_inq_varid(ncid_era, "ssrd",SSRDVarId)
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    status_era = nf90_get_var(ncid_era,SSRDVarId,scr,start=(/1, 1,itime /))
    if(status_era /= nf90_NoErr) call handle_err(status_era)
! Fill in masked-out values 
    scrp = 0.
    scrp = scr
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
    scr = scrp
! Smooth zonally
    CALL smooth_polar(scrp, scr, nlon, nlat)
! Write out in netCDF file
    status_era = nf90_put_var(ncid_era,SSRDVarId,scr,(/1,1,itime/))
    if(status_era /= nf90_NoErr) call handle_err(status_era)   
! _________________________________________________________________
! Read in STRD values
! Get STRD var id
    status_era = nf90_inq_varid(ncid_era, "strd",STRDVarId)
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    status_era = nf90_get_var(ncid_era,STRDVarId,scr,start=(/1, 1,itime /))
    if(status_era /= nf90_NoErr) call handle_err(status_era)
! Fill in masked-out values 
    scrp = 0.
    scrp = scr
    WHERE(emask>0.01) scrp = undef
    CALL fill(nlon,nlat,1,nlon,1,nlat,scrp,tx,critx,cor,mxs,work,error,nvalue)
    scr = scrp
! Smooth zonally
    CALL smooth_polar(scrp, scr, nlon, nlat)
! Write out in netCDF file
    status_era = nf90_put_var(ncid_era,STRDVarId,scr,(/1,1,itime/))
    if(status_era /= nf90_NoErr) call handle_err(status_era)   
ENDDO

! Close netcdf file
    status_era = nf90_close(ncid_era)
    if(status_era /= nf90_NoErr) call handle_err(status_era)
    write(*,*) 'Laplacian interpolation and Polar smoothing done for file: ', TRIM(ryear)// '0101FC.nc'
!
ENDDO
! ______________________________________________
END PROGRAM smooth_FC_fields



