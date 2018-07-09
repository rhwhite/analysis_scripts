PROGRAM smooth_v
! ______________________________________________
use netcdf
USE, INTRINSIC :: IEEE_ARITHMETIC
implicit none

INTEGER I,istart,istop
INTEGER Rcond, nchar, nv, k, nvalue, ict
INTEGER nc, status, ncout, status_out
INTEGER XDimid, YDimid, NDimid, Lp, Mp, N, VarId
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: mask3d
REAL, ALLOCATABLE, DIMENSION(:,:) :: scrin, scrout, scr2d
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: scr3d
REAL, ALLOCATABLE, DIMENSION(:,:) :: rmask, error
REAL :: background, txp
CHARACTER(len=160) pathin, pathout,filein, fileout
CHARACTER(len=80) varname,num_temp
CHARACTER(len=4)  num_record
REAL, PARAMETER :: tx=1.0e+20
REAL, PARAMETER :: critx=1.0e-4
REAL, PARAMETER :: cor=1.5
INTEGER, PARAMETER :: mxs=200

! ______________________________________________
varname = 'v'
background = 0.0
txp = 1.01*tx
! ______________________________________________

OPEN(unit=10,file='smooth_v.in')
READ(10,'(a)') pathin
READ(10,'(a)') pathout
READ(10,*) istart,istop
CLOSE(10)

OPEN(unit=11,file='v_files_corrupted.out')

DO I=istart,istop

write(num_temp,*) I
num_temp=ADJUSTL(num_temp)
if (len_trim(num_temp).EQ.1) num_record= '000'//TRIM(num_temp)
if (len_trim(num_temp).EQ.2) num_record= '00'//TRIM(num_temp)
if (len_trim(num_temp).EQ.3) num_record= '0'//TRIM(num_temp)
if (len_trim(num_temp).EQ.4) num_record=  TRIM(num_temp)

filein=TRIM(pathin)//'vvel_opendap_'//num_record//'.nc'
fileout=TRIM(pathout)//'vvel_opendap_'//num_record//'.nc'

CALL system('mv ' // TRIM(filein) // ' ' // TRIM(fileout))



! Open output file
   status = nf90_open(TRIM(fileout),nf90_write,nc)
   if(status /= nf90_NoErr) call handle_err(status)
! Get dimensions
   status = nf90_inq_dimid(nc, "X", XDimID)
   if(status /= nf90_NoErr) call handle_err(status)
   status = nf90_Inquire_Dimension(nc,XDimID,len = Lp)
   if(status /= nf90_NoErr) call handle_err(status) 

   status = nf90_inq_dimid(nc, "Y", YDimID)
   if(status /= nf90_NoErr) call handle_err(status)
   status = nf90_Inquire_Dimension(nc,YDimID,len = Mp)
   if(status /= nf90_NoErr) call handle_err(status) 

   status = nf90_inq_dimid(nc, "depth", NDimID)
   if(status /= nf90_NoErr) call handle_err(status)
   status = nf90_Inquire_Dimension(nc,NDimID,len = N)
   if(status /= nf90_NoErr) call handle_err(status) 

   IF(I.eq.istart) THEN
      write(*,*) 'smooth_v: Lp, Mp, N = ',Lp,Mp,N

!    Allocate arrays
      ALLOCATE(scrin(Lp,Mp))
      ALLOCATE(scrout(Lp,Mp))
      ALLOCATE(scr3d(N,Lp,Mp))
      ALLOCATE(mask3d(N,Lp,Mp))
      ALLOCATE(scr2d(Lp,Mp))
      ALLOCATE(rmask(Lp,Mp))
      ALLOCATE(error(Lp,Mp))
   ENDIF
! Get variable
   status = nf90_inq_varid(nc, TRIM(varname), VarId)
   if(status /= nf90_NoErr) call handle_err(status)
   status = nf90_get_var(nc,VarId,scr3d)
   if(status /= nf90_NoErr) call handle_err(status)

! If the file contains only zero's, move the file back into the input folder and end the program
   if ((MAXVAL(scr3d).eq. 0.).and.(MINVAL(scr3d).eq.0.))then
   write(11,*)'File : vvel_opendap_'//num_record//'.nc does not contain any data '
   CALL system('rm -f '//TRIM(fileout))
   write(*,*) 'Must re download file: vvel_opendap_'//num_record//'.nc'
   GO TO 100
   endif

! Do masking
    mask3d(:,:,:) = 1
    WHERE(isnan(scr3d)) mask3d=0
    DO k=1,N
      nv = SUM(mask3d(k,:,:))
      IF(nv.EQ.0) THEN
         scr2d(:,:) = background
      ELSE
         scr2d = scr3d(k,:,:)
         WHERE(isnan(scr2d)) scr2d=txp
         CALL fill(Lp,Mp,1,Lp,1,Mp,scr2d,tx,critx,cor,mxs,rmask,error,nvalue)
      ENDIF
      scr3d(k,:,:) = scr2d(:,:)
    ENDDO
!
! Write array
    status = nf90_put_var(nc,VarId,scr3d,(/1,1,1/))
    if(status /= nf90_NoErr) call handle_err(status)        

    write(*,*) 'Laplacian interpolation done for file: vvel_opendap_'//num_record//'.nc'

! Close netcdf file
100 status = nf90_close(nc)
    if(status /= nf90_NoErr) call handle_err(status)

ENDDO
CLOSE(11)
! ______________________________________________
END PROGRAM smooth_v



