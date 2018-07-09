SUBROUTINE spec_hum(td, ta, p, evapres, q, nlonp, nlats,IFL)
  IMPLICIT NONE
!
!     Laget av Morten Odegaard Februar 1999  
!     Beregner spesifikk fuktighet Q og vanndampstrykk 
!     Beregner forst metningsvanndamptrykket (E_AVTD) ut i fra Td, 
!     et problem er  her om man skal beregne i forhold til over vann eller is. 
!     TD fra ECMWF er beregnet for is hvis T<0, og ellers for vann (OBS! dette
!     gjelder frem til 1997) Spesifikk fuktighet beregnes saa fra 
!     metningsvanndamptrykket. Bruker at Q(T) = Q(TD) = QS(E_AVTD(TD)),
!     formlene er funnet i Rogers&Yau (1994).  
!     fajd Har ogs� puttet inn formlene som finnes i Gill 1982: Atmosphere-Ocean Dyn
!
!
!   INPUT PARAMETERS:
!
  INTEGER, INTENT(IN) :: nlats, nlonp
  INTEGER, INTENT(IN) :: IFL  ! Switch = 0 => TD= Specific humidity OR = 1 => TD= Dewpoint temp.
  REAL                :: td(nlonp,nlats) ! Specific humidity (-) or dewpoint temperature (deg C) (see IFL)
  REAL                :: ta(nlonp,nlats) ! Atmosphere temperature at 2 m (deg C)
  REAL                :: p(nlonp,nlats)  ! Atmospheric mean sea level pressure (Pa)
  

!
!   OUTPUT PARAMETER:
!
  REAL                :: evapres(nlonp,nlats) ! Water-vapour pressure in air (Pa)
  REAL                :: q(nlonp,nlats)     ! Atmosphere specific humidity (-)
!
! Local variables
!
  REAL                :: arg(nlonp,nlats)
  REAL                :: f(nlonp,nlats)
 
  IF (IFL.EQ.1) THEN    
     arg =(0.7859 + 0.03477 * td)/(1.0 + 0.00412*td)
     f = 1.0 + 1e-8*p*(4.5+0.0006*td*td)    
     ! Over ice or not ?
     WHERE(ta < 0.0) arg = arg + 0.00422*td
     ! pressure in Pa
     evapres = 1.e+2 * f * 10.**arg 
     q = 0.622*evapres/(p-0.378*evapres)
  ELSE
     q = td
     evapres = p * q/( 0.622 + 0.378*q)
  END IF
  RETURN
END SUBROUTINE spec_hum
