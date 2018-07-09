SUBROUTINE GR (GREG,JD)
!
!---COMPUTES THE GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY)
!   GIVEN THE JULIAN DATE (JD).
!
integer, parameter :: DefInt8 = selected_int_kind(8)

  INTEGER(kind = DefInt8) :: GREG(6),A,B,C,D,E,alpha,Z
  REAL F,Y,W,temp
  DOUBLE PRECISION JD
!

  JD = JD                                                     ! begin algorithm
  Z = INT(JD)
  F = JD - Z

  IF (Z .LT. 2299161) THEN
     A = Z
  ELSE
     ALPHA = INT((Z-1867216.25D0)/36524.25D0)
     A = Z + 1 + ALPHA - ALPHA/4
  END IF
  
  B = A + 1524
  C = INT((B-122.1D0)/365.25D0)

  D = INT(365.25D0*C)
  E = INT((B-D)/30.6001D0)
  
  temp = B - D - INT(30.6001D0*E) + F
  GREG(3)  = INT(temp)

  W = temp - GREG(3)
  HH = INT(W*24)
  IF(HH.ge.24)then
     GREG(4)= HH-24
  else
     GREG(4)= HH
  endif
  Y= W*24 - HH
  GREG(5) = INT(Y*60)
  W= Y*60 - GREG(5)
  GREG(6) = INT(W*60)
  
  IF (E .LT. 14) THEN
     GREG(2)  = E - 1
  ELSE
     GREG(2)  = E - 13
  END IF
  
  IF (GREG(2) .GT. 2) THEN
     GREG(1)  = C - 4716
  ELSE
     GREG(1)  = C - 4715
  END IF

  RETURN
END SUBROUTINE GR
