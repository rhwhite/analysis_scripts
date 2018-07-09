SUBROUTINE SSORT (X, Y, N, KFLAG)
!***BEGIN PROLOGUE  SSORT
!***PURPOSE  Sort an array and optionally make the same interchanges in
!            an auxiliary array.  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A2B
!***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
!***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   SSORT sorts array X and optionally makes the same interchanges in
!   array Y.  The array X may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!   Description of Parameters
!      X - array of values to be sorted   (usually abscissas)
!      Y - array to be (optionally) carried along
!      N - number of values in array X to be sorted
!      KFLAG - control parameter
!            =  2  means sort X in increasing order and carry Y along.
!            =  1  means sort X in increasing order (ignoring Y)
!            = -1  means sort X in decreasing order (ignoring Y)
!            = -2  means sort X in decreasing order and carry Y along.
!
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891024  Changed category.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain)
!   920501  Reformatted the REFERENCES section.  (DWL, WRB)
!   920519  Clarified error messages.  (DWL)
!   920801  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
!***END PROLOGUE  SSORT
IMPLICIT NONE
!     .. Intrinsic Functions ..
INTRINSIC ABS, INT
!     .. Scalar Arguments ..
INTEGER, INTENT(IN) :: KFLAG, N
!     .. Array Arguments ..
REAL :: X(N), Y(N)
!     .. Local Scalars ..
REAL :: R, T, TT, TTY, TY
INTEGER :: I, IJ, J, K, KK, L, M, NN
!     .. Local Arrays ..
INTEGER :: IL(21), IU(21)

!***FIRST EXECUTABLE STATEMENT  SSORT
NN = N
IF (NN .LT. 1) THEN
   PRINT *,'The number of values to be sorted is not positive.'
   RETURN
ENDIF
!
KK = ABS(KFLAG)
IF (KK.NE.1 .AND. KK.NE.2) THEN
   PRINT *,'The sort control parameter, K, is not 2, 1, -1, or -2.'
   RETURN
ENDIF
!
!     Alter array X to get decreasing order if needed
!
IF (KFLAG .LE. -1) THEN
   DO  I=1,NN
      X(I) = -X(I)
   END DO
END IF
!
IF (KK .EQ. 2) GO TO 100
!
!     Sort X only
!
M = 1
I = 1
J = NN
R = 0.375E0
!
20 IF (I .EQ. J) GO TO 60
IF (R .LE. 0.5898437E0) THEN
   R = R+3.90625E-2
ELSE
   R = R-0.21875E0
END IF
!
30 K = I
!
!     Select a central element of the array and save it in location T
!
IJ = I + INT((J-I)*R)
T = X(IJ)
!
!     If first element of array is greater than T, interchange with T
!
IF (X(I) .GT. T) THEN
   X(IJ) = X(I)
   X(I) = T
   T = X(IJ)
END IF
L = J
!
!     If last element of array is less than than T, interchange with T
!
IF (X(J) .LT. T) THEN
   X(IJ) = X(J)
   X(J) = T
   T = X(IJ)
!
!        If first element of array is greater than T, interchange with T
!
   IF (X(I) .GT. T) THEN
      X(IJ) = X(I)
      X(I) = T
      T = X(IJ)
   END IF
END IF
!
!     Find an element in the second half of the array which is smaller
!     than T
!
40 L = L-1
IF (X(L) .GT. T) GO TO 40
!
!     Find an element in the first half of the array which is greater
!     than T
!
50 K = K+1
IF (X(K) .LT. T) GO TO 50
!
!     Interchange these elements
!
IF (K .LE. L) THEN
   TT = X(L)
   X(L) = X(K)
   X(K) = TT
   GO TO 40
END IF
!
!     Save upper and lower subscripts of the array yet to be sorted
!
IF (L-I .GT. J-K) THEN
   IL(M) = I
   IU(M) = L
   I = K
   M = M+1
ELSE
   IL(M) = K
   IU(M) = J
   J = L
   M = M+1
END IF
GO TO 70
!
!     Begin again on another portion of the unsorted array
!
60 M = M-1
IF (M .EQ. 0) GO TO 190
I = IL(M)
J = IU(M)
!
70 IF (J-I .GE. 1) GO TO 30
IF (I .EQ. 1) GO TO 20
I = I-1
!
80 I = I+1
IF (I .EQ. J) GO TO 60
T = X(I+1)
IF (X(I) .LE. T) GO TO 80
K = I
!
90 X(K+1) = X(K)
K = K-1
IF (T .LT. X(K)) GO TO 90
X(K+1) = T
GO TO 80
!
!     Sort X and carry Y along
!
100 M = 1
I = 1
J = NN
R = 0.375E0
!
110 IF (I .EQ. J) GO TO 150
IF (R .LE. 0.5898437E0) THEN
   R = R+3.90625E-2
ELSE
   R = R-0.21875E0
END IF
!
120 K = I
!
!     Select a central element of the array and save it in location T
!
IJ = I + INT((J-I)*R)
T = X(IJ)
TY = Y(IJ)
!
!     If first element of array is greater than T, interchange with T
!
IF (X(I) .GT. T) THEN
   X(IJ) = X(I)
   X(I) = T
   T = X(IJ)
   Y(IJ) = Y(I)
   Y(I) = TY
   TY = Y(IJ)
END IF
L = J
!
!     If last element of array is less than T, interchange with T
!
IF (X(J) .LT. T) THEN
   X(IJ) = X(J)
   X(J) = T
   T = X(IJ)
   Y(IJ) = Y(J)
   Y(J) = TY
   TY = Y(IJ)
!
!        If first element of array is greater than T, interchange with T
!
   IF (X(I) .GT. T) THEN
      X(IJ) = X(I)
      X(I) = T
      T = X(IJ)
      Y(IJ) = Y(I)
      Y(I) = TY
      TY = Y(IJ)
   END IF
END IF
!
!     Find an element in the second half of the array which is smaller
!     than T
!
130 L = L-1
IF (X(L) .GT. T) GO TO 130
!
!     Find an element in the first half of the array which is greater
!     than T
!
140 K = K+1
IF (X(K) .LT. T) GO TO 140
!
!     Interchange these elements
!
IF (K .LE. L) THEN
   TT = X(L)
   X(L) = X(K)
   X(K) = TT
   TTY = Y(L)
   Y(L) = Y(K)
   Y(K) = TTY
   GO TO 130
END IF
!
!     Save upper and lower subscripts of the array yet to be sorted
!
IF (L-I .GT. J-K) THEN
   IL(M) = I
   IU(M) = L
   I = K
   M = M+1
ELSE
   IL(M) = K
   IU(M) = J
   J = L
   M = M+1
END IF
GO TO 160
!
!     Begin again on another portion of the unsorted array
!
150 M = M-1
IF (M .EQ. 0) GO TO 190
I = IL(M)
J = IU(M)
!
160 IF (J-I .GE. 1) GO TO 120
IF (I .EQ. 1) GO TO 110
I = I-1
!
170 I = I+1
IF (I .EQ. J) GO TO 150
T = X(I+1)
TY = Y(I+1)
IF (X(I) .LE. T) GO TO 170
K = I
!
180 X(K+1) = X(K)
Y(K+1) = Y(K)
K = K-1
IF (T .LT. X(K)) GO TO 180
X(K+1) = T
Y(K+1) = TY
GO TO 170
!
!     Clean up
!
190 IF (KFLAG .LE. -1) THEN
   DO  I=1,NN
      X(I) = -X(I)
   END DO
END IF
RETURN
END SUBROUTINE SSORT
