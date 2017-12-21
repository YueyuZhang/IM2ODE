
MODULE EIGENVAL_mod
  USE KINDS
  IMPLICIT NONE
  TYPE eigen_t
     INTEGER:: NIONS
     INTEGER:: ISPIN,NELECT,NKPTS,NB_TOT
     REAL(dp),POINTER:: VKPT(:,:)=>NULL(),WTKPT(:)=>NULL()
     REAL(dp),POINTER:: CELTOT(:,:,:)=>NULL()      !(NB_TOT,NKPTS,ISPIN)
  END TYPE eigen_t
CONTAINS
  SUBROUTINE read_EIGENVAL(iu,eigen_data,lgood)
    IMPLICIT NONE
    INTEGER,INTENT(in):: iu
    TYPE(eigen_t),INTENT(out):: eigen_data
    LOGICAL,INTENT(out):: lgood 
    !local var
    INTEGER:: it,NK,N
    !=======================
    lgood=.TRUE.
    
    READ(iu,'(4I5)',ERR=100,END=100)eigen_data%NIONS,it,it,eigen_data%ISPIN

    READ(iu,'(5E15.7)',ERR=100,END=100)
    READ(iu,*,ERR=100,END=100)
    READ(iu,*,ERR=100,END=100)
    READ(iu,*,ERR=100,END=100)
    READ(iu,'(3I5)',ERR=100,END=100) eigen_data%NELECT,eigen_data%NKPTS,eigen_data%NB_TOT
    ALLOCATE(eigen_data%VKPT(3,eigen_data%NKPTS))
    ALLOCATE(eigen_data%WTKPT(eigen_data%NKPTS))
    ALLOCATE(eigen_data%CELTOT(eigen_data%NB_TOT,eigen_data%NKPTS,eigen_data%ISPIN))

    DO NK=1,eigen_data%NKPTS
       READ(iu,*,ERR=100,END=100)
       READ(iu,'(4E15.7)',ERR=100,END=100) eigen_data%VKPT(1,NK),eigen_data%VKPT(2,NK),eigen_data%VKPT(3,NK),&
            & eigen_data%WTKPT(NK)
       DO N=1,eigen_data%NB_TOT
          IF (eigen_data%ISPIN==1) READ(iu,852,ERR=100,END=100)eigen_data%CELTOT(N,NK,1) 
          IF (eigen_data%ISPIN==2) &
               READ(iu,8852,ERR=100,END=100)eigen_data%CELTOT(N,NK,1),eigen_data%CELTOT(N,NK,2)
       ENDDO
    ENDDO
    
!852 FORMAT(1X,I3,4X,F10.4)
!8852 FORMAT(1X,I3,4X,F10.4,2X,F10.4)
852 FORMAT(1X,3X,4X,F10.4)
8852 FORMAT(1X,3X,4X,F10.4,2X,F10.4)

    RETURN

100 CONTINUE
    lgood=.FALSE.
  END SUBROUTINE read_EIGENVAL

  SUBROUTINE del_eigen_t(eigen_data)
    IMPLICIT NONE
    TYPE(eigen_t),INTENT(inout):: eigen_data

    IF(ASSOCIATED(eigen_data%VKPT))DEALLOCATE(eigen_data%VKPT)
    
    IF(ASSOCIATED(eigen_data%WTKPT))DEALLOCATE(eigen_data%WTKPT)

    IF(ASSOCIATED(eigen_data%CELTOT))DEALLOCATE(eigen_data%CELTOT)

  END SUBROUTINE del_eigen_t


END MODULE EIGENVAL_mod
MODULE sort_mod
  USE KINDS
  IMPLICIT NONE
  INTERFACE SSORT
     MODULE PROCEDURE SSORT_I,SSORT_F
  END INTERFACE
CONTAINS
  SUBROUTINE SSORT_I (X, Y, N, KFLAG)
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
    !     .. Scalar Arguments ..
    INTEGER KFLAG, N
    !     .. Array Arguments ..
    ! xhj change here
    !      REAL X(*), Y(*)
    INTEGER :: X(*)
    REAL(dp):: Y(*)
    !xhj change end here
    !     .. Local Scalars ..
    REAL(dp):: R, T, TT, TTY, TY
    INTEGER I, IJ, J, K, KK, L, M, NN
    !     .. Local Arrays ..
    INTEGER IL(21), IU(21)
    !     .. External Subroutines ..
    !     None
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, INT
    !***FIRST EXECUTABLE STATEMENT  SSORT
    NN = N
    IF (NN .LT. 1) THEN
       PRINT *, &
            &      'The number of values to be sorted is not positive.'
       RETURN
    ENDIF
    !
    KK = ABS(KFLAG)
    IF (KK.NE.1 .AND. KK.NE.2) THEN
       PRINT *, &
            &      'The sort control parameter, K, is not 2, 1, -1, or -2.'
       RETURN
    ENDIF
    !
    !     Alter array X to get decreasing order if needed
    !
    IF (KFLAG .LE. -1) THEN
       DO 10 I=1,NN
          X(I) = -X(I)
10     END DO
       !   10    CONTINUE
    ENDIF
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
20  IF (I .EQ. J) GO TO 60
    IF (R .LE. 0.5898437E0) THEN
       R = R+3.90625E-2
    ELSE
       R = R-0.21875E0
    ENDIF
    !
30  K = I
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
    ENDIF
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
       ENDIF
    ENDIF
    !
    !     Find an element in the second half of the array which is smaller
    !     than T
    !
40  L = L-1
    IF (X(L) .GT. T) GO TO 40
    !
    !     Find an element in the first half of the array which is greater
    !     than T
    !
50  K = K+1
    IF (X(K) .LT. T) GO TO 50
    !
    !     Interchange these elements
    !
    IF (K .LE. L) THEN
       TT = X(L)
       X(L) = X(K)
       X(K) = TT
       GO TO 40
    ENDIF
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
    ENDIF
    GO TO 70
    !
    !     Begin again on another portion of the unsorted array
    !
60  M = M-1
    IF (M .EQ. 0) GO TO 190
    I = IL(M)
    J = IU(M)
    !
70  IF (J-I .GE. 1) GO TO 30
    IF (I .EQ. 1) GO TO 20
    I = I-1
    !
80  I = I+1
    IF (I .EQ. J) GO TO 60
    T = X(I+1)
    IF (X(I) .LE. T) GO TO 80
    K = I
    !
90  X(K+1) = X(K)
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
    ENDIF
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
    ENDIF
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
       ENDIF
    ENDIF
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
    ENDIF
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
    ENDIF
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
       DO 200 I=1,NN
          X(I) = -X(I)
200    END DO
       !  200    CONTINUE
    ENDIF
    RETURN
  END SUBROUTINE SSORT_I

  SUBROUTINE SSORT_F (X, Y, N, KFLAG)
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
    !     .. Scalar Arguments ..
    INTEGER KFLAG, N
    !     .. Array Arguments ..
    ! xhj change here
    !      REAL X(*), Y(*)
    !      INTEGER :: X(*)
    !change back to real
    REAL(dp):: X(*)
    REAL(dp):: Y(*)
    !xhj change end here
    !     .. Local Scalars ..
    REAL(dp):: R, T, TT, TTY, TY
    INTEGER I, IJ, J, K, KK, L, M, NN
    !     .. Local Arrays ..
    INTEGER IL(21), IU(21)
    !     .. External Subroutines ..
    !     None
    !     .. Intrinsic Functions ..
    INTRINSIC ABS, INT
    !***FIRST EXECUTABLE STATEMENT  SSORT
    NN = N
    IF (NN .LT. 1) THEN
       PRINT *, &
            &      'The number of values to be sorted is not positive.'
       RETURN
    ENDIF
    !
    KK = ABS(KFLAG)
    IF (KK.NE.1 .AND. KK.NE.2) THEN
       PRINT *, &
            &      'The sort control parameter, K, is not 2, 1, -1, or -2.'
       RETURN
    ENDIF
    !
    !     Alter array X to get decreasing order if needed
    !
    IF (KFLAG .LE. -1) THEN
       DO 10 I=1,NN
          X(I) = -X(I)
10     END DO
       !   10    CONTINUE
    ENDIF
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
20  IF (I .EQ. J) GO TO 60
    IF (R .LE. 0.5898437E0) THEN
       R = R+3.90625E-2
    ELSE
       R = R-0.21875E0
    ENDIF
    !
30  K = I
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
    ENDIF
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
       ENDIF
    ENDIF
    !
    !     Find an element in the second half of the array which is smaller
    !     than T
    !
40  L = L-1
    IF (X(L) .GT. T) GO TO 40
    !
    !     Find an element in the first half of the array which is greater
    !     than T
    !
50  K = K+1
    IF (X(K) .LT. T) GO TO 50
    !
    !     Interchange these elements
    !
    IF (K .LE. L) THEN
       TT = X(L)
       X(L) = X(K)
       X(K) = TT
       GO TO 40
    ENDIF
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
    ENDIF
    GO TO 70
    !
    !     Begin again on another portion of the unsorted array
    !
60  M = M-1
    IF (M .EQ. 0) GO TO 190
    I = IL(M)
    J = IU(M)
    !
70  IF (J-I .GE. 1) GO TO 30
    IF (I .EQ. 1) GO TO 20
    I = I-1
    !
80  I = I+1
    IF (I .EQ. J) GO TO 60
    T = X(I+1)
    IF (X(I) .LE. T) GO TO 80
    K = I
    !
90  X(K+1) = X(K)
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
    ENDIF
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
    ENDIF
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
       ENDIF
    ENDIF
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
    ENDIF
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
    ENDIF
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
       DO 200 I=1,NN
          X(I) = -X(I)
200    ENDDO
       !200  CONTINUE
    ENDIF
    RETURN
  END SUBROUTINE SSORT_F
END MODULE sort_mod
MODULE OUTCAR_mod
  USE KINDS
  IMPLICIT NONE 
  TYPE outcar_t
     INTEGER:: na
     !Startparameter for this run
     INTEGER:: ISPIN 
     LOGICAL::  LNONCOLLINEAR
     LOGICAL:: LSORBIT 
     REAL(dp):: NELECT
     !output
     REAL(dp):: EFERMI
  END TYPE outcar_t
  
CONTAINS
  SUBROUTINE read_OUTCAR(iu,out_data,lgood)
    !read na
    INTEGER,INTENT(in):: iu
    TYPE(outcar_t),INTENT(out):: out_data
    LOGICAL,INTENT(out):: lgood 
    !local var
    CHARACTER(LEN=100):: tag,ch_read
    LOGICAL:: lmatch 
    !==================================================
    lgood=.TRUE.

    tag='NIONS ='
    CALL read_val(iu,tag,ch_read,0,lmatch)    
    IF(lmatch)THEN
       READ(ch_read,*)out_data%na
    ELSE
       PRINT*,"read_OUTCAR:",tag, " not found"
       !STOP 
       lgood=.FALSE.
       RETURN
    END IF

    tag='ISPIN  ='
    CALL read_val(iu,tag,ch_read,0,lmatch)
    IF(lmatch)THEN
       READ(ch_read,*)out_data%ISPIN
    ELSE
       PRINT*,"read_OUTCAR:",tag, " not found"
       lgood=.FALSE.
       RETURN
    END IF

    tag='LNONCOLLINEAR ='
    CALL read_val(iu,tag,ch_read,0,lmatch)
    IF(lmatch)THEN
       READ(ch_read,*)out_data%LNONCOLLINEAR
    ELSE
       PRINT*,"read_OUTCAR:",tag, " not found"
       lgood=.FALSE.
       RETURN
    END IF    

    tag='LSORBIT ='
    CALL read_val(iu,tag,ch_read,0,lmatch)
    IF(lmatch)THEN
       READ(ch_read,*)out_data%LSORBIT
    ELSE
       PRINT*,"read_OUTCAR:",tag, " not found"
       lgood=.FALSE.
       RETURN
    END IF

    tag='NELECT ='
    CALL read_val(iu,tag,ch_read,0,lmatch)
    IF(lmatch)THEN
       READ(ch_read,*)out_data%NELECT
    ELSE
       PRINT*,"read_OUTCAR:",tag, " not found"
       lgood=.FALSE.
       RETURN
    END IF
    
  END SUBROUTINE read_OUTCAR

  SUBROUTINE read_val(iu,tag,ch_read,ind,lmatch)
    IMPLICIT NONE 
    INTEGER,INTENT(in):: iu,ind
    !ind means which one you need
    !sometimes, in a file, there are several occurences of the same tag
    CHARACTER(LEN=*),INTENT(in):: tag 
    CHARACTER(LEN=*),INTENT(out):: ch_read    
    LOGICAL,INTENT(out):: lmatch   
    !local var
    CHARACTER(LEN=200)::linestring
    LOGICAL:: found,I_OPENED
    !==========================================
    INQUIRE (iu, OPENED=I_OPENED)
    
    IF( .NOT. I_OPENED)THEN
       STOP "read_val: file not open"
    END IF
    
    IF(ind .NE. 0)THEN
       IF(ind .EQ. -1)THEN
          !the last occurence of the tag in this file
          found=.FALSE.
          DO 
             READ(iu,'(A)',END=100)linestring

             CALL match_line(linestring,tag,ch_read,lmatch)
             IF(lmatch)THEN
                found=.TRUE.
             END IF
          END DO
100       CONTINUE
          
          IF(.NOT. found)THEN
             REWIND(iu)
             DO 
                READ(iu,'(A)',END=200)linestring
                
                CALL match_line(linestring,tag,ch_read,lmatch)
                IF(lmatch)THEN
                   found=.TRUE.
                END IF
             END DO
200          CONTINUE
          END IF

          IF(found)lmatch=.TRUE.
       END IF
    ELSE
       !the default behavior is to find one place with the tag
       DO 
          READ(iu,'(A)',END=300)linestring
          
          CALL match_line(linestring,tag,ch_read,lmatch)          
          IF(lmatch)RETURN
       END DO
300    CONTINUE
       IF(.NOT. lmatch)THEN
          REWIND(iu)
          !another try
          DO 
             READ(iu,'(A)',END=400)linestring
             
             CALL match_line(linestring,tag,ch_read,lmatch)             
             IF(lmatch)RETURN
          END DO          
400       CONTINUE
       END IF
    END IF
  END SUBROUTINE read_val

  SUBROUTINE match_line(linestring,tag,ch_read,lmatch)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(in):: linestring,tag
    CHARACTER(LEN=*),INTENT(out):: ch_read
    LOGICAL,INTENT(out):: lmatch
    !local var
    INTEGER:: j,k,istart,n 
    !===============================================    
    ch_read=''
       
    j = INDEX(linestring, TRIM(tag))
    IF (j /= 0) THEN
       DO k=j+LEN_TRIM(tag),LEN_TRIM(linestring)
          IF(linestring(k:k) .NE. ' ')THEN
             istart=k
             EXIT
          END IF
       END DO
       
       n=0
       DO k=istart,LEN_TRIM(linestring)
          IF(linestring(k:k) .NE. ' ')THEN
             n=n+1
             ch_read(n:n)=linestring(k:k)
          ELSE
             EXIT
          END IF
       END DO

       lmatch=.TRUE.
    ELSE
       lmatch=.FALSE.
    END IF
  END SUBROUTINE match_line

END MODULE OUTCAR_mod


MODULE Elec_Str_mod
  IMPLICIT NONE
CONTAINS
  SUBROUTINE ES_gap(tag)
    USE EIGENVAL_mod,ONLY: read_EIGENVAL,eigen_t,del_eigen_t
    USE sort_mod,ONLY: SSORT
    USE OUTCAR_mod,ONLY: read_OUTCAR,outcar_t
    USE parameters, ONLY : pstruct
    USE KINDS
    !USE parameters ,      ONLY : Icode
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: tag
    TYPE(eigen_t):: eigen_data
    TYPE(outcar_t):: out_data
    INTEGER:: NK,is,iu,NB_OCC
    LOGICAL:: lmetal,lgoodeig,lgoodOUTCAR,lexist
    REAL(dp):: Eg_d,Eg_id,VBM,CBM,tran_direct_gap(3)
    !=================================================
    !IF(icode .EQ. 5)THEN
    !   RETURN
    !END IF

    lmetal=.FALSE.
    
    !CALL io_assign(iu)
    iu = 4411
    OPEN(iu,file='EIGENVAL')
    CALL read_EIGENVAL(iu,eigen_data,lgoodeig)
    !CALL io_close(iu)    
    close(iu)

    IF(lgoodeig)THEN
       !sort energy in eigen_data%CELTOT                                                   
       DO NK=1,eigen_data%NKPTS
          DO is=1,eigen_data%ISPIN
             CALL SSORT(eigen_data%CELTOT(:,NK,is), eigen_data%CELTOT(:,NK,is), &
                  & eigen_data%NB_TOT, 1)
          END DO
       END DO
    END IF
    
    !now read OUTCAR
    !CALL io_assign(iu)
    iu = 4412
    OPEN(iu,file='OUTCAR')    
    CALL read_OUTCAR(iu,out_data,lgoodOUTCAR)
    !CALL io_close(iu)    
    close(iu)

    !now calculate the band gap
    !for the global band gap (including direct and indirect gap)
    IF(lgoodeig .AND. lgoodOUTCAR)THEN
       IF(out_data%LNONCOLLINEAR)THEN
          NB_OCC=eigen_data%NELECT
       ELSEIF(out_data%ISPIN == 1)THEN
          IF(MOD(eigen_data%NELECT,2) .EQ. 1)THEN
             !we can not have a gap in the case of spin unpolarized calculation with odd number of electrons
             lmetal=.TRUE.
          ELSE
             NB_OCC=eigen_data%NELECT/2
          END IF
       ELSEIF(out_data%ISPIN == 2)THEN
          STOP "ES_gap: spin polarized case not implemented yet"
       END IF
    END IF
    
    IF(lgoodeig .AND. lgoodOUTCAR)THEN
       IF(.NOT. lmetal)THEN
          Eg_d=1.0d9
          Eg_id=1.0d9
          !spin polarized system is not supported yet
          IF(eigen_data%ISPIN .EQ. 1)THEN
             DO NK=1,eigen_data%NKPTS
                Eg_d=MIN(eigen_data%CELTOT(NB_OCC+1,NK,1)-eigen_data%CELTOT(NB_OCC,NK,1), Eg_d)
             END DO
             VBM=MAXVAL(eigen_data%CELTOT(NB_OCC,1:eigen_data%NKPTS,1))
             CBM=MINVAL(eigen_data%CELTOT(NB_OCC+1,1:eigen_data%NKPTS,1))
             Eg_id=MINVAL(eigen_data%CELTOT(NB_OCC+1,:,1))-MAXVAL(eigen_data%CELTOT(NB_OCC,:,1))
             
!!$          IF(ABS(Eg_d-Eg_id) .LT. 1.0d-9)THEN
!!$             WRITE(6,'(a,f20.12,a)')"Direct band gap: ",Eg_d, " eV"
!!$          ELSE
!!$             WRITE(6,'(a,f20.12,a)')"Direct band gap: ",Eg_d, " eV"
!!$             WRITE(6,'(a,f20.12,a)')"Indirect band gap: ",Eg_id, " eV"
!!$          END IF
!!$          WRITE(6,'(a,f20.12)')"VBM:", VBM
!!$          WRITE(6,'(a,f20.12)')"CBM:", CBM
          END IF

          IF(Eg_d .LT. 0)THEN
             Eg_d=0.0d0
          END IF
          
          !IF(Eg_id .LT. 0)THEN
          !   Eg_id=0.0d0
          !END IF
       ELSE
          !This should not happen
          Eg_d=10000 !number of electrons is odd
          Eg_id=10000
       END IF
    ELSEIF( (.NOT. lgoodeig) .OR. (.NOT. lgoodOUTCAR) )THEN
       Eg_d=1.0d10 !vasp results are not complete
       Eg_id=1.0d10
    END IF    

    !read Matrix_elements if it exists    
    CALL read_dipole_mat(tran_direct_gap,lexist)
    
    IF(lexist)THEN
       WRITE(6,'(4x,a,3f10.3)')"Dipole_Matrix_Element: ",tran_direct_gap
       IF(SQRT(DOT_PRODUCT(tran_direct_gap,tran_direct_gap)) .LT. 1.0d-3)THEN
          !vanishingly small transition
          Eg_d=Eg_d+5000
          Eg_id=Eg_id+5000
       END IF
    ELSE
       Eg_d=7000
       Eg_id=7000
    END IF

    CALL del_eigen_t(eigen_data)
    
    !CALL io_assign(iu)
    !OPEN( iu, file="Eg.dat" )
    
    !WRITE(iu,'(2f24.8)')Eg_id,Eg_d
    !CALL io_close(iu)
    pstruct(tag) % Eg_id = Eg_id
    
    pstruct(tag) % Eg_d = Eg_d
    
  END SUBROUTINE ES_gap

  SUBROUTINE read_dipole_mat(tran_direct_gap,lexist)
    USE KINDS
    USE parameters, only : Es_opt
    IMPLICIT NONE
    REAL(dp),INTENT(out):: tran_direct_gap(3)
    LOGICAL,INTENT(out):: lexist 
    !local var
    INTEGER:: iu,iv,ic,ik
    REAL(dp):: VK(3),eigv,eigc,nabij(3),Eg_d
    CHARACTER(len=6) :: str1, str2
    !=====================================
    INQUIRE(file='Matrix_elements', exist= lexist )    

    IF(.NOT. lexist)THEN
       tran_direct_gap=0.0d0
       RETURN
    END IF

    !CALL io_assign(iu)
    iu = 4415
    OPEN(iu,file='Matrix_elements')
    READ(iu,*)
    
    tran_direct_gap=0.0d0
    
    eg_d=1.0d10

    DO 
       READ(iu,'(1X,2(i3,f9.4,1x), 3f10.4,1x,i4,1x,A1,3f6.2,A2)',END=100) &
            & iv, eigv, ic, eigc, nabij(1:3), ik, str1, VK(1:3), str2
       IF(eigc-eigv .LT. eg_d)THEN
          eg_d=eigc-eigv
          tran_direct_gap=nabij
       END IF
       ! zyy add here
       IF(eigc-eigv .LT. Es_opt) THEN
          IF(SQRT(DOT_PRODUCT(nabij,nabij))>SQRT(DOT_PRODUCT(tran_direct_gap,tran_direct_gap))) THEN
             tran_direct_gap = nabij
          END IF
       END IF
       ! zyy add end here
    END DO

100 CONTINUE    

    !CALL io_close(iu)
    close(iu)
  END SUBROUTINE read_dipole_mat


END MODULE Elec_Str_mod
