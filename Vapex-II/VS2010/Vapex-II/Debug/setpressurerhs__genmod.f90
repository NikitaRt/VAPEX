        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:39:59 2014
        MODULE SETPRESSURERHS__genmod
          INTERFACE 
            SUBROUTINE SETPRESSURERHS(QCOEFMATRIX,RHS,R1X,Q1X,RNX,QNX,  &
     &R1Y,Q1Y,RMY,QMY,NVARS,NRHS,N,M,N9,M9)
              INTEGER(KIND=4), INTENT(IN) :: M9
              INTEGER(KIND=4), INTENT(IN) :: N9
              INTEGER(KIND=4), INTENT(IN) :: NRHS
              INTEGER(KIND=4), INTENT(IN) :: NVARS
              REAL(KIND=8), INTENT(INOUT) :: QCOEFMATRIX(NVARS,NRHS,N9, &
     &M9)
              REAL(KIND=8), INTENT(INOUT) :: RHS(N9,M9)
              REAL(KIND=8), INTENT(IN) :: R1X(M9)
              REAL(KIND=8), INTENT(IN) :: Q1X(M9)
              REAL(KIND=8), INTENT(IN) :: RNX(M9)
              REAL(KIND=8), INTENT(IN) :: QNX(M9)
              REAL(KIND=8), INTENT(IN) :: R1Y(N9)
              REAL(KIND=8), INTENT(IN) :: Q1Y(N9)
              REAL(KIND=8), INTENT(IN) :: RMY(N9)
              REAL(KIND=8), INTENT(IN) :: QMY(N9)
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: M
            END SUBROUTINE SETPRESSURERHS
          END INTERFACE 
        END MODULE SETPRESSURERHS__genmod