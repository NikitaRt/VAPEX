        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:39:58 2014
        MODULE PRESSURECORRECTION__genmod
          INTERFACE 
            SUBROUTINE PRESSURECORRECTION(N,M,XX,YY,R1X,Q1X,RNX,QNX,R1Y,&
     &Q1Y,RMY,QMY,QCOEFMATRIX,NV,NR)
              INTEGER(KIND=4) :: NR
              INTEGER(KIND=4) :: NV
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: XX(N,M)
              REAL(KIND=8) :: YY(N,M)
              REAL(KIND=8), INTENT(IN) :: R1X(M)
              REAL(KIND=8), INTENT(IN) :: Q1X(M)
              REAL(KIND=8), INTENT(IN) :: RNX(M)
              REAL(KIND=8), INTENT(IN) :: QNX(M)
              REAL(KIND=8), INTENT(IN) :: R1Y(N)
              REAL(KIND=8), INTENT(IN) :: Q1Y(N)
              REAL(KIND=8), INTENT(IN) :: RMY(N)
              REAL(KIND=8), INTENT(IN) :: QMY(N)
              REAL(KIND=8), INTENT(IN) :: QCOEFMATRIX(NV,NR,N,M)
            END SUBROUTINE PRESSURECORRECTION
          END INTERFACE 
        END MODULE PRESSURECORRECTION__genmod
