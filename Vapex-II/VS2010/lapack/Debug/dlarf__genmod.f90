        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:42 2014
        MODULE DLARF__genmod
          INTERFACE 
            SUBROUTINE DLARF(SIDE,M,N,V,INCV,TAU,C,LDC,WORK)
              INTEGER(KIND=4) :: LDC
              CHARACTER(LEN=1) :: SIDE
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: V(*)
              INTEGER(KIND=4) :: INCV
              REAL(KIND=8) :: TAU
              REAL(KIND=8) :: C(LDC,*)
              REAL(KIND=8) :: WORK(*)
            END SUBROUTINE DLARF
          END INTERFACE 
        END MODULE DLARF__genmod
