        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:47 2014
        MODULE DLARFX__genmod
          INTERFACE 
            SUBROUTINE DLARFX(SIDE,M,N,V,TAU,C,LDC,WORK)
              INTEGER(KIND=4) :: LDC
              CHARACTER(LEN=1) :: SIDE
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: V(*)
              REAL(KIND=8) :: TAU
              REAL(KIND=8) :: C(LDC,*)
              REAL(KIND=8) :: WORK(*)
            END SUBROUTINE DLARFX
          END INTERFACE 
        END MODULE DLARFX__genmod
