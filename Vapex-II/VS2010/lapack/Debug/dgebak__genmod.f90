        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:49 2014
        MODULE DGEBAK__genmod
          INTERFACE 
            SUBROUTINE DGEBAK(JOB,SIDE,N,ILO,IHI,SCALE,M,V,LDV,INFO)
              INTEGER(KIND=4) :: LDV
              CHARACTER(LEN=1) :: JOB
              CHARACTER(LEN=1) :: SIDE
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ILO
              INTEGER(KIND=4) :: IHI
              REAL(KIND=8) :: SCALE(*)
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: V(LDV,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGEBAK
          END INTERFACE 
        END MODULE DGEBAK__genmod
