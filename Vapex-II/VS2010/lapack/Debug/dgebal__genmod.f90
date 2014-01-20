        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:46 2014
        MODULE DGEBAL__genmod
          INTERFACE 
            SUBROUTINE DGEBAL(JOB,N,A,LDA,ILO,IHI,SCALE,INFO)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: JOB
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: ILO
              INTEGER(KIND=4) :: IHI
              REAL(KIND=8) :: SCALE(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGEBAL
          END INTERFACE 
        END MODULE DGEBAL__genmod
