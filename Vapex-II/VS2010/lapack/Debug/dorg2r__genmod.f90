        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:46 2014
        MODULE DORG2R__genmod
          INTERFACE 
            SUBROUTINE DORG2R(M,N,K,A,LDA,TAU,WORK,INFO)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: TAU(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DORG2R
          END INTERFACE 
        END MODULE DORG2R__genmod