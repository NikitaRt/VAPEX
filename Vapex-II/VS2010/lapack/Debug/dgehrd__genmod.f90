        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:38 2014
        MODULE DGEHRD__genmod
          INTERFACE 
            SUBROUTINE DGEHRD(N,ILO,IHI,A,LDA,TAU,WORK,LWORK,INFO)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ILO
              INTEGER(KIND=4) :: IHI
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: TAU(*)
              REAL(KIND=8) :: WORK(LWORK)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGEHRD
          END INTERFACE 
        END MODULE DGEHRD__genmod
