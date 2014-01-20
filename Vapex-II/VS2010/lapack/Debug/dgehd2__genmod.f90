        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:41 2014
        MODULE DGEHD2__genmod
          INTERFACE 
            SUBROUTINE DGEHD2(N,ILO,IHI,A,LDA,TAU,WORK,INFO)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ILO
              INTEGER(KIND=4) :: IHI
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: TAU(*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGEHD2
          END INTERFACE 
        END MODULE DGEHD2__genmod
