        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:45 2014
        MODULE DGETRF__genmod
          INTERFACE 
            SUBROUTINE DGETRF(M,N,A,LDA,IPIV,INFO)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IPIV(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGETRF
          END INTERFACE 
        END MODULE DGETRF__genmod
