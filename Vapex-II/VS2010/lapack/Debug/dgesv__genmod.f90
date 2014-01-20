        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:52 2014
        MODULE DGESV__genmod
          INTERFACE 
            SUBROUTINE DGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NRHS
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IPIV(*)
              REAL(KIND=8) :: B(LDB,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGESV
          END INTERFACE 
        END MODULE DGESV__genmod
