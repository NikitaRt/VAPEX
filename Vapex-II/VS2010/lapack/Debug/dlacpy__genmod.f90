        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:51 2014
        MODULE DLACPY__genmod
          INTERFACE 
            SUBROUTINE DLACPY(UPLO,M,N,A,LDA,B,LDB)
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: UPLO
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: B(LDB,*)
            END SUBROUTINE DLACPY
          END INTERFACE 
        END MODULE DLACPY__genmod
