        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:43 2014
        MODULE DLALN2__genmod
          INTERFACE 
            SUBROUTINE DLALN2(LTRANS,NA,NW,SMIN,CA,A,LDA,D1,D2,B,LDB,WR,&
     &WI,X,LDX,SCALE,XNORM,INFO)
              INTEGER(KIND=4) :: LDX
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDA
              LOGICAL(KIND=4) :: LTRANS
              INTEGER(KIND=4) :: NA
              INTEGER(KIND=4) :: NW
              REAL(KIND=8) :: SMIN
              REAL(KIND=8) :: CA
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: D1
              REAL(KIND=8) :: D2
              REAL(KIND=8) :: B(LDB,*)
              REAL(KIND=8) :: WR
              REAL(KIND=8) :: WI
              REAL(KIND=8) :: X(LDX,*)
              REAL(KIND=8) :: SCALE
              REAL(KIND=8) :: XNORM
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DLALN2
          END INTERFACE 
        END MODULE DLALN2__genmod
