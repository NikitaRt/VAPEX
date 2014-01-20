        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:50 2014
        MODULE DLASY2__genmod
          INTERFACE 
            SUBROUTINE DLASY2(LTRANL,LTRANR,ISGN,N1,N2,TL,LDTL,TR,LDTR,B&
     &,LDB,SCALE,X,LDX,XNORM,INFO)
              INTEGER(KIND=4) :: LDX
              INTEGER(KIND=4) :: LDB
              INTEGER(KIND=4) :: LDTR
              INTEGER(KIND=4) :: LDTL
              LOGICAL(KIND=4) :: LTRANL
              LOGICAL(KIND=4) :: LTRANR
              INTEGER(KIND=4) :: ISGN
              INTEGER(KIND=4) :: N1
              INTEGER(KIND=4) :: N2
              REAL(KIND=8) :: TL(LDTL,*)
              REAL(KIND=8) :: TR(LDTR,*)
              REAL(KIND=8) :: B(LDB,*)
              REAL(KIND=8) :: SCALE
              REAL(KIND=8) :: X(LDX,*)
              REAL(KIND=8) :: XNORM
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DLASY2
          END INTERFACE 
        END MODULE DLASY2__genmod
