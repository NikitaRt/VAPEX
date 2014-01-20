        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:49 2014
        MODULE DLAHQR__genmod
          INTERFACE 
            SUBROUTINE DLAHQR(WANTT,WANTZ,N,ILO,IHI,H,LDH,WR,WI,ILOZ,   &
     &IHIZ,Z,LDZ,INFO)
              INTEGER(KIND=4) :: LDZ
              INTEGER(KIND=4) :: LDH
              LOGICAL(KIND=4) :: WANTT
              LOGICAL(KIND=4) :: WANTZ
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ILO
              INTEGER(KIND=4) :: IHI
              REAL(KIND=8) :: H(LDH,*)
              REAL(KIND=8) :: WR(*)
              REAL(KIND=8) :: WI(*)
              INTEGER(KIND=4) :: ILOZ
              INTEGER(KIND=4) :: IHIZ
              REAL(KIND=8) :: Z(LDZ,*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DLAHQR
          END INTERFACE 
        END MODULE DLAHQR__genmod
