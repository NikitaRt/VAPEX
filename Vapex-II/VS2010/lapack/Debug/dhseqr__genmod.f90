        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:40 2014
        MODULE DHSEQR__genmod
          INTERFACE 
            SUBROUTINE DHSEQR(JOB,COMPZ,N,ILO,IHI,H,LDH,WR,WI,Z,LDZ,WORK&
     &,LWORK,INFO)
              INTEGER(KIND=4) :: LDZ
              INTEGER(KIND=4) :: LDH
              CHARACTER(LEN=1) :: JOB
              CHARACTER(LEN=1) :: COMPZ
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ILO
              INTEGER(KIND=4) :: IHI
              REAL(KIND=8) :: H(LDH,*)
              REAL(KIND=8) :: WR(*)
              REAL(KIND=8) :: WI(*)
              REAL(KIND=8) :: Z(LDZ,*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DHSEQR
          END INTERFACE 
        END MODULE DHSEQR__genmod
