        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:45 2014
        MODULE DLAEXC__genmod
          INTERFACE 
            SUBROUTINE DLAEXC(WANTQ,N,T,LDT,Q,LDQ,J1,N1,N2,WORK,INFO)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDT
              LOGICAL(KIND=4) :: WANTQ
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: T(LDT,*)
              REAL(KIND=8) :: Q(LDQ,*)
              INTEGER(KIND=4) :: J1
              INTEGER(KIND=4) :: N1
              INTEGER(KIND=4) :: N2
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DLAEXC
          END INTERFACE 
        END MODULE DLAEXC__genmod
