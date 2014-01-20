        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:36 2014
        MODULE DTREXC__genmod
          INTERFACE 
            SUBROUTINE DTREXC(COMPQ,N,T,LDT,Q,LDQ,IFST,ILST,WORK,INFO)
              INTEGER(KIND=4) :: LDQ
              INTEGER(KIND=4) :: LDT
              CHARACTER(LEN=1) :: COMPQ
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: T(LDT,*)
              REAL(KIND=8) :: Q(LDQ,*)
              INTEGER(KIND=4) :: IFST
              INTEGER(KIND=4) :: ILST
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DTREXC
          END INTERFACE 
        END MODULE DTREXC__genmod
