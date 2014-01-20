        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:39 2014
        MODULE DTREVC__genmod
          INTERFACE 
            SUBROUTINE DTREVC(JOB,HOWMNY,SELECT,N,T,LDT,VL,LDVL,VR,LDVR,&
     &MM,M,WORK,INFO)
              INTEGER(KIND=4) :: LDVR
              INTEGER(KIND=4) :: LDVL
              INTEGER(KIND=4) :: LDT
              CHARACTER(LEN=1) :: JOB
              CHARACTER(LEN=1) :: HOWMNY
              LOGICAL(KIND=4) :: SELECT(*)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: T(LDT,*)
              REAL(KIND=8) :: VL(LDVL,*)
              REAL(KIND=8) :: VR(LDVR,*)
              INTEGER(KIND=4) :: MM
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DTREVC
          END INTERFACE 
        END MODULE DTREVC__genmod
