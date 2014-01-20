        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:44 2014
        MODULE DGEEV__genmod
          INTERFACE 
            SUBROUTINE DGEEV(JOBVL,JOBVR,N,A,LDA,WR,WI,VL,LDVL,VR,LDVR, &
     &WORK,LWORK,INFO)
              INTEGER(KIND=4) :: LDVR
              INTEGER(KIND=4) :: LDVL
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: JOBVL
              CHARACTER(LEN=1) :: JOBVR
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: WR(*)
              REAL(KIND=8) :: WI(*)
              REAL(KIND=8) :: VL(LDVL,*)
              REAL(KIND=8) :: VR(LDVR,*)
              REAL(KIND=8) :: WORK(*)
              INTEGER(KIND=4) :: LWORK
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGEEV
          END INTERFACE 
        END MODULE DGEEV__genmod
