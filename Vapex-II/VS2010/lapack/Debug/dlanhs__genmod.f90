        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:51 2014
        MODULE DLANHS__genmod
          INTERFACE 
            FUNCTION DLANHS(NORM,N,A,LDA,WORK)
              INTEGER(KIND=4) :: LDA
              CHARACTER(LEN=1) :: NORM
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: WORK(*)
              REAL(KIND=8) :: DLANHS
            END FUNCTION DLANHS
          END INTERFACE 
        END MODULE DLANHS__genmod
