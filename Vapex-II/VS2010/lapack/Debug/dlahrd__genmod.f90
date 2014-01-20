        !COMPILER-GENERATED INTERFACE MODULE: Mon Jan 20 14:28:46 2014
        MODULE DLAHRD__genmod
          INTERFACE 
            SUBROUTINE DLAHRD(N,K,NB,A,LDA,TAU,T,LDT,Y,LDY)
              INTEGER(KIND=4) :: LDY
              INTEGER(KIND=4) :: LDT
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: NB
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: A(LDA,*)
              REAL(KIND=8) :: TAU(NB)
              REAL(KIND=8) :: T(LDT,NB)
              REAL(KIND=8) :: Y(LDY,NB)
            END SUBROUTINE DLAHRD
          END INTERFACE 
        END MODULE DLAHRD__genmod
