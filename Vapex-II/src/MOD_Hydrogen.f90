!------------------------------------------------------------------!
!  Hydrogen release                                                !
!------------------------------------------------------------------!
!  $Id: MOD_Hydrogen.f90 5 2014-01-13 11:44:50Z Sergey $
!------------------------------------------------------------------!
MODULE HYDROGEN
  REAL(kind=8):: Gamma_a
  REAL(kind=8):: YH2_Gamma_a, TH2_Gamma_a
  REAL(kind=8):: RoH2, YH2, XH2 ! Partial density, mass and molar fraction
  REAL(kind=8):: RoH2n
  REAL(kind=8):: W_a ! Molecular weight of non-condensible gas

  DATA YH2_Gamma_a/1.D0/
  REAL(kind=8),PARAMETER:: omegaH2 = 1.D0
  REAL(kind=8),PARAMETER:: eps_H2 = 1.D-6
  INTEGER,PARAMETER:: niter_H2 = 20
CONTAINS

  SUBROUTINE SetHydrogenRelease(G_a, Y_H2, T_H2, Y_H2_Ga)
    IMPLICIT NONE
    REAL(kind=8), INTENT(in):: G_a, Y_H2, T_H2
    REAL(kind=8),OPTIONAL, INTENT(in):: Y_H2_Ga

    Gamma_a = G_a
    YH2 = Y_H2
    TH2_Gamma_a = T_H2
    IF(PRESENT(Y_H2_Ga)) YH2_Gamma_a = Y_H2_Ga
  END SUBROUTINE SetHydrogenRelease

  SUBROUTINE HydrogenRelease0D(tstep)
    USE GLOBALS
    USE LOCAL_VARS
    IMPLICIT NONE

    REAL(kind=8), INTENT(in):: tstep

    RoH2 = (RoH2n*Alphan + tstep*Gamma_a*YH2_Gamma_a)/MAX(Alpha, allim)
    YH2 = RoH2/Roa

  END SUBROUTINE HydrogenRelease0D
  !
  SUBROUTINE HYDROGEN_2D(aYH2,  &
       arH2, brH2, &
       au2, av2, ara, aa2, ba2, &
       allim, &
       gammaH2, &
       r12, rad, t, cpr, cpz, n, m, n9, m9, &
       niter, okprint)
    IMPLICIT NONE

    INTEGER:: n, m, n9, m9
    REAL(kind=8), INTENT(inout):: aYH2(n9,m9)
    REAL(kind=8), INTENT(inout):: arH2(n9,m9) 
    REAL(kind=8), INTENT(in):: brH2(n9,m9)
    REAL(kind=8), INTENT(in):: au2(n,m9), av2(n9,m), ara(n9,m9), &
         aa2(n9,m9), ba2(n9,m9), gammaH2(n9,m9)
    REAL(kind=8), INTENT(in):: allim
    REAL(kind=8):: r12(n), rad(n)
    REAL(kind=8):: t, cpr, cpz
    INTEGER:: niter
    LOGICAL:: okprint
    !
    ! Local variables
    !
    INTEGER:: i,j,iter
    REAL(kind=8):: Flux_R, Flux_L, Flux_T, Flux_B
    REAL(kind=8):: C,RHS, Yold, eps_reached
    REAL(kind=8):: deltaY, RoA2_lim, gamH2

    iter = 0
    DO iter = 1,niter
       ! Forward pass
       DO j = 2,m
          DO i = 2,n
             gamH2 = t*gammaH2(i,j)*YH2_Gamma_a

             ! Time derivative 
             C = aa2(i,j)
             RHS = brH2(i,j)*ba2(i,j)

             ! Source terms
             RHS = RHS +  gamH2

             ! Convective terms
             Flux_R = r12(i)/rad(i)*au2(i,j)
             Flux_L = r12(i-1)/rad(i)*au2(i-1,j)
             Flux_T = av2(i,j)
             Flux_B = av2(i,j-1)

             IF(au2(i,j) > 0.d0) THEN
                C = C + cpr*Flux_R*aa2(i,j)
             ELSE
                RHS = RHS - cpr*Flux_R*aa2(i+1,j)*arH2(i+1,j)
             ENDIF
             IF(au2(i-1,j) > 0.d0) THEN
                RHS = RHS + cpr*Flux_L*aa2(i-1,j)*arH2(i-1,j)
             ELSE
                C = C - cpr*Flux_L*aa2(i,j)
             ENDIF
             IF(av2(i,j) > 0.d0) THEN
                C = C + cpz*Flux_T*aa2(i,j)
             ELSE
                RHS = RHS - cpz*Flux_T*aa2(i,j+1)*arH2(i,j+1)
             ENDIF
             IF(av2(i,j-1) > 0.d0) THEN
                RHS = RHS + Flux_B*cpz*aa2(i,j-1)*arH2(i,j-1)
             ELSE
                C = C - Flux_B*cpz*aa2(i,j)
             ENDIF

             ! Update hydrogen density
             IF(C>allim) THEN
                arH2(i,j) = RHS/C
             ELSE
                arH2(i,j) = brH2(i,j)
             ENDIF
          ENDDO
       ENDDO
       arH2(1,:) = arH2(2,:); arH2(n9,:) = arH2(n,:)
       arH2(:,1) = arH2(:,2); arH2(:,m9) = arH2(:,m)

       eps_reached = 0.D0
       ! Backward pass
       DO j = m,2,-1
          DO i = n,2,-1
             Yold = arH2(i,j)

             gamH2 = t*gammaH2(i,j)*YH2_Gamma_a

             ! Time derivative 
             C = aa2(i,j)
             RHS = brH2(i,j)*ba2(i,j)

             ! Source terms
             RHS = RHS +  gamH2

             ! Convective terms
             Flux_R = r12(i)/rad(i)*au2(i,j)
             Flux_L = r12(i-1)/rad(i)*au2(i-1,j)
             Flux_T = av2(i,j)
             Flux_B = av2(i,j-1)

             IF(au2(i,j) > 0.d0) THEN
                C = C + cpr*Flux_R*aa2(i,j)
             ELSE
                RHS = RHS - cpr*Flux_R*aa2(i+1,j)*arH2(i+1,j)
             ENDIF
             IF(au2(i-1,j) > 0.d0) THEN
                RHS = RHS + cpr*Flux_L*aa2(i-1,j)*arH2(i-1,j)
             ELSE
                C = C - cpr*Flux_L*aa2(i,j)
             ENDIF
             IF(av2(i,j) > 0.d0) THEN
                C = C + cpz*Flux_T*aa2(i,j)
             ELSE
                RHS = RHS - cpz*Flux_T*aa2(i,j+1)*arH2(i,j+1)
             ENDIF
             IF(av2(i,j-1) > 0.d0) THEN
                RHS = RHS + Flux_B*cpz*aa2(i,j-1)*arH2(i,j-1)
             ELSE
                C = C - Flux_B*cpz*aa2(i,j)
             ENDIF

             ! Update hydrogen density
             IF(C>allim) THEN
                arH2(i,j) = RHS/C
             ELSE
                arH2(i,j) = brH2(i,j)
             ENDIF

             ! Estimate the maximum correaction
             eps_reached = MAX(eps_reached, ABS(Yold-arH2(i,j)))
          ENDDO
       ENDDO

       arH2(1,:) = arH2(2,:); arH2(n9,:) = arH2(n,:)
       arH2(:,1) = arH2(:,2); arH2(:,m9) = arH2(:,m)

       aYH2 = arH2/ara; aYH2 = MIN(aYH2,1.d0); aYH2 = MAX(aYH2,0.d0)
       aYH2(1,:) = aYH2(2,:); aYH2(n9,:) = aYH2(n,:)
       aYH2(:,1) = aYH2(:,2); aYH2(:,m9) = aYH2(:,m)

       IF(okprint) THEN
          PRINT *,'            H2: iter=',iter,'(',niter,') => ',eps_reached
       ENDIF
       IF(eps_reached < eps_H2) RETURN
    ENDDO
  END SUBROUTINE hydrogen_2D


END MODULE HYDROGEN
