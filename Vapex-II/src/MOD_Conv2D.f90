!------------------------------------------------------------------!
!  Coefficients of finite-volume approximations                    !
!  Convective terms                                                !
!------------------------------------------------------------------!
!  $Id: MOD_Conv2D.f90 5 2014-01-13 11:44:50Z Sergey $
!------------------------------------------------------------------!
MODULE CONVECTION_2D
  USE GLOBALS_2D
  USE COEFFICIENTS
  REAL(KIND=8):: ConvJac(NVARS,NVARS), PressJac(NVARS,NRHS-1)
  INTEGER,PARAMETER,PRIVATE:: isAlphaDiagonalImpl = 1
CONTAINS
  SUBROUTINE ADD_CONVECTIVE_TERMS(tstep,i,j)
    USE GLOBALS
    USE LOCAL_VARS
    USE WATER_PROPS
    USE CORRELATIONS
    USE EVAPORATION
    USE HYDROGEN
    USE DISPERSED_PHASE
    USE VARIABLES_2D
    USE GRID_2D

    IMPLICIT NONE

    REAL(KIND=8), INTENT(in):: tstep
    INTEGER, INTENT(in):: i,j
    !
    ! Local variables
    !
    REAL(KIND=8):: A1_L,A1_R,A1_B,A1_T
    REAL(KIND=8):: A2_L,A2_R,A2_B,A2_T
    REAL(KIND=8):: R1_L,R1_R,R1_B,R1_T
    REAL(KIND=8):: R2_L,R2_R,R2_B,R2_T
    REAL(KIND=8):: E1_L,E1_R,E1_B,E1_T
    REAL(KIND=8):: E2_L,E2_R,E2_B,E2_T
    REAL(KIND=8):: Ro1_L,Ro1_R,Ro1_B,Ro1_T
    REAL(KIND=8):: Ro2_L,Ro2_R,Ro2_B,Ro2_T

    REAL(KIND=8):: DYa_L, DYa_R, DYa_B, DYa_T
    REAL(KIND=8):: DE1_L, DE1_R, DE1_B, DE1_T
    REAL(KIND=8):: DE2_L, DE2_R, DE2_B, DE2_T

    REAL(KIND=8):: A1L_diag, A1R_diag, A1B_diag, A1T_diag
    REAL(KIND=8):: A1L_off, A1R_off, A1B_off, A1T_off
    REAL(KIND=8):: A2L_diag, A2R_diag, A2B_diag, A2T_diag
    REAL(KIND=8):: A2L_off, A2R_off, A2B_off, A2T_off

    REAL(KIND=8):: C1L_diag, C1R_diag, C1B_diag, C1T_diag
    REAL(KIND=8):: C1L_off, C1R_off, C1B_off, C1T_off
    REAL(KIND=8):: C2L_diag, C2R_diag, C2B_diag, C2T_diag
    REAL(KIND=8):: C2L_off, C2R_off, C2B_off, C2T_off

    REAL(KIND=8):: E1L_diag, E1R_diag, E1B_diag, E1T_diag
    REAL(KIND=8):: E1L_off, E1R_off, E1B_off, E1T_off
    REAL(KIND=8):: E2L_diag, E2R_diag, E2B_diag, E2T_diag
    REAL(KIND=8):: E2L_off, E2R_off, E2B_off, E2T_off
    !
    REAL(KIND=8):: cpr, cpz
    REAL(KIND=8):: cvt,cvb,cul,cur,c_u,c_v
    REAL(KIND=8):: RoA1_lim, RoA2_lim, A1_lim, A2_lim

    REAL(KIND=8):: Geom_L, Geom_R

    REAL(KIND=8):: Pij, Ro1ij, Ro2ij
    !
    ! Initialise convective Jacobian
    !
    ConvJac = 0.D0
    PressJac = 0.D0
    cpr = tstep/hi; cpz = tstep/hj
    Geom_R = dxu(i)/dxt(i)
    Geom_L = dxu(i-1)/dxt(i)
    !===================================================Convective terms
    ! Interpolation to cell faces
    Ro1_R = 0.5D0*(ar1(i+1,j)+ar1(i,j))
    Ro1_L = 0.5D0*(ar1(i,j)+ar1(i-1,j))
    Ro1_T = 0.5D0*(ar1(i,j+1)+ar1(i,j))
    Ro1_B = 0.5D0*(ar1(i,j)+ar1(i,j-1))
    A1_R = 0.5D0*(aa1(i+1,j)+aa1(i,j))
    A1_L = 0.5D0*(aa1(i,j)+aa1(i-1,j))
    A1_T = 0.5D0*(aa1(i,j+1)+aa1(i,j))
    A1_B = 0.5D0*(aa1(i,j)+aa1(i,j-1))
    E1_R = 0.5D0*(ae1(i+1,j)+ae1(i,j))
    E1_L = 0.5D0*(ae1(i,j)+ae1(i-1,j))
    E1_T = 0.5D0*(ae1(i,j+1)+ae1(i,j))
    E1_B = 0.5D0*(ae1(i,j)+ae1(i,j-1))

    Ro2_R = 0.5D0*(ar2(i+1,j)+ar2(i,j))
    Ro2_L = 0.5D0*(ar2(i,j)+ar2(i-1,j))
    Ro2_T = 0.5D0*(ar2(i,j+1)+ar2(i,j))
    Ro2_B = 0.5D0*(ar2(i,j)+ar2(i,j-1))
    A2_R = 0.5D0*(aa2(i+1,j)+aa2(i,j))
    A2_L = 0.5D0*(aa2(i,j)+aa2(i-1,j))
    A2_T = 0.5D0*(aa2(i,j+1)+aa2(i,j))
    A2_B = 0.5D0*(aa2(i,j)+aa2(i,j-1))
    E2_R = 0.5D0*(ae2(i+1,j)+ae2(i,j))
    E2_L = 0.5D0*(ae2(i,j)+ae2(i-1,j))
    E2_T = 0.5D0*(ae2(i,j+1)+ae2(i,j))
    E2_B = 0.5D0*(ae2(i,j)+ae2(i,j-1))

    DYa_L = aYa(i,j)-aYa(i-1,j)
    DYa_R = aYa(i+1,j)-aYa(i,j)
    DYa_B = aYa(i,j)-aYa(i,j-1)
    DYa_T = aYa(i,j+1)-aYa(i,j)

    DE1_L = ae1(i,j)-ae1(i-1,j)
    DE1_R = ae1(i+1,j)-ae1(i,j)
    DE1_B = ae1(i,j)-ae1(i,j-1)
    DE1_T = ae1(i,j+1)-ae1(i,j)

    DE2_L = ae2(i,j)-ae2(i-1,j)
    DE2_R = ae2(i+1,j)-ae2(i,j)
    DE2_B = ae2(i,j)-ae2(i,j-1)
    DE2_T = ae2(i,j+1)-ae2(i,j)

    Pij = bp(i,j)
    Ro1ij = ar1(i,j)
    Ro2ij = ar2(i,j)

    ! Backward differences: A1R = A1R_diag+A1R_off; A1R_diag = c1r_diag*aa1(i,j) etc.
    IF(au1(i,j)>=0) THEN
       A1R_diag = aa1(i,j); c1r_diag = 1.D0; A1R_off = 0.D0; c1r_off = 0.D0
       E1R_diag = ae1(i,j); E1R_off = 0.D0
    ELSE
       A1R_diag = 0.D0; c1r_diag = 0.D0; A1R_off = aa1(i+1,j); c1r_off = 1.D0
       E1R_diag = 0.D0; E1R_off = ae1(i+1,j)
    ENDIF
    IF(au1(i-1,j)>=0) THEN
       A1L_diag = 0.D0; c1L_diag = 0.D0; A1L_off = aa1(i-1,j); c1l_off = 1.D0
       E1L_diag = 0.D0; E1L_off = ae1(i-1,j)
    ELSE
       A1L_diag = aa1(i,j); c1L_diag = 1.D0; A1L_off = 0.D0; c1l_off = 0.D0
       E1L_diag = ae1(i,j); E1L_off = 0.D0
    ENDIF
    IF(av1(i,j)>=0) THEN
       A1T_diag = aa1(i,j); c1T_diag = 1.D0; A1T_off = 0.D0; c1t_off = 0.D0
       E1T_diag = ae1(i,j); E1T_off = 0.D0
    ELSE
       A1T_diag = 0.D0; c1T_diag = 0.D0; A1T_off = aa1(i,j+1); C1T_off = 1.D0
       E1T_diag = 0.D0; E1T_off = ae1(i,j+1)
    ENDIF
    IF(av1(i,j-1)>=0) THEN
       A1B_diag = 0.D0; c1B_diag = 0.D0; A1B_off = aa1(i,j-1); C1B_off = 1.D0
       E1B_diag = 0.D0; E1B_off = ae1(i,j-1)
    ELSE
       A1B_diag = aa1(i,j); c1B_diag = 1.D0; A1B_off = 0.D0; C1B_off = 0.D0
       E1B_diag = ae1(i,j); E1B_off = 0.D0
    ENDIF
    IF(au2(i,j)>=0) THEN
       A2R_diag = aa2(i,j); c2r_diag = 1.D0; A2R_off = 0.D0; C2R_off = 0.D0;
       E2R_diag = ae2(i,j); E2R_off = 0.D0
    ELSE
       A2R_diag = 0.D0; c2r_diag = 0.D0; A2R_off = aa2(i+1,j); C2R_off = 1.D0
       E2R_diag = 0.D0; E2R_off = ae2(i+1,j)
    ENDIF
    IF(au2(i-1,j)>=0) THEN
       A2L_diag = 0.D0; c2L_diag = 0.D0; A2L_off = aa2(i-1,j); c2L_off = 1.D0
       E2L_diag = 0.D0; E2L_off = ae2(i-1,j)
    ELSE
       A2L_diag = aa2(i,j); c2L_diag = 1.D0; A2L_off = 0.D0; c2L_off = 0.D0
       E2L_diag = ae2(i,j); E2L_off = 0.D0
    ENDIF
    IF(av2(i,j)>=0) THEN
       A2T_diag = aa2(i,j); c2T_diag = 1.D0; A2T_off = 0.D0; c2T_off = 0.D0
       E2T_diag = ae2(i,j); E2T_off = 0.D0
    ELSE
       A2T_diag = 0.D0; c2T_diag = 0.D0; A2T_off = aa2(i,j+1); c2T_off = 1.D0
       E2T_diag = 0.D0; E2T_off = ae2(i,j+1)
    ENDIF
    IF(av2(i,j-1)>=0) THEN
       A2B_diag = 0.D0; c2B_diag = 0.D0; A2B_off = aa2(i,j-1); c2B_off = 1.D0
       E2B_diag = 0.D0; E2B_off = ae2(i,j-1)
    ELSE
       A2B_diag = aa2(i,j); c2B_diag = 1.D0; A2B_off = 0.D0; c2B_off = 0.D0
       E2B_diag = ae2(i,j); E2B_off = 0.D0
    ENDIF
    !
    ! Limited denominators
    !
    A1_lim = MAX(ba1(i,j),allim)
    A2_lim = MAX(ba2(i,j),allim)
    RoA1_lim = br1(i,j)*A1_lim
    RoA2_lim = br2(i,j)*A2_lim
    !====================================================================
    !
    ! Liquid continuity
    !
    RHS(1) = RHS(1) - &
         cpr*(au1(i,j)*Ro1_R*Geom_R*(A1R_diag+A1R_off)- &
         au1(i-1,j)*Ro1_L*Geom_L*(A1L_diag+A1L_off)) - &
         cpz*(av1(i,j)*Ro1_T*(A1T_diag+A1T_off) - &
         av1(i,j-1)*Ro1_B*(A1B_diag+A1B_off))
    IF(isAlphaDiagonalImpl /= 0) THEN
       ConvJac(1,1) = ConvJac(1,1) - &
            cpr*(au1(i,j)*Ro1_R*Geom_R*c1R_diag - &
            au1(i-1,j)*Ro1_L*Geom_L*C1L_diag)- &
            cpz*(av1(i,j)*Ro1_T*c1T_diag - &
            av1(i,j-1)*Ro1_B*c1B_diag)
    ENDIF
    ! Terms with pressure corrections
    PressJac(1,1) =  cpr*fu1l(i-1,j)*Ro1_L*Geom_L*A1_L
    PressJac(1,2) =  cpr*fu1r(i,j)*Ro1_R*Geom_R*A1_R
    PressJac(1,3) =  cpz*fv1l(i,j-1)*Ro1_B*A1_B
    PressJac(1,4) =  cpz*fv1r(i,j)*Ro1_T*A1_T

    ConvJac(1,4) = ConvJac(1,4) + &
         cpr*fu1r(i-1,j)*Ro1_L*Geom_L*A1_L + &
         cpr*fu1l(i,j)*Ro1_R*Geom_R*A1_R + &
         cpz*fv1r(i,j-1)*Ro1_B*A1_B + cpz*fv1l(i,j)*Ro1_T*A1_T
    !====================================================================
    !
    ! Vapour continuity
    !
    RHS(2) = RHS(2) - &
         cpr*(au2(i,j)*Ro2_R*Geom_R*(A2R_diag+A2R_off)- &
         au2(i-1,j)*Ro2_L*Geom_L*(A2L_diag+A2L_off)) - &
         cpz*(av2(i,j)*Ro2_T*(A2T_diag+A2T_off) - &
         av2(i,j-1)*Ro2_B*(A2B_diag+A2B_off))
    IF(isAlphaDiagonalImpl /= 0) THEN
       ConvJac(2,1) = ConvJac(2,1) + &
            cpr*(au2(i,j)*Ro2_R*Geom_R*c2R_diag - &
            au2(i-1,j)*Ro2_L*Geom_L*c2L_diag)+ &
            cpz*(av2(i,j)*Ro2_T*c2T_diag - &
            av2(i,j-1)*Ro2_B*c2B_diag)
    ENDIF
    ! Terms with pressure corrections
    PressJac(2,1) =  cpr*fu2l(i-1,j)*Ro2_L*Geom_L*A2_L
    PressJac(2,2) =  cpr*fu2r(i,j)*Ro2_R*Geom_R*A2_R
    PressJac(2,3) =  cpz*fv2l(i,j-1)*Ro2_B*A2_B
    PressJac(2,4) =  cpz*fv2r(i,j)*Ro2_T*A2_T

    ConvJac(2,4) = ConvJac(2,4) + &
         cpr*fu2r(i-1,j)*Ro2_L*Geom_L*A2_L + &
         cpr*fu2l(i,j)*Ro2_R*Geom_R*A2_R + &
         cpz*fv2r(i,j-1)*Ro2_B*A2_B + cpz*fv2l(i,j)*Ro2_T*A2_T

    !====================================================================
    !
    ! Liquid energy
    !
    IF(au1(i,j) >=0.D0) THEN
       cur = 0.D0
    ELSE
       cur = cpr*Ro1_R*Geom_R*A1R_off/RoA1_lim
    ENDIF
    IF(au1(i-1,j) < 0.D0 ) THEN
       cul = 0.D0
    ELSE
       cul = cpr*Ro1_L*Geom_L*A1L_off/RoA1_lim
    ENDIF

    c_u = cul*au1(i-1,j)-cur*au1(i,j)
    IF(c_u < 0.D0) THEN
       PRINT *,'c_u1: ',i,j,c_u,cur,cul
       STOP
    ENDIF

    IF(av1(i,j) >=0.D0 ) THEN
       cvt = 0.D0
    ELSE
       cvt = cpz*Ro1_T*A1T_off/RoA1_lim
    ENDIF
    IF(av1(i,j-1) < 0.D0) THEN
       cvb = 0.D0
    ELSE
       cvb = cpz*Ro1_B*A1B_off/RoA1_lim
    ENDIF
    c_v = cvb*av1(i,j-1)-cvt*av1(i,j)
    IF(c_v < 0.D0) THEN
       PRINT *,'c_v1: ',i,j,c_v,cvt,cvb
       STOP
    ENDIF

    RHS(3) = RHS(3) - &
         (cur*au1(i,j)*DE1_R+ cul*au1(i-1,j)*DE1_L + &
         cvt*av1(i,j)*DE1_T + cvb*av1(i,j-1)*DE1_B)
    ConvJac(3,2) = ConvJac(3,2) + dE1_dT1*(c_u+c_v)
    ConvJac(3,4) = ConvJac(3,4) + dE1_dP *(c_u+c_v)
    !
    ! Terms with pressure corrections
    !
    PressJac(3,1) = cul*fu1l(i-1,j)*DE1_L
    PressJac(3,2) = -cur*fu1r(i,j)*DE1_R
    PressJac(3,3) = cvb*fv1l(i,j-1)*DE1_B
    PressJac(3,4) = -cvt*fv1r(i,j)*DE1_T
    ConvJac(3,4) = ConvJac(3,4) + &
         cul*fu1r(i-1,j)*DE1_L - cur*fu1l(i,j)*DE1_R + &
         cvb*fv1r(i,j-1)*DE1_B - cvt*fv1l(i,j)*DE1_T
    !
    ! Convective terms in (P/Ro)*a*DRo/Dt:
    !
    cur = Pij/Ro1ij*cpr*Geom_R*A1R_off/RoA1_lim
    cul = Pij/Ro1ij*cpr*Geom_L*A1L_off/RoA1_lim
    cvt = Pij/Ro1ij*cpz*A1T_off/RoA1_lim
    cvb = Pij/Ro1ij*cpz*A1B_off/RoA1_lim

    RHS(3) = RHS(3) + &
         (cur*au1(i,j)*(Ro1_R-Ro1ij)+cul*au1(i-1,j)*(Ro1ij-Ro1_L)+&
         cvt*av1(i,j)*(Ro1_T-Ro1ij)+cvb*av1(i,j-1)*(Ro1ij-Ro1_B))
    !
    ! Terms with pressure correction
    !
    PressJac(3,1) = PressJac(3,1) + cul*(Ro1ij-Ro1_L)*fu1l(i-1,j)
    PressJac(3,2) = PressJac(3,2) + cur*(Ro1_R-Ro1ij)*fu1r(i,j)
    PressJac(3,3) = PressJac(3,3) + cvb*(Ro1ij-Ro1_B)*fv1l(i,j-1)
    PressJac(3,4) = PressJac(3,4) + cvt*(Ro1_T-Ro1ij)*fv1r(i,j)

    ConvJac(3,4) = ConvJac(3,4) + &
         cul*(Ro1ij-Ro1_L)*fu1r(i-1,j) + &
         cur*(Ro1_R-Ro1ij)*fu1l(i,j) + &
         cvb*(Ro1ij-Ro1_B)*fv1r(i,j-1) + &
         cvt*(Ro1_T-Ro1ij)*fv1l(i,j)
    !====================================================================
    !
    ! Combined gas energy
    !
    IF(au2(i,j) >=0.D0) THEN
       cur = 0.D0
    ELSE
       cur = cpr*Ro2_R*Geom_R*A2R_off/RoA2_lim
    ENDIF
    IF(au2(i-1,j) <= 0.D0) THEN
       cul = 0.D0
    ELSE
       cul = cpr*Ro2_L*Geom_L*A2L_off/RoA2_lim
    ENDIF

    c_u = cul*au2(i-1,j)-cur*au2(i,j)
    IF(c_u < 0.D0) THEN
       PRINT *,'c_u2: ',i,j,c_u,cur,cul
       STOP
    ENDIF

    IF(av2(i,j) >=0.D0) THEN
       cvt = 0.D0
    ELSE
       cvt = cpz*Ro2_T*A2T_off/RoA2_lim
    ENDIF
    IF(av2(i,j-1) <= 0.D0) THEN
       cvb = 0.D0
    ELSE
       cvb = cpz*Ro2_B*A2B_off/RoA2_lim
    ENDIF
    c_v = cvb*av2(i,j-1)-cvt*av2(i,j)
    IF(c_v < 0.D0) THEN
       PRINT *,'c_v2: ',i,j,c_v
       STOP
    ENDIF

    RHS(4) = RHS(4) - &
         (cur*au2(i,j)*DE2_R+cul*au2(i-1,j)*DE2_L + &
         cvt*av2(i,j)*DE2_T+cvb*av2(i,j-1)*DE2_B)

    ConvJac(4,3) = ConvJac(4,3) + dE2_dT2*(c_u + c_v)
    ConvJac(4,4) = ConvJac(4,4) + dE2_dP *(c_u + c_v)
    ConvJac(4,5) = ConvJac(4,5) + dE2_dPa*(c_u + c_v)
    !
    ! Terms with pressure corrections
    !
    PressJac(4,1) = cul*fu2l(i-1,j)*DE2_L
    PressJac(4,2) = -cur*fu2r(i,j)*DE2_R
    PressJac(4,3) = cvb*fv2l(i,j-1)*DE2_B
    PressJac(4,4) = -cvt*fv2r(i,j)*DE2_T
    ConvJac(4,4) = ConvJac(4,4) + &
         cul*fu2r(i-1,j)*DE2_L - cur*fu2l(i,j)*DE2_R + &
         cvb*fv2r(i,j-1)*DE2_B - cvt*fv2l(i,j)*DE2_T
    !
    ! Convective terms in (P/Ro)*a*DRo/Dt:
    !
    cur = Pij/Ro2ij*cpr*Geom_R*A2R_off/RoA2_lim
    cul = Pij/Ro2ij*cpr*Geom_L*A2L_off/RoA2_lim
    cvt = Pij/Ro2ij*cpz*A2T_off/RoA2_lim
    cvb = Pij/Ro2ij*cpz*A2B_off/RoA2_lim

    RHS(4) = RHS(4) + &
         (cur*au2(i,j)*(Ro2_R-Ro2ij)+cul*au2(i-1,j)*(Ro2ij-Ro2_L)+&
         cvt*av2(i,j)*(Ro2_T-Ro2ij)+cvb*av2(i,j-1)*(Ro2ij-Ro2_B))
    !
    ! Terms with pressure correction
    !
    PressJac(4,1) = PressJac(4,1) + cul*(Ro2ij-Ro2_L)*fu2l(i-1,j)
    PressJac(4,2) = PressJac(4,2) + cur*(Ro2_R-Ro2ij)*fu2r(i,j)
    PressJac(4,3) = PressJac(4,3) + cvb*(Ro2ij-Ro2_B)*fv2l(i,j-1)
    PressJac(4,4) = PressJac(4,4) + cvt*(Ro2_T-Ro2ij)*fv2r(i,j)

    ConvJac(4,4) = ConvJac(4,4) + &
         cul*(Ro2ij-Ro2_L)*fu2r(i-1,j) + &
         cur*(Ro2_R-Ro2ij)*fu2l(i,j) + &
         cvb*(Ro2ij-Ro2_B)*fv2r(i,j-1) + &
         cvt*(Ro2_T-Ro2ij)*fv2l(i,j)
    !====================================================================
    !
    ! Non-condensible
    !
    IF(au2(i,j) >=0.D0 ) THEN
       cur = 0.D0
    ELSE
       cur = cpr*Ro2_R*Geom_R*A2R_off/RoA2_lim
    ENDIF
    IF(au2(i-1,j) <= 0.D0 ) THEN
       cul = 0.D0
    ELSE
       cul = cpr*Ro2_L*Geom_L*A2L_off/RoA2_lim
    ENDIF
    c_u = cul*au2(i-1,j)-cur*au2(i,j)
    IF(c_u < 0.D0) THEN
       PRINT *,'c_u2: ',i,j,c_u,cur,cul
       STOP
    ENDIF

    IF(av2(i,j) >=0.D0 ) THEN
       cvt = 0.D0
    ELSE
       cvt = cpz*Ro2_T*A2T_off/RoA2_lim
    ENDIF
    IF(av2(i,j-1) <= 0.D0 ) THEN
       cvb = 0.D0
    ELSE
       cvb = cpz*Ro2_B*A2B_off/RoA2_lim
    ENDIF
    c_v = cvb*av2(i,j-1)-cvt*av2(i,j)
    IF(c_v < 0.D0) THEN
       PRINT *,'c_v2_a: ',i,j,c_v
    ENDIF

    RHS(5) = RHS(5) - &
         (cur*au2(i,j)*DYa_R + cul*au2(i-1,j)*DYa_L + &
         cvt*av2(i,j)*DYa_T + cvb*av2(i,j-1)*DYa_B)

    ConvJac(5,3) = ConvJac(5,3) + dYa_dT2*(c_u + c_v)
    ConvJac(5,4) = ConvJac(5,4) + dYa_dP *(c_u + c_v)
    ConvJac(5,5) = ConvJac(5,5) + dYa_dPa*(c_u + c_v)
    !
    ! Terms with pressure corrections
    !
    PressJac(5,1) = cul*fu2l(i-1,j)*DYa_L
    PressJac(5,2) = -cur*fu2r(i,j)*DYa_R
    PressJac(5,3) = cvb*fv2l(i,j-1)*DYa_B
    PressJac(5,4) = -cvt*fv2r(i,j)*DYa_T
    ConvJac(5,4) = ConvJac(5,4) + &
         cul*fu2r(i-1,j)*DYa_L - cur*fu2l(i,j)*DYa_R + &
         cvb*fv2r(i,j-1)*DYa_B - cvt*fv2l(i,j)*DYa_T

  END SUBROUTINE ADD_CONVECTIVE_TERMS
END MODULE CONVECTION_2D
