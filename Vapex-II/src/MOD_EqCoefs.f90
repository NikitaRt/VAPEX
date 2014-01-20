!------------------------------------------------------------------!
!  Coefficients of finite-difference approximations                !
!  Time derivatives, source terms                                  !
!------------------------------------------------------------------!
!  $Id: MOD_EqCoefs.f90 5 2014-01-13 11:44:50Z Sergey $
!------------------------------------------------------------------!
MODULE COEFFICIENTS
  USE GLOBALS_2D
  !
  ! Flags
  !
  INTEGER,PRIVATE:: isPhaseChange = 0
  INTEGER,PRIVATE:: isGasAddition = 0

  REAL(KIND=8):: TimeJac(NVARS,NVARS), SourceJac(NVARS,NVARS)
  REAL(KIND=8):: GMatrix(NVARS,NVARS)

  REAL(KIND=8):: RHS(5)
  REAL(KIND=8), PRIVATE:: Ro1Al_lim, Ro2Al_lim, Alpha_lim
  REAL(KIND=8), PRIVATE:: h_liq, h_gas, h_add

  REAL(KIND=8):: isUseDeltaAlpha3

CONTAINS
  SUBROUTINE EQUATION_COEFS(tstep)
    USE GLOBALS
    USE LOCAL_VARS
    USE WATER_PROPS
    USE CORRELATIONS
    USE EVAPORATION
    USE HYDROGEN
    USE DISPERSED_PHASE

    IMPLICIT NONE

    REAL(KIND=8), INTENT(in):: tstep
    !
    ! Local variables
    !
    REAL(KIND=8):: G_lim
    !
    ! Initialise Jacobian and Right-hand side matrices
    !
    TimeJac = 0.D0
    SourceJac = 0.D0
    RHS = 0.D0
    !
    ! Set the limited (Ro*Alpha)_n
    !
    Ro1Al_lim = 1.D0/(Ro1n*MAX(Alpha_ln,allim))
    Ro2Al_lim = 1.D0/(Ro2n*MAX(Alphan,allim))
    !
    ! Set G_lim for calculation of true Gamma_Evap in continuity equations
    !
    G_lim = MIN(1.D0,Alpha/allim,Alpha_l/allim)
    !
    ! Liquid continuity: Eq. No. 1
    ! 
    RHS(1) = -(Ro1*Alpha_l-Ro1n*Alpha_ln)
    IF(isThirdPhase /= 0) RHS(1) = RHS(1) + &
         (tstep/tstep_disp)*Ro1*(Alpha_3-Alpha_3n)*isUseDeltaAlpha3
    IF(isPhaseChange /= 0) RHS(1) = RHS(1) - tstep*Gamma_Evap*G_lim
    TimeJac(1,1) = -Ro1
    TimeJac(1,2) = Alpha_l*dRo1_dT1
    TimeJac(1,4) = Alpha_l*dRo1_dP
    IF(isPhaseChange /= 0) THEN
       SourceJac(1,1) = -dGamma_dAlpha*tstep
       SourceJac(1,2) = -dGamma_dT1*tstep
       SourceJac(1,3) = -dGamma_dT2*tstep
       SourceJac(1,4) = -dGamma_dP*tstep
       SourceJac(1,5) = -dGamma_dPa*tstep
    ENDIF
    !
    ! Combined gas continuity: Eq. No. 2
    !
    Alpha_lim = MAX(Alpha,allim)
    RHS(2) = -(Ro2*Alpha-Ro2n*Alphan)
    IF(isPhaseChange /= 0) RHS(2) = RHS(2) + tstep*Gamma_Evap*G_lim
    IF(isGasAddition /= 0) RHS(2) = RHS(2) + tstep*Gamma_a
    TimeJac(2,1) = Ro_2
    TimeJac(2,3) = Alpha*dRo2_dT2
    TimeJac(2,4) = Alpha_lim*dRo2_dP
    TimeJac(2,5) = Alpha_lim*dRo2_dPa
    IF(isPhaseChange /= 0) THEN
       SourceJac(2,1) = dGamma_dAlpha*tstep
       SourceJac(2,2) = dGamma_dT1*tstep
       SourceJac(2,3) = dGamma_dT2*tstep
       SourceJac(2,4) = dGamma_dP*tstep
       SourceJac(2,5) = dGamma_dPa*tstep
    ENDIF
    !
    ! Energy of liquid water: Eq. No. 3
    !
    RHS(3) = -(E1 - E1n) & ! Time derivative
         + P*(Ro1-Ro1n)/(Ro1*Ro1n) & ! Compressibility
         + tstep*R1s*(Tsv - T1)*Ro1Al_lim & ! Liquid-interface heat transfer
         + tstep*RF1s*(T_sat - T1)*Ro1Al_lim & ! Flash (works when T1>T_sat only)
         + tstep*RG*(Pa/P)*(T2 - T1)*Ro1Al_lim  ! Sensible heat flux
    IF(isThirdPhase /= 0) THEN
       RHS(3) = RHS(3) + tstep*R31*(T3-T1)*Ro1Al_lim ! Interaction with third phase
    ENDIF

    TimeJac(3,2) = dE1_dT1
    TimeJac(3,4) = dE1_dP

    SourceJac(3,2) = P/Ro1**2*dRo1_dT1 + &
         tstep*Ro1Al_lim*(-R1s - RG*Pa/P - RF1s)
    IF(isThirdPhase /= 0) THEN
       SourceJac(3,2) = SourceJac(3,2) - tstep*Ro1Al_lim*R31
    ENDIF
    SourceJac(3,3) = tstep*Ro1Al_lim*RG*Pa/P
    SourceJac(3,4) = (Ro1-Ro1n)/(Ro1*Ro1n) + &
         P/Ro1**2*dRo1_dP + &
         tstep*Ro1Al_lim*(R1s*dTsv_dP - RG*Pa/P**2*(T2-T1) + RF1s*dTsat_dP)
    SourceJac(3,5) = tstep*Ro1Al_lim*(R1s*dTsv_dPa + RG/P*(T2-T1))
    ! Condensation only !
    IF(isPhaseChange /= 0 .AND. Gamma_Evap < 0.D0) THEN
       h_liq = E1 + P/Ro1
       RHS(3) = RHS(3) + tstep*Gamma_Evap*(h_liq-hl_prime)*Ro1Al_lim
       SourceJac(3,2) = SourceJac(3,2)&
            + tstep*Gamma_Evap*Ro1Al_lim*(dE1_dT1 - P/Ro1**2*dRo1_dT1) &
            + tstep*Ro1Al_lim*(h_liq - hl_prime)*dGamma_dT1
       SourceJac(3,3) = SourceJac(3,3)&
            + tstep*Ro1Al_lim*(h_liq - hl_prime)*dGamma_dT2
       SourceJac(3,4) = SourceJac(3,4)&
            + tstep*Gamma_Evap*Ro1Al_lim*(dE1_dP + 1/Ro1 - P/Ro1**2*dRo1_dP - dhl_dP)&
            + tstep*Ro1Al_lim*(h_liq - hl_prime)*dGamma_dP
       SourceJac(3,5) = SourceJac(3,5)&
            - tstep*Gamma_Evap*Ro1Al_lim*dhl_dPa&
            + tstep*Ro1Al_lim*(h_liq - hl_prime)*dGamma_dPa
    ENDIF
    !
    ! Energy of combined gas (steam+non-condensable): Eq. No. 4
    !
    RHS(4) = -(E2 - E2n) & ! Time derivative
         + P*(Ro2-Ro2n)/(Ro2*Ro2n) & ! Compressibility
         + tstep*(Pv/P)*R2s*(Tsv - T2)*Ro2Al_lim & ! Gas-interface heat transfer
         + tstep*RG*(Pa/P)*(T1 - T2)*Ro2Al_lim  ! Sensible heat flux
    IF(isThirdPhase /= 0) THEN
       RHS(4) = RHS(4) + tstep*R32*(T3-T2)*Ro2Al_lim ! Interaction with third phase
    ENDIF

    TimeJac(4,3) = dE2_dT2
    TimeJac(4,4) = dE2_dP
    TimeJac(4,5) = dE2_dPa

    SourceJac(4,2) = tstep*Ro2Al_lim*RG*(Pa/P)
    SourceJac(4,3) = P/Ro2**2*dRo2_dT2 &
         + tstep*Ro2Al_lim*(-R2s*Pv/P - RG*Pa/P)
    IF(isThirdPhase /= 0) THEN
       SourceJac(4,3) = SourceJac(4,3) - tstep*Ro2Al_lim*R32
    ENDIF
    SourceJac(4,4) = (Ro2-Ro2n)/(Ro2*Ro2n) + P/Ro2**2*dRo2_dP + &
         tstep*Ro2Al_lim*R2s*((Pv/P)*dTsv_dP+(Pa/P**2)*(Tsv-T2)) + &
         tstep*Ro2Al_lim*RG*(-(T1-T2)*(Pa/P**2))
    SourceJac(4,5) = P/Ro2**2*dRo2_dPa + &
         tstep*Ro2Al_lim*R2s*((Pv/P)*dTsv_dPa - (Tsv-T2)/P) + &
         tstep*Ro2Al_lim*RG*(T1-T2)/P
    IF(isPhaseChange /= 0) THEN
       h_gas = E2 + P/Ro2
       RHS(4) = RHS(4) - tstep*Gamma_Evap*(h_gas-hv_prime)*Ro2Al_lim
       SourceJac(4,2) = SourceJac(4,2) - tstep*(h_gas-hv_prime)*Ro2Al_lim*dGamma_dT1
       SourceJac(4,3) = SourceJac(4,3) - tstep*(h_gas-hv_prime)*Ro2Al_lim*dGamma_dT2 &
            - tstep*Gamma_Evap*Ro2Al_lim*(dE2_dT2-P/Ro2**2*dRo2_dT2 - dhv_dT2)
       SourceJac(4,4) = SourceJac(4,4) - tstep*(h_gas-hv_prime)*Ro2Al_lim*dGamma_dP &
            - tstep*Gamma_Evap*Ro2Al_lim*(dE2_dP-P/Ro2**2*dRo2_dP + 1/Ro2 - dhv_dP)
       SourceJac(4,5) = SourceJac(4,5) - tstep*(h_gas-hv_prime)*Ro2Al_lim*dGamma_dPa &
            - tstep*Gamma_Evap*Ro2Al_lim*(dE2_dPa-P/Ro2**2*dRo2_dPa - dhv_dPa)
    ENDIF
    IF(isGasAddition /= 0) THEN
       h_gas = E2 + P/Ro2
       h_add = h_Ar_H2(TH2_Gamma_a,YH2_Gamma_a,YH2)
       RHS(4) = RHS(4) - tstep*Gamma_a*(h_gas-h_add)*Ro2Al_lim
       SourceJac(4,3) = SourceJac(4,3)   &
            - tstep*Gamma_a*Ro2Al_lim*(dE2_dT2-P/Ro2**2*dRo2_dT2)
       SourceJac(4,4) = SourceJac(4,4)   &
            - tstep*Gamma_a*Ro2Al_lim*(dE2_dP-P/Ro2**2*dRo2_dP + 1/Ro2)
       SourceJac(4,5) = SourceJac(4,5)   &
            - tstep*Gamma_a*Ro2Al_lim*(dE2_dPa-P/Ro2**2*dRo2_dPa)
    ENDIF
    !
    ! Non-condensable gas continuity: Eq. No. 5
    !
    TimeJac(5,3) = dYa_dT2
    TimeJac(5,4) = dYa_dP
    TimeJac(5,5) = dYa_dPa
    RHS(5) = -(Ya-Yan)
    IF(isGasAddition /= 0) THEN
       RHS(5) = RHS(5) + tstep*Gamma_a*(1.D0-Ya)*Ro2Al_lim
       SourceJac(5,3) = SourceJac(5,3) - tstep*Gamma_a*Ro2Al_lim*dYa_dT2 
       SourceJac(5,4) = SourceJac(5,4) - tstep*Gamma_a*Ro2Al_lim*dYa_dP
       SourceJac(5,5) = SourceJac(5,5) - tstep*Gamma_a*Ro2Al_lim*dYa_dPa 
    ENDIF
    IF(isPhaseChange /= 0) THEN
       RHS(5) = RHS(5) - tstep*Gamma_Evap*Ya*Ro2Al_lim

       SourceJac(5,3) = SourceJac(5,3) - tstep*Gamma_Evap*Ro2Al_lim*dYa_dT2 
       SourceJac(5,4) = SourceJac(5,4) - tstep*Gamma_Evap*Ro2Al_lim*dYa_dP
       SourceJac(5,5) = SourceJac(5,5) - tstep*Gamma_Evap*Ro2Al_lim*dYa_dPa 

       SourceJac(5,2) = SourceJac(5,2) - Ya*dGamma_dT1*tstep*Ro2Al_lim
       SourceJac(5,3) = SourceJac(5,3) - Ya*dGamma_dT2*tstep*Ro2Al_lim
       SourceJac(5,4) = SourceJac(5,4) - Ya*dGamma_dP*tstep*Ro2Al_lim
       SourceJac(5,5) = SourceJac(5,5) - Ya*dGamma_dPa*tstep*Ro2Al_lim
    ENDIF
  END SUBROUTINE EQUATION_COEFS

  SUBROUTINE SetGassAdditionFlag(isGasRelease)
    IMPLICIT NONE
    INTEGER, INTENT(in):: isGasRelease
    IF(isGasRelease /= 0) THEN
       isGasAddition = 1
    ELSE
       isGasAddition = 0
    ENDIF
  END SUBROUTINE SetGassAdditionFlag

  SUBROUTINE SetPhaseChangeFlag(isEvaporation)
    IMPLICIT NONE
    INTEGER, INTENT(in):: isEvaporation
    IF(isEvaporation /= 0) THEN
       isPhaseChange = 1
    ELSE
       isPhaseChange = 0
    ENDIF
  END SUBROUTINE SetPhaseChangeFlag

END MODULE COEFFICIENTS

