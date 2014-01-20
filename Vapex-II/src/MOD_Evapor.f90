!------------------------------------------------------------------!
!  Evaporation rate and its partial derivatives                    !
!------------------------------------------------------------------!
!  $Id: MOD_Evapor.f90 5 2014-01-13 11:44:50Z Sergey $
!------------------------------------------------------------------!
MODULE EVAPORATION
  !
  ! Evaporation rate and its partial derivatives
  !
  REAL(kind=8):: Gamma_Evap
  REAL(kind=8):: dGamma_dAlpha
  REAL(kind=8):: dGamma_dT1, dGamma_dT2, dGamma_dP, dGamma_dPa
  REAL(kind=8):: hl_prime
  REAL(kind=8):: dhl_dT1, dhl_dP, dhl_dPa
  REAL(kind=8):: hv_prime
  REAL(kind=8):: dhv_dT2, dhv_dP, dhv_dPa
CONTAINS
  SUBROUTINE Evap_Rate
    USE WATER_PROPS
    USE LOCAL_VARS
    USE CORRELATIONS
    IMPLICIT NONE
    !
    ! Local variables
    !
    REAL(kind=8):: hl_plus,hl_minus
    REAL(kind=8):: hv_plus,hv_minus
    REAL(kind=8):: Qlv
    REAL(kind=8):: qi1, dqi1dt1, dqi1dp, dqi1dpa, dqi1da
    REAL(kind=8):: qi2, dqi2dt2, dqi2dp, dqi2dpa, dqi2da
    REAL(kind=8):: gam
    REAL(kind=8):: flash
    !
    !
    !
    Pv = MAX(P - Pa,0.D0)
    hl_plus = E_1 + P/Ro_1
    hl_minus = hl_sat

    hv_plus = hv_sat
    hv_minus = (E_2*Ro_2-E_a*Ro_a)/Ro_v + Pv/Ro_v

    qi1=R1s*(Tsv-T1)
    qi2=R2s*(Pv/P)*(Tsv-T2)

    dqi1dt1 = -R1s
    dqi1dp = R1s*dTsv_dP
    dqi1dpa = R1s*dTsv_dPa
    dqi1da = dR1s_dAlpha*(Tsv-T1)

    dqi2dt2 = -R2s*(Pv/P)
    dqi2dp = R2s*((Pv/P)*dTsv_dP+(Tsv-T2)*Pa/P**2)
    dqi2dpa = R2s*((Pv/P)*dTsv_dPa-(Tsv-T2)/P)
    dqi2da = dR2s_dAlpha*(Pv/P)*(Tsv-T2)

    IF(T1 > T_sat) THEN
       flash = RF1s*(T_sat - T1)
       dqi1dT1 = dqi1dT1 - RF1s
       dqi1dP = dqi1dP + RF1s*dTsat_dP
       dqi1da = dqi1da + dRF1s_dAlpha*(T_sat - T1)
    ELSE
       flash = 0.D0
    ENDIF

    gam = - (qi1+qi2+flash) ! Just to get sign!

    IF(gam >= 0.D0) THEN
       hl_prime = hl_plus
       dhl_dT1= dE1_dT1 - dRo1_dT1*(P/Ro_1**2)
       dhl_dP = dE1_dP + (Ro_1 - P*dRo1_dP)/(Ro_1**2)
       dhl_dPa = 0

       hv_prime = hv_plus
       dhv_dT2 = 0.d0
       dhv_dP  = dhv_sat_dP
       dhv_dPa = dhv_sat_dPa
    ELSE
       hl_prime = hl_minus
       hv_prime = hv_minus

       dhl_dT1 = 0.d0
       dhl_dP  = dhl_sat_dP
       dhl_dPa = dhl_sat_dPa

       dhv_dT2 = dEv_dT2 - Pv*dRov_dT2/Ro_v**2
       dhv_dP  = dEv_dP  + (dpsdp-Pv/Ro_v*dRov_dP)/Ro_v
       dhv_dPa = dEv_dPa + (dpsdpa - Pv/Ro_v*dRov_dPa)/Ro_v
    ENDIF

    Qlv = hv_prime - hl_prime
    Qlv = MAX(Qlv,1.D3)

    Gamma_Evap = gam/Qlv
    dGamma_dT1 = -(dqi1dt1-Gamma_Evap*dhl_dT1)/Qlv
    dGamma_dT2 = -(dqi2dt2+Gamma_Evap*dhv_dT2)/Qlv
    dGamma_dP  = -(dqi1dp+dqi2dp+Gamma_Evap*(dhv_dP-dhl_dP))/Qlv
    dGamma_dPa = -(dqi1dpa+dqi2dpa+Gamma_Evap*(dhv_dPa-dhl_dPa))/Qlv
    dGamma_dAlpha = -(dqi1da+dqi2da)/Qlv
  END SUBROUTINE Evap_Rate
END MODULE EVAPORATION
