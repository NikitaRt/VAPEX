!------------------------------------------------------------------!
!  Correlation for drag and heat transfer coefficients             !
!------------------------------------------------------------------!
!  $Id: MOD_Correl.f90 5 2014-01-13 11:44:50Z Sergey $
!------------------------------------------------------------------!
MODULE CORRELATIONS
  USE GLOBALS
  REAL(kind=8), PRIVATE:: Phi

  REAL(kind=8):: DGasBubble, DWaterDrop
  REAL(kind=8):: VGas(NDIMS),VLiq(NDIMS)
  REAL(kind=8),PARAMETER:: DBubbleMin=1.d-4,DBubbleMax=1.d-1
  REAL(kind=8),PARAMETER:: DDropletMin=1.d-5,DDropletMax=1.d-2

  REAL(kind=8):: DDispPhase, VDisp(NDIMS)

  REAL(kind=8):: SurfaceTension
  REAL(kind=8):: Visc_Liq, Cp_Liq, Lambda_Liq
  REAL(kind=8):: Visc_Gas, Cp_Gas, Lambda_Gas

  REAL(kind=8):: R1s, dR1s_dAlpha
  REAL(kind=8):: R2s, dR2s_dAlpha
  REAL(kind=8):: RG
  REAL(kind=8):: RF1s, dRF1s_dAlpha

  REAL(kind=8):: CDrag12, CDrag12_a1, CDrag12_a2

  REAL(kind=8):: qpl, qpg ! Heat transfer coefficients of an individual droplet

  REAL(kind=8):: c13, c23 ! Drag coefficients for individual particle
  REAL(kind=8):: c13_a1, c23_a2

  REAL(kind=8), PRIVATE:: Vel, Re, Pr, Area
  REAL(kind=8), PRIVATE:: VRel(NDIMS)

  REAL(kind=8), PARAMETER, PRIVATE:: RFLASH = 1.D8

  REAL(kind=8), PARAMETER, PRIVATE:: al3lim = 1.D0 !0.62D0

CONTAINS
  SUBROUTINE HeatExchangeGasLiq
    USE GLOBALS
    USE LOCAL_VARS
    USE WATER_PROPS

    IMPLICIT NONE
    REAL(kind=8):: hl, hv
    REAL(kind=8):: R1s_a1, R2s_a1

    Phi = Alpha/MAX(Alpha_l+Alpha,allim) ! Void fraction
    VRel = VGas-VLiq
    Vel = SQRT(DOT_PRODUCT(VRel,VRel))
    hl = E1 + P/Ro1
    hv = E2 + P/Ro2

    Visc_Liq = viscl(hl,P)
    Cp_Liq = cpll(hl,P)

    IF(Phi < 0.7) THEN ! Water with dispersed bubbles
       DGasBubble = 8.d0*SurfaceTension/(Ro1*Vel**2 + 1.d-30)
       DGasBubble = MIN(DGasBubble,DBubbleMax)
       DGasBubble = MAX(DGasBubble,DBubbleMin)

       Re = Ro1*Vel*DGasBubble/Visc_Liq ! Reynolds number
       Pr = Visc_Liq*Cp_Liq/Lambda_Liq
       Area = Alpha_l/(Alpha_l+Alpha_3)

       dR1s_dAlpha = 6.D0*Area*Lambda_Liq/DGasBubble**2* &
            & (2.D0 + 0.6D0*SQRT(Re)*Pr**(1.D0/3.D0))

       dR2s_dAlpha = 6.D0*Area*Lambda_Gas/DGasBubble**2*2.D0

       R1s = dR1s_dAlpha*MAX(Alpha,allim)
       R2s = dR2s_dAlpha*MAX(Alpha,allim)

    ELSE ! Gas with dispersed water droplets
       DWaterDrop = 12.d0*SurfaceTension/(Ro2*Vel**2 + 1.d-30)
       DWaterDrop = MIN(DWaterDrop,DDropletMax)
       DWaterDrop = MAX(DWaterDrop,DDropletMin)

       Re = Ro2*Vel*DWaterDrop/Visc_Gas ! Reynolds number
       Pr = Visc_Gas*Cp_Gas/Lambda_Gas
       Area = Alpha/(Alpha+Alpha_3)

       R2s_a1 = 6.D0*Area*Lambda_Gas/DWaterDrop**2* &
            & (2.D0 + 0.6D0*SQRT(Re)*Pr**(1.D0/3.D0))
       R1s_a1 = 6.D0*Area*Lambda_Liq/DWaterDrop**2*2.D0

       R1s = R1s_a1*MAX(Alpha_l,allim)
       R2s = R2s_a1*MAX(Alpha_l,allim)

       dR1s_dAlpha = -R1s_a1
       dR2s_dAlpha = -R2s_a1
    ENDIF
    !
    ! Sensible heat flux: taking RG = R2s
    !
    RG = R2s
    !
    ! Flashing
    !
    IF(T1 >T_sat) THEN
       RF1s = RFLASH*MAX(Alpha_l,allim)
       dRF1s_dAlpha = -RFLASH
    ELSE
       RF1s = 0.D0
    ENDIF

  END SUBROUTINE HeatExchangeGasLiq
  !
  ! Heat exchange of individual droplet with liquid/gas
  !
  SUBROUTINE HeatExchangeDispParticle
    USE GLOBALS
    USE LOCAL_VARS
    USE WATER_PROPS
    USE DISP_PROPS

    IMPLICIT NONE

    REAL(kind=8):: hr, hc, DeltaT, Hevap

    Phi = Alpha/MAX(Alpha_l+Alpha,allim) ! Void fraction
    Hevap = HEvaporation(T1)

!!!DEBUG!!!
    Hevap = 2.25912676712702d06   ! AUXILIARY


    IF(Phi< 0.7d0) THEN ! Fuel-Liquid Heat Transfer 

       hr=sigem*(T3+T1)*(T3**2+T1**2) ! Radiative heat transfer

       VRel = VDisp-VLiq
       Vel = SQRT(DOT_PRODUCT(VRel,VRel))
       DeltaT = MAX(ABS(T3-T1),1.D-3)

       hc=2.98d0*SQRT(Ro2*Lambda_Gas* &
            & (Hevap + 0.68d0*Cp_Gas*DeltaT)*Vel/(DDispPhase*DeltaT))

       Area = Alpha_l/MAX((Alpha_l+Alpha),allim)

       qpl=1.D0*(hr+hc)*prc6/DDispPhase*Area
       qpg=0.d0 ! Heat exchange with liquid water only

    ELSE ! Fuel-Vapour-Liquid Heat Transfer

       VRel = VDisp-VGas
       Vel = SQRT(DOT_PRODUCT(VRel,VRel))

       Re = Ro2*Vel*DDispPhase/Visc_Gas ! Reynolds number
       Pr = Visc_Gas*Cp_Gas/Lambda_Gas

       Area = Alpha/MAX((Alpha_l+Alpha),allim)

       hc = Lambda_Gas/DDispPhase*(2.d0+0.6d0*SQRT(Re)*Pr**(1.d0/3.d0))

       qpl=0.d0 ! Heat exchange with gas phase only
       qpg= hc*prc6/DDispPhase*Area

    ENDIF

  END SUBROUTINE HeatExchangeDispParticle
  !
  ! Drag between liquid/gas phases
  !
  SUBROUTINE DragGasLiq
    USE GLOBALS
    USE LOCAL_VARS
    USE WATER_PROPS
    IMPLICIT NONE

    VRel = VGas-VLiq
    Vel = SQRT(DOT_PRODUCT(VRel,VRel))

    CALL dragiph (CDrag12,CDrag12_a1,CDrag12_a2,&
         Ro1,Ro2,Vel,SurfaceTension,Alpha,Alpha_l, &
         Visc_Liq,Visc_Gas)

  CONTAINS
    !****************************************************************
    ! Subroutine  dragiph    Theofanous
    ! Calculation of interphase drag

    ! cdragi - vapour-liquid drag coefficient
    ! cd12_a1 = cdragi/alliq (calculated to avoid 0/0 division for low alliq)
    ! cd12_a2 = cdragi/alvap (calculated to avoid 0/0 division for low alvap)
    ! fi - void fraction
    ! rol -liquid density
    ! rov -vapor  density
    ! vmod -
    ! sigma - surface tension
    ! alvap - vapour fraction
    ! alliq - liquid fraction
    !****************************************************************

    SUBROUTINE dragiph (cdragi,cd12_a1,cd12_a2,&
         rol,rov,vmod,sigma,alvap,alliq, &
         mul,muv)

      IMPLICIT NONE

      REAL(kind=8), INTENT(out):: cdragi, cd12_a1, cd12_a2
      REAL(kind=8), INTENT(in):: rol,rov,vmod,sigma,alvap,alliq, &
           mul,muv
      !
      ! Local variables
      !
      REAL(kind=8):: fi, f, ql, area, Cd, Re, al3
      !
      al3 = 1.D0 - (alvap+alliq)
      fi=alvap/(alvap+alliq)

      IF (fi <= 0.3d0) THEN ! Vapour is dispersed phase
         f=(1.d0 - fi)**1.5d0
         ql=8.d0*sigma/(rol*MAX(vmod,1.D-5)**2)
         Cd=(2.d0/3.d0) *ql*dsqrt(9.8d0*(rol-rov)/sigma)* &
              ( (1.d0 + 17.67d0*f**(6.d0/7.d0))/(18.67d0*f) )**2
         area=alliq/(alliq+al3)
         cd12_a2 = 0.75d0*area*rol*Cd/ql
         cdragi = alvap*cd12_a2
         cd12_a1 = cdragi/alliq
      ELSEIF (fi > 0.7d0) THEN ! Liquid is dispersed phase
         f=fi*fi*fi
         ql=12.d0*sigma/(rov*MAX(vmod,1.D-5)**2)
         Cd=2.d0/3.d0 *ql*dsqrt(9.8d0*(rol-rov)/sigma)* &
              ( (1.d0 + 17.67d0*f**(6.d0/7.d0))/(18.67d0*f) )**2
         area=alvap/(alvap+al3)
         cd12_a1 = 0.75d0*area*rov*Cd/ql
         cdragi = alliq*cd12_a1
         cd12_a2 = cdragi/alvap
      ELSE
         ql=4.d0*dsqrt(sigma/(9.8d0*(rol-rov)))
         Cd=8.d0/3.d0*(1.d0 - fi)**2
         area=alliq/(alliq+al3)
         cdragi=0.75d0*alvap*area*rol*Cd/ql
         cd12_a1 = cdragi/alliq
         cd12_a2 = cdragi/alvap
      ENDIF
    END SUBROUTINE dragiph
  END SUBROUTINE DragGasLiq
  !
  ! Drag coefficient of individual fuel particle
  !
  SUBROUTINE DragDispParticle
    USE GLOBALS
    USE LOCAL_VARS
    USE WATER_PROPS
    USE DISP_PROPS

    IMPLICIT NONE
    REAL(kind=8):: Vmod13, Vmod23

    VRel = VDisp-VLiq
    Vmod13 = SQRT(DOT_PRODUCT(VRel,VRel))

    VRel = VDisp-VGas
    Vmod23 = SQRT(DOT_PRODUCT(VRel,VRel))

    CALL dragf (Ro1,Ro2,Alpha,Alpha_l,DDispPhase, &
         Vmod13,Vmod23,Visc_Liq,Visc_Gas,c13,c23)

  CONTAINS
    !****************************************************************
    ! Subroutine  dragf    Theofanous
    ! Calculation of drag of fuel particles

    ! rol -liquid density
    ! rov -vapor  density
    ! alvap - vapour fraction
    ! alliq - liquid fraction
    ! d - diameter of fuel particle
    ! c13 - fuel-liquid drag coefficient
    ! c23 - fuel-vapour drag coefficient
    !****************************************************************
    SUBROUTINE dragf (rol,rov,alvap,alliq,d, &
         vmod13,vmod23,mul,muv,c13,c23)

      IMPLICIT NONE
      REAL(kind=8),INTENT(in):: rol,rov,alvap,alliq,d
      REAL(kind=8),INTENT(in):: vmod13,vmod23,mul,muv
      REAL(kind=8),INTENT(out):: c13,c23
      REAL(kind=8), PARAMETER:: al3lim = 0.62D0
      !
      ! Local variables
      !
      REAL(kind=8):: fi, alfuel, Cd, Re, f, area

      fi=alvap/(alvap+alliq)

      alfuel=1.d0 - alliq - alvap

      IF(alfuel < al3lim) THEN
         f=dsqrt(1.d0 - alfuel)*(1.d0 - alfuel/al3lim)**1.55d0

         Cd=0.d0
         IF(fi >= 0.7) THEN
            Re=d*rov*vmod23/muv
         ELSE
            Re=d*rol*vmod13/mul
         ENDIF

         IF(Re > 0) THEN
            Cd=(24.d0/Re) + (4.d0/dsqrt(Re)) + 0.4d0
            Cd=Cd*( (1.d0 + 17.67d0*f**(6.d0/7.d0))/(18.67d0*f) )**2
         ENDIF

         area=alliq/(alliq+alvap)
         c13=0.75d0*area*rol*Cd/d

         area=alvap/(alliq+alvap)
         c23=0.75d0*area*rov*Cd/d
      ELSE
         c13=0.d0
         c23=0.d0
      ENDIF

    END SUBROUTINE dragf

  END SUBROUTINE DragDispParticle

  SUBROUTINE DragGasLiqDispersed
    USE GLOBALS
    USE LOCAL_VARS
    USE WATER_PROPS
    IMPLICIT NONE
    !
    ! Local variables
    !
    REAL(kind=8):: al_3,f,fi,Cd
    REAL(kind=8):: VMod13, VMod23

    IF(Alpha_3 < 1.D-12) THEN
       c13=0.d0; c13_a1 = 0.D0
       c23=0.d0; c23_a2 = 0.D0
    ELSE
       al_3 = MIN(Alpha_3,al3lim)
       f=SQRT(1.d0 - al_3)*(1.d0 - al_3/al3lim)**1.55d0

       f = MAX(f,1.D-3)

       fi=Alpha/MAX(Alpha+Alpha_l,allim)
       VMod13 = SQRT(DOT_PRODUCT(VLiq-VDisp,VLiq-VDisp))
       VMod23 = SQRT(DOT_PRODUCT(VGas-VDisp,VGas-VDisp))

       IF(fi > 0.7) THEN
          Re=DDispPhase*Ro2*Vmod23/Visc_Gas
       ELSE
          Re=DDispPhase*Ro1*Vmod13/Visc_Liq
       ENDIF

       IF(Re > 0.d0) THEN
          ! Cd=24.d0/Re*(1.d0 + 0.179d0*dsqrt(Re) + 0.013d0*Re)
          ! Cd=24.d0/Re*(1.d0 + 0.167d0*dsqrt(Re) + 0.0167d0*Re)
          Cd=(24.d0/Re) + (4.d0/dsqrt(Re)) + 0.4d0
          Cd=Cd*( (1.d0 + 17.67d0*f**(6.d0/7.d0))/(18.67d0*f) )**2
       ELSE
          Cd = 0.D0
       ENDIF

       !     area=alliq/(alliq+alvap)
       c13_a1 = 0.75d0*Alpha_3*Ro1*Cd/DDispPhase/MAX(Alpha_l+Alpha,allim)
       c13 = c13_a1*Alpha_l

       !     area=alvap/(alliq+alvap)
       c23_a2 = 0.75d0*Alpha_3*Ro2*Cd/DDispPhase/MAX(Alpha_l+Alpha,allim)
       c23 = c23_a2*Alpha
    ENDIF
  END SUBROUTINE DragGasLiqDispersed

END MODULE CORRELATIONS
