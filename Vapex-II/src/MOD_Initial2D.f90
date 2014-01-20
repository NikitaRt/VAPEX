!--------------------------------------------------------------------!
!  Initialization/restart                                            !
!--------------------------------------------------------------------!
!  $Id: MOD_Initial2D.f90 9 2014-01-14 14:22:27Z Sergey $
!--------------------------------------------------------------------!
MODULE INITIAL_2D

CONTAINS
  SUBROUTINE INITIALIZE
    USE GLOBALS
    USE WATER_PROPS
    USE DISP_PROPS
    USE GLOBALS_2D
    USE INOUT_2D
    USE GRID_2D
    USE VARIABLES_2D
    USE PROBLEM_DATA
    USE SOLVER_2D
    USE COEFFICIENTS
    USE DISPERSED_PHASE
    USE ELLEQ
    USE OUTTRANS
    IMPLICIT NONE

    CALL DATE_AND_TIME(date=date,values=date_time)
    CALL ReadInputData
    !
    ! Allocate arrays
    !
    CALL SetGrid_2D
    CALL AllocateVariables_2D
    CALL AllocateGlobals_2D
    !
    ! Initialize elliptic solver
    !
    CALL InitEllipticSolver
    !
    ! Initialize equation of state
    !
    CALL SetNonCondensGasType(nonCondGasType,isMixtureWithH2)
    CALL SetGassAdditionFlag(isHydrogenRelease)
    CALL SetPhaseChangeFlag(isEvaporation)
    CALL InitWaterProps
    CALL InitDisp(isDispersedPhase)
    CALL SetAdaptiveTimeStep(2.D0,5)
    !
    ! Set initial values
    !
    IF(NREAD == 0) THEN
       CALL SetInitialValues

       tstep = tstepinit
    ELSE
       CALL ReadRestart
    ENDIF
    !
    ! Update phase densities, energies etc
    !
    CALL UpdatePhaseVars
    CALL SetBoundaryConditions
    CALL SetVelocityBoundaryConditions
    CALL UpdatePhaseViscosities
    !
    ! Initialize transducers
    !
    CALL InitTransducers
    IF(NREAD == 0)  CALL ResetTransducerFile
    CALL WriteTransducers
    !
  END SUBROUTINE INITIALIZE


  SUBROUTINE Finalize
    IMPLICIT NONE

    CALL SaveData

  END SUBROUTINE Finalize


  SUBROUTINE UpdatePhaseVars
    USE GLOBALS_2D
    USE VARIABLES_2D
    USE HYDROGEN
    USE WATER_PROPS
    USE LOCAL_VARS

    IMPLICIT NONE

    INTEGER:: i,j

    DO j = 2,m
       DO i = 2,n
          YH2 = aYH2(i,j)
          T1  = at1(i,j)
          T2  = at2(i,j)
          P   = bp(i,j)
          Pa  = bpa(i,j)
          CALL Thermo_Props
          ar1(i,j) = Ro_1
          ae1(i,j) = E_1
          ar2(i,j) = Ro_2
          ae2(i,j) = E_2
          bpv(i,j) = MAX(P-Pa,0.D0)
          ats(i,j) = sattmp(bp(i,j))
          atsv(i,j) = sattmp(bpv(i,j))
          ara(i,j) = Ro_a
          arH2(i,j) = aYH2(i,j)*Ro_a
          aa1(i,j) = (1.D0 - aa2(i,j) - al3(i,j))
          aa1(i,j) = MIN(MAX(aa1(i,j),0.D0),1.D0-al3(i,j))
          aYa(i,j) = Ro_a/Ro_2
       ENDDO
    ENDDO
  END SUBROUTINE UpdatePhaseVars

  SUBROUTINE InitialBalances
    USE GLOBALS
    USE PROBLEM_DATA
    USE DISPERSED_PHASE
    USE WATER_PROPS
    USE LOCAL_VARS
    USE HYDROGEN
    USE GRID_2D
    USE VARIABLES_2D
    USE GLOBALS_2D
    IMPLICIT NONE

    INTEGER:: i, j
    DOUBLE PRECISION:: dV
    !
    ! Initialize integral values
    !
    qmpart=0.d0
    qmdebris=0.d0
    qmjet=0.d0
    qmdrops=0.d0
    qmzero=0.d0

    ! Initial mass and energy balances
    amas10=0.d0;amas20=0.d0;amas30=0.d0; amasAr_0=0.d0;amasH2_0=0.d0;
    ener10=0.d0;ener20=0.d0;ener30=0.d0; FRtotal = 0.d0
    DO i=2,n
       dV = pi*(dxu(i)**2 - dxu(i-1)**2)*hj
       amas10=amas10 + &
            dV*SUM(aa1(i,2:m)*ar1(i,2:m))
       amas20=amas20 + &
            dV*SUM(aa2(i,2:m)*ar2(i,2:m))
       amas30=amas30 + &
            dV*rmelt*SUM(al3(i,2:m))
       amasAr_0=amasAr_0 + &
            dV*SUM(aa2(i,2:m)*ara(i,2:m)*(1.D0-aYH2(i,2:m)))
       amasH2_0=amasH2_0 + &
            dV*SUM(aa2(i,2:m)*arH2(i,2:m))
       ener10 = ener10 + &
            dV*SUM(aa1(i,2:m)*ar1(i,2:m)*ae1(i,2:m))
       ener20 = ener20 + &
            dV*SUM(aa2(i,2:m)*ar2(i,2:m)*ae2(i,2:m))
       ener30 = ener30 + &
            dV*rmelt*cmelt*SUM(al3(i,2:m)*at3(i,2:m))
       FRtotal = FRtotal + &
            dV*SUM(FRateLocal(i,2:m))
    ENDDO
    enermix0 = ener10+ener20+ar2FB*ae2FB*VolumeFB
    ener3disp = 0.d0 ! Energy of dispersed particles (calculated from Lagrangian particles)

    ! Mass balance of steam
    amasH2O_0=0.d0
    DO i = 2,n
       dV = pi*(dxu(i)**2 - dxu(i-1)**2)*hj
       DO j = 2, m
          YH2 = aYH2(i,j)
          CALL Thermo_Props(at1(1,j),at2(1,j),bp(i,j),bpa(i,j))
          amasH2O_0=amasH2O_0 + &
               dV*(Ro_2-Ro_a)*aa2(i,j)
       ENDDO
    ENDDO

    ! Freeboard mass balance
    amas2_fb = ar2FB*VolumeFB                 ! Mass of vapour (non-condens+steam)
    amasAr_fb = araFB*(1.D0-aYH2FB)*VolumeFB ! Mass of argon
    amasH2_fb = arH2FB*VolumeFB        ! Mass of hydrogen
    amasH2O_fb = arH2OFB*VolumeFB             ! Mass of steam
    !
    amasH2O_total0 = amas10 + amasH2O_0 + amasH2O_fb ! Total mass of water
    amasAr_total0  = amasAr_0 + amasAr_fb  ! Total mass of argone
    amasH2_total0  = amasH2_0 + amasH2_fb  ! Total mass of hydrogen
    ! 
    ! Fragmented mass
    FragMass = 0.d0
    ! 
    ! Min. and max. values
    !
    alphaAv = 0.d0
    ! Max. velocity of liquid
    vlmax=0.5*SQRT(MAXVAL((au1(2:n,2:m)+au1(1:n-1,2:m))**2+&
         (av1(2:n,2:m)+av1(2:n,1:m-1))**2))
    ! Max. velocity of vapour
    vvmax=0.5*SQRT(MAXVAL((au2(2:n,2:m)+au2(1:n-1,2:m))**2+&
         (av2(2:n,2:m)+av2(2:n,1:m-1))**2))
    ! Max. T liquid and Tliquid-Ts (superheat)
    tlmax=MAXVAL(at1(2:n,2:m))
    tl_ts=MAXVAL(at1(2:n,2:m)-ats(2:n,2:m))
    ! Min. T liquid
    tlmin=MINVAL(at1(2:n,2:m))
    ! Max. T vapour
    tvmax=MAXVAL(at2(2:n,2:m))
    ! Min. T vapour
    tvmin=MINVAL(at2(2:n,2:m))
    ! Max. vol. fraction Alpha_3
    al3max=MAXVAL(al3(2:n,2:m))
    ! Min. height of dispersed phase
    zmin = hrelease; zminj = hrelease; zmindr = hrelease
    ! Max and min total pressure
    ptotmax=MAXVAL(bp(2:n,2:m))
    ptotmin=MINVAL(bp(2:n,2:m))
    !
    ! Dispersed phase
    !
    IF(isThirdPhase /= 0) dMfrag = 0

    IF(isOpenTop /= 0) THEN
       topflux_w = 2.*PI*hi*SUM(ar1(2:n,m)*aa1(2:n,m)*dn(2:n)*av1(2:n,m))
       topflux_g = 2.*PI*hi*SUM(ar2(2:n,m)*aa2(2:n,m)*dn(2:n)*av2(2:n,m))
       topflux_a = 2.*PI*hi*SUM(ara(2:n,m)*aa2(2:n,m)*dn(2:n)*av2(2:n,m))
       topflux_v = topflux_g - topflux_a
    ENDIF

  END SUBROUTINE InitialBalances

  SUBROUTINE CurrentBalances
    USE GLOBALS
    USE PROBLEM_DATA
    USE DISPERSED_PHASE
    USE LOCAL_VARS
    USE WATER_PROPS
    USE GRID_2D
    USE VARIABLES_2D
    USE GLOBALS_2D
    USE SOLVER_2D
    USE HYDROGEN
    IMPLICIT NONE

    INTEGER:: i,j
    DOUBLE PRECISION:: dV
    !======================================= Balances =========================
    ! Total mass and heat transfer rates
    agam=0.d0; agamH2=0.d0; qq31=0.d0; 
    quenchRate=0.d0; qrWater=0.d0; qrVapour=0.d0
    alphaAv = 0.d0
    DO i=2,n
       dV = pi*(dxu(i)**2 - dxu(i-1)**2)*hj
       agam   = agam + dV*SUM(gam(i,2:m))
       agamH2 = agamH2 + dV*SUM(gamH2(i,2:m))
       qq31   = qq31 + dV*SUM(htc1(i,2:m))
       quenchRate = quenchRate +&
            dV*SUM(R31c(i,2:m)*(at3(i,2:m)-at1(i,2:m))+&
            R32c(i,2:m)*(at3(i,2:m)-at2(i,2:m)))
       qrWater = qrWater +&
            dV*SUM(R31c(i,2:m)*(at3(i,2:m)-at1(i,2:m)))
       qrVapour = qrVapour +&
            dV*SUM(R32c(i,2:m)*(at3(i,2:m)-at2(i,2:m)))
       alphaAv = alphaAv +&
            dV*SUM(aa2(i,2:m)/(aa1(i,2:m)+aa2(i,2:m)),mask=dm(2:m)<=WaterLevel)
    END DO
    alphaAv = alphaAv/(pi*VesselRad**2*WaterLevel)
    ! Time integrals:
    qgam = qgam + agam*tstep
    qgamH2 = qgamH2 + agamH2*tstep
    quenchEnergy = quenchEnergy + quenchRate*tstep
    qeWater = qeWater + qrWater*tstep
    qeVapour = qeVapour + qrVapour*tstep
    ! Mass balance in the vessel
    amas1=0.d0;amas2=0.d0;amas3=0.d0;amasAr=0.d0;amasH2=0.d0;
    ener1=0.d0;ener2=0.d0;ener3=0.d0;FRtotal = 0.d0
    DO i=2,n
       dV = pi*(dxu(i)**2 - dxu(i-1)**2)*hj
       amas1=amas1 + &
            dV*SUM(aa1(i,2:m)*ar1(i,2:m))
       amas2=amas2 + &
            dV*SUM(aa2(i,2:m)*ar2(i,2:m))
       amas3=amas3 + &
            dV*rmelt*SUM(al3(i,2:m))
       amasAr=amasAr + &
            dV*SUM(aa2(i,2:m)*ara(i,2:m)*(1.D0-aYH2(i,2:m)))
       amasH2=amasH2 + &
            dV*SUM(aa2(i,2:m)*arH2(i,2:m))
       ener1 = ener1 + &
            dV*SUM(aa1(i,2:m)*ar1(i,2:m)*ae1(i,2:m))
       ener2 = ener2 + &
            dV*SUM(aa2(i,2:m)*ar2(i,2:m)*ae2(i,2:m))
       ener3 = ener3 + &
            dV*rmelt*cmelt*SUM(al3(i,2:m)*at3(i,2:m))
       FRtotal = FRtotal + &
            dV*SUM(FRateLocal(i,2:m))
    ENDDO
    enermix = ener1+ener2+ar2FB*ae2FB*VolumeFB

    ! Mass balance of steam
    amasH2O=0.d0
    DO i = 2,n
       dV = pi*(dxu(i)**2 - dxu(i-1)**2)*hj
       DO j = 2, m
          YH2 = aYH2(i,j)
          CALL Thermo_Props(at1(1,j),at2(1,j),bp(i,j),bpa(i,j))
          amasH2O=amasH2O + &
               dV*(Ro_2-Ro_a)*aa2(i,j)
       ENDDO
    ENDDO
    ! Freeboard mass balance
    amas2_fb = ar2FB*VolumeFB                 ! Mass of vapour (non-condens+steam)
    amasAr_fb = araFB*(1.D0-aYH2FB)*VolumeFB ! Mass of argone
    amasH2_fb = arH2FB*VolumeFB        ! Mass of hydrogen
    amasH2O_fb = arH2OFB*VolumeFB             ! Mass of steam

    amasH2O_total = amas1 + amasH2O + amasH2O_fb ! Total mass of water
    amasAr_total = amasAr + amasAr_fb            ! Total mass of argone
    amasH2_total = amasH2 + amasH2_fb            ! Total mass of hydrogen
    !
    ! Total fragmented mass
    !
    FragMass = FragMass + FRTotal*tstep
    !
    !------- Calculation of maximum and minimum values ---------------------
    ! Max. velocity of liquid
    vlmax=0.5*SQRT(MAXVAL((au1(2:n,2:m)+au1(1:n-1,2:m))**2+&
         (av1(2:n,2:m)+av1(2:n,1:m-1))**2))
    ! Max. velocity of vapour
    vvmax=0.5*SQRT(MAXVAL((au2(2:n,2:m)+au2(1:n-1,2:m))**2+&
         (av2(2:n,2:m)+av2(2:n,1:m-1))**2))
    ! Max. T liquid and Tliquid-Ts (superheat)
    tlmax=MAXVAL(at1(2:n,2:m))
    tl_ts=MAXVAL(at1(2:n,2:m)-ats(2:n,2:m))
    ! Min. T liquid
    tlmin=MINVAL(at1(2:n,2:m))
    ! Max. T vapour
    tvmax=MAXVAL(at2(2:n,2:m))
    ! Min. T vapour
    tvmin=MINVAL(at2(2:n,2:m))
    ! Max. vol. fraction Alpha_3
    al3max=MAXVAL(al3(2:n,2:m))
    ! Max and min total pressure
    ptotmax=MAXVAL(bp(2:n,2:m))
    ptotmin=MINVAL(bp(2:n,2:m))

    IF(isDispersedPhase /= 0) THEN
       IF(ll>0) THEN
          zmin = MINVAL(z(1:ll),mask=ind<=3)
          zminj = MINVAL(z(1:ll),mask=ind==1)
          zmindr = MINVAL(z(1:ll),mask=ind==2)
       ELSE
          zmin = hrelease
          zminj = hrelease
          zmindr = hrelease
       ENDIF
    ENDIF

    IF(isOpenTop /= 0) THEN
       topflux_w = 2.*PI*hi*SUM(ar1(2:n,m)*aa1(2:n,m)*dn(2:n)*av1(2:n,m))
       topflux_g = 2.*PI*hi*SUM(ar2(2:n,m)*aa2(2:n,m)*dn(2:n)*av2(2:n,m))
       topflux_a = 2.*PI*hi*SUM(ara(2:n,m)*aa2(2:n,m)*dn(2:n)*av2(2:n,m))
       topflux_v = topflux_g - topflux_a
    ENDIF

  END SUBROUTINE CurrentBalances

END MODULE INITIAL_2D


