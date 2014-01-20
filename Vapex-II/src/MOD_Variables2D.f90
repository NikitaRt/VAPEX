!---------------------------------------------------------------------!
!  Variables for 2D calculation                                       !
!---------------------------------------------------------------------!
!  $Id: MOD_Variables2D.f90 10 2014-01-20 10:41:15Z Sergey $
!---------------------------------------------------------------------!
MODULE GLOBALS_2D
  INTEGER, PARAMETER:: NVARS = 5
  INTEGER, PARAMETER:: NRHS = 5

  INTEGER:: n,m,n9,m9 ! Problem size

  REAL(kind=8),ALLOCATABLE:: &
       QCoefMatrix(:,:,:,:)
  REAL(kind=8),ALLOCATABLE:: &
       dVars(:,:,:)
  REAL(kind=8),ALLOCATABLE:: &
       dVels(:,:,:)
CONTAINS
  SUBROUTINE AllocateGlobals_2D
    IMPLICIT NONE

    ALLOCATE(QCoefMatrix(NVARS,NRHS,n9,m9))
    ALLOCATE(dVars(NVARS,n9,m9))
    ALLOCATE(dVels(4,n9,m9))
  END SUBROUTINE AllocateGlobals_2D
END MODULE GLOBALS_2D

MODULE VARIABLES_2D
  !
  ! Liquid (phase 1)
  !
  REAL(kind=8),ALLOCATABLE:: &
       ar1(:,:), br1(:,:),  & ! Density of liquid (phase 1)
       aa1(:,:), ba1(:,:),  & ! Volume fraction alpha_1
       at1(:,:), bt1(:,:),  & ! Temperature, T1
       au1(:,:), bu1(:,:),  & ! U-Velocity, U1
       av1(:,:), bv1(:,:),  & ! V-Velocity, V1
       ae1(:,:), be1(:,:)     ! Internal energy, E1
  !
  ! Combined Gas (phase 2)
  !
  REAL(kind=8),ALLOCATABLE:: &
       ar2(:,:), br2(:,:),  & ! Density of gas (phase 2)
       aa2(:,:), ba2(:,:),  & ! Volume fraction alpha_2
       at2(:,:), bt2(:,:),  & ! Temperature, T2
       au2(:,:), bu2(:,:),  & ! U-Velocity, U2
       av2(:,:), bv2(:,:),  & ! V-Velocity, V2
       ae2(:,:), be2(:,:)     ! Internal energy, E2
  !
  ! Non-condensable gas
  !
  REAL(kind=8),ALLOCATABLE:: &
       ara(:,:),            & ! Density of non-condensible
       aYa(:,:), bYa(:,:),  & ! Mass fraction of non-condensable
       arH2(:,:),brH2(:,:), & ! Density of H2
       aYH2(:,:)              ! Mass fraction of H2
  !
  ! Pressure
  !
  REAL(kind=8),ALLOCATABLE:: &
       dp(:,:), & ! Pressure increment between iterations
       bp(:,:), & ! Total pressure
       bpa(:,:),& ! Partial pressure of non-condensible
       bpv(:,:)   ! Partial pressure of H2O steam
  REAL(kind=8),ALLOCATABLE:: &
       bpn(:,:), & ! Total pressure on previous time step
       bpan(:,:)   ! Partial pressure of non-condensible on n-th step
  !
  ! Phase velocities and pressure coefficients
  !
  REAL(kind=8),ALLOCATABLE:: &
       u1n(:,:), v1n(:,:),  & ! Velocities on previous
       u2n(:,:), v2n(:,:)   ! (n-th) time step

  REAL(kind=8),ALLOCATABLE:: &
       fu1l(:,:), fu1r(:,:),  & ! Velocity coefficients
       fv1l(:,:), fv1r(:,:),  & 
       fu2l(:,:), fu2r(:,:),  & ! Velocity coefficients
       fv2l(:,:), fv2r(:,:)

  REAL(kind=8),ALLOCATABLE:: &
       bmul(:,:), bmuv(:,:) ! Viscosities of liquid and gas phases
  !
  ! Dispersed phase
  !
  REAL(kind=8),ALLOCATABLE:: &
       al3(:,:),bl3(:,:),   & ! Volume fraction of dispersed phase
       al3j(:,:),           & ! -------//------ jet
       al3dr(:,:),          & ! -------//------ droplets
       al3de(:,:)             ! -------//------ debris

  REAL(kind=8),ALLOCATABLE:: &
       al3Molten(:,:),      & ! Volume fraction of particles with T>TMelting
       al3drMolten(:,:),    & ! -------//------ droplets
       at3Molten(:,:),      & ! T of particles with T>TMelting
       at3drMolten(:,:)       ! T of molten droplets

  REAL(kind=8),ALLOCATABLE:: &
       at3(:,:),            & ! Temperature of dispersed phase
       at3j(:,:),           & ! -------//------ jet
       at3dr(:,:),          & ! -------//------ droplets
       at3de(:,:)             ! -------//------ debris

  REAL(kind=8),ALLOCATABLE:: &
       dr3(:,:)               ! Diameter of dispersed phase

  REAL(kind=8),ALLOCATABLE:: &
       au3(:,:),av3(:,:)      ! Velocities of dispersed phase

  REAL(kind=8),ALLOCATABLE:: & 
       FRateLocal(:,:)        ! Local fragmentation rate 

  !
  ! Saturation temperature
  !
  REAL(kind=8),ALLOCATABLE:: ats(:,:), atsv(:,:)
  !
  ! Exchange terms
  !
  REAL(kind=8),ALLOCATABLE:: &
       gam(:,:) ! Mass exchange term
  REAL(kind=8),ALLOCATABLE:: &
       gamH2(:,:)           ! Hydrogen production source term
  REAL(kind=8),ALLOCATABLE:: &
       r31c(:,:),r32c(:,:),htc1(:,:),htc2(:,:)

CONTAINS
  SUBROUTINE AllocateVariables_2D
    USE GLOBALS
    USE GLOBALS_2D,ONLY: n,m,n9,m9
    IMPLICIT NONE

    ALLOCATE( &
         ar1(n9,m9), br1(n9,m9), & ! Density of liquid (phase 1)
         aa1(n9,m9), ba1(n9,m9), & ! Volume fraction alpha
         at1(n9,m9), bt1(n9,m9), & ! Temperature
         au1(n,m9) , bu1(n,m9) , & ! U-Velocity
         av1(n9,m) , bv1(n9,m) , & ! V-Velocity
         ae1(n9,m9), be1(n9,m9), & ! Internal energy
         ar2(n9,m9), br2(n9,m9), & ! Density of vapour (phase 2)
         aa2(n9,m9), ba2(n9,m9), & ! Volume fraction alpha
         at2(n9,m9), bt2(n9,m9), & ! Temperature
         au2(n,m9) , bu2(n,m9) , & ! U-Velocity
         av2(n9,m) , bv2(n9,m) , & ! V-Velocity
         ae2(n9,m9), be2(n9,m9), & ! Internal energy
         ats(n9,m9), atsv(n9,m9) & ! Saturation temperature at P
         )
    ALLOCATE(  &
         ara(n9,m9),             & ! Density of non-condensible
         aYa(n9,m9), bYa(n9,m9), & ! Mass fraction of non-condensable
         arH2(n9,m9),brH2(n9,m9),& ! Density of H2
         aYH2(n9,m9)             & ! Mass fraction of H2
         )
    ALLOCATE( &
         dp(n9,m9), & ! Pressure increment between iterations
         bp(n9,m9), & ! Total pressure
         bpa(n9,m9),& ! Partial pressure of non-condensible
         bpv(n9,m9),& ! Partial pressure of H2O steam
         bpn(n9,m9), & ! Total pressure on previous time step
         bpan(n9,m9) &
         )
    ALLOCATE( &
         u1n(n,m9), v1n(n9,m), & ! Velocities on previous
         u2n(n,m9), v2n(n9,m)  & ! time step
         )
    ALLOCATE( &
         bmul(n9,m9),bmuv(n9,m9) &
         )
    ALLOCATE( &
         fu1l(n,m9), fu1r(n,m9),  & ! Velocity coefficients
         fv1l(n9,m), fv1r(n9,m),  &
         fu2l(n,m9), fu2r(n,m9),  &
         fv2l(n9,m), fv2r(n9,m)   &
         )

    ALLOCATE(al3(n9,m9), bl3(n9,m9), al3j(n9,m9), al3dr(n9,m9), al3de(n9,m9))
    ALLOCATE(at3(n9,m9), at3j(n9,m9), at3dr(n9,m9), at3de(n9,m9),dr3(n9,m9))
    ALLOCATE(al3Molten(n9,m9),al3drMolten(n9,m9),at3Molten(n9,m9),at3drMolten(n9,m9))
    ALLOCATE(gam(n9,m9), gamH2(n9,m9))
    ALLOCATE(FRateLocal(n9,m9))
    ALLOCATE(r31c(n9,m9),r32c(n9,m9),htc1(n9,m9),htc2(n9,m9))
    ALLOCATE(au3(n,m9), av3(n9,m))

  END SUBROUTINE AllocateVariables_2D
END MODULE VARIABLES_2D



MODULE GRID_2D
  ! Geometry of vessel
  REAL(kind=8):: &
       VesselRad, VesselHeight,& ! Radius and height
       WaterLevel, &         ! Water level 
       HFreeSpace            ! Free space above the liquid level

  REAL(kind=8):: &
       VolumeFB   ! Freeboard volume

  INTEGER,PARAMETER:: Symmetry = 1
  INTEGER:: isOpenTop,isOpenSide

  REAL(kind=8):: HOpening0, HOpening
  INTEGER:: JOpening0,JOpening

  REAL(kind=8):: &
       hi, hj ! Cell sizen in radial and vertical directions
  REAL(kind=8):: &
       hi_inv, hj_inv ! 1/hi and 1/hj
  REAL(kind=8):: s_inv ! 1/(hi*hj)

  REAL(kind=8),ALLOCATABLE:: &
       dn(:), dm(:),        & ! Grid coordinates in radial and vertical directions
       dxu(:),dzu(:),       & ! Coordinates for U-Velocity
       dxv(:),dzv(:),       & ! -------//------ V_Velocity
       dxt(:),dzt(:),       & ! -------//------ Cell-centred vars (e.g., T)
       vol(:,:),            & ! Volume of cell (i,j)
       volc(:)                ! Volume of cell (reciprocal)

CONTAINS
  SUBROUTINE SetGrid_2D
    USE GLOBALS,ONLY:pi
    USE GLOBALS_2D, ONLY: n,m,n9,m9
    IMPLICIT NONE
    INTEGER:: i, j

    ALLOCATE( &
         dn(n9), dm(m9), & 
         dxu(n), dzu(m9),&
         dxv(n9),dzv(m), &
         dxt(n9),dzt(m9),volc(n9),vol(n9,m9) &
         )
    !------------------------------------------------------------------------
    ! Grid coordinates
    !------------------------------------------------------------------------
    hi=VesselRad/(n-1)
    hj=VesselHeight/(m-1)
    hi_inv = 1.D0/hi
    hj_inv = 1.D0/hj
    s_inv=1.D0/(hi*hj)

    dn(1)=-0.5d0*hi
    dm(1)=-0.5d0*hj
    DO i=2,n9
       dn(i)=0.5d0*hi + (i-2)*hi
    ENDDO
    DO j=2,m9
       dm(j)=0.5d0*hj + (j-2)*hj
    ENDDO
    DO i=1,n
       dxu(i)=hi*(i-1)
    ENDDO
    DO j=1,m9
       dzu(j)=-0.5d0*hj + (j-1)*hj
    ENDDO
    DO i=1,n9
       dxv(i)=-0.5d0*hi + (i-1)*hi
    ENDDO
    DO j=1,m
       dzv(j)=hj*(j-1)
    ENDDO
    DO i=1,n9
       dxt(i)=-0.5d0*hi + (i-1)*hi
    ENDDO
    DO j=1,m9
       dzt(j)=-0.5d0*hj + (j-1)*hj
    ENDDO
    DO i = 2,n
       vol(i,:) = pi*hj * (dxu(i)**2 - dxu(i-1)**2)
       volc(i)=1.d0/ (pi*hj * (dxu(i)**2 - dxu(i-1)**2) )
    ENDDO
    volc(1) = volc(2);volc(n9)=volc(n)
    vol(1,:) = vol(2,:);vol(n9,:)=vol(n,:)
    !
    ! Find indices corresponding to side opening
    !
    IF(isOpenSide /= 0) THEN
       DO j = 2,m
          IF(dzt(j) >= HOpening0) THEN
             JOpening0 = j
             EXIT
          ENDIF
       ENDDO
       DO j = m,2,-1
          IF(dzt(j) <= HOpening) THEN
             JOpening = j
             EXIT
          ENDIF
       ENDDO
       IF(JOpening == JOpening0) THEN
          IF(JOpening == m) THEN
             JOpening0 = JOpening0-1
          ELSE
             JOpening = JOpening + 1
          ENDIF
       ENDIF
    ENDIF
  END SUBROUTINE SetGrid_2D
END MODULE GRID_2D

MODULE PROBLEM_DATA
  !
  ! Liquid and gas
  !
  REAL(kind=8):: P0      ! Initial pressure, [Pa]
  REAL(kind=8):: TLiq0   ! Initial temperature of liquid
  REAL(kind=8):: TGas0   ! Initial temperature ov vapour
  REAL(kind=8):: Gravity ! Gravity acceleration
  REAL(kind=8):: Alpha_0 ! Initial gas contents in liquid
  !
  REAL(kind=8):: RoLiq0  ! Initial density of liquid
  REAL(kind=8):: RoGas0  ! Initial density of liquid
  REAL(kind=8):: PvTLiq0 ! Saturated vapour pressure at TLiq0
  REAL(kind=8):: PvTGas0 ! Saturated vapour pressure at TGas0
  !
  ! Freeboard
  !
  REAL(kind=8):: TempFB_0! Initial temperature of gas in freeboard
  REAL(kind=8)::&
       ar2FB, br2FB, &  ! Vapour density
       araFB, braFB, &  ! Non-condensible density
       apFB,  bpFB,  &  ! Pressure
       apaFB, &         ! Pressure of non-condensible
       ae2FB, be2FB, &  ! Vapour energy
       at1FB, at2FB, &  ! Temperatures
       aYH2FB, bYH2FB,& ! Hydrogen mass fractions
       arH2FB, brH2FB,& ! Partial densities of hydrogen
       aXH2FB, &        ! Molar (volume) fraction
       apH2FB, &        ! Hydrogen partial pressure
       arH2OFB          ! Steam density
  REAL(kind=8)::&
       PvFB_0, PaFB_0
  !
  ! Non-condensable gas
  !
  INTEGER:: nonCondGasType   ! 1 - Air, 2 - H2, 3 - He, 4 - Ar
  INTEGER:: isMixtureWithH2  ! 0 - pure gas, 1 - mixture with hydrogen
  INTEGER:: isHydrogenRelease! 0/1
  REAL(kind=8):: C_H2 ! Constant in hydrogen generation model
  REAL(kind=8):: T_H2 ! Temperature of the released hydrogen
  !
  ! Phase changes
  !
  INTEGER:: isEvaporation
  !
  ! Balances
  !
  REAL(kind=8):: amas10, amas20, amas30, amasAr_0, amasH2_0, amasH2O_0 ! Initial 
  REAL(kind=8):: amas1, amas2, amas3, amasAr, amasH2, amasH2O     ! Current
  REAL(kind=8):: amas2_fb, amasAr_fb, & ! Freeboard
       amasH2_fb, amasH2O_fb
  REAL(kind=8):: amasH2O_total0, amasAr_total0, & ! Total (vessel+freeboard)
       amasH2_total0                                  ! (initial values)
  REAL(kind=8):: amasH2O_total, amasAr_total, &   ! Total (vessel+freeboard)
       amasH2_total                                   ! (current values)
  REAL(kind=8):: agam, agamH2                     ! Vol. integrals of gam and gamH2
  REAL(kind=8):: ener10, ener20, ener30, enermix0 ! Initial energies
  REAL(kind=8):: ener1, ener2, ener3, enermix     ! Current energies
  REAL(kind=8):: quenchRate, quenchEnergy ! Heat transfer phase3->mixture and its time integral
  REAL(kind=8):: qrWater, qeWater ! Heat transfer phase3->water and its time integral
  REAL(kind=8):: qrVapour, qeVapour ! Heat transfer phase3->vapour and its time integral
  REAL(kind=8):: delm1, delm2, qq31,qqq31, apot1, apot2, qmp1, qmp2, pot1, pot2
  REAL(kind=8):: enerd, enerd1, enerd2
  REAL(kind=8):: ener3disp ! Energy of the released dispersed phase (calculated from Lagr. parts)
  REAL(kind=8):: qgam, qgamH2, gmax
  REAL(kind=8):: alphaAv
  !
  ! Flowrates through top boundary
  !
  REAL(kind=8):: topflux_w, topflux_g, topflux_a, topflux_v
  !
  ! Min. and max. values
  !
  REAL(kind=8):: vlmax,vvmax,tlmax,tlmin,tvmax,tvmin,al3max,tl_ts
  REAL(kind=8):: zmin, zminj,zmindr
  REAL(kind=8):: ptotmin,ptotmax
  !
  ! Dispersed phase on/off flag
  !
  INTEGER:: isDispersedPhase

END MODULE PROBLEM_DATA
