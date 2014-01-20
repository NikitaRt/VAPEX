!------------------------------------------------------------------!
!  Dispersed phase properties                                      !
!------------------------------------------------------------------!
!  $Id: MOD_Disp.f90 5 2014-01-13 11:44:50Z Sergey $
!------------------------------------------------------------------!
MODULE DISP_PROPS

  REAL(kind=8):: &
       rmelt,&    ! Melt density
       cmelt,&    ! Melt heat capacity
       TMelting   ! Melting temperature

  REAL(kind=8), PARAMETER, PRIVATE:: em0 = 0.7d0
  REAL(kind=8), PARAMETER, PRIVATE:: el0 = 0.3d0
  REAL(kind=8), PARAMETER:: sigem=5.67d-08*em0
  REAL(kind=8), PARAMETER:: sigemel=sigem*el0

  REAL(kind=8):: pr6, prc6

CONTAINS
  SUBROUTINE InitDispProps
    USE GLOBALS
    pr6=6.d0/(pi*rmelt)
    prc6=6.d0/(rmelt*cmelt)
    TMelting = 2830.D0
  END SUBROUTINE InitDispProps
END MODULE DISP_PROPS



!------------------------------------------------------------------!
! Lagrangian particles (melt jet, droplets, debris)                !
!------------------------------------------------------------------!
MODULE DISPERSED_PHASE
  USE DISP_PROPS
  !
  ! Flag
  !
  INTEGER:: isThirdPhase
  !
  ! Particles
  !
  REAL(kind=8):: &
       qmpart0, & ! Total mass of melt
       qmpart, &  ! Current mass of melt released
       qmjet, &   ! Mass of melt jet (ind = 1)
       qmdrops, & ! Mass of melt droplets (ind = 2)
       qmdebris,& ! Mass of debris (ind = 3)
       qmzero     ! "Zero" particles (ind = 4)

  REAL(kind=8):: &
       qmpart_in  ! Mass of macroparticles representing jet
  INTEGER:: nump_in
  INTEGER:: nump

  INTEGER:: &
       ll,lljet,lldrops,lldebris,llzero ! Current numbers of macroparticles

  INTEGER:: &
       lid, &     ! Index of the leading particle
       llj        ! Index of the particle last to leave nozzle

  REAL(kind=8):: &
       HRelease,& ! Height of melt supply 
       Diam, &    ! Nozzle diameter
       Spour, &   ! Nozzle area
       Djet, Ddrop ! Diameters of macroparticles representing jet and drops

  REAL(kind=8):: &
       dropConeAngle ! Tangent of drop cone angle

  REAL(kind=8), PARAMETER:: sigme=45.0d-2 ! Surface tension

  INTEGER:: isUseVin
  REAL(kind=8):: Vin ! Maximum vertical velocity of jet particles
  REAL(kind=8):: Tin, Din ! Temperature and diameter of macroparticles leaving nozzle

  !------Z_JetPart - distance between jet particles corresponding to qmpart_in kg of melt
  !------and condition that jet is incompressible - distance between the centers
  !------of the jet particles is the same
  REAL(kind=8):: Z_JetPart ! Height of cylindrical macroparticle representing melt jet

  REAL(kind=8):: debrisHeight ! Height of debris bed

  REAL(kind=8):: tstep_disp
  REAL(kind=8):: time_disp


  INTEGER:: ll3
  LOGICAL:: isSaito
  REAL(kind=8):: fragL_D
  INTEGER:: nfrag ! Number of macroparticles (ind=2) created by fragmentation of a jet particle
  INTEGER:: nskip_frag ! fragmentation occurs every nskip_frag time steps
  INTEGER:: ncount_frag ! step counter used to trigger fragmentation
  REAL(kind=8):: DropM

  REAL(kind=8):: VdropMin

  REAL(kind=8):: FRate    ! Total fragmentation rate
  REAL(kind=8):: FragMass ! Total fragmented mass
  REAL(kind=8):: FRtotal  ! Vol. integral of FRateLocal

  REAL(kind=8),ALLOCATABLE:: &
       x(:),z(:),           & ! Coordinates of macroparticles
       xu3(:),xv3(:),       & ! Velocities of macroparticles
       xt3(:),              & ! Temperature of macroparticles
       xd3(:),              & ! Diameters of macroparticles
       qmk(:),              & ! Masses of macroparticles
       we(:),sai(:),epf(:), & ! Criteria for fragmentation
       c31k(:),c32k(:),     & ! Drag coefficients
       r31k(:),r32k(:),     & ! Heat transfer coefficients
       dMfrag(:)              ! Mass residuals after fragmentaion

  INTEGER, ALLOCATABLE:: &
       ind(:) ! Indices of macroparticles (1-jet, 2-droplet, 3-debris)

  REAL(kind=8):: R31, R32
CONTAINS
  SUBROUTINE InitDisp(isDispPhase) ! Melt jet and droplet parameters
    USE GLOBALS
    USE DISP_PROPS
    IMPLICIT NONE

    INTEGER,INTENT(in):: isDispPhase

    IF(isDispPhase /= 0) THEN
       isThirdPhase = 1
    ELSE
       isThirdPhase = 0
    ENDIF

    CALL InitDispProps

    Djet = Diam
    Spour=pi*Diam**2/4.d0
    qmpart_in = qmpart0/nump_in
    Z_JetPart=qmpart_in/(rmelt*Spour) 
    DropM = rmelt*pi*Ddrop**3/6.D0

    IF(fragL_D < 0.d0) THEN
       isSaito = .TRUE.
    ELSE
       isSaito = .FALSE.
    ENDIF

    IF(isThirdPhase) THEN
       ALLOCATE( &
            x(nump),z(nump),xu3(nump),xv3(nump), &
            xt3(nump),xd3(nump),qmk(nump),ind(nump))
       ALLOCATE( &
            we(nump),sai(nump),epf(nump), &
            dMfrag(nump))
       ALLOCATE(c31k(nump),c32k(nump),r31k(nump),r32k(nump))
    ENDIF
  END SUBROUTINE InitDisp

  SUBROUTINE SetLocalHeatExchangeDisp(R_3_1,R_3_2)
    IMPLICIT NONE
    REAL(kind=8),INTENT(in):: R_3_1,R_3_2

    R31 = R_3_1
    R32 = R_3_2
  END SUBROUTINE SetLocalHeatExchangeDisp

  SUBROUTINE SetDispPhaseTimeStep(time_step)
    IMPLICIT NONE
    REAL(kind=8),INTENT(in):: time_step

    tstep_disp = time_step
  END SUBROUTINE SetDispPhaseTimeStep

END MODULE DISPERSED_PHASE
