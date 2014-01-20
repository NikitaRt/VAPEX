!------------------------------------------------------------------!
!  Solver for 2D calculation                                       !
!------------------------------------------------------------------!
!  $Id: MOD_Solver2D.f90 9 2014-01-14 14:22:27Z Sergey $
!------------------------------------------------------------------!
MODULE SOLVER_2D
  USE GLOBALS_2D
  USE GLOBALS
  USE PROBLEM_DATA
  IMPLICIT NONE
  REAL(kind=8):: tstep         ! Current time step
  REAL(kind=8):: tstepmin      ! Minimum possible tstep
  REAL(kind=8):: tstepmax      ! Maximum possible tstep
  REAL(kind=8):: tstepinit     ! Initial time tstep
  INTEGER:: isVariableTStep        ! if 0, tstep = tstepmax
  REAL(kind=8):: tstep_factor  ! tstep is divided or multiplied 
  REAL(kind=8):: n_tsteps      ! by tstep_factor, after n_tsteps
  ! attempt is made to increase tstep
  REAL(kind=8):: i_tsteps      ! Current time step counter
  REAL(kind=8):: CFLmax        ! Maximum CFL number
  REAL(kind=8):: CFL           ! Actual CFL number

  INTEGER:: MaxNonLinear           ! Max. no. of non-linear iterations
  REAL(kind=8):: TolNonLinear  ! Termination criterion

  INTEGER:: NOuter, NInner
  INTEGER:: iouter, iinner
  INTEGER:: IterNL

  REAL(kind=8):: ResidualNL(NVARS)
  REAL(kind=8):: ResidualUV(2*NDIMS)
  REAL(kind=8):: ResidualBCGS

  REAL(kind=8):: TimeMax ! Calculations stopped when either all
  ! iterations were carried out, or
  ! time > TimeMax

  INTEGER,PARAMETER:: GoOn = 1
  INTEGER,PARAMETER:: Converged = 0
  INTEGER,PARAMETER:: BadIteration = -1
  INTEGER,PARAMETER:: BadTimeStep = -2

CONTAINS
  SUBROUTINE SOLVER2D
    USE DISPERSED_PHASE_2D
    USE VARIABLES_2D
    USE LOCAL_VARS
    USE GLOBALS_2D
    USE LOCAL_VARS
    IMPLICIT NONE

    INTEGER:: info
    !--------------------- UNSTEADY CALCULATIONS -----------------
    CALL BeforeFirstStep
    OUTER_LOOP:  DO iouter=1,NOuter       ! Outer Loop

       INNER_LOOP: DO iinner=1,NInner     ! Inner Loop
          CALL BeforeEachStep

          CALL SolveForDispersedPhase

          CALL Current2PreviousTimeStep    ! Assign f^n
1000      CONTINUE

          CALL InitialGuess

          NEWTON_ITERATION:  DO IterNL = 1,MaxNonLinear

             CALL Current2PreviousIteration 

             CALL UpdatePhaseViscosities

             CALL SolveForVelocities

             CALL SetEquationCoefficients

             CALL SolveForVarIncrements

             CALL CheckIncrements(info)

             SELECT CASE(info)
             CASE(Converged)
                EXIT
             CASE(BadTimeStep)
                CALL Previous2Current
                CALL ReduceTimeStep
                GOTO 1000
             CASE(BadIteration)
                CALL WaterPack
                CYCLE
             CASE default
                CALL UpdateVariables
             END SELECT

             CALL SolveForHydrogen

          ENDDO NEWTON_ITERATION
          time = time + tstep
          CALL AfterOneStep
          CALL AdaptTimeStep
          CALL ScreenOutput(time,iinner,iouter)
          IF(isFinished(time) /= 0) GOTO 2000
       ENDDO INNER_LOOP
       CALL SaveData
    ENDDO OUTER_LOOP
2000 CALL AfterLastStep
  END SUBROUTINE SOLVER2D

  SUBROUTINE SetAdaptiveTimeStep(tfactor,n_steps)
    IMPLICIT NONE
    REAL(kind=8), INTENT(in):: tfactor
    INTEGER, INTENT(in):: n_steps

    tstep_factor = tfactor ! tstep is divided or multiplied 
    n_tsteps = n_steps        
  END SUBROUTINE SetAdaptiveTimeStep

  INTEGER FUNCTION isFinished(time)
    IMPLICIT NONE
    REAL(kind=8), INTENT(in):: time

    isFinished = 0
    IF(time >= TimeMax) isFinished = 1
  END FUNCTION isFinished


  SUBROUTINE ScreenOutput(time,iinner,iouter)
    REAL(8):: time
    INTEGER:: iinner, iouter


!    PRINT *,iinner,iouter,time

  END  SUBROUTINE ScreenOutput


END MODULE SOLVER_2D
