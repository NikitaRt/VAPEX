!----------------------------------------------------------------!
!  Vapex-II: main program                                         !
!----------------------------------------------------------------!
!  $Id: Vapex-II.f90 5 2014-01-13 11:44:50Z Sergey $
!----------------------------------------------------------------!
PROGRAM VAPEX_II
  USE SOLVER_2D
  USE INITIAL_2D
  implicit none

  CALL Initialize
  CALL Solver2D
  CALL Finalize

END PROGRAM VAPEX_II
