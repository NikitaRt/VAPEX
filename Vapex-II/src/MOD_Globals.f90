!------------------------------------------------------------------!
!  Global and local variables                                      !
!------------------------------------------------------------------!
!  $Id: MOD_Globals.f90 5 2014-01-13 11:44:50Z Sergey $
!------------------------------------------------------------------!
MODULE GLOBALS
  REAL(kind=8), PARAMETER:: pi=3.14159265358979D0
  REAL(kind=8) :: allim
  INTEGER, PARAMETER:: NDIMS = 2

  CHARACTER(len=40) TestName
  CHARACTER(len=10) date
  INTEGER:: date_time(8)

  REAL(kind=8):: time

!!!Nikita
  CHARACTER (LEN=10) :: FileRep = "report.lst"
  CHARACTER (LEN=10) :: FilePAR = "vapex.inp"
  INTEGER, PARAMETER :: UNITPAR = 1
  INTEGER, PARAMETER :: UNITREP = 2
  CHARACTER (LEN=15) :: TRACE_LINE  
  REAL(8), PARAMETER :: UNASSIGNED = -1.1D30
  REAL(8), PARAMETER :: IFASSIGNED = -1.D30
  INTEGER(4), PARAMETER :: UNASSIGNED_INT = -10000000
  CHARACTER(LEN=75) :: ERR_READING
  CHARACTER(LEN=75) :: SCFL_READING
  INTEGER :: isSodium, isWater

END MODULE GLOBALS


MODULE LOCAL_VARS
  ! Current-iteration values:
  REAL(kind=8):: Ro1, Ro2, Roa
  REAL(kind=8):: Ya
  REAL(kind=8):: T1, T2, T3, Tsat
  REAL(kind=8):: P, Pa, Pv
  REAL(kind=8):: E1, E2
  REAL(kind=8):: Alpha, Alpha_l,Alpha_3
  ! Previous-iteration values:
  REAL(kind=8):: Ro1n, Ro2n
  REAL(kind=8):: Yan
  REAL(kind=8):: T1n, T2n
  REAL(kind=8):: E1n, E2n
  REAL(kind=8):: Alphan, Alpha_ln,Alpha_3n
CONTAINS
  SUBROUTINE SetLocalVars(Alpha1,Alpha2,Alpha3,TLiq,TGas,TDisp,Ptot,Pncon,Y_a)
    IMPLICIT NONE
    REAL(kind=8):: Alpha1,Alpha2,Alpha3,TLiq,TGas,TDisp,Ptot,Pncon,Y_a

    Alpha_l = Alpha1; Alpha = Alpha2; Alpha_3 = Alpha3
    T1 = TLiq; T2 = TGas; T3 = TDisp
    P = Ptot; Pa = Pncon
    Ya = Y_a
  END SUBROUTINE SetLocalVars

  SUBROUTINE SetPreviousValues(Alpha1_n,Alpha2_n,TLiq_n,TGas_n,ELiq_n,EGas_n,&
       Ro1_n,Ro2_n,Ya_n)
    IMPLICIT NONE
    REAL(kind=8):: Alpha1_n,Alpha2_n,TLiq_n,TGas_n,ELiq_n,EGas_n,&
         Ro1_n,Ro2_n,Ya_n

    Alpha_ln = Alpha1_n; Alphan = Alpha2_n
    T1n  = TLiq_n; T2n = TGas_n
    E1n  = ELiq_n; E2n = EGas_n
    Ro1n = Ro1_n; Ro2n = Ro2_n
    Yan  = Ya_n
  END SUBROUTINE SetPreviousValues


END MODULE LOCAL_VARS

