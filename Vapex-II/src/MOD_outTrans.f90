!------------------------------------------------------------------!
!  Transducers                                                     !
!------------------------------------------------------------------!
!  $Id: MOD_outTrans.f90 10 2014-01-20 10:41:15Z Sergey $
!------------------------------------------------------------------!
MODULE OUTTRANS

  TYPE TRANS_TYPE
     CHARACTER(len=256) :: Name = ""
     CHARACTER(len=256) :: ControlParameter = ""
     REAL(8) :: Coord(2) = (/0.D0,0.D0/)
     INTEGER :: Ind(2) = (/-1,-1/)
     INTEGER :: TimeDep = 1
     LOGICAL :: Active = .FALSE.
  endtype TRANS_TYPE

  INTEGER,PARAMETER:: TRANSDUCER_POINT = 0
  INTEGER,PARAMETER:: TRANSDUCER_INTEGRAL = 1
  INTEGER,PARAMETER:: TRANSDUCER_STORAGE = 2

  TYPE(TRANS_TYPE), POINTER :: Trans(:)

  INTEGER:: NTrans = 0
  INTEGER, PARAMETER:: UNITTRANS = 45
  CHARACTER(len=256):: FileTRANS = 'Trans.dia'
  INTEGER,PARAMETER:: lenTRANS = 14
  CHARACTER(len=20):: formTRANS = '(1x,ES12.5)'
  LOGICAL:: IS_WRITE_TRANS = .TRUE.

  CHARACTER(len=256):: FileINFO = 'VAPEX.info'
  INTEGER,PARAMETER:: UNITINFO = 95

CONTAINS

  REAL(kind=8) FUNCTION MassOfWater()
    USE Globals
    USE GLOBALS_2D
    USE VARIABLES_2D
    USE GRID_2D

    IMPLICIT NONE

    MassOfWater = SUM(ar1(2:n,2:m)*aa1(2:n,2:m)*vol(2:n,2:m))
  END FUNCTION MassOfWater

  REAL(kind=8) FUNCTION MassOfGas()
    USE Globals
    USE GLOBALS_2D
    USE VARIABLES_2D
    USE GRID_2D

    IMPLICIT NONE

    MassOfGas = SUM(ar2(2:n,2:m)*aa2(2:n,2:m)*vol(2:n,2:m))
  END FUNCTION MassOfGas

  SUBROUTINE WriteTransducers
    USE Globals
    USE GLOBALS_2D
    USE VARIABLES_2D
    USE SOLVER_2D
    USE DISPERSED_PHASE

    IMPLICIT NONE
    CHARACTER(len=2048):: buf
    CHARACTER(len=lenTRANS):: valbuf
    REAL(kind=4):: val
    INTEGER:: iTR
    LOGICAL:: EXISTS

    IF(.NOT.IS_WRITE_TRANS) RETURN
    IF(LEN_TRIM(FileTRANS)<= 0 .OR. NTrans <= 0) RETURN

    INQUIRE(file=TRIM(ADJUSTL(FileTRANS)),exist = EXISTS)
    IF(.NOT.EXISTS) CALL ResetTransducerFile

    OPEN(UNITTRANS,file=FileTRANS,position='append')
    buf = ""
    WRITE(valbuf,formTRANS) Time
    buf = ADJUSTR(valbuf)
    DO iTR = 1,NTrans
       IF(.NOT.Trans(iTR)%Active) CYCLE
       SELECT CASE(Trans(iTR)%ControlParameter)
       CASE("MLiq")
          val = MassOfWater()
       CASE("MGas")
          val = MassOfGas()
       CASE("MDisp")
          val = qmpart
       CASE("MJet")
          val = qmJet 
       CASE("MDrops")
          val = qmDrops
       CASE("MDebr")
          val = qmDebris
       CASE("TStep")
          val = TStep
       CASE("CFL")
          val = CFL
       CASE("TLiq")
          val = at1(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("TVap")
          val = at2(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("TSat")
          val = ats(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("DDisp")
          val = dr3(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("P")
          val = bp(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("ULiq")
          val = au1(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("VLiq")
          val = av1(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("UVap")
          val = au2(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("VVap")
          val = av2(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("eLiq")
          val = ae1(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("eVap")
          val = ae2(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("ALiq")
          val = aa1(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("AGas")
          val = aa2(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("ADisp")
          val = al3(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("TDisp")
          val = at3(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("TJet")
          val = at3j(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("TDrops")
          val = at3dr(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE("TDebris")
          val = at3de(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       CASE default
          CYCLE
       END SELECT
       WRITE(valbuf,formTRANS) val
       buf = TRIM(buf)//ADJUSTR(valbuf)
    ENDDO
    WRITE(UNITTRANS,'(a)') TRIM(buf)
    CLOSE(UNITTRANS)
  END SUBROUTINE WriteTransducers


  INTEGER FUNCTION TransducerType(TR)
    IMPLICIT NONE
    TYPE(TRANS_TYPE), POINTER:: TR

    TransducerType = TRANSDUCER_POINT
    SELECT CASE(TR%ControlParameter)
    CASE("MLiq","MGas","MDisp","MJet","MDrops","MDebris","TStep","CFL") ! <=== Только интегральные!
       TransducerType = TRANSDUCER_INTEGRAL
    END SELECT
  END FUNCTION TransducerType


  SUBROUTINE CutTransducerFile
    USE Globals
    USE GLOBALS_2D
    USE VARIABLES_2D

    IMPLICIT NONE

    CHARACTER(len=2048):: buf
    LOGICAL:: EXISTS
    REAL(kind=4):: t

    IF(.NOT.IS_WRITE_TRANS) RETURN
    IF(LEN_TRIM(FileTRANS)<= 0) RETURN
    INQUIRE(file=TRIM(ADJUSTL(FileTRANS)),exist = EXISTS)
    IF(.NOT.EXISTS) RETURN

    OPEN(UNITTRANS,file=FileTRANS)
    READ(UNITTRANS,'(a)',END=200) buf ! Skip header
200 CONTINUE
    DO
       READ(UNITTRANS,'(a)',END=300) buf ! Skip header
       READ(buf,*) t
       IF(t > Time) THEN
          BACKSPACE UNITTRANS
          ENDFILE UNITTRANS
          GOTO 300
       ENDIF
    ENDDO
300 CONTINUE
    CLOSE(UNITTRANS)
  END SUBROUTINE CutTransducerFile

  SUBROUTINE ResetTransducerFile
    USE Globals
    USE GLOBALS_2D
    USE VARIABLES_2D

    IMPLICIT NONE
    CHARACTER(len=2048):: buf
    CHARACTER(len=lenTRANS):: valbuf
    INTEGER:: iTR

    IF(.NOT.IS_WRITE_TRANS) RETURN
    IF(LEN_TRIM(FileTRANS)<= 0) RETURN
    OPEN(UNITTRANS,file=FileTRANS)
    REWIND(UNITTRANS)
    ENDFILE UNITTRANS
    IF(NTrans > 0) THEN
       REWIND(UNITTRANS)
       buf = ""
       valbuf = "Time"
       buf = ADJUSTR(valbuf)
       DO iTR = 1,NTrans
          IF(.NOT.Trans(iTR)%Active) CYCLE
          SELECT CASE(Trans(iTR)%ControlParameter)
          CASE("MLiq","MGas","MDisp","MJet","MDrops","MDebris","TStep","CFL","TLiq","TVap",&
               "TSat","DDisp","P","ULiq","UVap","VLiq","VVap","eLiq","eVap",&
               "ALiq","AGas","ADisp","TDisp","TJet","TDrops","TDebris")
             WRITE(valbuf,'(a)') TRIM(ADJUSTL(Trans(iTR)%Name))
             buf = TRIM(buf)//ADJUSTR(valbuf)
          END SELECT
       ENDDO
       WRITE(UNITTRANS,'(a)') TRIM(buf)
    ENDIF
    CLOSE(UNITTRANS)
  END SUBROUTINE ResetTransducerFile

  SUBROUTINE InitTransducers
    USE Globals
    USE GLOBALS_2D
    USE VARIABLES_2D
    USE GRID_2D

    IMPLICIT NONE
    INTEGER:: iTR,i,j, TType
    REAL(kind=8):: C(NDIMS)
    TYPE(TRANS_TYPE), POINTER:: ptr

    IF(.NOT.IS_WRITE_TRANS) RETURN
    TRANS_LOOP: DO iTR = 1,NTrans
       IF(MINVAL(Trans(iTR)%Ind) <= 0) THEN ! Define indices from coordinates
          C = Trans(iTR)%Coord
          Trans(iTR)%Ind = -1
          j = 2
          DO i = 2,n
             IF(C(1) <= dxu(i) .AND. C(1) >= dxu(i-1)) THEN
                Trans(iTR)%Ind(1) = i
                EXIT
             ENDIF
          ENDDO
          i = 2
          DO j = 2,m
             IF(C(2) <= dzv(j) .AND. C(2) >= dzv(j-1)) THEN
                Trans(iTR)%Ind(2) = j
                EXIT
             ENDIF
          ENDDO
       ELSE ! Define coordinates from indices
          Trans(iTR)%Coord = UNASSIGNED
          i = Trans(iTR)%Ind(1)
          j = Trans(iTR)%Ind(2)
          IF(i>=2 .AND. i<= N .AND. j>= 2 .AND. j<= m) THEN
             Trans(iTR)%Coord(1) = dn(i)
             Trans(iTR)%Coord(2) = dm(j)
          ENDIF
       ENDIF
       ptr => Trans(iTR)
       TType = TransducerType(ptr)
       SELECT CASE(TType)
       CASE(TRANSDUCER_POINT)
          Trans(iTR)%Active = MINVAL(Trans(iTR)%Ind) > 0 .AND. MINVAL(Trans(iTR)%Coord) > IFASSIGNED
       CASE(TRANSDUCER_INTEGRAL)
          Trans(iTR)%Active = .TRUE.
       END SELECT
    ENDDO TRANS_LOOP
  END SUBROUTINE InitTransducers

END MODULE OUTTRANS
