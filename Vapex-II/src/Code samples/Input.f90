!###################################################################!
!                                                                   !
! DECOSIM:  DEbris COolability SIMulator                            !
!                                                                   !
!###################################################################!
! $Id: Input.f90 5 2014-01-13 11:44:50Z Sergey $
!-------------------------------------------------------------------!
! Input of initial data                                             !
!-------------------------------------------------------------------!
subroutine INPUT
  use Globals
  use SolutionVector
  use Elleq
  use Water
  use DebrisBed
  use Correlations
  use PVolFrac
  use PhaseFlowStagger
  use Geometry,only:SetSymmetry
  use Multiphase
  use Turb
  use Material
  use SolidMaterial
  use Particle
  use Output
  use InitVal
!  use IFPORT,ONLY:makedirqq
  implicit none
  !--------------------------------------------------------------------------!
  character(len=256) buf
  logical EXISTS
  integer:: NDB
  integer:: iMat, pos
  integer:: size_seed, iFeed
  integer:: res
  logical:: resMkDir,resDirExists
  !--------------------------------------------------------------------------!

  FilePAR='Initial.par'
  call get_command_argument(1,buf)
  if(len_trim(buf) > 0) then
     FilePAR = trim(buf)
  endif
  !---------------------------- Check if the file exists -------------------!
  inquire(file = trim(adjustl(FilePAR)), exist = EXISTS)
  if(.not.EXISTS) then
     call Salute
     write(*,444) trim(adjustl(FilePAR))
     stop
  endif
  !-------------------------------- Open file ------------------------------!
  open(unit = UNITPAR,file=FilePAR,form='formatted', status ='old')
  TRACE_LINE = ""
  call SetupRegisters(10)
  !
  if(.not.ReadProblem()) goto 2000
  !
  ! Options
  !
  if(.not.ReadOptions()) goto 2000
  !
  ! Gap/Tooth
  !
  if(.not.ReadGapTooth()) goto 2000
  !
  ! Problem scales for non-dimensionalization
  !
  if(.not. ReadScales())  goto 2000
  !
  ! Phase materials (liquid, gas)
  !
  if(.not.ReadPhaseMaterials()) goto 2000
  !
  ! Multiphase model
  !
  if(.not.ReadMultiphaseModel()) goto 2000
  !
  ! Numerical scheme
  !
  if(.not.ReadNumerics()) goto 2000
  !
  ! Correlations
  !
  if(.not.ReadCorrelations()) goto 2000
  !
  ! Energy model
  !
  if(.not.ReadEnergyModel()) goto 2000
  !
  ! Turbulence model
  !
  if(.not.ReadTurbulenceModel()) goto 2000
  !
  ! Solid materials
  !
  if(.not.ReadSolidMaterials()) goto 2000
  !
  ! Storages
  !
  if(.not.ReadStorages()) goto 2000
  !
  ! Time stepping
  !
  if(.not.ReadTimeStepping()) goto 2000
  !
  ! Scheme for flow equations
  !
  if(.not.ReadSchemeFlow()) goto 2000
  !
  ! Scheme for Alpha-Pressure equations
  !
  if(.not.ReadSchemeAlpha()) goto 2000
  !
  ! Iterations
  !
  if(.not.ReadIterations()) goto 2000
  !
  ! Boundary conditions
  !
  if(.not.ReadBoundaryConditions()) goto 2000
  !
  ! Initial conditions on rectangular patches
  !
  if(.not.ReadPatches()) goto 2000
  !
  ! Output frequency
  !
  if(.not.ReadOutput()) goto 2000
  !
  ! Output files
  !
  if(.not.ReadOutputFiles()) goto 2000
  !
  ! Transducers
  !
  if(.not.ReadTransducers()) goto 2000
  !
  ! Output: profiles
  !
  call ReadProfiles
  !
  ! Restart/Save controls, including Gap/Tooth
  !
  if(.not.ReadRestartControls()) goto 2000
  !
  ! Porous regions
  !
  if(.not.ReadPorousRegions()) goto 2000
  !
  ! Particles
  !
  if(.not.ReadParticles()) goto 2000
  !
  ! Coolability task - always after regions read
  !
  if(.not.ReadCoolability()) goto 2000
  ! 
  ! Make all output go to OutputDir
  !
!!$ INQUIRE(DIRECTORY=trim(OutputDirPrefix),EXIST=resDirExists)
!!$  if(.not.resDirExists) then
!!$      ! Try to create the output directory
!!$      resMkDir=makedirqq(trim(OutputDirPrefix))
!!$      if(.not.resMkDir) then
!!$          print *,'Unable to create the output directory ',trim(OutputDirPrefix)
!!$        stop
!!$      endif
!!$    endif
      
  if(len_trim(OutputDirPrefix) > 0 .and. trim(adjustl(OutputDirPrefix)) /= './') then
     pos = len_trim(OutputDirPrefix)
     if(OutputDirPrefix(pos:pos) /= '/' .or. OutputDirPrefix(pos:pos) /= '\\') then
        OutputDirPrefix = trim(adjustl(OutputDirPrefix))//'/'
     endif
     if(.not.IS_MAKE_GLOBAL_STEP) &
          FileRST = trim(adjustl(OutputDirPrefix))//trim(adjustl(FileRST))
     FileSAVE = trim(adjustl(OutputDirPrefix))//trim(adjustl(FileSAVE))
     do iFeed = 1,NFeeders
        if(len_trim(Feeders(iFeed)%FeederFile) > 0) then
           Feeders(iFeed)%FeederFile = trim(adjustl(OutputDirPrefix))//&
                trim(adjustl(Feeders(iFeed)%FeederFile))
        endif
        if(len_trim(Feeders(iFeed)%CatcherFile) > 0) then
           Feeders(iFeed)%CatcherFile = trim(adjustl(OutputDirPrefix))//&
                trim(adjustl(Feeders(iFeed)%CatcherFile))
        endif
        if(len_trim(Feeders(iFeed)%CatcherCDFFile) > 0) then
           Feeders(iFeed)%CatcherCDFFile = trim(adjustl(OutputDirPrefix))//&
                trim(adjustl(Feeders(iFeed)%CatcherCDFFile))
        endif
     enddo
     if(len_trim(FileMassBal) > 0) then
        FileMassBal = trim(adjustl(OutputDirPrefix))//trim(adjustl(FileMassBal))
     endif
     if(len_trim(FileMinMax) > 0) then
        FileMinMax = trim(adjustl(OutputDirPrefix))//trim(adjustl(FileMinMax))
     endif
     if(len_trim(FileCheckAV) > 0) then
        FileCheckAV = trim(adjustl(OutputDirPrefix))//trim(adjustl(FileCheckAV))
     endif
     if(len_trim(FileSTOP) > 0) then
        FileSTOP = trim(adjustl(OutputDirPrefix))//trim(adjustl(FileSTOP))
     endif
     if(len_trim(FileTRANS) > 0) then
        FileTRANS = trim(adjustl(OutputDirPrefix))//trim(adjustl(FileTRANS))
     endif
     if(len_trim(FileINFO) > 0) then
        FileINFO = trim(adjustl(OutputDirPrefix))//trim(adjustl(FileINFO))
     endif
  endif

  !----------------------------------------------------------------------------!
  NNonCondensable = 0
  !----------------------------------------------------------------------------!
  close(UNITPAR)
  return
  !----------------------------------------------------------------------------!
2000 continue
  write(*,560) trim(adjustl(FilePAR)),trim(adjustl(TRACE_LINE))
  close(UNITPAR)
  stop

2200 continue
  write(*,777) trim(adjustl(FilePAR)),NDB,trim(adjustl(FileDB)),NDBParts
  close(UNITPAR)
  stop


2500 continue
  write(*,810) SolidMaterials(iMat)%NTable, SolidMaterials(iMat)%Name 
  close(UNITPAR)
  stop


5000 continue
  close(UNITPAR)
  stop
  !----------------------------------------------------------------------------!
444 format(1x,'Unable to open file ',a)
555 format(1x,'End of file or error reading ',a)
560 format(1x,'End of file or error reading ',a,' in ',a)

777 format(1x,'Different number of porous regions in files'//&
       1x,a,' (',i0,') and ',a,' (',i0,')')


888 format(1x,'Bad data for parameter ',a,' = ',i0)
810 format(1x,'Insufficient no. of entries (',i0,' in heat release rate table for material ',a)

contains
  logical function GetLine(buf)
    implicit none
    character(len=*) buf
    integer:: pos

    GetLine = .true.
    do
       read(UNITPAR,'(a)',end=2000,err=2000) buf
       buf = adjustl(trim(buf));  pos = scan(buf,'!%')
       if(pos > 0) then
          if(pos == 1) cycle
          if(pos > 1) buf = buf(1:pos-1)
       endif
       if(len_trim(buf) == 0) cycle
       return
    enddo
2000 continue
    GetLine = .false.
  end function GetLine

  logical function ReadPatches()
    real(kind=8):: XBnd(2) 
    real(kind=8):: YBnd(2)
    real(kind=8):: Void
    real(kind=8):: K, Eps
    real(kind=8):: UPH(NPHASE), VPH(NPHASE)
    real(kind=8):: Temp(NPHASE)
    integer:: Saturated
    integer:: iP, iPH
    character(len=256):: Name

    namelist/PATCH/ Name, XBnd, YBnd, Void, UPH, VPH, K, Eps, Temp, Saturated

    ReadPatches = .true.
    !
    ! Count patches
    !
    NPatches = 0
    rewind UNITPAR
    do
       read(UNITPAR,NML=PATCH,end=7000,err=2000)
       NPatches = NPatches + 1
    enddo
7000 continue
    if(NPatches <= 0) return
    !
    ! Do actual reading
    !
    allocate(Patches(NPatches))
    do iP = 1,NPatches
       ! Initialize by values which are used to detect if a parameter was set by user
       Patches(iP)%XBnd(1) = -1.d30;  Patches(iP)%XBnd(2) = 1.d30
       Patches(iP)%YBnd(1) = -1.d30;  Patches(iP)%YBnd(2) = 1.d30
       Patches(iP)%Void = UNASSIGNED
       Patches(iP)%UPH = UNASSIGNED; Patches(iP)%VPH = UNASSIGNED
       Patches(iP)%K = UNASSIGNED; Patches(iP)%Eps = UNASSIGNED
       Patches(iP)%Temp = UNASSIGNED
    enddo
    rewind UNITPAR
    do iP = 1, NPatches
       XBnd = (/-1.d30, 1.d30/)
       YBnd = (/-1.d30, 1.d30/)
       Void = UNASSIGNED
       K = UNASSIGNED; Eps = UNASSIGNED
       Temp = UNASSIGNED; Saturated = UNASSIGNED_INT
       read(UNITPAR,NML=PATCH,end=2000,err=2000)
       Patches(iP)%Name = Name
       ! Set coordinate ranges
       Patches(iP)%XBnd(1) = max(Patches(iP)%XBnd(1),XBnd(1))
       Patches(iP)%XBnd(2) = min(Patches(iP)%XBnd(2),XBnd(2))
       Patches(iP)%YBnd(1) = max(Patches(iP)%YBnd(1),YBnd(1))
       Patches(iP)%YBnd(2) = min(Patches(iP)%YBnd(2),YBnd(2))
       ! Set values (if any)
       if(Void >= 0.d0) Patches(iP)%Void = Void
       do iPh = 1,NPHASE
          Patches(iP)%UPH(iPh) = max(Patches(iP)%UPH(iPh),UPH(iPh))
          Patches(iP)%VPH(iPh) = max(Patches(iP)%VPH(iPh),VPH(iPh))
       enddo
       if(K > 0.d0) Patches(iP)%K = K
       if(Eps > 0.d0) Patches(iP)%Eps = Eps
       if(minval(Temp) > 0.d0) Patches(iP)%Temp = Temp
       Patches(iP)%Saturated = (Saturated /= 0)
    enddo
    return
2000 continue
    ReadPatches = .false.      
  end function ReadPatches

  subroutine ReadProfiles()
    use Globals
    use Output
    implicit none

    integer:: iP
    real(kind=8):: Coord
    integer:: type
    character(len=256):: File
    character(len=256):: Name
    integer:: Dir

    namelist/PROFILE_WN/ file, Name, type, Coord
    !
    ! Wall normal profiles
    !
    NWallNormalProfiles = 0
    rewind UNITPAR
    do
       read(UNITPAR,NML=PROFILE_WN,end=1000,err=2000)
       NWallNormalProfiles = NWallNormalProfiles+1
    enddo
1000 continue

    if(NWallNormalProfiles <= 0) return
    !
    ! Do actual reading
    !
    allocate(WallNormalProfiles(NWallNormalProfiles))
    rewind UNITPAR
    do iP = 1,NWallNormalProfiles
       Dir = 0; Coord = -1.d30; type = -1
       read(UNITPAR,NML=PROFILE_WN,end=2000,err=2000)

       if(len_trim(File) <= 0) goto 3000
       WallNormalProfiles(iP)%FileName = File
       WallNormalProfiles(iP)%Name = Name
       if(Dir /= WALL_NORMAL_Y_POS .and. &
            Dir /= WALL_NORMAL_Y_NEG) goto 3500
       WallNormalProfiles(iP)%Dir = Dir

       if(type /=  HOR_WALL_PROFILE .and. &
            type /= VERT_WALL_PROFILE) goto 3500
       WallNormalProfiles(iP)%type = type

       if(Coord < -1.d29) goto 3500
       WallNormalProfiles(iP)%Coord = Coord

    enddo
    return
3000 continue
    write(*,666)
    stop
3500 continue
    write(*,777)
    stop
2000 continue
    write(*,555) 
    stop
555 format(1x,'Error reading profiles')
666 format(1x,'Bad data for profile')
777 format(1x,'Missing data for profile')
  end subroutine ReadProfiles

  logical function ReadPhaseMaterials()
    use Globals
    use Material
    implicit none

    integer::  DensModel = UNASSIGNED_INT
    integer::  ViscModel = UNASSIGNED_INT
    integer::  CondModel = UNASSIGNED_INT
    integer::  Phase = UNASSIGNED_INT
    integer:: Standard = UNASSIGNED_INT
    real(kind=8):: Density
    real(kind=8):: Viscosity
    real(kind=8):: Conductivity
    real(kind=8):: MolMass
    real(kind=8):: AlphaP
    real(kind=8):: BetaT
    real(kind=8):: PRef
    real(kind=8):: TRef
    real(kind=8):: Visc0
    real(kind=8):: bVisc
    real(kind=8):: TRefV
    real(kind=8):: CSuth(3)
    real(kind=8):: Cond0
    real(kind=8):: cCond
    real(kind=8):: TRefC
    real(kind=8):: PrGas
    character(len=128):: Name

    namelist/PHASEMAT/ Name, Phase, Standard, &
         Density, DensModel, &
         Viscosity, ViscModel, &
         Conductivity, CondModel,&
         MolMass, AlphaP, BetaT, PRef, TRef,&
         Visc0, bVisc, TRefV, CSuth, PrGas, Cond0, cCond, TRefC

    integer:: iM
    type(MAT_TYPE), pointer:: Mat

    ReadPhaseMaterials = .true.

    rewind UNITPAR
    iM = 0
    do
       read(UNITPAR,NML=PHASEMAT,end=100,err=2000)
       iM = iM+1
    enddo
100 continue
    if(iM <= 0) then
       NMaterials = 0
       return
    endif

    NMaterials = iM
    allocate(Materials(NMaterials))

    rewind UNITPAR
    do iM = 1,NMaterials
       Standard = UNASSIGNED_INT; Phase = UNASSIGNED_INT; Name =""
       Density = UNASSIGNED; Viscosity = UNASSIGNED; MolMass = UNASSIGNED
       Conductivity = UNASSIGNED
       PRef = UNASSIGNED; TRef = UNASSIGNED; AlphaP = UNASSIGNED;
       BetaT = UNASSIGNED;  Visc0 = UNASSIGNED; bVisc = UNASSIGNED;
       TRefV =  UNASSIGNED; CSuth = UNASSIGNED
       PrGas =  UNASSIGNED; Cond0 = UNASSIGNED
       cCond = UNASSIGNED; TRefC =  UNASSIGNED;
       read(UNITPAR,NML=PHASEMAT,end=2000,err=2000)
       if(len_trim(Name) <= 0) goto 2000

       Mat => Materials(iM)
       call INIT_MATERIAL(Name, Mat)
       if(Standard /= Mat%Standard) goto 2000

       if(Standard == 1) then
          if(DensModel == UNASSIGNED_INT) DensModel = Mat%DensModel
          if(ViscModel == UNASSIGNED_INT) ViscModel = Mat%ViscModel
          if(CondModel == UNASSIGNED_INT) CondModel = Mat%CondModel
          if(Phase /= Mat%Phase) goto 2000
       else ! User-defined materials
          if(min(DensModel,ViscModel,CondModel) == UNASSIGNED_INT) goto 2000
          if(Phase == UNASSIGNED_INT) goto 2000
          Mat%Phase = Phase
       endif
       !
       ! Assign the density, viscosity, and conductivity models
       !
       if(Phase == MATERIAL_PHASE_LIQUID) then ! Liquid
          if(DensModel == VARIABLE_DENSITY) then
             if(Standard == 0) then
                if(min(PRef,TRef,AlphaP,BetaT) < IFASSIGNED) goto 2000
                Mat%PRef = PRef*1.d5 ! [bar] -> [Pa]
                Mat%TRef = TRef; Mat%AlphaP= AlphaP; Mat%BetaT = BetaT
             endif
          else
             if(Density < IFASSIGNED) goto 2000
          endif
          if(ViscModel == VARIABLE_VISC) then
             if(Standard == 0) then
                if(min(Visc0,bVisc,TRefV) < IFASSIGNED) goto 2000
                Mat%Visc0 = Visc0; Mat%bVisc = bVisc; Mat%TRefV = TRefV
             endif
          else
             if(Viscosity  < IFASSIGNED) goto 2000
          endif
          if(CondModel == VARIABLE_COND) then
             if(Standard == 0) then
                if(min(Cond0, cCond, TRefC) < IFASSIGNED) goto 2000
                Mat%Cond0 = Cond0
                Mat%cCond = cCond
                Mat%TRefC = TRefC
             endif
          else
             if(Conductivity  < IFASSIGNED) goto 2000
          endif
       else
          if(MolMass < IFASSIGNED) goto 2000
          Mat%MolMass = MolMass
          if(DensModel == VARIABLE_DENSITY) then
             if(Standard == 0) then
                ! Ideal gas law
             endif
          else
             if(Density < IFASSIGNED) goto 2000
          endif
          if(ViscModel == VARIABLE_VISC) then
             if(Standard == 0) then
                if( minval(CSuth) < IFASSIGNED) goto 2000
                Mat%CSuth = CSuth
             endif
          else
             if(Viscosity  < IFASSIGNED) goto 2000
          endif
          if(CondModel == VARIABLE_COND) then
             if(Standard == 0) then
                if(PrGas < IFASSIGNED) goto 2000
                Mat%PrGas = PrGas
             endif
          else
             if(Conductivity  < IFASSIGNED) goto 2000
          endif
       endif
       call SetMatDensityModel(Mat,Model = DensModel,Density = Density)
       call SetMatViscosityModel(Mat,Model = ViscModel,Viscosity = Viscosity)
       call SetMatConductivityModel(Mat,Model = CondModel,Conductivity = Conductivity)
    enddo
    return
2000 continue
    write(*,555) 
    ReadPhaseMaterials = .false.
555 format(1x,'Error reading phase materials')
  end function ReadPhaseMaterials

  logical function ReadSolidMaterials()
    use Globals
    use SOLIDMATERIAL
    implicit none
    character(len=128):: Name
    integer:: HRR_Model
    integer:: type
    real(kind=8):: HRR_Value
    real(kind=8):: Time_HRR(2,MAX_HRR_TABLE_SIZE)    
    real(kind=8):: TimeShutdown
    real(kind=8):: Density
    real(kind=8):: DensityLiq
    real(kind=8):: SpecificHeat
    real(kind=8):: SpecificHeatLiq
    real(kind=8):: Conductivity
    real(kind=8):: ConductivityLiq
    real(kind=8):: TMelt
    real(kind=8):: HMelt

    namelist/SOLIDMAT/ Name, type, Density, DensityLiq, HRR_Model, HRR_Value, Time_HRR, &
         SpecificHeat, Conductivity, SpecificHeatLiq, ConductivityLiq, &
         TimeShutdown, TMelt, HMelt

    integer:: iM, nT, iT

    ReadSolidMaterials = .true.
    rewind UNITPAR
    iM = 0
    do
       read(UNITPAR,NML=SOLIDMAT,end=100,err=2000)
       iM = iM+1
    enddo
100 continue
    if(iM <= 0) then
       NSolidMaterials = 0
       return
    endif

    NSolidMaterials = iM
    allocate(SolidMaterials(NSolidMaterials))
    !
    ! Do actual reading
    !
    rewind UNITPAR
    do iM = 1,NSolidMaterials
       Name = ""
       HRR_Model = -1
       type = -1
       HRR_Value = UNASSIGNED
       Time_HRR =  UNASSIGNED
       Density = UNASSIGNED
       SpecificHeat = UNASSIGNED
       Conductivity = UNASSIGNED
       TimeShutdown = UNASSIGNED
       TMelt = UNASSIGNED
       HMelt = UNASSIGNED
       DensityLiq = UNASSIGNED
       SpecificHeatLiq = UNASSIGNED
       ConductivityLiq = UNASSIGNED

       read(UNITPAR,NML=SOLIDMAT,end=2000,err=2000)
       if(len_trim(Name) <= 0) goto 3010
       SolidMaterials(iM)%Name = Name

       if(type < 0) goto 3005
       SolidMaterials(iM)%type = type


       if(Density <= 0.d0) goto 3020
       SolidMaterials(iM)%Density = Density
       if(DensityLiq < 0.d0) DensityLiq = Density
       SolidMaterials(iM)%DensityLiq = DensityLiq

       if(SpecificHeat <= 0.d0) goto 3030
       SolidMaterials(iM)%SpecificHeat = SpecificHeat
       if(SpecificHeatLiq < 0.d0) SpecificHeatLiq = SpecificHeat
       SolidMaterials(iM)%SpecificHeatLiq = SpecificHeatLiq

       if(Conductivity < 0.d0) goto 3040
       SolidMaterials(iM)%Conductivity = Conductivity
       if(ConductivityLiq < 0.d0) ConductivityLiq = Conductivity
       SolidMaterials(iM)%ConductivityLiq = ConductivityLiq

       select case(type)
       case(SOLID_TYPE_CORIUM)
          if(TMelt <= 0.d0) goto 3042
          SolidMaterials(iM)%TMelt = TMelt
          if(HMelt <= 0.d0) goto 3044
          SolidMaterials(iM)%HMelt = HMelt
       case(SOLID_TYPE_PASSIVE)
          SolidMaterials(iM)%TMelt = -1
          SolidMaterials(iM)%HMelt = -1
       case default
          goto 3006
       end select
       !
       ! Heat release model
       !
       if(HRR_Model <= 0) goto 3080
       SolidMaterials(iM)%HRR_Model = HRR_Model
       select case(HRR_Model)
       case(CONSTANT_HRR)
          if(HRR_Value < IFASSIGNED) goto 3050
          SolidMaterials(iM)%HRRPUM = HRR_Value
       case(TRANSIENT_TABLE_HRR)
          nT = 0
          do iT = 1,MAX_HRR_TABLE_SIZE
             if(minval(Time_HRR(:,iT)) < IFASSIGNED) exit
             nT = nT+1
          enddo
          if(nT <= 0) goto 3060
          SolidMaterials(iM)%NTable = nT
          allocate(SolidMaterials(iM)%Table_Time(nT))
          allocate(SolidMaterials(iM)%Table_HRRPUM(nT))
          SolidMaterials(iM)%Table_Time(1:nT) = Time_HRR(1,1:nT)
          SolidMaterials(iM)%Table_HRRPUM(1:nT) = Time_HRR(2,1:nT)
       case(TRANSIENT_CORIUM_HRR)
          if(TimeShutdown < IFASSIGNED) goto 3070
          SolidMaterials(iM)%TimeShutdown = TimeShutdown*3600.d0 ! [h] -> [sec]
       case default
          goto 3090
       end select
    enddo
    return
3000 continue
    write(*,777) iM
    stop
3005 continue
    write(*,805) iM
    stop
3006 continue
    write(*,806) iM
    stop
3010 continue
    write(*,810) iM
    stop
3020 continue
    write(*,820) iM,trim(adjustl(SolidMaterials(iM)%Name))
    stop
3030 continue
    write(*,830) iM,trim(adjustl(SolidMaterials(iM)%Name))
    stop
3040 continue
    write(*,840) iM,trim(adjustl(SolidMaterials(iM)%Name))
    stop
3042 continue
    write(*,842) iM,trim(adjustl(SolidMaterials(iM)%Name))
    stop
3044 continue
    write(*,844) iM,trim(adjustl(SolidMaterials(iM)%Name))
    stop
3050 continue
    write(*,850) iM,trim(adjustl(SolidMaterials(iM)%Name)), &
         HRR_Model,trim(adjustl(KnownHeatReleaseRateModelNames(HRR_Model)))
    stop
3060 continue
    write(*,860) iM,trim(adjustl(SolidMaterials(iM)%Name)), &
         HRR_Model,trim(adjustl(KnownHeatReleaseRateModelNames(HRR_Model)))
    stop
3070 continue
    write(*,870) iM,trim(adjustl(SolidMaterials(iM)%Name)), &
         HRR_Model,trim(adjustl(KnownHeatReleaseRateModelNames(HRR_Model)))
    stop
3080 continue
    write(*,880) iM,trim(adjustl(SolidMaterials(iM)%Name))
    stop
3090 continue
    write(*,890) HRR_Model, iM,trim(adjustl(SolidMaterials(iM)%Name))
    stop
2000 continue
    write(*,555) 
    ReadSolidMaterials = .false.
    return
555 format(1x,'Error reading solid materials')
777 format(1x,'Missing data for solid material #',i0)
810 format(1x,'Missing Name for solid material #',i0)
805 format(1x,'Missing Type for solid material #',i0)
806 format(1x,'Bad Type for solid material #',i0)
820 format(1x,'Missing Density for solid material #',i0,' (',a,')')
830 format(1x,'Missing SpecificHeat for solid material #',i0,' (',a,')')
840 format(1x,'Missing Conductivity for solid material #',i0,' (',a,')')
842 format(1x,'Missing TMelt for solid material #',i0,' (',a,')')
844 format(1x,'Missing HMelt for solid material #',i0,' (',a,')')
850 format(1x,'Missing HRR_Value for solid material #',i0,' (',a,'), HRR_Model #',i0,' (',a,')')
860 format(1x,'Missing Time_HRR table for solid material #',i0,' (',a,'), HRR_Model #',i0,' (',a,')')
870 format(1x,'Missing TimeShutdown for solid material #',i0,' (',a,'), HRR_Model #',i0,' (',a,')')!
880 format(1x,'Missing HRR_Model for solid material #',i0,' (',a,')')
890 format(1x,'Bad HRR_Model #',i0,' for solid material #',i0,' (',a,')')
  end function ReadSolidMaterials

  logical function ReadEnergyModel()
    use Globals
    use Heat
    use Material
    implicit none

    integer:: EnergyModel,Evaporation,Saturated,PorousEnergyModel
    integer:: iM
    type(MAT_TYPE), pointer:: Mat

    namelist/HEATMODEL/ EnergyModel,Evaporation,Saturated,PorousEnergyModel

    TRACE_LINE = "ReadEnergyModel"
    ReadEnergyModel = .true.
    EnergyModel = UNASSIGNED_INT
    Evaporation = UNASSIGNED_INT
    Saturated  = UNASSIGNED_INT
    PorousEnergyModel = UNASSIGNED_INT

    rewind(UNITPAR)
    read(UNITPAR,NML=HEATMODEL,end=100,err = 2000)
100 continue
    IS_ENERGY = EnergyModel > 0
    IS_EVAPORATION = Evaporation > 0
    IS_SATURATED = Saturated > 0
    IS_POROUS_ENERGY = PorousEnergyModel > 0

    if(IS_ENERGY) then
       if(IS_SATURATED) goto 3000
    else
       if(.not.IS_SATURATED .and. IS_EVAPORATION) goto 3500
    endif
    call SET_EQUATION_SET
    if(IS_SATURATED) then
       do iM = 1,NMaterials
          Mat => Materials(iM)
          call SET_EOS_SATURATED(Mat,IS_SATURATED)
       enddo
    endif
    return
3000 continue
    write(*,300)
    goto 2000
3500 continue
    write(*,350)
    goto 2000
2000 continue
    ReadEnergyModel = .false.
300 format("Saturated conditions cannot be used with EnergyModel > 0")
350 format("Evaporation model cannot be used for EnergyModel=0 and Saturated=0")
  end function ReadEnergyModel

  logical function ReadTurbulenceModel()
    use Globals
    use Multiphase
    use Turb
    implicit none

    integer:: Model, TurbDisp, UseViscEff
    integer:: RampNearPorous
    real(kind=8):: ViscEff
    namelist/TURBMODEL/ Model, TurbDisp, UseViscEff, ViscEff, RampNearPorous,&
         TinyKEps_K, TinyKEps_KG, TinyKEps_Eps, OmegaTurb, InterpTypeTurb, FLTypeTurb

    TRACE_LINE = "ReadTurbulenceModel"
    ReadTurbulenceModel = .true.

    Model = 0
    TurbDisp = 1
    ViscEff = UNASSIGNED
    UseViscEff = UNASSIGNED_INT
    RampNearPorous = UNASSIGNED_INT
    TinyKEps_KG = UNASSIGNED
    InterpTypeTurb = UNASSIGNED_INT
    FLTypeTurb = UNASSIGNED_INT

    rewind(UNITPAR)
    read(UNITPAR,NML=TURBMODEL,end=100,err = 2000)
100 continue
    if(RampNearPorous > UNASSIGNED_INT) IS_RAMP_NEAR_POROUS = RampNearPorous > 0
    call SetupTurbulenceModel(Model)

    IS_PHASE_DISPERSION = TurbDisp /= 0
    IS_EFFECTIVE_VISCOSITY = UseViscEff /= UNASSIGNED_INT .and. UseViscEff /= 0
    if(IS_EFFECTIVE_VISCOSITY) then
       if(ViscEff < IFASSIGNED) goto 2000
       Visc_Eff = ViscEff
    endif

    if(TinyKEps_KG <= 0.d0) TinyKEps_KG = TinyKEps_K
    !
    ! Numerical scheme
    !
    if(HighOrderSpace > 0) then
       if(InterpTypeTurb == UNASSIGNED_INT) InterpTypeTurb = INTERP_TYPE_TURB_DEFAULT
       if(FLTypeTurb == UNASSIGNED_INT) FLTypeTurb = FL_TYPE_TURB_DEFAULT
    else
       InterpTypeTurb = 0
       FLTypeTurb = 0
    endif
    return
2000 continue
    ReadTurbulenceModel = .false.
  end function ReadTurbulenceModel

  logical function ReadMultiphaseModel()
    use Globals
    use Multiphase
    implicit none

!    integer:: SinglePhase = 0
    character(len=256):: Phase1Material
    character(len=256):: Phase2Material
    namelist /MPHASEMODEL/ SinglePhase, Phase1Material, Phase2Material

    integer:: iMat

    TRACE_LINE = "ReadMultiphaseModel"
    ReadMultiphaseModel = .true.

    rewind(UNITPAR)
    Phase1Material = ""; Phase2Material = ""
    read(UNITPAR,NML=MPHASEMODEL,end=100,err=2000)
100 continue
    if(len_trim(Phase1Material) <= 0) goto 2000
    if(len_trim(Phase2Material) <= 0) goto 2000
    IS_SINGLE_PHASE = SinglePhase /= 0
    NPhaseActive = NPHASE
    if(IS_SINGLE_PHASE) NPhaseActive = 1
    PhaseMatName = ""

    iMat = FindMaterial(Phase1Material)
    if(iMat < 0) goto 2000
    MatPhase(IndLiq) = Materials(iMat)
    PhaseMatName(IndLiq) = Phase1Material

    if(.not.IS_SINGLE_PHASE) then
       iMat = FindMaterial(Phase2Material)
       if(iMat < 0) goto 2000
       MatPhase(IndGas) = Materials(iMat)
       PhaseMatName(IndGas) = Phase2Material
    endif
    return
2000 continue
    ReadMultiphaseModel = .false.
  end function ReadMultiphaseModel

  logical function ReadCorrelations()
    use Globals
    use Correlations
    implicit none

    integer:: PorousDragModel
    namelist/CORREL/ PorousDragModel

    TRACE_LINE = "ReadCorrelations"
    ReadCorrelations = .true.

    rewind(UNITPAR)
    read(UNITPAR,NML=CORREL,end=100,err=2000)
100 continue
    PorousModel = min(max(PorousDragModel,0),N_MOD_POROUS)
    return
2000 continue
    ReadCorrelations = .false.
  end function ReadCorrelations

  logical function ReadProblem()
    use Globals
    implicit none

    namelist/DECOSIM/  FileGEO, FileGRD, PSystem000, TSystem000

    TRACE_LINE = "ReadProblem"

    ReadProblem = .true.
    rewind(UNITPAR)
    FileGEO = ""; FileGRD = ""
    read(UNITPAR, NML=DECOSIM, end = 2100, ERR = 2000)
    PSystem000 = PSystem000*1.d5  ! Convert bar->Pa 
    call SET_SYSTEM_P_T(PSystem000,TSystem000)
    FileGEO = trim(adjustl(FileGEO))
    if(len_trim(FileGEO) <= 0) goto 2105
    FileDB = FileGEO
    call SetSymmetry
    call ReadDebrisBedGeometry

    FileGRD = trim(adjustl(FileGRD))
    if(len_trim(FileGRD) <= 0) goto 2106
    return
2100 continue
    write(*,666) trim(adjustl(FilePAR))
    close(UNITPAR)
    stop
2105 continue
    write(*,680) trim(adjustl(FilePAR))
    close(UNITPAR)
    stop
2106 continue
    write(*,681) trim(adjustl(FilePAR))
    close(UNITPAR)
    stop
2000 continue
    ReadProblem = .false.
666 format(1x,'Input file ',a,' is not for program DECOSIM')
680 format(1x,'Missing geometry file (FileGEO) in ',a)
681 format(1x,'Missing grid file (FileGRD) in ',a)
  end function ReadProblem

  logical function ReadOptions()
    use Globals
    implicit none

    integer:: UseGravity
    namelist/OPTIONS/  Gravity, UseGravity, SEED, ForceDDisp

    TRACE_LINE = "ReadOptions"
    ReadOptions = .true.
    rewind(UNITPAR)
    SEED = -1
    UseGravity = 1
    ForceDDisp = UNASSIGNED
    read(UNITPAR,NML=OPTIONS, end = 7000, err= 2000)
7000 continue
    IS_GRAVITY = UseGravity /= 0
    if(maxval(SEED) > 0) then
       ! Check if enough values were supplied
       call random_seed(size=size_seed)
       if(minval(SEED(1:size_seed)) < 0) then
          write(*,890) size_seed
       endif
    endif
    IS_FORCE_DDISP = ForceDDisp > 0.d0
    return
2000 continue
    ReadOptions = .false.
890 format(1x,'Insufficient data for random generator seed: '/&
         1x,i0,' integer numbers must be given on SEED= line in OPTIONS namelist'/&
         1x,'Using system-dependent seed')
  end function ReadOptions

  logical function ReadGapTooth()
    use Globals
    implicit none

    integer:: Active, Master
    namelist/GAPTOOTH/ Active, Master

    TRACE_LINE = "ReadGapTooth"
    ReadGapTooth = .true.

    Active = 0
    Master = UNASSIGNED_INT ! Unused by decosim
    read(UNITPAR,NML=GAPTOOTH, end = 7000, err= 2000)
7000 continue
    IS_GAP_TOOTH = Active > 0
    return
2000 continue
    ReadGapTooth = .false.
  end function ReadGapTooth

  logical function ReadScales()
    use Globals
    implicit none
    integer:: UseDimensional


    namelist/SCALES/ UseDimensional, XScale, VScale, RoScale,TempScale

    TRACE_LINE = "ReadScales"
    ReadScales = .true.
    UseDimensional = -1
    rewind(UNITPAR)
    read(UNITPAR,NML=SCALES, end = 7010, err= 2000)
7010 continue
    IS_DIMENSIONAL = UseDimensional /= 0
    return
2000 continue
    ReadScales = .false.
  end function ReadScales

  logical function ReadTimeStepping()
    use Globals
    implicit none

    integer:: Transient
    integer:: StepUpdateFreq

    integer,parameter:: TABLELEN = 1000
    real(kind=8):: MaxTimeStep_Table(2,TABLELEN)
    real(kind=8):: CFL_Max_Table(2,TABLELEN)
    integer:: i,N

    namelist/STEPPING/ Transient, NTotalSteps, TimeMax, &
         FluxUpdateTime, &
         MaxTimeStep, MaxTimeStep_Table,&
         MinTimeStep, TimeStepInit,&
         AdaptiveTimeStep, &
         CFL_Max, CFL_Max_Table, &
         CFL_Frac, StepUpdateFreq


    TRACE_LINE = "ReadTimeStepping"
    ReadTimeStepping = .true.

    FluxUpdateTime = -1
    StepUpdateFreq = UNASSIGNED_INT
    MaxTimeStep_Table = UNASSIGNED
    CFL_Max_Table = UNASSIGNED
    rewind(UNITPAR)
    read(UNITPAR,NML=STEPPING,end=7001,err=2000)
    if(StepUpdateFreq > 0) NStepsAfterStepChange = StepUpdateFreq

    ! Setup registers for time step and CFL
    N = 0
    do i = 1,TABLELEN
       if(MaxTimeStep_Table(1,i) > IFASSIGNED .and. MaxTimeStep_Table(2,i) > IFASSIGNED) then
          N = N+1
          cycle
       endif
       exit
    enddo
    if(N > 0) then
       call SetupRegister("MaxTimeStep",N,t_v = MaxTimeStep_Table)
    else
       call SetupRegister("MaxTimeStep",1,Val=MaxTimeStep)
    endif
    MaxTimeStep_Reg => GetRegisterByName("MaxTimeStep")

    N = 0
    do i = 1,TABLELEN
       if(CFL_Max_Table(1,i) > IFASSIGNED .and. CFL_Max_Table(2,i) > IFASSIGNED) then
          N = N+1
          cycle
       endif
       exit
    enddo
    if(N > 0) then
       call SetupRegister("CFL_Max",N,t_v = CFL_Max_Table)
    else
       call SetupRegister("CFL_Max",1,Val=CFL_Max)
    endif
    CFL_Max_Reg => GetRegisterByName("CFL_Max")

7001 continue
    IS_TRANSIENT = Transient /= 0
    IS_ADAPTIVE_TIME_STEP = AdaptiveTimeStep /= 0
    CFL = CFL_Max
    IS_DUMP_MODE = NTotalSteps == 0
    return
2000 continue
    ReadTimeStepping = .false.
  end function ReadTimeStepping

  logical function ReadNumerics()
    use Globals
    implicit none

    namelist/NUMERICS/ KindTimeScheme, HighOrderSpace, InterpType, FLType

    TRACE_LINE = "ReadNumerics"
    ReadNumerics = .true.

    rewind(UNITPAR)
    HighOrderSpace = UNASSIGNED_INT
    InterpType = UNASSIGNED_INT; FLType = UNASSIGNED_INT
    read(UNITPAR,NML=NUMERICS, end = 7020, err= 2000)
7020 continue
    if(HighOrderSpace <= UNASSIGNED_INT) goto 2000
    KindTimeScheme = min(KindTimeScheme,KIND_TIME_2ND_ORDER)
    KindTimeScheme = max(KindTimeScheme,KIND_TIME_1ST_ORDER)
    !
    ! Numerical scheme
    !
    if(HighOrderSpace > 0) then
       if(InterpType == UNASSIGNED_INT) InterpType = 3
       if(FLType == UNASSIGNED_INT) FLType = 4
    else
       InterpType = 0
       FLType = 0
    endif
    if(KindTimeScheme == KIND_TIME_2ND_ORDER) then
       Theta = 0.5d0
    else
       Theta = 1.d0
    endif
    return
2000 continue
    ReadNumerics = .false.
  end function ReadNumerics

  logical function ReadSchemeFlow()
    use Globals
    use Level
    implicit none

    integer:: PoolRegime

    namelist/SCHEME_FLOW/ PoolRegime

    TRACE_LINE = "ReadSchemeFlow"
    ReadSchemeFlow = .true.
    rewind(UNITPAR)
    PoolRegime = UNASSIGNED_INT

    read(UNITPAR,NML=SCHEME_FLOW, end = 100, err= 2000)
100 continue
    if(PoolRegime /= UNASSIGNED_INT) then
       IS_POOL = PoolRegime /= 0
       IS_POOL_MODIFY_EQ = PoolRegime > 1
    endif
    return
2000 continue
    ReadSchemeFlow = .false.
  end function ReadSchemeFlow

  logical function ReadSchemeAlpha()
    use Globals
    use PVolFrac
    implicit none

    integer:: UpwindDragdA, UpwindDragdAD, DragdA, DragdAD, SummupPhases

    namelist/SCHEME_ALPHA/ UpwindDragdA, UpwindDragdAD, DragdA, DragdAD, SummupPhases

    TRACE_LINE = "ReadSchemeAlpha"
    ReadSchemeAlpha = .true.
    rewind(UNITPAR)
    UpwindDragdA = UNASSIGNED_INT
    UpwindDragdAD = UNASSIGNED_INT
    DragdA = UNASSIGNED_INT
    DragdAD = UNASSIGNED_INT
    SummupPhases = UNASSIGNED_INT
    read(UNITPAR,NML=SCHEME_ALPHA, end = 100, err= 2000)
100 continue
    if(UpwindDragdA /= UNASSIGNED_INT) IS_UPWIND_DRAG_dA = UpwindDragdA /= 0
    if(UpwindDragdAD /= UNASSIGNED_INT) IS_UPWIND_DRAG_dAD = UpwindDragdAD /= 0
    if(DragdA /= UNASSIGNED_INT) IS_DRAG_dA = DragdA /= 0
    if(DragdAD /= UNASSIGNED_INT) IS_DRAG_dAD = DragdAD /= 0
    if(SummupPhases /= UNASSIGNED_INT) IS_SUMMUP_PHASES = SummupPhases /= 0
    return
2000 continue
    ReadSchemeAlpha = .false.
  end function ReadSchemeAlpha

  logical function ReadIterations()
    use Globals
    use PVolFrac
    use Heat
    implicit none

    namelist/ITERS/ NIterPreRun, MaxItersNonLinear,&
         EpsItersVelocity, EpsItersVelocityUP,&
         EpsItersPressure, EpsItersPressureUP,&
         EpsItersAlpha, EpsItersAlphaUP,&
         EpsItersEnergy, EpsItersEnergyUP,&
         MaxItersBSGSOuter, EpsItersStep,&
         TolUpLowFactor,&
         OmegaP, OmegaAlpha, OmegaVelocity, OmegaEnergy,&
         dAmaxRelStep, dPmaxRelStep,&
         dUMaxRelStep, dTmaxRelStep

    TRACE_LINE = "ReadIterations" 
    ReadIterations = .true.
    rewind(UNITPAR)
    EpsItersVelocityUP = UNASSIGNED
    EpsItersPressureUP = UNASSIGNED
    EpsItersAlphaUP = UNASSIGNED
    EpsItersEnergyUP = UNASSIGNED    
    read(UNITPAR,NML=ITERS, end = 100, err= 2000)
100 continue
    if(EpsItersVelocityUP < IFASSIGNED) EpsItersVelocityUP = EpsItersVelocity*TolUpLowFactor
    if(EpsItersPressureUP < IFASSIGNED) EpsItersPressureUP = EpsItersPressure*TolUpLowFactor
    if(EpsItersAlphaUP < IFASSIGNED) EpsItersAlphaUP = EpsItersAlpha*TolUpLowFactor
    if(EpsItersEnergyUP < IFASSIGNED) EpsItersEnergyUP = EpsItersEnergy*TolUpLowFactor
    return
2000 continue
    ReadIterations = .false.
  end function ReadIterations

  logical function ReadOutput()
    use Globals
    use Output
    implicit none

    namelist/OUT/ NOutScreen, NOutFields, NOutRestart, NOutInfo, NOutTrans, MsgLevel,&
         DTOutScreen, DTOutFields, DTOutRestart, DTOutInfo, DTOutTrans,&
         OutputDirPrefix

    TRACE_LINE = "ReadOutput" 
    ReadOutput = .true.
    rewind(UNITPAR)
    OutputDirPrefix=""
    DTOutScreen = -1.d0
    DTOutFields = -1.d0
    DTOutRestart = -1.d0
    DTOutInfo = -1.d0
    DTOutTrans = 1.d-1

    read(UNITPAR,NML=OUT, end = 100, err= 2000)
100 continue
    return
2000 continue
    ReadOutput = .false.
  end function ReadOutput

  logical function ReadTransducers()
    use Globals
    use Output
    implicit none
    character(len=256):: Name
    character(len=256):: Label
    real(kind=8):: Coord(NDIMS)
    integer:: Ind(NDIMS)
    integer:: iTR
    integer:: Storage

    namelist/TRANSDUCER/ Name, Label, Coord, Ind, Storage


    TRACE_LINE = "ReadTransducers"
    ReadTransducers = .true.
    if(.not.IS_WRITE_TRANS) return 
    NTrans = 0
    rewind(UNITPAR)
    do
       read(UNITPAR,NML=TRANSDUCER, end = 100, err= 2000)
       NTrans = NTrans+1
    enddo
100 continue
    if(NTrans > 0) then
       !
       ! Do actual reading
       !
       allocate(Trans(NTrans))
       rewind(UNITPAR)
       do iTR = 1, NTrans
          Name = ""; Label = ""; Coord = UNASSIGNED
          Ind = -1
          Storage = -1
          read(UNITPAR,NML=TRANSDUCER, end = 2000, err= 2000)
          if(len_trim(Name) <= 0 .or. len_trim(Label) <= 0) goto 2000
          !    if(minval(Coord) < IFASSIGNED .and. minval(Ind) <= 0) goto 2000
          Trans(iTR)%Name = Name
          Trans(iTR)%Label = Label
          Trans(iTR)%Coord = Coord
          Trans(iTR)%Ind = Ind
          Trans(ITR)%Storage = Storage
       enddo
       Trans(1:NTrans)%Active = .false.   
    endif
    IS_WRITE_TRANS = IS_WRITE_TRANS .and. NTrans > 0
    return
2000 continue
    ReadTransducers = .false.
    IS_WRITE_TRANS = .false.
  end function ReadTransducers

  logical function ReadRestartControls()
    use Globals
    implicit none
    integer:: StartNewTooth, StartFrom

    namelist /START/ StartNewTooth, GapTime, NREAD, NWRITE, &
         FileRST, FileSAVE, FileRSTGLOBAL, FileGAPTH, StartFrom

    TRACE_LINE = "ReadRestartControls"
    ReadRestartControls = .true.
    rewind(UNITPAR)
    StartNewTooth = 0
    StartFrom = 0
    GapTime = UNASSIGNED
    FileRSTGLOBAL = ""
    NREAD = 0; NWRITE = 0
    FileRST = ""; FileSAVE = ""
    FileGAPTH = ""
    read(UNITPAR,NML=START, end = 100, err= 2000)
100 continue
    IS_MAKE_GLOBAL_STEP = StartNewTooth /= 0
    IS_START_FROM_PREV_GAP = StartFrom /= 0
    !
    ! Check consistency of flags
    !
    if(IS_GAP_TOOTH) then
       if(IS_MAKE_GLOBAL_STEP) then
          if(NREAD /= 0) then
             print *,'NREAD cannot be active if StartNewTooth = 1'
             stop
          endif
          if(IS_START_FROM_PREV_GAP) then
             if(len_trim(FileRSTGLOBAL) <= 0) then
                write(*,200) "FileRSTGLOBAL","&START"
                goto 2000
             endif
             FileRST = FileRSTGLOBAL         
             if(GapTime < IFASSIGNED) then
                write(*,200) "GapTime","&START"
                goto 2000
             endif
          else
             GapTime = 0.d0
          endif
          GapTime = GapTime*60.d0 ! [min] -> [sec]
       else
          if(NREAD == 0) then
             print *,'NREAD must be active if StartNewTooth = 0'
             stop
          endif
       endif
    else
       if(IS_MAKE_GLOBAL_STEP) then
          print *,'StartNewTooth cannot be active if Active = 0 in GAPTOOTH'
          stop
       endif
       if(IS_START_FROM_PREV_GAP) then
          print *,'StartFrom cannot be active if Active = 0 in GAPTOOTH'
          stop
       endif
    endif
    if(NREAD > 0 .and. len_trim(FileRST) <= 0) then
       write(*,200) "FileRST","&START"
       goto 2000
    endif
    if(NWRITE > 0 .and. len_trim(FileSAVE) <= 0) then
       write(*,200) "FileSAVE","&START"
       goto 2000
    endif
    if(len_trim(FileGAPTH) <= 0) FileGAPTH = "GapTooth.dat"
    IS_FIRST_GLOBAL_STEP = IS_GAP_TOOTH .and. IS_MAKE_GLOBAL_STEP .and. &
         (.not. IS_START_FROM_PREV_GAP)
    if(IS_FIRST_GLOBAL_STEP) GlobalTime = 0.d0         
    return
200 format(1x,'Missing parameter ',a,' in namelist ',a)
2000 continue
    ReadRestartControls = .false.
  end function ReadRestartControls

  logical function ReadPorousRegions()
    use InOut
    use DebrisBed
    implicit none

    integer:: NDB
    logical:: found
    character(len=256) Name
    real(kind=8):: Porosity, DDisp, PorInternal, PermFactor, T0, RhoC
    integer:: isT0, isRhoC
    integer:: SolidMat, HeatModel
    logical, allocatable:: RegFound(:)
    logical:: set

    namelist/POROUSREG/ Name, Porosity, DDisp, PorInternal, PermFactor, SolidMat,&
         isT0, T0, isRhoC, RhoC, HeatModel



    if(NDBParts <= 0) then
       ReadPorousRegions = .true.
       return
    endif

    TRACE_LINE = "ReadPorousRegions"
    allocate(RegFound(NDBParts))
    RegFound = .false.
    ReadPorousRegions = .false.
    !
    ! Read Namelist to get debris bed parameters
    !
    rewind(UNITPAR)
    REG_LOOP: do
       Name=""
       isT0 = 0
       T0 = -1.D0
       isRhoC = 0
       RhoC = -1.D0
       HeatModel = 0
       read(UNITPAR,NML=POROUSREG,end=1000,err=2000)
       !
       ! Find the porous region by name
       !
       found = .false.
       PART_LOOP: do NDB = 1,NDBParts
          if(trim(adjustl(name)) == trim(adjustl(DBPart(NDB)%Name))) then
             found = .true.
             RegFound(NDB) = .true.
             exit PART_LOOP
          endif
       enddo PART_LOOP
       if(.not.found) cycle REG_LOOP

       DBPart(NDB)%Porosity = Porosity
       DBPart(NDB)%PorInternal = PorInternal*1.d-2 ! [%] -> [-]
       DBPart(NDB)%Mat = SolidMat
       DBPart(NDB)%HeatModel = HeatModel
       if(DBPart(NDB)%type == POR_TYPE_NORMAL) then
          DBPart(NDB)%DDisp = DDisp*1.d-3 ! [mm] -> [m]       
          DBPart(NDB)%PermFactor = PermFactor       
       else
          DBPart(NDB)%DDisp = -1.D0       
          DBPart(NDB)%PermFactor = -1.D0       
       endif
       DBPart(NDB)%isT0 = isT0
       if(isT0 > 0) then
          if(T0 < 0.D0) then
             write(*,900) trim(adjustl(DBPart(NDB)%Name))
             goto 2000
          endif
          DBPart(NDB)%T0 = T0
          isRhoC = 0  ! If temperature is given, heat capacity of structures is not taken into account
       endif
       DBPart(NDB)%isRhoC = isRhoC
       if(isT0 <= 0) then            
          if(isRhoC > 0) then
             if(RhoC < 0.D0) then
                write(*,901) trim(adjustl(DBPart(NDB)%Name))
                goto 2000
             endif
             DBPart(NDB)%RhoC = RhoC
          endif
       endif

    enddo REG_LOOP
1000 continue

    set = .true.
    do NDB = 1,NDBParts
       set = set .and. RegFound(NDB)
    enddo
    ReadPorousRegions =  set
    return
2000 continue
    ReadPorousRegions = .false.

900 format(1x,'Initial temperature not provided for porous region ',a)
901 format(1x,'Volumetric heat capacity RhoC not provided for porous region ',a)
  end function ReadPorousRegions

  logical function ReadOutputFiles()
    use Globals
    use Output
    implicit none

    integer:: WriteDBShape, WriteMassBal, WriteMinMax, WriteTrans
    namelist/OUTFILES/ WriteDBShape, FileDBShape,&
         WriteMassBal, FileMassBal,&
         WriteMinMax, FileMinMax,&
         WriteTrans, FileTRANS,&
         FileINFO,&
         tecPrefix, tecSuffix


    TRACE_LINE = "ReadOutputFiles"
    ReadOutputFiles = .true.
    rewind(UNITPAR)
    WriteDBShape = 0; FileDBShape = ""
    WriteMassBal = 0; FileMassBal = ""
    WriteMinMax = 0; FileMinMax = ""
    WriteTrans = 0; FileTRANS = ""
    FileINFO = "DECOSIM.info"
    read(UNITPAR,NML=OUTFILES, end = 100, err= 2000)
    IS_WRITE_DB_SHAPE = WriteDBShape /= 0 .and. len_trim(FileDBShape) > 0
    IS_WRITE_MASS_BALANCE = WriteMassBal /= 0 .and. len_trim(FileMassBal) > 0
    IS_WRITE_MINMAX = WriteMinMax /= 0 .and. len_trim(FileMinMax) > 0
    IS_WRITE_TRANS = WriteTrans /= 0 .and. len_trim(FileTRANS) > 0
    if(len_trim(tecPrefix) <= 0) tecPrefix = "fn"
    if(len_trim(tecSuffix) <= 0) tecSuffix = "tec"
100 continue
    return
2000 continue
    ReadOutputFiles = .false. 
  end function ReadOutputFiles

  logical function ReadCoolability()
    use Globals
    use Output
    use Heat
    implicit none

    integer:: IsCoolability                           ! flag for coolability
    integer:: IsManualHRR
    integer:: i
    character(len=256):: CoriumName
    namelist/COOLABILITY/ IsCoolability, FileCheckAV, NCheckAV, DTCheckAV, AGasMinHeat, &
         CoriumName, HRRInitial, HRRInitialWet, HRRStep, FunctionThreshold, AGasDryout, AGasRewet,&
         IsManualHRR

    TRACE_LINE = "ReadCoolability"
    ReadCoolability = .true.
    IsManualHRR = 0
    rewind(UNITPAR)
    HRRInitialWet = -1.d0
    read(UNITPAR,NML=COOLABILITY, end = 200, err= 2000)
    IS_COOLABILITY = IsCoolability > 0 .and. len_trim(FileCheckAV) > 0
    if(.not.IS_COOLABILITY) return

    IS_MANUAL_HRR_ON_RESTART = IsManualHRR /= 0
    CoriumiMat = -1
    do i=1,NSolidMaterials
       if(SolidMaterials(i)%Name == CoriumName) then
          CoriumiMat = i
          exit
       endif
    enddo
    if(CoriumiMat<0) goto 2000
    SolidMaterials(CoriumiMat)%HRR_Model = CONSTANT_HRR
    HRRCurrent = HRRInitial
    SolidMaterials(CoriumiMat)%HRRPUM = HRRCurrent
    if(HRRInitialWet >= 0.d0) then
       HRRWet = HRRInitialWet
    endif
100 continue
    return
200 continue
    IS_COOLABILITY = .false.
    return
2000 continue
    ReadCoolability = .false.
  end function ReadCoolability

  logical function ReadParticles()
    use Globals
    use Particle
    use DebrisBed
    implicit none

    integer:: iFeed, iSp, iDB, N, i

    character(len=256):: Name
    character(len=256):: FeederFile, CatcherFile, CatcherCDFFile, Prefix
    integer:: Distr, SolidMat, NFeedPart, NFeedCell, WriteFeeder, WriteCatcher
    integer:: NCatchCell, WriteCatcherCDF, NTrace
    real(kind=8):: RMin, RMax, RadGauss, RCenterGauss, MFeed, &
         TimeStart, TimeFeed, ZFeed, VFeed, CatcherCDFTime, &
         CatcherCDFdT, TFeed

    integer,parameter:: TABLELEN = 1000
    real(kind=8):: MFeed_Table(2,TABLELEN)

    namelist/FEEDS/ Name, RMin, RMax, Distr, RadGauss, RCenterGauss, MFeed, TimeStart, TimeFeed, &
         MFeed_Table, &
         ZFeed, VFeed, TFeed, &
         SolidMat, NFeedPart, NFeedCell, WriteFeeder, WriteCatcher, FeederFile, CatcherFile,&
         NCatchCell, WriteCatcherCDF, CatcherCDFFile, CatcherCDFTime, CatcherCDFdT,&
         NTrace, Prefix

    real(kind=8):: D, PorInternal, MassFrac
    integer:: NSort, WritePath, ForceSediment
    character(len=256):: FeedName
    namelist/PART/ FeedName, D, PorInternal, MassFrac, WritePath, ForceSediment

    integer:: RandomWalk
    integer:: Avalanche
    namelist/SCHEME_PART/ TStepPart, TimePartMax, MaxStepWritePart, RandomWalk,&
         Avalanche, AngleOfRepose

    TRACE_LINE = "ReadParticles"
    ReadParticles = .true.
    !
    ! Count feeders
    !
    NFeeders = 0
    rewind(UNITPAR)
    do
       read(UNITPAR,NML=FEEDS,end=100,err=2000)
       NFeeders = NFeeders + 1
    enddo
100 continue
    !    if(NFeeders <= 0) return
    rewind(UNITPAR)

    if(IS_DEBRIS_BED) then
       allocate(Feeders(NFeeders+NDBParts)) ! Add fictitious feeder, just to setup DebrisBedTop
       NTotalFeeders = NFeeders+NDBParts
    else
       allocate(Feeders(NFeeders))
       NTotalFeeders = NFeeders
    endif
    FEEDER_LOOP: do iFeed = 1,NFeeders
       Feed => Feeders(iFeed)
       FeederFile = ""; CatcherFile=""
       ZFeed = UNASSIGNED
       VFeed = UNASSIGNED
       TFeed = UNASSIGNED
       TimeStart = UNASSIGNED
       MFeed = UNASSIGNED
       MFeed_Table = UNASSIGNED
       read(UNITPAR,NML=FEEDS,end=2000,err=2000)
       Feed%ACTIVE = .true.
       Feed%Name = Name
       Feed%XMin = RMin
       Feed%XMax = RMax
       if(VFeed > IFASSIGNED) then
          Feed%VFeed = VFeed
       else
          Feed%VFeed = 0.d0
       endif
       Feed%ZFeed = ZFeed
       Feed%TFeed = TFeed
       Feed%DISTR = Distr
       Feed%RadGauss = RadGauss
       Feed%RCenterGauss = RCenterGauss
       Feed%MFeed = MFeed       ![tonnes]
       if(TimeStart > IFASSIGNED) then
          Feed%TimeStart = TimeStart ! [hours]
       else
          Feed%TimeStart = 0.d0
       endif
       Feed%TimeFeed = TimeFeed ! [hours]       
       N = 0
       do i = 1,TABLELEN
          if(MFeed_Table(1,i) > IFASSIGNED .and. MFeed_Table(2,i) > IFASSIGNED) then
             N = N+1
             cycle
          endif
          exit
       enddo
       if(N <= 0) then
          ! Create artificial register
          MFeed_Table(1,1) = -1.d30; MFeed_Table(2,1) = 0.d0
          MFeed_Table(1,2) = Feed%TimeStart; MFeed_Table(2,2) = 0.d0
          MFeed_Table(1,3) = Feed%TimeFeed; MFeed_Table(2,3) = Feed%MFeed
          N = 3
       
       do i = 1,N
          MFeed_Table(1,i) = MFeed_Table(1,i)*3600.d0 ! [hours] -> [sec]
          MFeed_Table(2,i) = MFeed_Table(2,i)*1000.d0 ! [t] -> [kg] 
       enddo
       endif
       Name = trim(Name)//'_MFeed_Reg'
       call SetupRegister(Name,N,t_v = MFeed_Table)
       Feed%MFeed_Reg => GetRegisterByName(Name)

       Feed%Mat = SolidMat
       if(.not.IsValidSolidMaterial(Feed%Mat)) goto 2300
       Feed%NP = NFeedPart
       Feed%NCell = NFeedCell
       Feed%WriteFeeder = WriteFeeder /= 0  
       Feed%FeederFile = trim(adjustl(FeederFile))
       Feed%WriteCatcher = WriteCatcher /= 0 
       Feed%CatcherFile = trim(adjustl(CatcherFile))
       Feed%NCatchCell = NCatchCell
       Feed%WriteCatcherCDF = WriteCatcherCDF /= 0
       Feed%CatcherCDFFile = trim(adjustl(CatcherCDFFile))
       Feed%CatcherCDFTime = CatcherCDFTime
       Feed%CatcherCDFdT = CatcherCDFdT
       ! Write particle paths
       Feed%NTrace = NTrace
       Feed%Prefix = trim(adjustl(Prefix))
    enddo FEEDER_LOOP

    FEED_PART_LOOP: do iFeed = 1,NFeeders
       Feed => Feeders(iFeed)
       NSort = 0
       rewind(UNITPAR)
       PART_COUNT: do
          read(UNITPAR,NML=PART,end=300,err=2000)
          if(FeedName == Feed%Name) NSort = NSort+1
       enddo PART_COUNT
300    continue
       if(NSort <= 0) goto 2400
       Feed%NSort = NSort
       allocate(Feed%DPart(Feed%NSort),Feed%EpsPart(Feed%NSort),&
            Feed%MassFracPart(Feed%NSort),Feed%WritePath(Feed%NSort),&
            Feed%ForceSediment(Feed%NSort))
       rewind(UNITPAR)
       iSp = 0
       PART_SORT_LOOP: do 
          read(UNITPAR,NML=PART,end=350,err=2000)
          if(FeedName /= Feed%Name) cycle PART_SORT_LOOP
          iSp = iSp+1
          Feed%DPart(iSP) = D*1.d-3 ! [mm]->[m]
          Feed%EpsPart(iSp) = PorInternal*1.d-2   ! [%] -> [-]
          Feed%MassFracPart(iSp) = MassFrac*1.d-2 ! [%] -> [-]
          Feed%WritePath(iSp) = WritePath /= 0
          Feed%ForceSediment(iSp) = ForceSediment /= 0
       enddo PART_SORT_LOOP
350    continue

    enddo FEED_PART_LOOP
    !
    ! Setup fictitious feeder to represent the initial debris bed
    !
    if(IS_DEBRIS_BED) then
       do iDB = 1,NDBParts          
          if(DBPart(iDB)%type /= POR_TYPE_NORMAL) cycle
          NFeeders = NFeeders+1
          Feed => Feeders(NFeeders)
          Feed%ACTIVE = .false.
          Feed%Name = DBPart(iDB)%Name
          Feed%Mat = DBPart(iDB)%Mat
          if(.not.IsValidSolidMaterial(Feed%Mat)) goto 2300
          Feed%NSort = 1
          allocate(Feed%DPart(Feed%NSort),Feed%EpsPart(Feed%NSort),&
               Feed%MassFracPart(Feed%NSort))
          Feed%DPart(1) = DBPart(iDB)%DDisp
          Feed%EpsPart(1) = DBPart(iDB)%PorInternal
          Feed%MassFracPart(1) = 1.d0
          Feed%NCatchCell = 1
          Feed%WriteCatcherCDF = .false.
          Feed%MLastFed = UNASSIGNED
          Feed%dtLastFed = UNASSIGNED
          Feed%NP = UNASSIGNED_INT
       enddo
       do iFeed = NFeeders+1,NTotalFeeders
          Feed => Feeders(iFeed)
          Feed%NSort = -1
       enddo
    endif
    !
    ! Numerical scheme for particles
    !
    rewind(UNITPAR)
    RandomWalk = 0
    Avalanche = UNASSIGNED_INT
    read(UNITPAR,NML=SCHEME_PART,end=400,err=2000)
    IS_RANDOM_WALK = RandomWalk /= 0
    IS_AVALANCHE = Avalanche > 0

400 continue

    return
2300 continue
    write(*,788) Feed%Mat,iFeed,trim(adjustl(FilePAR))
    close(UNITPAR)
    stop
2400 continue
    write(*,799) Feed%NSort,iFeed,trim(adjustl(FilePAR))
    close(UNITPAR)
    stop
788 format(1x,'Bad material ',i0', for Feeder #',i0,' in file ',a)
799 format(1x,'Bad number of particle sorts ',i0', for Feeder #',i0,' in file ',a)
2000 continue
    ReadParticles = .false.
  end function ReadParticles

  logical function ReadBoundaryConditions()
    use Globals
    use InOut
    use Boundary
    use Heat
    use Storage
    implicit none
    !
    logical:: found
    integer:: NB
    character(len=256) buf
    integer:: iB

    character(len=256):: Name, type
    real(kind=8):: UPH(NPHASE)
    real(kind=8):: VPH(NPHASE)
    real(kind=8):: Temp(NPHASE)
    real(kind=8):: Void, TurbLevel, TurbScale
    real(kind=8):: Prel, K, Eps, DHydr
    integer:: TurbMode, TempMode
    integer:: lName, iName
    character(len=128):: StorageName

    namelist/BND/ Name, type, UPH, VPH, Prel, Void, &
         TurbLevel, TurbScale, K, Eps, DHydr,&
         TurbMode, TempMode, Temp, StorageName

    TRACE_LINE = "ReadBoundaryConditions"
    ReadBoundaryConditions = .true.
    !
    ! Count boundary parts
    !
    rewind UNITPAR
    NBndCond = 0
    do
       read(UNITPAR,NML=BND,end=200,err=2000)
       NBndCond = NBndCond+1
    enddo
200 continue
    if(NBndCond == 0) goto 2000
    !
    ! Allocate boundary descriptors
    !
    allocate(BoundaryCondition(NBndCond))
    do iB = 1,NBndCond
       allocate(BoundaryCondition(iB)%U000(NPHASE))
       allocate(BoundaryCondition(iB)%V000(NPHASE))
       allocate(BoundaryCondition(iB)%U(NPHASE))
       allocate(BoundaryCondition(iB)%V(NPHASE))
       BoundaryCondition(iB)%U000 = UNASSIGNED
       BoundaryCondition(iB)%V000 = UNASSIGNED
       BoundaryCondition(iB)%U    = UNASSIGNED
       BoundaryCondition(iB)%V    = UNASSIGNED
       BoundaryCondition(iB)%P000 = UNASSIGNED
       BoundaryCondition(iB)%P = UNASSIGNED
       if(IS_ENERGY) then
          allocate(BoundaryCondition(iB)%Temp(NPHASE))
          allocate(BoundaryCondition(iB)%T(NPHASE))
          BoundaryCondition(iB)%Temp  = UNASSIGNED
       endif
       BoundaryCondition(iB)%TempMode = -1
       BoundaryCondition(iB)%StorageName = "none"
       BoundaryCondition(iB)%nStor = -1
    enddo
    !
    ! Read all boundary conditions
    !
    rewind UNITPAR
    NB=0
    BND_LOOP: do

       UPH = UNASSIGNED; VPH = UNASSIGNED; Prel = UNASSIGNED
       Void = UNASSIGNED; Temp = UNASSIGNED
       TurbLevel = UNASSIGNED; TurbScale = UNASSIGNED
       K = UNASSIGNED; Eps = UNASSIGNED; DHydr = UNASSIGNED
       TurbMode = -1; TempMode = -1; 
       StorageName = ""
       read(UNITPAR,NML=BND,end=300,err=2000)
       NB = NB+1
       !
       ! Read boundary condition
       !
       BoundaryCondition(NB)%Name = trim(adjustl(Name))
       found = .false.
       do iB = 1,KNOWNBOUNDARYTYPE
          if(trim(adjustl(type)) == trim(adjustl(KnownBoundaryTypeNames(iB)))) then
             found = .true.
             BoundaryCondition(NB)%type = KnownBoundaryTypes(iB)
             exit
          endif
       enddo
       if(.not.found) goto 3000
       lName = len_trim(BoundaryCondition(NB)%Name)
       do iB = 1,NB-1
          iName  = len_trim(BoundaryCondition(iB)%Name)
          if(iName /= lName) cycle
          if(trim(BoundaryCondition(NB)%Name) == trim(BoundaryCondition(iB)%Name)) then
             write(*,1777) BoundaryCondition(iB)%Name
             stop
          endif
       enddo
       !
       ! Initialize all boundary variables
       !
       if(IS_SINGLE_PHASE) Void = 0.d0
       BoundaryCondition(NB)%P000 = Prel
       BoundaryCondition(NB)%Alpha = Void
       BoundaryCondition(NB)%U000 = UPH
       BoundaryCondition(NB)%V000 = VPH
       BoundaryCondition(NB)%modeK = -1
       BoundaryCondition(NB)%K000 = K
       BoundaryCondition(NB)%TurbLevel = TurbLevel
       BoundaryCondition(NB)%modeEps = -1
       BoundaryCondition(NB)%Eps000 = Eps
       BoundaryCondition(NB)%TurbScale = TurbScale
       BoundaryCondition(NB)%DHydr = DHydr
       if(IS_ENERGY) then
          BoundaryCondition(NB)%TempMode = TempMode
          BoundaryCondition(NB)%Temp = Temp
       endif
       !
       ! Assign values specific to boundary types
       !
       select case(BoundaryCondition(NB)%type)
       case(BND_SOLID_WALL)
          if(minval(UPH) < IFASSIGNED) BoundaryCondition(NB)%U000 = 0.d0
          if(minval(VPH) < IFASSIGNED) BoundaryCondition(NB)%V000 = 0.d0
          if(IS_ENERGY .and. TempMode == WALL_TEMP) then
             if(minval(Temp) < IFASSIGNED ) goto 3510
          endif
       case(BND_VELOCITY_INLET,BND_PRESSURE_OUTLET)
          if(BoundaryCondition(NB)%type == BND_VELOCITY_INLET .and. &
               ( minval(UPH) < IFASSIGNED .or. minval(VPH) < IFASSIGNED)) goto 3500
          if(BoundaryCondition(NB)%type == BND_PRESSURE_OUTLET .and. &
               Prel < IFASSIGNED) goto 3700
          if(Void < 0.d0 .or. Void > 1.d0) then
             if(BoundaryCondition(NB)%type == BND_VELOCITY_INLET) goto 3600
             if(BoundaryCondition(NB)%type == BND_PRESSURE_OUTLET) goto 3605
          endif
          if(IS_ENERGY) then
             if(minval(Temp) < IFASSIGNED ) goto 3510
          endif
          if(IS_K_Eps_TurbModel) then
             if(TurbMode < 0) goto 3610
             !
             ! Set turbulent quantities
             !
             select case(TurbMode)
             case(TURB_BND_K_AND_EPS)
                BoundaryCondition(NB)%modeK = BND_K_VALUE
                if(BoundaryCondition(NB)%K000 < IFASSIGNED) goto 3620
                BoundaryCondition(NB)%modeEps = BND_EPS_VALUE
                if(BoundaryCondition(NB)%Eps000 < IFASSIGNED) goto 3630
             case(TURB_BND_INTENS_SCALE_DHYDR)
                BoundaryCondition(NB)%modeK = BND_K_INTENSITY
                if(BoundaryCondition(NB)%TurbLevel < IFASSIGNED) goto 3625
                BoundaryCondition(NB)%TurbLevel = 1.d-2*BoundaryCondition(NB)%TurbLevel ! [%] -> [-]
                BoundaryCondition(NB)%modeEps = BND_EPS_SCALE
                if(BoundaryCondition(NB)%TurbScale < IFASSIGNED) goto 3635
                if(BoundaryCondition(NB)%DHydr < IFASSIGNED) goto 3640
             case default
                goto 3610
             end select
          endif
       case(BND_STORAGE_CONNECT)
          if(len_trim(StorageName)<= 0) goto 3800
          BoundaryCondition(NB)%StorageName = StorageName
          BoundaryCondition(NB)%nStor = FindStorage(StorageName)
          if(BoundaryCondition(NB)%nStor <= 0) goto 3810
          if(IS_K_Eps_TurbModel) then
             if(TurbMode < 0) goto 3610
             !
             ! Set turbulent quantities
             !
             select case(TurbMode)
             case(TURB_BND_K_AND_EPS)
                BoundaryCondition(NB)%modeK = BND_K_VALUE
                if(BoundaryCondition(NB)%K000 < IFASSIGNED) goto 3620
                BoundaryCondition(NB)%modeEps = BND_EPS_VALUE
                if(BoundaryCondition(NB)%Eps000 < IFASSIGNED) goto 3630
             case(TURB_BND_INTENS_SCALE_DHYDR)
                BoundaryCondition(NB)%modeK = BND_K_INTENSITY
                if(BoundaryCondition(NB)%TurbLevel < IFASSIGNED) goto 3625
                BoundaryCondition(NB)%TurbLevel = 1.d-2*BoundaryCondition(NB)%TurbLevel ! [%] -> [-]
                BoundaryCondition(NB)%modeEps = BND_EPS_SCALE
                if(BoundaryCondition(NB)%TurbScale < IFASSIGNED) goto 3635
                if(BoundaryCondition(NB)%DHydr < IFASSIGNED) goto 3640
             case default
                goto 3610
             end select
          endif
       end select
    enddo BND_LOOP
    !
    ! Normal exit
    !
300 continue
    ReadBoundaryConditions = .true.
    return
    !
    ! Error
    !
2000 continue
    write(*,555) trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3000 continue
    write(*,777)   trim(adjustl(FilePAR)),&
         trim(adjustl(buf)),&
         KnownBoundaryTypeNames(1:KNOWNBOUNDARYTYPE)
    ReadBoundaryConditions = .false.
    return
3500 continue
    write(*,888) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3510 continue
    write(*,868) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3600 continue
    write(*,889) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3607 continue
    write(*,858) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3605 continue
    write(*,899) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3610 continue
    write(*,910) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3620 continue
    write(*,920) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3625 continue
    write(*,925) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3630 continue
    write(*,930) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3635 continue
    write(*,935) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3640 continue
    write(*,940) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3650 continue
    write(*,890) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3700 continue
    write(*,891) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
3800 continue
    write(*,1800) trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
3810 continue
    write(*,1810) trim(adjustl(StorageName)),trim(adjustl(Name)), trim(adjustl(FilePAR))
    ReadBoundaryConditions = .false.
    return
444 format(1x,'Unable to open file ',a)
555 format(1x,'End of file or error reading ',a)
777 format(1x,'Bad boundary type found in file ',a,': ',a/&
         1x,'Known boundary types:'/(5x,a24))
888 format(1x,'Missing inlet velocities (UPH or VPH) for boundary ',a,' in file ',a)
858 format(1x,'Missing wall heat flag (TempMode) for boundary ',a,' in file ',a)
868 format(1x,'Missing temperatures (Temp) for boundary ',a,' in file ',a)
889 format(1x,'Missing or bad inlet void fraction for boundary ',a,' in file ',a)
899 format(1x,'Missing or bad backflow void fraction for boundary ',a,' in file ',a)
890 format(1x,'Missing hydraulic diameter for boundary ',a,' in file ',a)
891 format(1x,'Missing outlet pressure for boundary ',a,' in file ',a)
910 format(1x,'Missing or bad TurbMode for boundary ',a,' in file ',a)
920 format(1x,'Missing or bad K for boundary ',a,' in file ',a)
925 format(1x,'Missing or bad TurbLevel for boundary ',a,' in file ',a)
930 format(1x,'Missing or bad Eps for boundary ',a,' in file ',a)
935 format(1x,'Missing or bad TurbScale for boundary ',a,' in file ',a)
940 format(1x,'Missing or bad DHydr for boundary ',a,' in file ',a)
1800 format(1x,'Missing or bad StorageName for boundary ',a,' in file ',a)
1810 format(1x,'Storage ',a,' not found for boundary ',a,' in file ',a)
1777 format(1x,'Duplicate boundary names are not allowed: ',a)
  end function ReadBoundaryConditions

  logical function ReadStorages()
    use Globals
    use InOut
    use Multiphase
    use Boundary
    use Heat
    use Storage
    implicit none  

    integer:: iStor

    real(kind=8):: ZBottom
    real(kind=8):: Area
    real(kind=8):: Level
    real(kind=8):: Vol
    real(kind=8):: Temp
    real(kind=8):: PTop

    integer:: type
    integer:: LevType
    real(kind=8):: fracCondens

    character(len=128):: Name
    character(len=128):: Parm

    namelist/STOR/ Name,ZBottom,Area,Level,Vol,type,LevType,fracCondens, Temp, PTop

    TRACE_LINE="ReadStorages"

    ReadStorages = .true.
    !
    ! Count storages
    !
    rewind UNITPAR
    NStor = 0
    do
       read(UNITPAR,NML=STOR,end=200,err=2000)
       NStor = NStor+1
    enddo
200 continue
    if(NStor <= 0) return
    !
    ! Read storage parameters
    !
    allocate(Storages(NStor))
    rewind UNITPAR
    STOR_LOOP: do iStor = 1,NStor
       Str => Storages(iStor)
       Name = ""
       ZBottom = UNASSIGNED
       Area = UNASSIGNED
       Level = UNASSIGNED
       Temp = UNASSIGNED
       PTop = UNASSIGNED
       type = UNASSIGNED_INT
       Vol = UNASSIGNED
       LevType = UNASSIGNED_INT
       fracCondens = UNASSIGNED
       read(UNITPAR,NML=STOR,end=2000,err=2000)
       !
       if(len_trim(Name) <= 0) goto 2100
       Str%Name = Name
       if(ZBottom < IFASSIGNED) then
          Parm = "ZBottom"; goto 2200
       endif
       Str%ZBottom = ZBottom
       if(Area < IFASSIGNED) then
          Parm = "Area"; goto 2200
       endif
       Str%Area = Area

       if(Vol < IFASSIGNED) then
          if(Level < IFASSIGNED) then 
             Parm = "Level"; goto 2200
          endif
          Str%Level = Level
          Str%Vol = Str%Area*Str%Level
       else
          Str%Vol = Vol
          Str%Level = Str%Vol/Str%Area
       endif

       if(type == UNASSIGNED_INT) then
          Str%type = STOR_CONDENSER
       else
          Str%type = type
       endif

       if(LevType == UNASSIGNED_INT) then
          Str%LevType = STOR_CONST_MASS
       else
          Str%LevType = LevType
       endif
       !
       ! Operating conditions
       !
       if(Temp < IFASSIGNED) then
          Str%Temp = TSystem000
       else
          Str%Temp = Temp
       endif
       if(PTop < IFASSIGNED) then
          Str%PTop = PSystem000
       else
          Str%PTop = PTop
       endif

       if(fracCondens < IFASSIGNED) then
          Str%fracCondens = 1.d0
       else
          Str%fracCondens = fracCondens
       endif
    enddo STOR_LOOP
    return
    !
    ! Error
    !
2000 continue
    write(*,555) trim(adjustl(FilePAR))
    ReadStorages = .false.
    return
2100 continue
    write(*,600) trim(adjustl(FilePAR))
    ReadStorages = .false.
    return
2200 continue
    write(*,610) trim(adjustl(Parm)),trim(adjustl(Name)),trim(adjustl(FilePAR))
    ReadStorages = .false.
    return
555 format(1x,'End of file or error reading ',a)
600 format(1x,'Missing or bad Storage name in file ',a)
610 format(1x,'Missing or bad ',a,' for Storage ',a,' in file ',a)
  end function ReadStorages
end subroutine INPUT

