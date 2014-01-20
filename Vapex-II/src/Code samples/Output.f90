!###################################################################!
!                                                                   !
! DECOSIM:  DEbris COolability SIMulator                            !
!                                                                   !
!###################################################################!
! $Id: Output.f90 5 2014-01-13 11:44:50Z Sergey $
!-------------------------------------------------------------------!
! Output of results for post-processing and restart                 !
!-------------------------------------------------------------------!
module OUTPUT
  use Dimensions

  ! Near-wall profiles
  integer,parameter:: WALL_NORMAL_Y_POS = 1
  integer,parameter:: WALL_NORMAL_Y_NEG = -1

  integer,parameter:: HOR_WALL_PROFILE = 1
  integer,parameter:: VERT_WALL_PROFILE = 2

  
  ! 2D fields in Tecplot format
  logical:: IS_TRANSIENT_TECPLOT
  logical:: IS_WRITE_DP_TECPLOT = .false.
  character(len=20):: tecPrefix = "fn"
  character(len=20):: tecSuffix = "tec"

  type WALL_NORMAL_PROFILE
     character(len=256):: Name
     character(len=256):: filename
     real(kind=8):: Coord
     integer:: Dir
     integer:: type
  end type WALL_NORMAL_PROFILE

  type(WALL_NORMAL_PROFILE), allocatable:: WallNormalProfiles(:)
  integer::  NWallNormalProfiles = 0

  character(len=256):: FileDRAG = 'drag.dat'
  integer,parameter:: UNITDRAG = 35

  type TRANSDUCER_TYPE
     character(len=256):: Name = ""
     character(len=256):: Label = ""
     real(kind=8):: Coord(NDIMS) = (/0.d0,0.d0/)
     integer:: Ind(NDIMS) = (/-1,-1/)
     logical:: Active = .false.
     integer:: Storage = 1
  end type TRANSDUCER_TYPE

  integer,parameter:: TRANSDUCER_POINT = 0
  integer,parameter:: TRANSDUCER_INTEGRAL = 1
  integer,parameter:: TRANSDUCER_STORAGE = 2

  type(TRANSDUCER_TYPE), pointer:: Trans(:)
  integer:: NTrans = 0
  integer, parameter:: UNITTRANS = 45
  character(len=256):: FileTRANS = 'Trans.dat'
  integer,parameter:: lenTRANS = 14
  character(len=20):: formTRANS = '(1x,ES12.5)'
  logical:: IS_WRITE_TRANS = .false.

  character(len=256):: FileINFO = 'DECOSIM.info'
  integer,parameter:: UNITINFO = 95
  


contains
 

  subroutine WriteTECPLOT(NPict,isGrid)
    use Globals
    use GlobalArrays
    use Dimensions
    use SolutionVector
    use MultiPhase
    use Geometry
    use Turb
    use Level
    use Heat
    use DebrisBed
    use SolidMaterial
    use Particle

    implicit none
    integer:: NPict ! Number of picture
    logical, intent(in), optional:: isGrid ! if .true., grid file is written separately
    integer::i,j
    character(len=256) filename
    character(len=12) title
    character(len=2048) line
    real(kind=4):: Vals(NARRAYS),  ValTLiq, ValTGas, ValPor, ValCong, ValK, ValEps
    real(kind=4):: ValVisc,ValDp, ValGamma, ValTSup, ValTGSup,ValQM
    real(kind=4):: ValsUPH(NPHASE), ValsVPH(NPHASE), ValsAPH(NPHASE), ValsJUPH(NPHASE), ValsJVPH(NPHASE)
    real(kind=4):: ValTs, ValTc, ValTsat, ValMeltFrac
    real(kind=8):: APH(NPHASE,2,2) 
    integer:: SumMask, iV, iPh, IndVar, IndBnd
    real(kind=8):: QBnd,X,Y, DpOut
    logical:: isNeumann
    integer:: ii,jj,iMat,SumN


    write(title,'(f8.3,a1)') Time*TimeScale,'s'

999 format('VARIABLES= "R","Z"')
991 format('VARIABLES= "R","Z","Phi","mask"')
992 format('VARIABLES= "R","Z",',10(a6))
998 format('ZONE')
996 format(2X,'T= "',a12,'"',2X,'I=',I0,2X,'J=',I0)
800 format(1X,'SOLUTIONTIME = ',G15.6)

995 format('VARIABLES= "R[m]","Z[m]","U[m/s]","V[m/s]","Ug[m/s]","Vg[m/s]","JU[m/s]",'&
         '"JV[m/s]","JUg[m/s]","JVg[m/s]","P[Pa]","a_l[-]","a_v[-]","a_s[-]","Eps_c[-]","Dp[mm]","K","Eps","Visc"')
980 format('VARIABLES= "X[m]","Y[m]","U[m/s]","V[m/s]","Ug[m/s]","Vg[m/s]","JU[m/s]",'&
         '"JV[m/s]","JUg[m/s]","JVg[m/s]","P[Pa]","a_l[-]","a_v[-]","a_s[-]","Eps_c[-]","Dp[mm]"')
985 format('VARIABLES= "R[m]","Z[m]","U[m/s]","V[m/s]","Ug[m/s]","Vg[m/s]","JU[m/s]",'&
         '"JV[m/s]","JUg[m/s]","JVg[m/s]","P[Pa]","a_l[-]","a_v[-]","a_s[-]","Eps_c[-]","Dp[mm]"')
988 format(',"K","Eps","Visc"')


880 format('VARIABLES= "X","Y","Dp[mm]"')
890 format('VARIABLES= "R","Z","Dp[mm]"')
891 format('DATAPACKING=BLOCK, VARLOCATION=([3]=CELLCENTERED)')

    if(present(isGrid)) then
       if(isGrid) then
          !
          ! Writing grid into a separate file
          !
          filename = 'grid.dat'
          if(len_trim(OutputDirPrefix) > 0 .and. trim(adjustl(OutputDirPrefix)) /= './') then
             filename = trim(adjustl(OutputDirPrefix))//trim(adjustl(filename))
          endif
          open(50,file=filename)
          write(50,999)
          write(50,998)
          write(50,996) title,NX1,NY1
          do J = 1,NY1
             do I = 1,NX1
                write (50,'(18(g15.6,1x))') &
                     XK(i,j)*XScale,YK(i,j)*XScale
             end do
          end do
          close(50)
       endif
    endif

    if(IS_WRITE_DP_TECPLOT) then
       if(NDBParts > 0) then
          filename='Dp0000.dat'
          write(filename(3:6),'(i4.4)') NPict
          if(len_trim(OutputDirPrefix) > 0 .and. trim(adjustl(OutputDirPrefix)) /= './') then
             filename = trim(adjustl(OutputDirPrefix))//trim(adjustl(filename))
          endif
          open(50,file=filename)
          if(AXIAL_SYMMETRY_X) then
             write(line,890)
          else
             write(line,880)
          endif
          write(50,'(a)') trim(line)
          write(50,998)
          write(50,996) title,NX1,NY1
          write(50,891) 
          do J = 1,NY1
             do I = 1,NX1
                write (50,'(18(g15.6,1x))') XK(i,j)*XScale
             enddo
          enddo
          do J = 1,NY1
             do I = 1,NX1
                write (50,'(18(g15.6,1x))') YK(i,j)*XScale
             enddo
          enddo

          do J = 2,NY1
             do I = 2,NX1
                DpOut = D_Disp(i,j)*1.d3 
                if(Porosity(i,j) > 0.95) DpOut = 0.d0
                write (50,'(18(g15.6,1x))') DpOut
             enddo
          enddo

          close(50)
       endif
    endif


    filename='0000'
    write(filename(1:4),'(i4.4)') NPict
    filename = trim(adjustl(tecPrefix))//trim(filename)//'.'//trim(adjustl(tecSuffix))
    if(len_trim(OutputDirPrefix) > 0 .and. trim(adjustl(OutputDirPrefix)) /= './') then
       filename = trim(adjustl(OutputDirPrefix))//trim(adjustl(filename))
    endif
    open(50,file=filename)
    line='VARIABLES='
    if(AXIAL_SYMMETRY_X) then
       call AddName(line,"R[m]")
       call AddName(line,"Z[m]")
    else
       call AddName(line,"X[m]")
       call AddName(line,"Y[m]")
    endif
    call AddName(line,"Ul[m/s]")
    call AddName(line,"Vl[m/s]")
    call AddName(line,"Ug[m/s]")
    call AddName(line,"Vg[m/s]")
    call AddName(line,"JUl[m/s]")
    call AddName(line,"JVl[m/s]")
    call AddName(line,"JUg[m/s]")
    call AddName(line,"JVg[m/s]")

    call AddName(line,"P[Pa]")
    call AddName(line,"a_l[-]")
    call AddName(line,"a_g[-]")

    if(IS_ENERGY) then
       call AddName(line,"Tf[K]")
       call AddName(line,"TfSup[K]")
       call AddName(line,"Tg[K]")
       call AddName(line,"TgSup[K]")
    endif

    if(IS_EVAPORATION) then
       call AddName(line,"Gamma")
    endif

    if(IS_ACTIVE_FEEDERS) then
       call AddName(line,"QMELT")
    endif

    if(NDBParts > 0) then
       call AddName(line,"a_s[-]")
       call AddName(line,"Dp[mm]")
    endif

    if(NUM_POROUS > 0) then
       call AddName(line,"Ts[K]")
       call AddName(line,"Ts-Tsat[K]")
    endif

    if(NUM_CONGEST > 0) then
       call AddName(line,"Eps_c")
       call AddName(line,"Tc[K]")
    endif

    if(IS_MELT_MODEL .and. NUM_POROUS > 0) then
       call AddName(line,"MeltFrac[-]")
    endif

    if(IS_K_Eps_TurbModel) then
       call AddName(line,"K")
       call AddName(line,"Eps")
       call AddName(line,"Visc")
    endif

    write(50,'(a)') trim(line)
    write(50,998)
    write(50,996) title,NX1,NY1

    do  J = 1,NY1
       do I = 1,NX1

          SumMask = sum(abs(mask(i:i+1,j:j+1)))
          if(SumMask >0) then
             APH(1:NPHASE,1:2,1:2) = AlphaPH(1:NPHASE,i:i+1,j:j+1)
             if(IS_POOL) then
                if(PoolF(i,j) > 0) APH(:,1,1) = AlphaPH(:,i,j+1)
                if(PoolF(i+1,j) > 0) APH(:,2,1) = AlphaPH(:,i+1,j+1)
                if(PoolF(i,j+1) > 0) APH(:,1,2) = AlphaPH(:,i,j+1)
                if(PoolF(i+1,j+1) > 0) APH(:,2,2) = AlphaPH(:,i+1,j+1)
             endif
             do iV = 1,NARRAYS
                Vals(iV) = &
                                !      sum(Q(iV,i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/&
                                !      SumMask
                     sum(Q(iV,i:i+1,j:j+1))*0.25d0
             enddo
             do iV = 1,NPHASE

                ValsUPH(iV) = 0.5d0*(UFPH(iV,i,j)+UFPH(iV,i,j+1))
                ValsVPH(iV) = 0.5d0*(VFPH(iV,i,j)+VFPH(iV,i+1,j))

                if(IS_POOL) then
                   if(PoolF(i,j) /= 0) then
                      ValsUPH(iV) = 0.5d0*(UFPHstar(iV,i,j)+UFPHstar(iV,i,j+1))
                      ValsVPH(iV) = 0.5d0*(VFPHstar(iV,i,j)+VFPHstar(iV,i+1,j))
                   endif
                endif

                ValsAPH(iV) = &
                     sum(APH(iV,1:2,1:2)*abs(mask(i:i+1,j:j+1)))/SumMask

                ValsJUPH(iV) =0.5d0*(UJPH(iV,i,j)+UJPH(iV,i,j+1))
                ValsJVPH(iV) = 0.5d0*(VJPH(iV,i,j)+VJPH(iV,i+1,j))

             enddo

!!$             ValTGas = sum(Temp_Gas(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/&
!!$                  SumMask*TempScale
!!$             ValTLiq = sum(Temp_Liq(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/&
!!$                  SumMask*TempScale
!!$
!!$             ValGamma = sum(Gamma_Evap(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/&
!!$                  SumMask

             ValTGas = sum(Temp_Gas(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/SumMask*TempScale
             ValTLiq = sum(Temp_Liq(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/SumMask*TempScale
             ValGamma = sum(Gamma_Evap(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/SumMask
             if(abs(ValGamma) < 1.d-10) ValGamma = 0.d0
             ValTSup = sum(TSubSup(IndLiq,i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/SumMask*TempScale
             if(abs(ValTSup) < 1.d-8) ValTSup = 0.D0
             ValTGSup = sum(TSubSup(IndGas,i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/SumMask*TempScale
             if(abs(ValTSup) < 1.d-8) ValTGSup = 0.D0
             ValTSat = sum(TempSat(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/SumMask*TempScale
             ValQM = sum(QMELT(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/SumMask

!!$             if(ValsAPH(IndLiq) < 1.d-4) then
!!$                ValTGas = TSystem000
!!$                ValTLiq = TSystem000
!!$                ValTSup = TSystem000 - ValTSat
!!$             endif

             ValPor = 1.d0-sum(PorReal(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/&
                  SumMask
             ValCong = sum(PorCong(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/&
                  SumMask
             if(IS_MELT_MODEL) then
                ValMeltFrac = 0.d0
                SumN = 0
                do jj = j,j+1
                   do ii = i,i+1
                      iMat = maskMaterial(ii,jj)
                      if(iMat <= 0) cycle
                      if(SolidMaterials(iMat)%Type == SOLID_TYPE_PASSIVE) cycle
                      ValMeltFrac = ValMeltFrac + MeltFrac(iMat,ii,jj)
                      SumN = SumN + 1
                   enddo
                enddo
                if(SumN > 0) ValMeltFrac = ValMeltFrac/SumN

             endif

             ValDp = sum(D_Disp(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/&
                  SumMask*1.d3 ! [m] -> [mm]
             if(IS_K_Eps_TurbModel) then
                ValK = sum(k_Eps_K(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/&
                     SumMask
                ValEps = sum(k_Eps_Eps(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/&
                     SumMask
             else
                ValK = 0.d0; ValEps = 0.d0
             endif

             ValVisc = sum(ViscPH(IndLiq,i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/&
                  SumMask

             if(NUM_POROUS > 0) then
                !           ValTs = sum(TempPor(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/&
                !                SumMask     
                ValTs = sum(TempPor(i:i+1,j:j+1))*0.25D0
             endif
             if(NUM_CONGEST > 0) then
                ValTc = sum(TempCong(i:i+1,j:j+1)*abs(mask(i:i+1,j:j+1)))/&
                     SumMask
             endif
          else
             Vals = 0
             ValTLiq = TSystem000
             ValTGas = TSystem000
             ValTSup = 0.d0
             ValPor = 0.d0
             ValCong = 0.d0
             ValDp = 0.d0
             ValsUPH = 0.d0
             ValsVPH = 0.d0
             do iPh = 1,NPHASE
                ValsAPH(iPh) = sum(AlphaPH(iPh,i:i+1,j:j+1))*0.25d0
             enddo
             ValTs = TSystem000
             ValTSat = TSystem000
             ValTc = TSystem000
          endif

          where(Vals(IndAlpha_l:IndAlpha_g) < 1.d-8) Vals(IndAlpha_l:IndAlpha_g) = 0.D0
          where(Vals(IndAlpha_l:IndAlpha_g) > 1.d0)  Vals(IndAlpha_l:IndAlpha_g) = 1.D0

          if(i==1) ValsUPH = 0.d0 ! Axis
          if(i==NX1) ValsUPH = 0.d0 ! Wall


          Vals(IndUVelocity:IndUVelocity+NDIMS-1) = &
               Vals(IndUVelocity:IndUVelocity+NDIMS-1)*VScale

          if(j==1 .and. mask(i,j+1) /= 0) then
             IndBnd = -min(SS(LEFT_Y,i,j+1),SS(LEFT_Y,i+1,j+1))
             X = 0.5*(XC(i,j+1)+XC(i+1,j+1))
             Y = YB(i,j)
             do iPh = 1,NPHASE
                IndVar = IndUPH(iPh)
                call BOUNDARY_VALUE(IndVar, IndBnd, LEFT_Y,i,j+1,X, Y, isNeumann, QBnd)
                if(.not.isNeumann) ValsUPH(iPh) = QBnd
                IndVar = IndVPH(iPh)
                call BOUNDARY_VALUE(IndVar, IndBnd, LEFT_Y,i,j+1,X, Y, isNeumann, QBnd)
                if(.not.isNeumann) ValsVPH(iPh) = QBnd
                IndVar = IndAPH(iPh)
                call BOUNDARY_VALUE(IndVar, IndBnd, LEFT_Y,i,j+1,X, Y, isNeumann, QBnd)
                if(.not.isNeumann) ValsAPH(iPh) = QBnd
             enddo
          endif

          ValsUPH = ValsUPH*VScale; ValsVPH = ValsVPH*VScale
	  ValsJUPH = ValsJUPH*VScale; ValsJVPH = ValsJVPH*VScale

          line = ''
          call AddValue(line,real(XK(i,j)*XScale,kind=4))
          call AddValue(line,real(YK(i,j)*XScale,kind=4))
          call AddValue(line,ValsUPH(IndLiq))
          call AddValue(line,ValsVPH(IndLiq))
          call AddValue(line,ValsUPH(IndGas))
          call AddValue(line,ValsVPH(IndGas))
          call AddValue(line,ValsJUPH(IndLiq))
          call AddValue(line,ValsJVPH(IndLiq))
          call AddValue(line,ValsJUPH(IndGas))
          call AddValue(line,ValsJVPH(IndGas))
          call AddValue(line,Vals(IndPressure))
          call AddValue(line,ValsAPH(IndLiq))
          call AddValue(line,ValsAPH(IndGas))

          if(IS_ENERGY) then
             call AddValue(line,ValTLiq)            
             call AddValue(line,ValTSup)
             call AddValue(line,ValTGas)
             call AddValue(line,ValTGSup)
          endif

          if(IS_EVAPORATION) then
             call AddValue(line,ValGamma)
          endif

          if(IS_ACTIVE_FEEDERS) then
             call AddValue(line,ValQM)
          endif

          if(NDBParts > 0) then
             call AddValue(line,ValPor)
             call AddValue(line,ValDp)
          endif

          if(NUM_POROUS > 0) then
             call AddValue(line,ValTs)
             call AddValue(line,ValTs-ValTsat)
          endif

          if(NUM_CONGEST > 0) then
             call AddValue(line,ValCong)
             call AddValue(line,ValTc)
          endif

          if(IS_MELT_MODEL .and. NUM_POROUS > 0) then
             call AddValue(line,ValMeltFrac)
          endif

          if(IS_K_Eps_TurbModel) then
             call AddValue(line,ValK)
             call AddValue(line,ValEps)
             call AddValue(line,ValVisc)
          endif

          write(50,'(a)') trim(line)
       enddo
    enddo

    close(50)
    if(MsgLevel > 0) then
       write(*,100) trim(filename)
    endif
    call WriteVertProfile(NPict)

100 format(96('=')/' Writing graphics file ',a/96('='))

  contains

    subroutine AddName(line,name)
      character(*):: line
      character(*), intent(in):: name
      integer:: len

      len = len_trim(line)
      if(len <= 0) then
         line = '"'//trim(adjustl(name))//'"'
      else
         line = trim(adjustl(line))//',"'//trim(adjustl(name))//'"'
      endif
    end subroutine AddName

    subroutine AddValue(line,val)
      character(*):: line
      real(kind=4), intent(in):: val
      character(20):: buf
      integer:: len

      write(buf,'(E16.6E2,1x)') val
      line = trim(line)//trim(buf)
    end subroutine AddValue

  end subroutine WriteTECPLOT

 

  subroutine Salute(iunit)
    use BUILD
    implicit none
    integer,intent(in),optional:: iunit
    integer:: unit
    if(present(iunit)) then
       unit = iunit
    else
       unit = 6
    endif
    write(unit,100) 
    write(unit,200)'######            ######           ####    ##          '
    write(unit,200)'##    ##  ####   ##    ##   ###   ##   ##              '
    write(unit,200)'##    ## ##   ## ##       ##   ##  ###     ##  ##    ##'
    write(unit,200)'##    ## ####### ##       ##   ##     ##   ##  ###  ###'
    write(unit,200)'##    ## ##      ##    ## ##   ## ##   ##  ##  ## ## ##'
    write(unit,200)'######    #####   ######    ###    #####   ##  ## ## ##'
    write(unit,100)

    !    call UpdateRevision
    !    write(UNIT,300) trim(COMPILE_DATE),SVN_REVISION_NUMBER
    !    write(unit,100)
100 format(96('='))
200 format(17(' '),a60,4(' '))
    !300 format(27(' '),'BUILD: ',a,'  SVN Version: ',i0)
    call Print_Version(unit)

  end subroutine Salute

  subroutine Print_Version(iunit)
    use Globals
    use Build
    implicit none
    integer,intent(in),optional:: iunit
    integer:: unit
    if(present(iunit)) then
       unit = iunit
    else
       unit = 6
    endif
    call UpdateRevision
    if(SVN_REVISION_RANGE(1) == SVN_REVISION_RANGE(2)) then
       write(unit,50) SVN_REVISION_NUMBER,&
            trim(adjustl(REVISION_DATE)),trim(adjustl(COMPILE_DATE))
    else
       write(unit,60) SVN_REVISION_RANGE(1:2),&
            trim(adjustl(REVISION_DATE)),trim(adjustl(COMPILE_DATE))
    endif
    write(unit,100)
50  format(26x,'SVN ',i0/26x,'Revision Date:  ',a/26x,&
         'Build Date:     ',a)
60  format(26x,'SVN RANGE ',i0,':',i0/26x,'Revision Date:  ',a/26x,&
         'Build Date:     ',a)
100 format(79('='))
  end subroutine Print_Version

 

 

 


  real(kind=8) function MassOfWater()
    use Globals
    use GlobalArrays
    use Multiphase
    use Geometry
    use Water
    use PVolFrac
    use Turb    
    implicit none

    integer:: i,j
    MassOfWater = 0.d0
    do j = 2,NY1
       do i = 2,NX1
          MassOfWater = MassOfWater + &
               DensPH(IndLiq,i,j)*AlphaPH(IndLiq,i,j)*Porosity(i,j)*&
               2.d0*PI*XC(i,j)*(XB(i,j)-XB(i-1,j))*(YB(i,j)-YB(i,j-1))
       enddo
    enddo
    MassOfWater = MassOfWater*XScale**3*RoScale
    if(IS_INITIAL_BALANCE_WATER) then
       BalIncrMassOfWater = 0.d0
       IS_INITIAL_BALANCE_WATER = .false.
    else
       BalIncrMassOfWater = MassOfWater - BalMassOfWater
    endif
    BalMassOfWater = MassOfWater
  end function MassOfWater

  real(kind=8) function MassOfVap()
    use Globals
    use GlobalArrays
    use Multiphase
    use Geometry
    use Water
    use PVolFrac
    use Turb    
    implicit none

    integer:: i,j
    MassOfVap = 0.d0
    do j = 2,NY1
       do i = 2,NX1
          MassOfVap = MassOfVap + &
               DensPH(IndGas,i,j)*AlphaPH(IndGas,i,j)*Porosity(i,j)*&
               2.d0*PI*XC(i,j)*(XB(i,j)-XB(i-1,j))*(YB(i,j)-YB(i,j-1))
       enddo
    enddo
    MassOfVap = MassOfVap*XScale**3
    if(IS_INITIAL_BALANCE_VAP) then
       BalIncrMassOfVap = 0.d0
       IS_INITIAL_BALANCE_VAP = .false.
    else
       BalIncrMassOfVap = MassOfVap - BalMassOfVap
    endif
    BalMassOfVap = MassOfVap
  end function MassOfVap

  real(kind=8) function WaterFlowRate()
    use Globals
    use GlobalArrays
    use Multiphase
    WaterFlowRate = sum(VFPH(IndLiq,2:NX1,NY1)*AlphaY(IndLiq,2:NX1,NY1)*Por_Y(2:NX1,NY1)*&
         (XB(2:NX1,NY1)-XB(1:NX,NY1))*XC(2:NX1,NY1)*0.5d0*(DensPH(IndLiq,2:NX1,NY1)+DensPH(IndLiq,2:NX1,NY2)))*2.d0*PI
    WaterFlowRate = WaterFlowRate*RoScale*VScale*XScale**2
  end function WaterFlowRate

  real(kind=8) function  VaporFlowRate()
    use Globals
    use GlobalArrays
    use Multiphase
    VaporFlowRate = sum(VFPH(IndGas,2:NX1,NY1)*AlphaY(IndGas,2:NX1,NY1)*Por_Y(2:NX1,NY1)*&
         (XB(2:NX1,NY1)-XB(1:NX,NY1))*XC(2:NX1,NY1)*0.5d0*(DensPH(IndGas,2:NX1,NY1)+DensPH(IndGas,2:NX1,NY2)))*2.d0*PI
    VaporFlowRate = VaporFlowRate*RoScale*VScale*XScale**2
  end function VaporFlowRate

  real(kind=8) function EvapRateTot()
    use Globals
    use GlobalArrays
    use Multiphase
    EvapRateTot = sum(Gamma_Evap(2:NX1,2:NY1)*&
         (XB(2:NX1,2:NY1)-XB(1:NX,2:NY1))*XC(2:NX1,2:NY1)*&
         (YB(2:NX1,2:NY1)-YB(2:NX1,1:NY)))*2.d0*PI
    EvapRateTot = EvapRateTot*XScale**3
  end function EvapRateTot

  real(kind=8) function EvapRate()
    use Globals
    use GlobalArrays
    use Multiphase
    EvapRate = sum(Gamma_Evap(2:NX1,2:NY1)*&
         (XB(2:NX1,2:NY1)-XB(1:NX,2:NY1))*XC(2:NX1,2:NY1)*&
         (YB(2:NX1,2:NY1)-YB(2:NX1,1:NY)),mask=Gamma_Evap(2:NX1,2:NY1) > 0.d0)*2.d0*PI
    EvapRate = EvapRate*XScale**3
  end function EvapRate

  real(kind=8) function CondRate()
    use Globals
    use GlobalArrays
    use Multiphase
    CondRate = sum(Gamma_Evap(2:NX1,2:NY1)*&
         (XB(2:NX1,2:NY1)-XB(1:NX,2:NY1))*XC(2:NX1,2:NY1)*&
         (YB(2:NX1,2:NY1)-YB(2:NX1,1:NY)),mask=Gamma_Evap(2:NX1,2:NY1) < 0.d0)*2.d0*PI
    CondRate = CondRate*XScale**3
  end function CondRate

  real(kind=8) function UmaxPhase(iPh)
    use Globals
    use GlobalArrays
    use Multiphase
    use Geometry
    use Water
    use PVolFrac
    use CORRELATIONS
    implicit none
    integer, intent(in):: iPh

    UMaxPhase = maxval(abs(UPH(iPh,2:NX1,2:NY1)), mask=AlphaPH(iPh,2:NX1,2:NY1) > Alpha_lim)*VScale
  end function UmaxPhase

  real(kind=8) function VmaxPhase(iPh)
    use Globals
    use GlobalArrays
    use Multiphase
    use Geometry
    use Water
    use PVolFrac
    use CORRELATIONS
    implicit none
    integer, intent(in):: iPh

    VMaxPhase = maxval(abs(VPH(iPh,2:NX1,2:NY1)), mask=AlphaPH(iPh,2:NX1,2:NY1) > Alpha_lim)*VScale
  end function VmaxPhase

  real(kind=8) function VoidMax(isPor)
    use Globals
    use GlobalArrays
    use Multiphase
    implicit none
    logical, intent(in), optional:: isPor

    logical IS_POR
    if(present(isPor)) then
       IS_POR = isPor
    else
       IS_POR = .false.
    endif

    VoidMax = maxval(AlphaPH(IndGas,2:NX1,2:NY1), mask= (.not.IS_POR .or.&
         maskPor(2:NX1,2:NY1) /= 0))
  end function VoidMax

  real(kind=8) function MomentumPhase(iPh)
    use Globals
    use GlobalArrays
    use Multiphase
    use Geometry
    use Water
    use PVolFrac
    implicit none
    integer, intent(in):: iPh

    integer:: i,j
    real(kind=8):: Uabs
    MomentumPhase = 0.d0
    do j = 2,NY1
       do i = 2,NX1
          Uabs = sqrt(UPH(iPh,i,j)**2+VPH(iPh,i,j)**2)
          MomentumPhase = MomentumPhase + &
               2.d0*PI*XC(i,j)*(XB(i,j)-XB(i-1,j))*(YB(i,j)-YB(i,j-1))*&
               DensPH(iPh,i,j)*AlphaPH(iPH,i,j)*Uabs*Porosity(i,j)
       enddo
    enddo
    MomentumPhase = MomentumPhase*XScale**3*RoScale*VScale
  end function MomentumPhase


  real(kind=8) function MassOfEvap()
    use Globals
    use GlobalArrays
    use Multiphase
    use Geometry
    use Water
    use PVolFrac
    use Turb    
    implicit none

    integer:: i,j
    MassOfEvap = 0.d0
    do j = 2,NY1
       do i = 2,NX1
          MassOfEvap = MassOfEvap + &
               Gamma_Evap(i,j)*TimeStep*&
               2.d0*PI*XC(i,j)*(XB(i,j)-XB(i-1,j))*(YB(i,j)-YB(i,j-1))
       enddo
    enddo
    MassOfEvap = MassOfEvap*XScale**3*RoScale*TimeScale
  end function MassOfEvap

  subroutine WritePhaseDrag(P, T, Ui)
    use Globals
    use Correlations
    use Multiphase
    implicit none

    real(kind=8),intent(in):: P            ! Pressure, [Pa]
    real(kind=8),intent(in):: T            ! Temperature, [K]
    real(kind=8),intent(in):: Ui(NPHASE)   ! Phase velocities, [m/s]

    integer:: i, NDRAG
    real(kind=8):: ADi(NPHASE),Ai(NPHASE)
    real(kind=8):: CdP(NPHASE) ! Liquid-porous and Gas-porous
    real(kind=8):: CdPQ(NPHASE) ! Liquid-porous and Gas-porous
    real(kind=8):: Cd12(NPHASE)! Drag coefficient Gas-Liquid
    real(kind=8):: Cd12Q(NPHASE)! Drag coefficient Gas-Liquid
    real(kind=8):: dCd12_dA(NPHASE)! dCd12/dA (A is the void fraction)
    real(kind=8):: dCd12Q_dA(NPHASE) ! dCd12Q/dA
    real(kind=8):: dCdP_dA(NPHASE)! dCdP/dA (see above)
    real(kind=8):: dCdPQ_dA(NPHASE) ! dCdPQ/dA
    logical:: isDeriv = .true.

    real(kind=8):: Por = -1.d0
    real(kind=8):: Kp = -1.d0, Kps = -1.d0
    real(kind=8):: Dp = -1.d0

    real(kind=8):: UVi(NPHASE,NDIMS), Ti(NPHASE)


    UVi(:,1) = Ui
    UVi(:,2:NDIMS) = 0.d0
    Ti = T

    NDRAG = 101
    open(UNITDRAG,file=FileDRAG)
    write(UNITDRAG,'(20(a15,2x))') 'AG[-]','AL[-]','Cd12','dCd12_dA'

    call SET_LOCAL_PARMS(P, Ti)

    do i = 1,NDRAG
       Ai(IndGas) = (i-1.d0)/(NDRAG-1.d0)
       Ai(IndLiq) = 1.d0-Ai(IndGas)
       ADi = Ai

       call PHASE_DRAG_COEFS(Ai, ADi, Por, Dp, Kp, Kps, &
            CdP, CdPQ, Cd12, Cd12Q, isDeriv, dCdP_dA, dCdPQ_dA, dCd12_dA, dCd12Q_dA)

       write(UNITDRAG,'(20(es15.6,2x))') Ai(IndGas),Ai(IndLiq),Cd12(IndLiq),dCd12_dA(IndLiq)
    enddo
    close(UNITDRAG)
  end subroutine WritePhaseDrag


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ TRANSDUCERS

  subroutine WriteTransducers
    use Globals
    use GlobalArrays
    use Multiphase
    use Level
    use Particle
    use Storage
    implicit none
    character(len=2048):: buf
    character(len=lenTRANS):: valbuf
    real(kind=4):: val
    integer:: iTR
    logical:: EXISTS

    if(.not.IS_WRITE_TRANS) return
    if(len_trim(FileTRANS)<= 0 .or. NTrans <= 0) return

    inquire(file=trim(adjustl(FileTRANS)),exist = EXISTS)
    if(.not.EXISTS) call ResetTransducerFile
    open(UNITTRANS,file=FileTRANS,position='append')
    buf = ""
    write(valbuf,formTRANS) Time*TimeScale
    buf = adjustr(valbuf)
    do iTR = 1,NTrans
       if(.not.Trans(iTR)%Active) cycle
       select case(Trans(iTR)%Name)
       case("P")
          val = P(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))*PScale
       case("Al")
          val = AlphaPH(IndLiq,Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))*PScale
       case("Ag")
          val = AlphaPH(IndGas,Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))*PScale
       case("Ul")
          val = 0.5d0*(UFPH(IndLiq,Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))+&
               UFPH(IndLiq,Trans(iTR)%Ind(1)-1,Trans(iTR)%Ind(2)))*VScale
       case("Vl")
          val = 0.5d0*(VFPH(IndLiq,Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))+&
               VFPH(IndLiq,Trans(iTR)%Ind(1),Trans(iTR)%Ind(2)-1))*VScale
       case("Ug")
          val = 0.5d0*(UFPH(IndGas,Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))+&
               UFPH(IndGas,Trans(iTR)%Ind(1)-1,Trans(iTR)%Ind(2)))*VScale
       case("Vg")
          val = 0.5d0*(VFPH(IndGas,Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))+&
               VFPH(IndGas,Trans(iTR)%Ind(1),Trans(iTR)%Ind(2)-1))*VScale
       case("Tl")
          val = TempPH(IndLiq,Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       case("Tg")
          val = TempPH(IndGas,Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))
       case("Ts")
          val = TempPor(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))*TempScale
       case("TsSup")
          val = (TempPor(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))-&
               TempSat(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2)))*TempScale
       case("Lev")
          if(IS_POOL) then
             val = zLev(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2))+YB(Trans(iTR)%Ind(1),Trans(iTR)%Ind(2)-1)
          else
             val = 0.d0
          endif
       case("MLiq")
          val = TotalMassPhase(IndLiq)
       case("MGas")
          val = TotalMassPhase(IndGas)
       case("HLiq")
          val = TotalEnthPhase(IndLiq)
       case("HGas")
          val = TotalEnthPhase(IndGas)
       case("PhChg")
          val = EvapRateTot()
       case("Evap")
          val = EvapRate()
       case("Cond")
          val = CondRate()
       case("UmaxLiq")
          val = UmaxPhase(IndLiq)
       case("UmaxGas")
          val = UmaxPhase(IndGas)
       case("VmaxLiq")
          val = VMaxPhase(IndLiq)
       case("VmaxGas")
          val = VmaxPhase(IndGas)
       case("UavLiq")
          val = MomentumPhase(IndLiq)/max(TotalMassPhase(IndLiq),1.d-5)
       case("UavGas")
          val = MomentumPhase(IndGas)/max(TotalMassPhase(IndGas),1.d-5)
       case("Amax")
          val = VoidMax()
       case("AmaxPor")
          val = VoidMax(isPor = .true.)
       case("TStep")
          val = TimeStep*TimeScale
       case("CFL")
          val = CFL
       case("TGlob")
          val = GlobalTime/60.d0
       case("TStor") ! Temperature of fluid in storage
          val = StorageTemp(Trans(iTR)%Storage)
       case("HStor") ! Enthalpy of fluid in storage
          val = StorageEnth(Trans(iTR)%Storage)
       case("LevStor") ! Level of fluid in storage
          val = StorageLev(Trans(iTR)%Storage)
       case("MStor") ! Mass of fluid in storage
          val = StorageMass(Trans(iTR)%Storage)
       case("VStor") ! Volume of fluid in storage
          val = StorageVol(Trans(iTR)%Storage)
       case("TfSupMax") ! Maximum superheat of fluid
          val = maxval((TempPH(IndLiq,2:NX1,2:NY1)-TempSat(2:NX1,2:NY1))*TempScale,mask=maskPor /= 0)
       case("TfMax") ! Maximum temperature of fluid
          val = maxval(TempPH(IndLiq,2:NX1,2:NY1)*TempScale,mask=maskPor /= 0)
       case("TgSupMax") ! Maximum superheat of gas 
          val = maxval((TempPH(IndGas,2:NX1,2:NY1)-TempSat)*TempScale,mask=maskPor /= 0)
       case("TgMax") ! Maximum temperature of gas
          val = maxval(TempPH(IndGas,2:NX1,2:NY1)*TempScale,mask=maskPor /= 0)
       case("TsMax") ! Maximum temperature of porous material
          val = maxval(TempPor*TempScale,mask=maskPor /= 0)
       case("TsSupMax") ! Maximum superheat of porous material
          val = maxval((TempPor-TempSat)*TempScale,mask=maskPor /= 0)
       case("QsLiq") ! Total power from solid to liquid
          val = TotalQsPH(IndLiq)*1.d-6
       case("QsGas") ! Total power from solid to gas
          val = TotalQsPH(IndGas)*1.d-6
       case("QsCong") ! Total power from solid to structures
          val = TotalQsCong()*1.d-6
       case("QsTotal") ! Total power from solid to phases and structures
          val = (TotalQsPH(IndLiq)+TotalQsPH(IndGas)+TotalQsCong())*1.d-6
       case("QsDecay") ! Total decay heat power
          val = TotalQSolid()*1.d-6
       case("QsBalance") ! Heat gain - Heat sink
          val = (TotalQSolid()-(TotalQsPH(IndLiq)+TotalQsPH(IndGas)+TotalQsCong()))*1.d-6
       case default
          cycle
       end select
       write(valbuf,formTRANS) val
       buf = trim(buf)//adjustr(valbuf)
    enddo
    write(UNITTRANS,'(a)') trim(buf)
    close(UNITTRANS)
  end subroutine WriteTransducers

  integer function TransducerType(TR)
    implicit none
    type(TRANSDUCER_TYPE), pointer:: TR

    TransducerType = TRANSDUCER_POINT
    select case(TR%Name)
    case("MLiq","MGas","HLiq","HGas","PhChg","Evap","Cond",&
         "Amax","AmaxPor",&
         "UmaxLiq","UmaxGas",&
         "VmaxLiq","VmaxGas",&
         "UavLiq","UavGas","TStep","TGlob","TsMax","TsSupMax",&
         "TfMax","TfSupMax","TgMax","TgSupMax",&
         "QsLiq","QsGas","QsCong","QsTotal","QsDecay","QsBalance","CFL")
       TransducerType = TRANSDUCER_INTEGRAL
    case("TStor","HStor","LevStor","MStor","VStor")
       TransducerType = TRANSDUCER_STORAGE
    end select
  end function TransducerType


  subroutine CutTransducerFile
    use Globals
    use GlobalArrays
    use Multiphase
    implicit none

    character(len=2048):: buf
    logical:: EXISTS
    real(kind=4):: t

    if(.not.IS_WRITE_TRANS) return
    if(len_trim(FileTRANS)<= 0) return
    inquire(file=trim(adjustl(FileTRANS)),exist = EXISTS)
    if(.not.EXISTS) return

    open(UNITTRANS,file=FileTRANS)
    read(UNITTRANS,'(a)',end=200) buf ! Skip header
200 continue
    do
       read(UNITTRANS,'(a)',end=300) buf ! Skip header
       read(buf,*) t
       if(t > Time*TimeScale) then
          backspace UNITTRANS
          endfile UNITTRANS
          goto 300
       endif
    enddo
300 continue
    close(UNITTRANS)
  end subroutine CutTransducerFile

  subroutine ResetTransducerFile
    use Globals
    use Multiphase
    implicit none
    character(len=2048):: buf
    character(len=lenTRANS):: valbuf
    integer:: iTR

    if(.not.IS_WRITE_TRANS) return
    if(len_trim(FileTRANS)<= 0) return
    open(UNITTRANS,file=FileTRANS)
    rewind(UNITTRANS)
    endfile UNITTRANS
    if(NTrans > 0) then
       rewind(UNITTRANS)
       buf = ""
       valbuf = "Time"
       buf = adjustr(valbuf)
       do iTR = 1,NTrans
          if(.not.Trans(iTR)%Active) cycle
          select case(Trans(iTR)%Name)
          case("P","Al","Ag","Ul","Ug","Vl","Vg","Tl","Tg","Lev",&
               "MLiq","MGas","HLiq","HGas","PhChg","Evap","Cond",&
               "Amax","AmaxPor",&
               "UmaxLiq","UmaxGas",&
               "VmaxLiq","VmaxGas",&
               "UavLiq","UavGas", &
               "TStep","CFL","TGlob",&
               "TStor","HStor","VStor","MStor","LevStor",&
               "Ts","TsSup","TsMax","TsSupMax",&
               "TfMax","TfSupMax","TgMax","TgSupMax",&
               "QsLiq","QsGas","QsCong","QsTotal","QsDecay","QsBalance")
             write(valbuf,'(a)') trim(adjustl(Trans(iTR)%Label))
             buf = trim(buf)//adjustr(valbuf)
          end select
       enddo
       write(UNITTRANS,'(a)') trim(buf)
    endif
    close(UNITTRANS)
  end subroutine ResetTransducerFile

  subroutine InitTransducers
    use Globals
    use GlobalArrays
    use Storage
    implicit none
    integer:: iTR,i,j, TType
    real(kind=8):: C(NDIMS)
    type(TRANSDUCER_TYPE), pointer:: ptr

    if(.not.IS_WRITE_TRANS) return
    TRANS_LOOP: do iTR = 1,NTrans
       if(minval(Trans(iTR)%Ind) <= 0) then ! Define indices from coordinates
          C = Trans(iTR)%Coord
          Trans(iTR)%Ind = -1
          j = 2
          do i = 2,NX1
             if(C(1) <= XB(i,j) .and. C(1) >= XB(i-1,j)) then
                Trans(iTR)%Ind(1) = i
                exit
             endif
          enddo
          i = 2
          do j = 2,NY1
             if(C(2) <= YB(i,j) .and. C(2) >= YB(i,j-1)) then
                Trans(iTR)%Ind(2) = j
                exit
             endif
          enddo
       else ! Define coordinates from indices
          Trans(iTR)%Coord = UNASSIGNED
          i = Trans(iTR)%Ind(1)
          j = Trans(iTR)%Ind(2)
          if(i>=2 .and. i<= NX1 .and. j>= 2 .and. j<= NY1) then
             Trans(iTR)%Coord(1) = XC(i,j)
             Trans(iTR)%Coord(2) = YC(i,j)
          endif
       endif
       ptr => Trans(iTR)
       TType = TransducerType(ptr)
       select case(TType)
       case(TRANSDUCER_POINT)
          Trans(iTR)%Active = minval(Trans(iTR)%Ind) > 0 .and. minval(Trans(iTR)%Coord) > IFASSIGNED
       case(TRANSDUCER_INTEGRAL)
          Trans(iTR)%Active = .true.
       case(TRANSDUCER_STORAGE)
          Trans(iTR)%Active = NStor > 0 .and. Trans(iTR)%Storage >= 0 .and. Trans(iTR)%Storage <= NStor
       end select
    enddo TRANS_LOOP
  end subroutine InitTransducers


end module OUTPUT
