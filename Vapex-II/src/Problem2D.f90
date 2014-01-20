!--------------------------------------------------------------------!
!  Problem-specific subroutines                                      !
!--------------------------------------------------------------------!
!  $Id: Problem2D.f90 9 2014-01-14 14:22:27Z Sergey $
!--------------------------------------------------------------------!
SUBROUTINE SetBoundaryConditions
  USE VARIABLES_2D
  USE GLOBALS_2D
  USE GRID_2D
  USE INOUT_2D
  IMPLICIT NONE
  !
  ! Axis and outer boundary
  !
  au1(1,:)=0.d0; au1(n,:)=0.d0
  av1(1,:)=av1(2,:); av1(n9,:)=av1(n,:)
  aa1(1,:)=aa1(2,:); aa1(n9,:)=aa1(n,:)
  ar1(1,:)=ar1(2,:); ar1(n9,:)=ar1(n,:)
  at1(1,:)=at1(2,:); at1(n9,:)=at1(n,:)
  ae1(1,:)=ae1(2,:); ae1(n9,:)=ae1(n,:)

  au2(1,:)=0.d0; au2(n,:)=0.d0
  av2(1,:)=av2(2,:); av2(n9,:)=av2(n,:)
  aa2(1,:)=aa2(2,:); aa2(n9,:)=aa2(n,:)
  ar2(1,:)=ar2(2,:); ar2(n9,:)=ar2(n,:)

  at2(1,:)=at2(2,:); at2(n9,:)=at2(n,:)
  ae2(1,:)=ae2(2,:); ae2(n9,:)=ae2(n,:)

  ats(1,:)=ats(2,:); ats(n9,:)=ats(n,:)
  atsv(1,:)=atsv(2,:); atsv(n9,:)=atsv(n,:)

  ara(1,:)=ara(2,:); ara(n9,:)=ara(n,:)
  aYa(1,:)=aYa(2,:); aYa(n9,:)=aYa(n,:)

  arH2(1,:)=arH2(2,:); arH2(n9,:)=arH2(n,:)
  aYH2(1,:)=aYH2(2,:); aYH2(n9,:)=aYH2(n,:)

  gam(1,:)=gam(2,:); gam(n9,:)=gam(n,:)

  bp(1,:)=bp(2,:); bp(n9,:)=bp(n,:)
  bpa(1,:)=bpa(2,:); bpa(n9,:)=bpa(n,:)
  bpv(1,:)=bpv(2,:); bpv(n9,:)=bpv(n,:)    

  al3(1,:)=al3(2,:); al3(n9,:)=al3(n,:)
  !
  ! Top and bottom boundaries
  !
  av1(:,1)=0.d0; 
  au1(:,1)=au1(:,2); au1(:,m9)=au1(:,m)
  aa1(:,1)=aa1(:,2); aa1(:,m9)=aa1(:,m)
  ar1(:,1)=ar1(:,2); ar1(:,m9)=ar1(:,m)
  at1(:,1)=at1(:,2); at1(:,m9)=at1(:,m)
  ae1(:,1)=ae1(:,2); ae1(:,m9)=ae1(:,m)

  av2(:,1)=0.d0; 
  au2(:,1)=au2(:,2); au2(:,m9)=au2(:,m)
  aa2(:,1)=aa2(:,2); aa2(:,m9)=aa2(:,m)
  ar2(:,1)=ar2(:,2); ar2(:,m9)=ar2(:,m)

  IF(isOpenTop == 0) THEN
     av1(:,m)=0.d0
     av2(:,m)=0.d0
  ENDIF

  at2(:,1)=at2(:,2); at2(:,m9)=at2(:,m)
  ae2(:,1)=ae2(:,2); ae2(:,m9)=ae2(:,m)

  ats(:,1)=ats(:,2); ats(:,m9)=ats(:,m)
  atsv(:,1)=atsv(:,2); atsv(:,m9)=atsv(:,m)

  ara(:,1)=ara(:,2); ara(:,m9)=ara(:,m)
  aYa(:,1)=aYa(:,2); aYa(:,m9)=aYa(:,m)

  arH2(:,1)=arH2(:,2); arH2(:,m9)=arH2(:,m)
  aYH2(:,1)=aYH2(:,2); aYH2(:,m9)=aYH2(:,m)

  gam(:,1)=gam(:,2); gam(:,m9)=gam(:,m)

  bp(:,1)=bp(:,2); 
  bpa(:,1)=bpa(:,2); bpa(:,m9)=bpa(:,m)
  bpv(:,1)=bpv(:,2); bpv(:,m9)=bpv(:,m)    

  al3(:,1)=al3(:,2); al3(:,m9)=al3(:,m)
  !
  ! Connection to freeboard
  !
  IF(isOpenSide /= 0) THEN
  ENDIF
  !
  ! Open top boundary
  !
  IF(isOpenTop == 0) THEN
     bp(:,m9)=bp(:,m)
  ENDIF
END SUBROUTINE SetBoundaryConditions

SUBROUTINE BeforeFirstStep
  USE GLOBALS
  USE GLOBALS_2D
  USE VARIABLES_2D
  USE INOUT_2D
  USE INITIAL_2D
  USE SOLVER_2D
  USE GRID_2D
  USE GRAPHICS_TEC2D
  USE GRAPHICS_PDB2D
  USE PROBLEM_DATA
  USE DISPERSED_PHASE
  IMPLICIT NONE
  !
  ! Initial balances
  !
  IF(NREAD == 0) THEN
     CALL InitialBalances

     IF(isOUT /= 0) THEN

        ! Write the initial data
        OPEN(11,file='#out/zmin.dat',position='rewind')
        OPEN(14,file='#out/level.dat',position='rewind');
        OPEN(27,file='#out/qmpart.dat',position='rewind')
        WRITE(27,100) 'Time','qmpart','qmjet','qmdrops','qmdebris','qmzero',&
             'qmpart0','Balance qm','amas3'
        WRITE(27,5107) time,qmpart,qmjet,qmdrops,qmdebris,qmzero,qmpart0,&
             qmpart-(qmjet+qmdrops+qmdebris+qmzero),amas30

        OPEN(28,file='#out/energy.dat',position='rewind')
        WRITE(28,5107) time, ener10/EnergyScale,&
             (ener20+ar2FB*ae2FB*VolumeFB)/EnergyScale,&
             enermix0/EnergyScale,ener30/EnergyScale,&
             0.D0,0.D0,0.D0,&
             0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0


        OPEN(30,file='#out/npart.dat',position='rewind')
        OPEN(31,file='#out/freeboard.dat',position='rewind')

        OPEN(32,file='#out/mass.dat',position='rewind')
        WRITE(32,100) 'Time','MassLiq','MassGas','MassVapH2O',&
             'Mass_Ar','Mass_H2',&
             'MassTot_H2O','MassTot_Ar','MassTot_H2',&
             'DM_H2O','DM_Ar','DM_H2'
        WRITE(32,5107) time,amas10,amas20,amasH2O_0,amasAr_0,amasH2_0,& ! Initial masses
             amasH2O_total0,amasAr_total0,amasH2_total0,& ! Total masses
             0.D0,0.D0,0.D0 ! Disbalances (initially zero)


        OPEN(33,file='#out/mtransfer.dat',position='rewind')
        OPEN(34,file='#out/minmax.dat',position='rewind')
        WRITE(34,100) 'Time','U1max','U2max', &
             'T1max','T1min','T2max','T2min',&
             'Al3max','AlphaAv','(tl-ts)max','Pmax','Pmin'
        WRITE(34,5107) time, vlmax,vvmax,tlmax,tlmin, &
             tvmax,tvmin,al3max, alphaAv,tl_ts,ptotmax,ptotmin


        OPEN(35,file='#out/frate.dat',position='rewind')
        WRITE(11,5107) time,zmin,zminj,zmindr
        WRITE(14,5107) time,WaterLevel,WaterLevel,0.d0,&
             WaterLevel*1.d3,WaterLevel*1.d3
        WRITE(30,5108) time,ll,lljet,lldrops,lldebris,llzero,nump,&
             ll-(lljet+lldrops+lldebris+llzero)
        WRITE(31,5107) time,apFB/PScale,apaFB/PScale,(apFB-apaFB)/PScale, &
             apH2FB/PScale,at2FB-273.15d0,ar2FB,aYH2FB,aXH2FB,&
             amasAr_fb,&  ! Mass of argone
             amasH2_fb,&  ! Mass of hydrogen
             amasH2O_fb   ! Mass of steam
        WRITE(33,5107) time, agam, agamH2, qgam, qgamH2
        WRITE(35,5107) time, frate, FragMass, FragMass/qmpart0,FRTotal,fragL_D

        IF(isOpenTop /= 0) THEN
           OPEN(150,file='#out/topflux.dat',position='rewind')
           WRITE(150,5107) time, topflux_w, topflux_g, topflux_a, topflux_v
           CLOSE(150,status='keep')
        ENDIF

        CLOSE(11,status='keep')
        CLOSE(14,status='keep')
        CLOSE(27,status='keep')
        CLOSE(28,status='keep')
        CLOSE(30,status='keep')
        CLOSE(31,status='keep')
        CLOSE(32,status='keep')
        CLOSE(33,status='keep')
        CLOSE(34,status='keep')
        CLOSE(35,status='keep')

        timeOUT = timeOUT + dtOUT
     ENDIF

!!!DEBUG!!!
!!$     open(77,file='H2.dat')
!!$     write(77,*) time,aYH2(2,2),arH2(2,2)
!!$     close(77)
  ENDIF
  !
  !
  IF(NREAD == 0) THEN
     FileRST = 'none'; FileRST(5:LEN(FileRST)) = ' '
  ENDIF
  IF(FileMsgLevel > 0) THEN
     OPEN(MLOGFILE,file='#out/vapex.out',status='unknown', &
          access='sequential',position='rewind', &
          form='formatted')
     CALL Salute(MLOGFILE)
  ENDIF
  IF(MsgLevel>0)  THEN
     CALL Salute(6)
     WRITE(*,2003) date_time(3),&
          date_time(2),date_time(1),date_time(5),date_time(6),date_time(7),&
          TestName(1:LEN(TRIM(TestName))),&
          FileRST(1:LEN(TRIM(FileRst)))
  ENDIF
  IF(FileMsgLevel>0)  THEN
     WRITE(MLOGFILE,2003) date_time(3),&
          date_time(2),date_time(1),date_time(5),date_time(6),date_time(7),&
          TestName(1:LEN(TRIM(TestName))),&
          FileRST(1:LEN(TRIM(FileRst)))
  ENDIF
  !
  ! Graphics output
  !
  IF(isTEC /= 0 .AND. time >= timeTEC) THEN
     CALL WriteTECPLOT(NPictTEC, Time, aa2, al3,al3dr,al3j,&
          au1, av1, au2, av2, &
          bp, at1,at2,at3,gam,aYH2,ara,gamH2, &
          dn, dm, N, M, N9, M9,MsgLevel,FileMsgLevel,MLOGFILE)
     CALL WriteTECPLOTMesh(dxu,dzv,n,m)
     timeTEC = timeTEC + dtTEC
     NPictTEC = NPictTEC + 1
  ENDIF
  IF(isTECDisp /= 0 .AND. time >= timeTECDisp) THEN
     !     call WriteTECPLOTDisp(NPictTECDisp, Time,x,z,xu3,xv3,xt3,ll)
     timeTECDisp = timeTECDisp + dtTECDisp
     NPictTECDisp = NPictTECDisp + 1
  ENDIF
  IF(isPDB /= 0 .AND. time >= timePDB) THEN
     CALL WritePDB(NPictPDB, Time, aa1, aa2, al3,al3j,al3dr,al3de,&
          au1, av1, au2, av2, &
          ae1, ae2, aYa, &
          bp, bpv, at1, at2, at3, ats, gam, aYH2, ar2, gamH2, & !r31c, r32c,&
          dn, dm, N, M, N9, M9, MsgLevel, FileMsgLevel,MLOGFILE)
     timePDB = timePDB + dtPDB
     NPictPDB = NPictPDB + 1
  ENDIF
  IF(isPDBDisp /= 0 .AND. time >= timePDB) THEN
!!$     call  WritePDBDisp(NPictPDBDisp, Time, dn,dm,n9,m9,x, z, xu3, xv3, xt3, ind, ll,&
!!$          MsgLevel, FileMsgLevel)
     timePDBDisp = timePDBDisp + dtPDBDisp
     NPictPDBDisp = NPictPDBDisp + 1
  ENDIF
  !
  ! Adaptive time step counter
  !
  i_tsteps = 0


1 FORMAT(80('-'))

2003 FORMAT(' Date: ',i2.2,'/',i2.2,'/',i4.4,&
       '  | Time: ',i2.2,':',i2.2,':',i2.2,3x,'| Test: ',a,t61,' |   Restart: ',a)
5107 FORMAT(24(d12.5,' '))
5108 FORMAT(d12.5,10(1x,i8))

100 FORMAT(24(a12,' '))

END SUBROUTINE BeforeFirstStep


SUBROUTINE BeforeEachStep
  USE GLOBALS
  USE GLOBALS_2D
  USE VARIABLES_2D
  USE INOUT_2D
  USE INITIAL_2D
  USE SOLVER_2D
  USE GRID_2D
  USE GRAPHICS_TEC2D
  USE GRAPHICS_PDB2D
  USE PROBLEM_DATA
  USE DISPERSED_PHASE

  IMPLICIT NONE
100 FORMAT(80('#'))  
200 FORMAT('  TIME = ',ES10.3,1x,'|', &
       ' T = ',ES10.3,1x,'|',1x,'Int: ',i4,' (',i4,')',1x, &
       '|',1x,'Ext: ',i4,' (',i4,')')
300 FORMAT(' IT   Del-Ui',3X, &
       'Del-Vi',3X, &
       'Del-Al',3X,'Del-T1',3X, &
       'Del-T2',3X,'Del-DP',3X,'DelDPA',3x,'Tol-BCG'/80('-'))

  IF(MsgLevel > 0) THEN
     WRITE(*,100)
     WRITE(*,200) time,tstep, iinner, NInner, iouter,NOuter
     IF(MsgLevel>1) THEN
        WRITE(*,100)
        WRITE(*,300)
     ENDIF
  ENDIF
  IF(FileMsgLevel > 0) THEN
     WRITE(MLOGFILE,100)
     WRITE(MLOGFILE,200) time,tstep, iinner, NInner, iouter,NOuter
     IF(FileMsgLevel>1) THEN
        WRITE(MLOGFILE,100)
        WRITE(MLOGFILE,300)
     ENDIF
  ENDIF

END SUBROUTINE BeforeEachStep


SUBROUTINE AfterOneStep
  USE GLOBALS
  USE GLOBALS_2D
  USE VARIABLES_2D
  USE INITIAL_2D
  USE INOUT_2D
  USE SOLVER_2D
  USE GRID_2D
  USE GRAPHICS_TEC2D
  USE GRAPHICS_PDB2D
  USE PROBLEM_DATA
  USE DISPERSED_PHASE
  USE OUTTRANS
  IMPLICIT NONE
  INTEGER::i
  !
  ! Balances
  !
  CALL CurrentBalances

  IF(isOUT /= 0 .AND. time > timeOUT) THEN
     !
     ! Write balances
     !
     OPEN(11,file='#out/zmin.dat',position='append')
     OPEN(27,file='#out/qmpart.dat',position='append')
     OPEN(28,file='#out/energy.dat',position='append')
     OPEN(30,file='#out/npart.dat',position='append')
     OPEN(31,file='#out/freeboard.dat',position='append')
     OPEN(32,file='#out/mass.dat',position='append')
     OPEN(33,file='#out/mtransfer.dat',position='append')
     OPEN(34,file='#out/minmax.dat',position='append')
     OPEN(35,file='#out/frate.dat',position='append')
     WRITE(11,5107) time,zmin,zminj,zmindr
     WRITE(27,5107) time,qmpart,qmjet,qmdrops,qmdebris,qmzero,qmpart0,&
          qmpart-(qmjet+qmdrops+qmdebris+qmzero),amas3
     WRITE(28,5107) time, &     ! (1)
          ener1/EnergyScale,&   ! (2)
          (ener2+ar2FB*ae2FB*VolumeFB)/EnergyScale,&   ! (3)
          enermix/EnergyScale,& ! (4)
          ener3/EnergyScale,&   ! (5)
          quenchRate/EnergyScale,& ! (6)
          quenchEnergy/EnergyScale,& !(7)
          (enermix-quenchEnergy-enermix0)/enermix0,& !(8) - balance
          (ener1-ener10)/EnergyScale,& ! Energy gained by water (9)
          (ener2-ener20)/EnergyScale,& ! Energy gained by gas (10)
          ener3disp/EnergyScale,& ! Should coincide with ener3/EnergyScale (11)
          (qmpart*cmelt*Tin-ener3disp)/EnergyScale,& !Quenched energy<=>quenchEnergy/EnergyScale (12)
          quenchEnergy/(qmpart0*cmelt*(Tin-Tliq0)),& !Quenched energy as a % of total melt energy (13)
          qeWater/EnergyScale,qeVapour/EnergyScale,& ! (14),(15)
          enerd/EnergyScale,enerd1/EnergyScale,-enerd2/EnergyScale !(16),(17),(18)
     WRITE(30,5108) time,ll,lljet,lldrops,lldebris,llzero,nump,&
          ll-(lljet+lldrops+lldebris+llzero)
     WRITE(31,5107) time,apFB/PScale,apaFB/PScale,(apFB-apaFB)/PScale, &
          apH2FB/PScale,at2FB-273.15,ar2FB,aYH2FB,aXH2FB,&
          amasAr_fb,&  ! Mass of argone
          amasH2_fb,&  ! Mass of hydrogen
          amasH2O_fb   ! Mass of steam
     WRITE(32,5107) time,amas1,amas2,amasH2O,amasAr,amasH2,& ! Current masses
          amasH2O_total,amasAr_total,amasH2_total,& ! Total masses
          (amasH2O_total-amasH2O_total0)/amasH2O_total0,& ! Disbalances of H2O and Ar
          (amasAr_total-amasAr_total0)/amasAr_total0,&
          amasH2_total-qgamH2
     WRITE(33,5107) time, agam, agamH2, qgam, qgamH2
     WRITE(34,5107) time, vlmax,vvmax,tlmax,tlmin, &
          tvmax,tvmin,al3max,alphaAv,tl_ts,ptotmax,ptotmin
     WRITE(35,5107) time, FRate, FragMass, FragMass/qmpart0,FRtotal,fragL_D

     IF(isOpenTop /= 0) THEN
        OPEN(150,file='#out/topflux.dat',position='append')
        WRITE(150,5107) time, topflux_w, topflux_g, topflux_a, topflux_v
        CLOSE(150,status='keep')
     ENDIF

     CLOSE(11,status='keep')
     CLOSE(27,status='keep')
     CLOSE(28,status='keep')
     CLOSE(30,status='keep')
     CLOSE(31,status='keep')
     CLOSE(32,status='keep')
     CLOSE(33,status='keep')
     CLOSE(34,status='keep')
     CLOSE(35,status='keep')
5107 FORMAT(24(d12.5,' '))
5108 FORMAT(d12.5,10(1x,i8))

     IF(isThirdPhase /= 0) THEN
        OPEN(35,file='#out/disp.dat')
        WRITE(35,5107) time
        DO i = 1,ll
           WRITE(35,5107) float(i),float(ind(i)),&
                x(i),z(i),xu3(i),xv3(i),xd3(i),xt3(i)
        ENDDO
        CLOSE(35)
     ENDIF
     CALL WriteTransducers
     timeOUT = timeOUT + dtOUT
  ENDIF
  !
  ! Graphics output
  !
  IF(isTEC /= 0 .AND. time >= timeTEC) THEN
     CALL WriteTECPLOT(NPictTEC, Time, aa2, al3,al3dr,al3j,&
          au1, av1, au2, av2, &
          bp, at1,at2,at3,gam,aYH2,ara,gamH2, &
          dn, dm, N, M, N9, M9, MsgLevel, FileMsgLevel,MLOGFILE)
     timeTEC = timeTEC + dtTEC
     NPictTEC = NPictTEC + 1
  ENDIF
  IF(isTECDisp /= 0 .AND. time >= timeTECDisp) THEN
     !     call WriteTECPLOTDisp(NPictTECDisp, Time,x,z,xu1,xu2,xt3,ll)
     timeTECDisp = timeTECDisp + dtTECDisp
     NPictTECDisp = NPictTECDisp + 1
  ENDIF
  IF(isPDB /= 0.AND. time >= timePDB) THEN
     CALL WritePDB(NPictPDB, Time, aa1, aa2, al3,al3j,al3dr,al3de,&
          au1, av1, au2, av2, &
          ae1, ae2, aYa, &
          bp, bpv, at1, at2, at3, ats, gam, aYH2, ar2, gamH2, & 
          dn, dm, N, M, N9, M9, MsgLevel, FileMsgLevel,MLOGFILE)
     timePDB = timePDB + dtPDB
     NPictPDB = NPictPDB + 1
  ENDIF
  IF(isPDBDisp /= 0.AND. time >= timePDB) THEN
!!$     call  WritePDBDisp(NPictPDBDisp, Time, dn,dm,n9,m9,x, z, xu1, xu2, xt3, ind, ll,&
!!$          MsgLevel, FileMsgLevel)
     timePDBDisp = timePDBDisp + dtPDBDisp
     NPictPDBDisp = NPictPDBDisp + 1
  ENDIF
END SUBROUTINE AfterOneStep

SUBROUTINE AfterLastStep
  USE INOUT_2D
  IMPLICIT NONE
  IF(FileMsgLevel>0) CLOSE(MLOGFILE)
  IF(MsgLevel>0) PRINT *,'Bye...'
END SUBROUTINE AfterLastStep

SUBROUTINE SetEquationCoefficients
  USE VARIABLES_2D
  USE GLOBALS_2D
  USE CONVECTION_2D
  USE LOCAL_VARS
  USE SOLVER_2D
  USE COEFFICIENTS
  USE WATER_PROPS
  USE EVAPORATION
  USE CORRELATIONS
  USE HYDROGEN
  USE PROBLEM_DATA
  USE DISPERSED_PHASE
  IMPLICIT NONE

  INTEGER:: i,j
  INTEGER:: info, ipiv(NVARS)
  !
  !
  !
  DO j = 2,m
     DO i = 2,n

        CALL SetPreviousValues(Alpha1_n = ba1(i,j), &
             Alpha2_n = ba2(i,j), &
             TLiq_n = bt1(i,j), &
             TGas_n = bt2(i,j), &
             ELiq_n = be1(i,j), &
             EGas_n = be2(i,j), &
             Ro1_n  = br1(i,j), &
             Ro2_n  = br2(i,j), &
             Ya_n   = bYa(i,j))

        CALL SetLocalVars(Alpha1 = aa1(i,j),&
             Alpha2 = aa2(i,j), &
             Alpha3 = al3(i,j), &
             TLiq = at1(i,j), &
             TGas = at2(i,j), &
             TDisp = at3(i,j),&
             Ptot = bp(i,j), &
             Pncon = bpa(i,j),&
             Y_a = aYa(i,j))

        IF(isThirdPhase /= 0) THEN
           Alpha_3n = bl3(i,j)
        ELSE
           Alpha_3n = 0.D0
        ENDIF

        YH2 = aYH2(i,j)

        VLiq(1) = 0.5*(au1(i,j)+au1(i-1,j))
        VLiq(2) = 0.5*(av1(i,j)+av1(i,j-1))
        VGas(1) = 0.5*(au2(i,j)+au2(i-1,j))
        VGas(2) = 0.5*(av2(i,j)+av2(i,j-1))

        CALL Thermo_Props
        CALL SetLocalProps
        CALL SetLiqGasParms

        CALL HeatExchangeGasLiq

        IF(isEvaporation /= 0) THEN
           CALL Evap_Rate
           gam(i,j) = Gamma_Evap
        ENDIF

        IF(isDispersedPhase /= 0) THEN
           CALL SetLocalHeatExchangeDisp(R31c(i,j),R32c(i,j))
        ENDIF

        CALL SetHydrogenRelease(G_a = gamH2(i,j), &
             Y_H2 = aYH2(i,j), T_H2 = T_H2)

        CALL EQUATION_COEFS(tstep)
        CALL ADD_CONVECTIVE_TERMS(tstep,i,j)

        GMatrix = TimeJac + ConvJac - SourceJac
        QCoefMatrix(:,1,i,j) = RHS(:)
        QCoefMatrix(:,2:5,i,j) = PressJac

        ! Solve system of linear equations to find equation coefficients
        CALL DGESV(NVARS,NRHS,GMatrix,NVARS,ipiv,QCoefMatrix(1,1,i,j), &
             NVARS,info)
        IF(info /= 0) THEN
           PRINT *,'Error in DGESV: info = ',info,' at (',i,',',j,')'
           STOP
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE SetEquationCoefficients

SUBROUTINE SolveForVarIncrements
  USE GLOBALS_2D
  USE VARIABLES_2D
  USE ELLEQ
  USE INOUT_2D
  USE SOLVER_2D
  USE GRID_2D
  IMPLICIT NONE

  INTEGER:: iVar, mxmv
  REAL(kind=8):: tol
  !  external matvec_bcg2

  INTEGER i,j
  REAL(kind=8):: res

  dp = 0.d0
  !
  ! Pressure boundary conditions
  !
  r1x=1.d0; q1x=0.d0; rnx=1.d0; qnx=0.d0 ! Impermeable walls
  r1y=1.d0; q1y=0.d0; 

  IF(isOpenTop == 0) THEN
     rmy=1.d0; qmy=0.d0
  ELSE
     rmy=-1.D0; qmy = 0.D0
  ENDIF

  CALL SetPressureRHS(&
       QCoefMatrix,dpRHS,&
       r1x,q1x,rnx,qnx,&
       r1y,q1y,rmy,qmy,&
       NVARS,NRHS,n,m,n9,m9)
  work_bcg2 = 0
  tol = tol_bcg2
  mxmv = mxmv_bcg2
  nonzero_bcg2 = .FALSE.

  CALL BiCGSTAB2 (okprint_bcg2,L_bcg2, n_bcg2, &
       dp, dpRHS,           & ! x and rhs
       matvec_bcg2, nonzero_bcg2, tol, &
       typestop,mxmv, work_bcg2, ldw_bcg2, info_bcg2)

  ResidualBCGS = tol
  IF(MsgLevel > 2) THEN
     IF(info_bcg2 /= 0) THEN
        WRITE(*,*) 'BiCGSTAB2 : info = ',info_bcg2, &
             ' Tol(reached) = ',tol
     ENDIF
  ENDIF
  IF(FileMsgLevel > 2) THEN
     IF(info_bcg2 /= 0) THEN
        WRITE(MLOGFILE,*),'BiCGSTAB2 : info = ',info_bcg2, &
             ' Tol(reached) = ',tol
     ENDIF
  ENDIF

  CALL SetPressureBoundaryConditions(dp,&
       r1x,q1x,rnx,qnx,&
       r1y,q1y,rmy,qmy,&
       n,m,n9,m9)
  !
  ! Linear elliptic solver should go here
  !
  dVars = 0.D0
  DO iVar = 1,NVARS
     dVars(iVar,2:n,2:m) = QCoefMatrix(iVar,1,2:n,2:m) + &
          QCoefMatrix(iVar,2,2:n,2:m)*dp(1:n-1,2:m) + &
          QCoefMatrix(iVar,3,2:n,2:m)*dp(3:n9,2:m)  + &
          QCoefMatrix(iVar,4,2:n,2:m)*dp(2:n,1:m-1) + &
          QCoefMatrix(iVar,5,2:n,2:m)*dp(2:n,3:m9)
  ENDDO
!!!DEBUG!!! Here dVars(4,:,:) should coincide with dp!!!
  !
  ! Calculate the increments of phase velocities
  !
  dVels = 0.D0
  dVels(1,1:n,2:m) = au1(1:n,2:m) - bu1(1:n,2:m)&
       -(fu1r(1:n,2:m)*dP(2:n9,2:m)-fu1l(1:n,2:m)*dP(1:n,2:m))
  dVels(3,1:n,2:m) = au2(1:n,2:m) - bu2(1:n,2:m)&
       -(fu2r(1:n,2:m)*dP(2:n9,2:m)-fu2l(1:n,2:m)*dP(1:n,2:m))

  dVels(2,2:n,1:m) = av1(2:n,1:m) - bv1(2:n,1:m)&
       -(fv1r(2:n,1:m)*dP(2:n,2:m9)-fv1l(2:n,1:m)*dP(:,1:m))
  dVels(4,2:n,1:m) = av2(2:n,1:m) - bv2(2:n,1:m)&
       -(fv2r(2:n,1:m)*dP(2:n,2:m9)-fv2l(2:n,1:m)*dP(:,1:m))

END SUBROUTINE SolveForVarIncrements

SUBROUTINE CheckIncrements(info)
  USE SOLVER_2D
  USE GLOBALS_2D
  USE VARIABLES_2D
  USE WATER_PROPS
  USE INOUT_2D

  IMPLICIT NONE
  INTEGER,INTENT(out):: info
  INTEGER:: i,j
  REAL(kind=8):: P,T1,T2,deltaP,Pnew
  REAL(kind=8):: dV(NVARS)

  INTEGER:: iVar

  info = GoOn

  DO iVar = 1,NVARS
     ResidualNL(iVar) = MAXVAL(ABS(dVars(iVar,2:n,2:m)))
  ENDDO
  ResidualNL(2) = &
       ResidualNL(2)/MAXVAL(at1(2:n,2:m))
  ResidualNL(3) = &
       ResidualNL(3)/MAXVAL(at2(2:n,2:m))
  ResidualNL(4) = &
       ResidualNL(4)/MAXVAL(bp(2:n,2:m))
  ResidualNL(5) = &
       ResidualNL(5)/MAXVAL(bpa(2:n,2:m))

  DO iVar = 1,3,2
     ResidualUV(iVar) = MAXVAL(ABS(dVels(iVar,1:n,2:m)))
  ENDDO
  DO iVar = 2,4,2
     ResidualUV(iVar) = MAXVAL(ABS(dVels(iVar,2:n,1:m)))
  ENDDO

  IF(MsgLevel >= 2) THEN
     WRITE(*,400) IterNL,MAXVAL(ResidualUV(1:3:2)), &
          MAXVAL(ResidualUV(2:4:2)), &
          ResidualNL(1:NVARS),&
          ResidualBCGS
  ENDIF

400 FORMAT(I3,1x,12(ES9.2))


  IF(MAXVAL(ResidualNL) <  TolNonLinear .AND.&
       MAXVAL(ResidualUV) <  TolNonLinear) THEN
     info = Converged
     RETURN
  ENDIF

  DO j = 2,m
     DO i = 2,n
        dV(:) = dVars(:,i,j)
        P = bp(i,j)
        T1 = at1(i,j)
        T2 = at2(i,j)
        deltaP = dV(4)
        Pnew = P + deltaP

        IF(Pnew >= EOSPMax() .OR. Pnew <= EOSPMin()) info = BadTimeStep
        IF(T1+dV(2) >= EOSTLiqMax() .OR. &
             T1+dV(2) <= EOSTLiqMin()) info = BadTimeStep
        IF(T2+dV(3) >= EOSTVapMax() .OR. &
             T2+dV(3) <= EOSTVapMin()) info = BadTimeStep
     ENDDO
  ENDDO

END SUBROUTINE CheckIncrements


SUBROUTINE Current2PreviousTimeStep
  USE VARIABLES_2D
  USE HYDROGEN
  USE PROBLEM_DATA
  IMPLICIT NONE
  !
  ! Assign values from previous time step u1n ...  
  !
  u1n=au1; u2n=au2
  v1n=av1; v2n=av2

  !  br1 = ar1; br2 = ar2
  !  bt1 = at1; bt2 = at2

  !  be1 = ae1; be2 = ae2
  !  ba1 = aa1; ba2 = aa2
  !  bl3 = al3
  !
  !  bpn = bp; bpan = bpa
  !  
  !  brH2 = arH2
  !  bYa = aYa
  !
  !  br2FB = ar2FB
  !  braFB = araFB
  !  bpFB = apFB
  !  be2FB = ae2FB
  !  bYH2FB = aYH2FB
  !  brH2FB = arH2FB
  !
END SUBROUTINE Current2PreviousTimeStep


SUBROUTINE Current2PreviousIteration
  USE VARIABLES_2D
  USE HYDROGEN
  USE PROBLEM_DATA
  IMPLICIT NONE
  !
  ! Assign values from previous iteration to bu1....
  !
  bu1=au1; bu2=au2
  bv1=av1; bv2=av2

  br1 = ar1; br2 = ar2
  bt1 = at1; bt2 = at2

  be1 = ae1; be2 = ae2
  ba1 = aa1; ba2 = aa2
  bl3 = al3

  bpn = bp; bpan = bpa

  brH2 = arH2
  bYa = aYa

  br2FB = ar2FB
  braFB = araFB
  bpFB = apFB
  be2FB = ae2FB
  bYH2FB = aYH2FB
  brH2FB = arH2FB

END SUBROUTINE Current2PreviousIteration


SUBROUTINE Previous2Current
  USE COEFFICIENTS
  USE VARIABLES_2D
  USE HYDROGEN
  USE PROBLEM_DATA
  IMPLICIT NONE
  !
  ! Set n-th time step values equal current
  !
  au1=bu1; au2=bu2
  av1=bv1; av2=bv2

  ar1 = br1; ar2 = br2
  at1 = bt1; at2 = bt2

  ae1 = be1; ae2 = be2
  aa1 = ba1; aa2 = ba2
  al3 = bl3

  arH2 = brH2
  aYa =  bYa

  bp =bpn; bpa = bpan

  ar2FB = br2FB
  araFB = braFB
  apFB = bpFB
  ae2FB = be2FB
  aYH2FB = bYH2FB
  arH2FB = brH2FB

  isUseDeltaAlpha3 = 1.D0

END SUBROUTINE Previous2Current

SUBROUTINE InitialGuess
  USE VARIABLES_2D
  USE GLOBALS_2D
  USE GLOBALS
  USE HYDROGEN
  USE PROBLEM_DATA
  USE WATER_PROPS
  USE INITIAL_2D
  USE COEFFICIENTS
  IMPLICIT NONE

  INTEGER:: i, j
  IF(isEvaporation /= 0) THEN
     ! Set phase equilibrium (saturated pressure)
     DO j = 2,m
        DO i = 2,n
           IF(aa2(i,j)/MAX(aa1(i,j)+aa2(i,j),allim)<0.7D0) THEN
              at2(i,j) = at1(i,j)
              bpa(i,j) = bp(i,j) - SaturatedPressure(at2(i,j))
              bpa(i,j) = MAX(bpa(i,j),0.D0)
           ENDIF
        ENDDO
     ENDDO
     CALL UpdatePhaseVars
     CALL SetBoundaryConditions
  ENDIF
  !
  ! Set flag for dispersed phase
  !
  isUseDeltaAlpha3 = 1.D0
  !
END SUBROUTINE InitialGuess

SUBROUTINE UpdateVariables
  USE GLOBALS_2D
  USE VARIABLES_2D
  USE INITIAL_2D
  USE SOLVER_2D
  USE COEFFICIENTS
  IMPLICIT NONE

  INTEGER:: iVar,i,j
  REAL(kind=8):: Omega,v1,v2

  ! Underrelaxation
  Omega = MIN(1.D0,MINVAL(0.1*bp(2:n,2:m)/ABS(dVars(4,2:n,2:m))))
  dVars = dVars*Omega

!!$  do j = 2, m
!!$     do i = 1,n
!!$        if(dVars(1,i,j) > 0) then
!!$           aa1(i,j) = aa1(i,j)-dVars(1,i,j)
!!$           aa1(i,j) = min(max(aa1(i,j),0.D0),1.D0)
!!$           aa2(i,j) = 1.D0 - al3(i,j) - aa1(i,j)
!!$           aa2(i,j) = min(max(aa2(i,j),0.D0),1.D0)
!!$        else
!!$           aa2(i,j) = aa2(i,j)+dVars(1,i,j)
!!$           aa2(i,j) = min(max(aa2(i,j),0.D0),1.D0)
!!$           aa1(i,j) = 1.D0 - al3(i,j) - aa2(i,j)
!!$           aa1(i,j) = min(max(aa1(i,j),0.D0),1.D0)
!!$        endif
!!$     enddo
!!$  enddo
!!$
!!$  aa2(:,:) = aa2(:,:) + dVars(1,:,:)
!!$  aa2 = min(max(aa2,0.D0),1.D0)

!!$  aa2 = aa2 + dVars(1,:,:)
!!$  aa1 = aa1 - dVars(1,:,:)

  aa2(:,:)=aa2(:,:)+dVars(1,:,:)
  aa2(:,:) = MAX(0.D0,MIN(aa2(:,:),1.D0-al3(:,:)))

  aa1(:,:)=1.d0 - al3(:,:) - aa2(:,:)
  aa1(:,:) = MAX(0.D0,MIN(aa1(:,:),1.D0-al3(:,:)))


  at1(:,:) = at1(:,:) + dVars(2,:,:)
  at2(:,:) = at2(:,:) + dVars(3,:,:)
  bp(:,:)  = bp(:,:)  + dVars(4,:,:)
  bpa(:,:) = bpa(:,:) + dVars(5,:,:)

  WHERE(bpa<0.D0) bpa=0.D0

  au1(1:n,2:m) = bu1(1:n,2:m) + dVels(1,1:n,2:m)
  av1(2:n,1:m) = bv1(2:n,1:m) + dVels(2,2:n,1:m)
  au2(1:n,2:m) = bu2(1:n,2:m) + dVels(3,1:n,2:m)
  av2(2:n,1:m) = bv2(2:n,1:m) + dVels(4,2:n,1:m)

  CALL UpdatePhaseVars
  CALL SetBoundaryConditions
  !
  ! Update velocities on current iteration
  !
  bu1 = au1; bv1 = av1
  bu2 = au2; bv2 = av2
  !
  isUseDeltaAlpha3 = 0.D0 ! Increment of Alpha3 is accounted for 
  ! on the 1st iteration only!

END SUBROUTINE UpdateVariables


SUBROUTINE SetLiqGasParms
  USE LOCAL_VARS
  USE CORRELATIONS
  IMPLICIT NONE

  SurfaceTension = 22.5d-03 ! P=50 bar T=Tsat

  Lambda_Liq = 64.4d-2

  Visc_Gas = 2.27d-05
  Lambda_Gas = 1.77d-01
  Cp_Gas = 0.521d4
END SUBROUTINE SetLiqGasParms

SUBROUTINE SetInitialValues
  USE LOCAL_VARS
  USE WATER_PROPS
  USE HYDROGEN
  USE VARIABLES_2D
  USE SOLVER_2D
  USE GLOBALS_2D
  USE GRID_2D
  USE PROBLEM_DATA
  USE INOUT_2D
  USE EVAPORATION
  USE CORRELATIONS
  USE DISPERSED_PHASE
  IMPLICIT NONE

  INTEGER:: i,j
  REAL(kind=8):: P_top,P_j,P_a, P_j0, Ro1_j, Ro2_j
  REAL(kind=8):: TLiq, TGas

  REAL(kind=8):: check_hydro(m)

  !
  ! Time and time integrals
  !
  time = 0.D0
  time_disp = 0.D0
  ncount_frag = 1
  !
  ! Counters for graphics files
  !
  timePDB = 0.D0
  timePDBDisp = 0.D0
  NPictPDB = 0
  NPictPDBDisp = 0

  timeTEC = 0.D0
  timeTECDisp = 0.D0
  NPictTEC = 0
  NPictTECDisp = 0

  timeOUT= 0.D0
  !-------------------- Velocities ----------------------------
  au1 = 0.d0; bu1 = 0.d0; av1 = 0.d0; bv1 = 0.d0
  au2 = 0.d0; bu2 = 0.d0; av2 = 0.d0; bv2 = 0.d0

  !-------------------- Liquid properties ---------------------
  YH2 = 0.D0
  CALL Thermo_Props(TLiq0, TGas0, P0, 0.D0)
  RoLiq0 = Ro_1
  PvTLiq0 = satprs(TLiq0)   ! Saturated vapour pressure at Tliq0
  PvTGas0 = satprs(TGas0)   ! Saturated vapour pressure at Tvap0

  IF(PvTGas0 > p0) THEN
     PRINT *,' Saturated vapour pressure Pv = ',PvTGas0, &
          '  exceeds P0 = ',P0
     STOP
  ENDIF

  IF(PvTliq0 > P0) THEN
     PRINT *,' Saturated pressure at Tliq0 exceeds P0: ', &
          PvTliq0, P0
     STOP
  ENDIF
  !
  ! Volume fractions
  !
  DO j = 2, m
     IF(dzt(j) <= WaterLevel) THEN
        aa2(:,j) = Alpha_0
     ELSE
        aa2(:,j) = 1.D0
     ENDIF
  ENDDO



  aa2(:,1)=aa2(:,2); aa2(:,m9) = aa2(:,m)
  aa1 = 1.D0 - aa2
  al3 = 0.D0; al3dr = 0.D0; al3j = 0.D0; al3de = 0.D0
  !
  ! Densities and pressures: hydrostatic equilibrium
  !
  IF(ABS(Gravity)<1.D-8) THEN
     DO j = 1,m9
        IF(dzt(j) <= WaterLevel) THEN
           CALL Thermo_Props(TLiq0,TLiq0,P0,P0-PvTLiq0)
           ar1(:,j) = Ro_1;  ar2(:,j) = Ro_2; ara(:,j) = Ro_a
           at1(:,j) = TLiq0; at2(:,j) = TLiq0
           bp(:,j) = P0; bpa(:,j) = P0-PvTLiq0; bpv(:,j) = PvTLiq0
        ELSE
           TLiq = MIN(TGas0,sattmp(P0))
           CALL Thermo_Props(TLiq,TGas0,P0,P0-PvTGas0)
           ar1(:,j) = Ro_1;  ar2(:,j) = Ro_2; ara(:,j) = Ro_a
           at1(:,j) = TLiq; at2(:,j) = TGas0
           bp(:,j) = P0; bpa(:,j) = P0-PvTGas0; bpv(:,j) = PvTGas0
        ENDIF
     ENDDO
  ELSE ! Find self-consistent pressure and density fields
     ! Find pressure in the top cells
     IF(WaterLevel >= VesselHeight) THEN
        P_top = P0 - RoLiq0*Gravity*(dm(m9)-VesselHeight);
        CALL Thermo_Props(TLiq0,TLiq0,P_top,P_top-PvTLiq0)
        bp(:,m9) = P_top; bpa(:,m9) = P_top - PvTLiq0;
        bpv(:,m9) = PvTLiq0
        at1(:,m9) = TLiq0; at2(:,m9) = TLiq0
        ar1(:,m9) = Ro_1;  ar2(:,m9) = Ro_2; ara(:,m9) = Ro_a
     ELSE
        TLiq = MIN(TGas0,sattmp(P0))
        CALL Thermo_Props(TLiq,TGas0,P0,P0-PvTGas0)
        RoGas0 = Ro_2
        P_top = P0 - RoGas0*Gravity*(dm(m9)-VesselHeight)
        CALL Thermo_Props(TLiq,TGas0,P_top,P_top-PvTGas0)
        bp(:,m9) = P_top; bpa(:,m9) = P_top - PvTGas0;
        bpv(:,m9) = PvTGas0
        at1(:,m9) = TLiq; at2(:,m9) = TGas0
        ar1(:,m9) = Ro_1;  ar2(:,m9) = Ro_2; ara(:,m9) = Ro_a
     ENDIF
     Ro1 = ar1(1,m9); Ro2 = ar2(1,m9)
     DO j = m,1,-1
        Ro1_j = Ro_1; Ro2_j = Ro_2
        P_j0 = P_top + 0.5*Gravity*(Ro1_j*aa1(1,j)+Ro2_j*aa2(1,j)+ &
             Ro1*aa1(1,j+1)+Ro2*aa2(1,j+1))*(dm(j+1)-dm(j))

        DO ! Loop in which hydrostatic equilibrium is enforced
           IF(dzt(j) <= WaterLevel) THEN
              TLiq = TLiq0
              TGas = TLiq0
           ELSE
              TLiq = MIN(TGas0,sattmp(P_j0))
              TGas = TGas0
           ENDIF
           P_a = P_j0 - satprs(TGas)
           CALL Thermo_Props(TLiq,TGas0,P_j0,P_a)
           Ro1_j = Ro_1; Ro2_j = Ro_2
           P_j = P_top + 0.5*Gravity*(Ro1_j*aa1(1,j)+Ro2_j*aa2(1,j)+ &
                Ro1*aa1(1,j+1)+Ro2*aa2(1,j+1))*(dm(j+1)-dm(j))
           IF(ABS(P_j - P_j0)<1.D-8) EXIT
           P_j0 = P_j
        ENDDO
        ar1(:,j) = Ro_1;  ar2(:,j) = Ro_2; ara(:,j) = Ro_a
        at1(:,j) = TLiq;  at2(:,j) = TGas
        bp(:,j) = P_j; bpa(:,j) = P_a; bpv(:,j) = P_j - P_a
        P_top = P_j
        Ro2 = Ro_2; Ro1 = Ro_1
     ENDDO
  ENDIF

  check_hydro(1:m) = (bp(1,1:m)-bp(1,2:m9))/(dm(2:m9)-dm(1:m))-&
       0.5*Gravity*(aa2(1,1:m)*ar2(1,1:m)+aa2(1,2:m9)*ar2(1,2:m9)+&
       aa1(1,1:m)*ar1(1,1:m)+aa1(1,2:m9)*ar1(1,2:m9))

!!!DEBUG!!!
  !  TLiq0 = 350.D0
  !  TGas0 = 500.D0
  !  P0 = 2.d5
  !  Pa = 1.d5
  !  YH2 = 0.3D0
!!! Pure Liquid!
  !  TLiq0 = 350.D0
  !  TGas0 = 350.D0
  !  P0 = 2.D5
  !  Pa = P0-PvTLiq0
  !  YH2 = 0.D0

  ! Bubble !
  !  aa2(n9/2-1:n9/2+1,m/2-1:m/2+1) = 0.25
  !  at2(n9/2-1:n9/2+1,m/2-1:m/2+1) = 350
  !  aa1=1-aa2


  ! Heating by dispersed phase
!!$  r31c = 0;r32c = 0;T3=0
!!$  r31c(n9/2,m/2) = 5.D5
!!$  aT3(n9/2,m/2) = 3.D3



!!$
!!$  call Thermo_Props(TLiq0, TGas0, P0, Pa)
!!$  at1 = T1
!!$  at2 = T2
!!$  ar1 = Ro_1
!!$  ar2 = Ro_2
!!$  ae1 = E_1
!!$  ae2 = E_2
!!$  ara = Ro_a
!!$  aYa = Ro_a/Ro_2
!!$  bp = P0
!!$  bpa=Pa
!!$  bpv=bp-bpa
!!$
!!$  if(isHydrogenRelease /= 0) then
!!$     gamH2 = 1.D0
!!$  else
!!$     gamH2 = 0.D0
!!$  endif
!!$  aYH2 = YH2
!!$  YH2_Gamma_a = 0.2

  IF(isEvaporation /= 0) THEN
     DO j = 2, m
        DO i = 2, n
           CALL SetLocalVars(Alpha1 = aa1(i,j),&
                Alpha2 = aa2(i,j), &
                Alpha3 = al3(i,j), &
                TLiq = at1(i,j), &
                TGas = at2(i,j), &
                TDisp = at3(i,j),&
                Ptot = bp(i,j), &
                Pncon = bpa(i,j),&
                Y_a = aYa(i,j))

           YH2 = aYH2(i,j)

           VLiq(1) = 0.5*(au1(i,j)+au1(i-1,j))
           VLiq(2) = 0.5*(av1(i,j)+av1(i,j-1))
           VGas(1) = 0.5*(au2(i,j)+au2(i-1,j))
           VGas(2) = 0.5*(av2(i,j)+av2(i,j-1))

           CALL Thermo_Props
           CALL SetLocalProps
           CALL SetLiqGasParms

           CALL HeatExchangeGasLiq

           IF(isEvaporation /= 0) THEN
              CALL Evap_Rate
              gam(i,j) = Gamma_Evap
           ENDIF
        ENDDO
     ENDDO
  ELSE
     gam=0.D0
  ENDIF
END SUBROUTINE SetInitialValues


SUBROUTINE ReduceTimeStep
  USE PROBLEM_DATA
  USE SOLVER_2D
  USE INOUT_2D
  IMPLICIT NONE

  !  i_tsteps = 0


  tstep = tstep/tstep_factor

  IF(MsgLevel>=2) THEN
     WRITE(*,500), tstep
  ENDIF
500 FORMAT(24('>'),' Reducing time step: ',g9.2,t57,24('<'))


  IF(tstep < tstepmin) THEN
     PRINT *,'*** Minimum possible time step reached: ***'
     STOP
  ENDIF
END SUBROUTINE ReduceTimeStep

SUBROUTINE AdaptTimeStep
  USE PROBLEM_DATA
  USE SOLVER_2D
  USE VARIABLES_2D
  USE GRID_2D
  IMPLICIT NONE
  
  REAL(kind=8):: UMax, VMax

  i_tsteps = i_tsteps + 1
  IF(i_tsteps > n_tsteps .AND. tstep < tstepmax) THEN ! Try increasing time step
     tstep = MIN(tstepmax, tstep*tstep_factor)
     i_tsteps = 0
  ENDIF
  !
  ! Check the CFL condition
  !
  UMax = MAX(MAXVAL(ABS(au1)),MAXVAL(ABS(au2)))
  VMax = MAX(MAXVAL(ABS(av1)),MAXVAL(ABS(av2)))
  CFL = MAX(UMAX/hi,VMax/hj)*TStep
  if(CFL > CFLmax) then
    TStep = TStep*CFLMax/CFL
  endif
END SUBROUTINE AdaptTimeStep


SUBROUTINE WaterPack
END SUBROUTINE WaterPack


SUBROUTINE SolveForHydrogen
  USE HYDROGEN
  USE VARIABLES_2D
  USE GLOBALS_2D
  USE INOUT_2D
  USE GRID_2D
  USE GLOBALS
  USE SOLVER_2D
  IMPLICIT NONE

  LOGICAL:: okprint

  IF(isHydrogenRelease <= 0) RETURN

  okprint = (MsgLevel > 2)

  CALL HYDROGEN_2D(aYH2, &
       arH2, brH2, &
       au2, av2, ara,  aa2, ba2, &
       allim, &
       gamH2, &
       dxu, dxt, &
       tstep, tstep/hi, tstep/hj, &
       n, m, n9, m9, &
       niter_H2, okprint)

END SUBROUTINE SolveForHydrogen


SUBROUTINE SolveForVelocities
  USE VARIABLES_2D
  USE GLOBALS_2D
  USE INOUT_2D
  USE GRID_2D
  USE GLOBALS
  USE SOLVER_2D
  USE VELOCITY_2D
  USE PROBLEM_DATA

  IMPLICIT NONE  

  CALL velocity_coefs( &
       u1n, v1n, & ! u- and v-velocity of 1st phase on n-th time step
       u2n, v2n, & ! u- and v-velocity of 2nd phase on n-th time step
       au3, av3, & ! u- and v-velocity of 3rd phase
       bu1,bv1,  & ! u- and v-velocity of 1st phase on previous iteration 
       bu2,bv2,  & ! u- and v-velocity of 2nd phase on previous iteration
       ba1,ba2,  & ! phase volume fractions
       al3,dr3,  & ! 3rd phase volume fraction and diameter
       br1,br2,  & ! phase densities
       bp,       & ! pressure
       gam,      & ! phase change rate
       bmul,bmuv,& ! viscosities of liquid and gas phases
       au1, av1, & ! new (uncorrected) velocities on current iteration
       au2, av2, & !
       fu1l, fu1r, & ! u1(i,j) = au1(i,j) - (fu1l(i,j)*dp(i,j)+fu1r(i,j)*dp(i+1,j))
       fv1l, fv1r, & ! v1(i,j) = av1(i,j) - (fv1l(i,j)*dp(i,j)+fv1r(i,j)*dp(i,j+1))
       fu2l, fu2r, & ! u2(i,j) = au2(i,j) - (fu2l(i,j)*dp(i,j)+fu2r(i,j)*dp(i+1,j))
       fv2l, fv2r, & ! v2(i,j) = av2(i,j) - (fv2l(i,j)*dp(i,j)+fv2r(i,j)*dp(i,j+1))
       tstep, hi, hj, Gravity, &
       n,m,n9,m9)

  CALL SetVelocityBoundaryConditions

END SUBROUTINE SolveForVelocities


SUBROUTINE UpdatePhaseViscosities
  USE VARIABLES_2D
  USE GLOBALS_2D
  USE LOCAL_VARS
  USE SOLVER_2D
  USE COEFFICIENTS
  USE WATER_PROPS
  USE EVAPORATION
  USE CORRELATIONS
  USE HYDROGEN
  USE PROBLEM_DATA
  IMPLICIT NONE

  INTEGER:: i,j

  DO j = 2,m
     DO i = 2,n
        CALL SetLocalVars(Alpha1 = aa1(i,j),&
             Alpha2 = aa2(i,j), &
             Alpha3 = al3(i,j), &
             TLiq = at1(i,j), &
             TGas = at2(i,j), &
             TDisp = at3(i,j),&
             Ptot = bp(i,j), &
             Pncon = bpa(i,j),&
             Y_a = aYa(i,j))

        CALL Thermo_Props
        CALL SetLocalProps
        CALL SetLiqGasParms

        bmuv(i,j) = Visc_Gas
        bmul(i,j) = Visc_Liq
     ENDDO
  ENDDO

  bmul(1,:)=bmul(2,:); bmul(n9,:)=bmul(n,:)
  bmul(:,1)=bmul(:,2); bmul(:,m9)=bmul(:,m)
  bmuv(1,:)=bmuv(2,:); bmuv(n9,:)=bmuv(n,:)
  bmuv(:,1)=bmuv(:,2); bmuv(:,m9)=bmuv(:,m)

END SUBROUTINE UpdatePhaseViscosities


SUBROUTINE SetVelocityBoundaryConditions
  USE PROBLEM_DATA
  USE GLOBALS_2D
  USE VARIABLES_2D
  USE GRID_2D
  IMPLICIT NONE
  !
  ! Axis
  !
  au1(1,:) = 0.D0
  av1(1,:) = av1(2,:)
  au2(1,:) = 0.D0
  av2(1,:) = av2(2,:)
  !
  ! Bottom
  !
  av1(:,1) = 0.D0
  au1(:,1) = au1(:,2)
  av2(:,1) = 0.D0
  au2(:,1) = au2(:,2)
  !
  ! Side 
  !
  IF(isOpenSide /= 0) THEN
     au1(n,JOpening0:JOpening) = au1(n-1,JOpening0:JOpening)
     au2(n,JOpening0:JOpening) = au2(n-1,JOpening0:JOpening)
     au1(n,:JOpening0-1) = 0.D0
     au2(n,:JOpening0-1) = 0.D0
     au1(n,JOpening+1:m9) = 0.D0
     au2(n,JOpening+1:m9) = 0.D0
  ELSE
     au1(n,:) = 0.D0
     au2(n,:) = 0.D0
  ENDIF
  av1(n9,:) = av1(n,:)
  av2(n9,:) = av2(n,:)
  !
  ! Top
  !
  IF(isOpenTop /= 0) THEN
     !     av1(:,m) = av1(:,m-1)
     !     av2(:,m) = av2(:,m-1)
  ELSE
     av1(:,m) = 0.D0
     av2(:,m) = 0.D0
  ENDIF
  au1(:,m9) = au1(:,m)
  au2(:,m9) = au2(:,m)

END SUBROUTINE SetVelocityBoundaryConditions

!-------------------------------------- BiCGSTAB2
! Calculate y := A*x, where A is the iteration matrix
SUBROUTINE MATVEC_BCG2(nn,xx,yy)
  USE ELLEQ
  USE VARIABLES_2D
  USE GLOBALS_2D
  USE GRID_2D

  IMPLICIT NONE
  INTEGER:: nn
  REAL(kind=8):: xx(nn),yy(nn)

  CALL PressureCorrection(NX_ELL,NY_ELL,xx,yy,&
       r1x,q1x,rnx,qnx,r1y,q1y,rmy,qmy,QCoefMatrix,NVARS,NRHS)

END SUBROUTINE MATVEC_BCG2

!-------------------------------------- Problem-specific
SUBROUTINE PressureCorrection(n,m,xx,yy,&
     r1x,q1x,rnx,qnx,r1y,q1y,rmy,qmy,QCoefMatrix,NV,NR)
  USE GLOBALS
  IMPLICIT NONE
  INTEGER:: n,m,NV,NR
  REAL(kind=8):: xx(n,m),yy(n,m)
  REAL(kind=8),INTENT(in):: r1x(m),q1x(m),rnx(m),qnx(m)
  REAL(kind=8),INTENT(in):: r1y(n),q1y(n),rmy(n),qmy(n)
  REAL(kind=8),INTENT(in):: QCoefMatrix(NV,NR,n,m)
  INTEGER:: i, j

  ! Enforce boundary conditions
!!$  xx(1,:)=xx(2,:)*r1x(:) + q1x(:)
!!$  xx(n,:)=xx(n-1,:)*rnx(:) + qnx(:)
!!$  xx(:,1)=xx(:,2)*r1y(:) + q1y(:)
!!$  xx(:,m)=xx(:,m-1)*rmy(:) + qmy(:)

  yy = 0
  DO j = 3,m-2
     DO i = 3,n-2
        yy(i,j) = xx(i,j) - ( &
             QCoefMatrix(4,2,i,j)*xx(i-1,j) + &
             QCoefMatrix(4,3,i,j)*xx(i+1,j) + &
             QCoefMatrix(4,4,i,j)*xx(i,j-1) + &
             QCoefMatrix(4,5,i,j)*xx(i,j+1))
     ENDDO
  ENDDO
  ! Account for boundary conditions:
  DO i = 3,n-2
     j = 2
     yy(i,j) = xx(i,j) - ( &
          QCoefMatrix(4,2,i,j)*xx(i-1,j) + &
          QCoefMatrix(4,3,i,j)*xx(i+1,j) + &
          QCoefMatrix(4,4,i,j)*r1y(i)*xx(i,j) + &
          QCoefMatrix(4,5,i,j)*xx(i,j+1))
     j = m-1
     yy(i,j) = xx(i,j) - ( &
          QCoefMatrix(4,2,i,j)*xx(i-1,j) + &
          QCoefMatrix(4,3,i,j)*xx(i+1,j) + &
          QCoefMatrix(4,4,i,j)*xx(i,j-1) + &
          QCoefMatrix(4,5,i,j)*rmy(i)*xx(i,j))
  ENDDO
  DO j = 3, m-2
     i = 2
     yy(i,j) = xx(i,j) - ( &
          QCoefMatrix(4,2,i,j)*r1x(j)*xx(i,j)+&
          QCoefMatrix(4,3,i,j)*xx(i+1,j) + &
          QCoefMatrix(4,4,i,j)*xx(i,j-1) + &
          QCoefMatrix(4,5,i,j)*xx(i,j+1))
     i = n-1
     yy(i,j) = xx(i,j) - ( &
          QCoefMatrix(4,2,i,j)*xx(i-1,j) + &
          QCoefMatrix(4,3,i,j)*rnx(j)*xx(i,j) + &
          QCoefMatrix(4,4,i,j)*xx(i,j-1) + &
          QCoefMatrix(4,5,i,j)*xx(i,j+1))
  ENDDO
  !
  ! Corner points
  !
  i = 2; j = 2
  yy(i,j) = xx(i,j) - ( &
       QCoefMatrix(4,2,i,j)*r1x(j)*xx(i,j) + &
       QCoefMatrix(4,3,i,j)*xx(i+1,j) + &
       QCoefMatrix(4,4,i,j)*r1y(i)*xx(i,j) + &
       QCoefMatrix(4,5,i,j)*xx(i,j+1))
  j = m-1
  yy(i,j) = xx(i,j) - ( &
       QCoefMatrix(4,2,i,j)*r1x(j)*xx(i,j) + &
       QCoefMatrix(4,3,i,j)*xx(i+1,j) + &
       QCoefMatrix(4,4,i,j)*xx(i,j-1) + &
       QCoefMatrix(4,5,i,j)*rmy(i)*xx(i,j))

  i = n-1; j=2
  yy(i,j) = xx(i,j) - ( &
       QCoefMatrix(4,2,i,j)*xx(i-1,j) + &
       QCoefMatrix(4,3,i,j)*rnx(j)*xx(i,j) + &
       QCoefMatrix(4,4,i,j)*r1y(i)*xx(i,j) + &
       QCoefMatrix(4,5,i,j)*xx(i,j+1))

  j = m-1
  yy(i,j) = xx(i,j) - ( &
       QCoefMatrix(4,2,i,j)*xx(i-1,j) + &
       QCoefMatrix(4,3,i,j)*rnx(j)*xx(i,j) + &
       QCoefMatrix(4,4,i,j)*xx(i,j-1) + &
       QCoefMatrix(4,5,i,j)*rmy(i)*xx(i,j))

END SUBROUTINE PressureCorrection

SUBROUTINE SetPressureRHS(&
     QCoefMatrix,RHS,&
     r1x,q1x,rnx,qnx,&
     r1y,q1y,rmy,qmy,&
     NVARS,NRHS,n,m,n9,m9)

  IMPLICIT NONE
  INTEGER, INTENT(in):: NVARS,NRHS,n,m,n9,m9
  REAL(kind=8), INTENT(inout):: QCoefMatrix(NVARS,NRHS,n9,m9)
  REAL(kind=8), INTENT(inout):: RHS(n9,m9)
  REAL(kind=8), INTENT(in):: r1x(m9),q1x(m9),rnx(m9),qnx(m9)
  REAL(kind=8), INTENT(in):: r1y(n9),q1y(n9),rmy(n9),qmy(n9)

  INTEGER:: n1, m1

  ! Initialize the RHS
  RHS = 0.D0
  RHS(2:n,2:m) = QCoefMAtrix(4,1,2:n,2:m)

  n1 = n-1; m1 = m-1
  ! Correct RHS array in near-boundary points
  RHS(2,2:m1)  = RHS(2,2:m1)  + QCoefMatrix(4,2,2,2:m1)*q1x(2:m1)
  RHS(n1,2:m1) = RHS(n1,2:m1) + QCoefMatrix(4,3,n1,2:m1)*qnx(2:m1)
  RHS(2:n1,2)  = RHS(2:n1,2)  + QCoefMatrix(4,4,2:n1,2)*q1y(2:n1)
  RHS(2:n1,2)  = RHS(2:n1,2)  + QCoefMatrix(4,5,2:n1,2)*qmy(2:n1)

END SUBROUTINE SetPressureRHS

SUBROUTINE SetPressureBoundaryConditions(dp,&
     r1x,q1x,rnx,qnx,&
     r1y,q1y,rmy,qmy,&
     n,m,n9,m9)

  IMPLICIT NONE
  INTEGER, INTENT(in):: n,m,n9,m9

  REAL(kind=8), INTENT(inout):: dp(n9,m9)
  REAL(kind=8), INTENT(in):: r1x(m9),q1x(m9),rnx(m9),qnx(m9)
  REAL(kind=8), INTENT(in):: r1y(n9),q1y(n9),rmy(n9),qmy(n9)

  dp(1,:)=dp(2,:)*r1x(:) + q1x(:)
  dp(n9,:)=dp(n9-1,:)*rnx(:) + qnx(:)
  dp(:,1)=dp(:,2)*r1y(:) + q1y(:)
  dp(:,m9)=dp(:,m9-1)*rmy(:) + qmy(:)

  dp(1,1)=0.5*(dp(1,2)+dp(2,1))
  dp(n9,1) = 0.5*(dp(n,1)+dp(n9,2))
  dp(1,m9)=0.5*(dp(1,m)+dp(2,m9))
  dp(n9,m9) = 0.5*(dp(n9,m)+dp(n,m9))

END SUBROUTINE SetPressureBoundaryConditions

SUBROUTINE SaveData
  USE INOUT_2D
  IMPLICIT NONE

  IF(NWRITE /= 0) CALL WriteRestart
END SUBROUTINE SaveData

