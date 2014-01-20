program water_test
  USE WATER_PROPS
  USE DISP_PROPS
  USE LOCAL_VARS
  USE COEFFICIENTS
  USE CORRELATIONS
  USE EVAPORATION
  USE HYDROGEN

  implicit none
  double precision:: QVar(5), QVarn(5), DQVar(5)
  double precision:: dR31dA, dR32dA

  double precision:: TMatr(5,5),AMatr(5)
  integer:: ipiv(5)
  double precision:: tstep, tstepmax, tstepmin, ts_factor
  integer:: nsteps_ts, isteps_ts

  double precision:: Alpha_30

  integer:: nsteps, istep
  integer:: NIter, iter
  double precision:: Tol_NL
  integer:: info
  integer,parameter:: NRHS_WT=1

  logical BadIter
  double precision:: dP, Pnew, Pn, Pan

  ! Zeroth iteration values
  double precision:: Alpha_k0, T1_k0, T2_k0, P_k0, Pa_k0

  
  double precision:: Omega
  integer:: i,j

  integer:: igastype,isMixture,isEvaporation,isHydrogenRelease
  integer:: isDispersedPhase

  ! Setup model parameters
  !  DBubbleMin = 1.D-4;DBubbleMax=1.D-2
  !  DDropletMin=1.D-4;DDropletMax=1.D-2
  VGas = 1.5D-0
  VLiq = 0.D0

  SurfaceTension = 22.5d-03
  Lambda_Liq = 64.4d-2
  Visc_Gas = 2.27d-05
  Cp_Gas = 0.521d4  
  Lambda_Gas = 1.77d-01         ! Argon T = 300 K p=1   atm
  !
  open(90,file='water.inp')
  read(90,*) Alpha
  read(90,*) T1
  read(90,*) T2
  read(90,*) P
  read(90,*) Pa
  read(90,*) igastype
  read(90,*) isMixture
  read(90,*) YH2
  read(90,*) tstepmax
  read(90,*) tstepmin
  read(90,*) ts_factor
  read(90,*) nsteps_ts
  read(90,*) nsteps
  read(90,*) NIter
  read(90,*) Tol_NL
  read(90,*) isEvaporation
  read(90,*) isHydrogenRelease 
  read(90,*) Gamma_a
  read(90,*) YH2_Gamma_a
  read(90,*) TH2_Gamma_a
  read(90,*) isDispersedPhase
  read(90,*) Alpha_3
  read(90,*) T3
  read(90,*) dR31dA
  read(90,*) dR32dA
  close(90)


  call SetPhaseChangeFlag(isEvaporation)
  call SetGassAdditionFlag(isHydrogenRelease)  
  if(isHydrogenRelease == 0) Gamma_a = 0.D0

  if(isDispersedPhase == 0) Alpha_3 = 0
  Alpha_30 = Alpha_3

  call SetNonCondensGasType(igastype,isMixture)
  call InitWaterProps
  call InitDispProps


  CALL Thermo_Props()
  Pv = P-Pa
  E1 = E_1; E2 = E_2
  Ro1 = Ro_1; Ro2 = Ro_2; Roa = Ro_a
  Ya = Roa/Ro2
  if(isDispersedPhase /= 0) then
     Alpha_l = 1.D0-Alpha-Alpha_3
  else
     Alpha_l = 1.D0-Alpha
  endif
  RoH2 = Roa*YH2

  CALL HeatExchangeGasLiq
  CALL Evap_Rate

  QVar(1) = Alpha; QVar(2) = T1; QVar(3) = T2; QVar(4) = P; QVar(5) = Pa
  DQVar = 0.D0

  time = 0.D0
  tstep = tstepmax
  isteps_ts = 0

  !-------------------------------------- Output files ------------------------
  open(55,file="TestRes.dat")
  write(55,'(10(G14.8E3,1x))') time,QVar(1:5),tstep

  open(60,file="Temp.dat")
  write(60,'(10(G14.8,1x))') time,QVar(2:3),T_sat,Tsv

  open(65,file="Press.dat")
  write(65,'(10(G14.8,1x))') time,P,Pa,Pv,Pv/satprs(T1)*100,Pv/satprs(T2)*100

  open(70,file="Mass.dat")
  write(70,'(10(G14.8,1x))') time,Ro1*Alpha_l,Ro2*Alpha,Ro_v*Alpha,Roa*Alpha, &
       Ro1*Alpha_l+Ro2*Alpha,Ro1*Alpha_l+Ro_v*Alpha,Ya

  open(75,file="MTransfer.dat")
  write(75,'(10(G14.8,1x))') time,Gamma_Evap,Gamma_a,R1s,R2s

  open(80,file="Heat.dat")
  write(80,'(10(G14.8,1x))') time,Ro1*Alpha_l*E1,Ro2*Alpha*E2,&
       Ro1*Alpha_l*E1+Ro2*Alpha*E2

  open(85,file="Alpha.dat")
  write(85,'(10(G14.8,1x))') time,Alpha_l,Alpha,Alpha_3,&
       & Alpha_l+Alpha+Alpha_3

  open(90,file="Hydrogen.dat")
  write(90,'(10(G14.8,1x))') time,YH2,XH2,W_a,RoH2*Alpha,RoH2


  !----------------------------------------------------------------------------
  do istep = 1, nsteps

     QVarn = QVar
     Alphan = QVarn(1)

     if(isDispersedPhase /= 0) then
        Alpha_ln = 1.D0-Alphan-Alpha_3
     else
        Alpha_ln = 1.D0-Alphan
     endif

     E1n = E1; E2n = E2; Ro1n = Ro1; Ro2n = Ro2
     T1n = T1; T2n = T2
     Yan = Ya
     Pn = P; Pan = Pa
     RoH2n = RoH2

     if(isDispersedPhase /= 0) then
        ! Changes to Alpha_3 and T3 should be defined here
!        Alpha_3 = Alpha_30*cos(10*time)**2
     endif

3000 continue

!!! Trying to come closer to equilibrium state at the 1st iteration:
     T2 = T1; Pa = P - satprs(T2)
     QVar(3) = T2; QVar(5) = Pa
     CALL Thermo_Props
     Pv = P - Pa
     E1 = E_1; E2 = E_2
     Ro1 = Ro_1; Ro2 = Ro_2; Roa = Ro_a
     Ya = Y_a


     NewtonIter:   do iter = 1, NIter

        ! Fill Jacobian of time derivatives and source terms
        if(isHydrogenRelease /= 0) then
           if(Alpha < allim .and. Gamma_a > 0) YH2 = 1
        endif
        CALL Thermo_Props
        CALL HeatExchangeGasLiq
        CALL Evap_Rate
        CALL EQUATION_COEFS(tstep)

        Tmatr = TimeJac - SourceJac
        Amatr = RHS

        !if(istep==1 .and. iter==1) then
!!!DEBUG!!!

        open(81,file='TimeJac.dat')
        write(81,*) time,tstep
        do j = 5,1,-1
           write(81,*) (TimeJac(i,j),i=1,5)
        enddo
        close(81)
        open(81,file='SourceJac.dat')
        do j = 5,1,-1
           write(81,*) (SourceJac(i,j),i=1,5)
        enddo
        close(81)
        open(81,file='RHS.dat')
        do j = 5,1,-1
           write(81,*) RHS(j)
        enddo
        close(81)
        open(81,file='Vars.dat')
        write(81,*) Alpha_l
        write(81,*) Alpha
        write(81,*) T1
        write(81,*) T2
        write(81,*) P
        write(81,*) Pa
        write(81,*) Gamma_Evap
        close(81)

        !endif
        call DGESV(NVARS,NRHS_WT,Tmatr,NVARS,ipiv,Amatr, &
             NVARS,info)
        if(info /= 0) then
           print *,'Error in DGESV: info = ',info
           stop
        endif

        DQVar = Amatr

        ! Check the resulting pressure and temperature fields
        BadIter = .false.
        dP = DQVar(4)
        Pnew = P + dP
        if(Pnew >= EOSPMax() .OR. Pnew <= EOSPMin()) BadIter = .true.
        if(T1+DQVar(2) >= EOSTLiqMax() .OR. T1+DQVar(2) <= EOSTLiqMin()) BadIter = .true.
        if(T2+DQVar(3) >= EOSTVapMax() .OR. T2+DQVar(3) <= EOSTVapMin()) BadIter = .true.
!        if(dP > 0.D0 .and. dP > 0.2*P) BadIter = .true.
!        if(dP < 0.D0 .and. abs(dP) > 0.1*P) BadIter = .true.

        if(BadIter) then ! Skip assignment, repeat current time step with smaller tstep
           Alpha = Alphan; Alpha_l = Alpha_ln;
           T1 = T1n; T2 = T2n; P = Pn; Pa = Pan;
           Ya = Yan; RoH2 = RoH2n
           tstep = tstep/ts_factor
           print *,'****************************> tstep = ',tstep

           if(tstep < tstepmin) then
              print *,'*** Minimum possible time step reached: ***'
              print *,'DQVar = ',dQVar(1:5)
              print *,'Pnew = ',Pnew
              stop
           endif
           goto 3000
        endif

        ! Update variables
        Omega = min(0.1*P/max(abs(dP),1.D-2),1.D0)
!        Omega=1

        QVar = QVar + DQVar*Omega
        Alpha = QVar(1); T1 = QVar(2) ; T2 = QVar(3);  P = QVar(4); Pa = QVar(5)

        Alpha = min(max(Alpha,0.D0),1.D0)

        if(isDispersedPhase /= 0) then
           Alpha_l = 1.D0-Alpha-Alpha_3
        else
           Alpha_l = 1.D0-Alpha
        endif
        Pa = min(max(Pa,0.D0),P)

        CALL HydrogenRelease0D(tstep)

        if(isHydrogenRelease /= 0) then
           if(Alpha < allim .and. Gamma_a > 0) YH2 = 1
        endif
        CALL Thermo_Props
        Pv = P - Pa
        E1 = E_1; E2 = E_2
        Ro1 = Ro_1; Ro2 = Ro_2; Roa = Ro_a
        Ya = Y_a
        print *,'=========>',maxval(abs(DQVAr))

        if(maxval(abs(DQVar)) < Tol_NL) exit
     enddo NewtonIter

     time = time + tstep
     write(55,'(10(G14.8,1x))') time,QVar(1:5),tstep
     write(60,'(10(G14.8,1x))') time,QVar(2:3),T_sat,Tsv
     write(65,'(10(G14.8,1x))') time,P,Pa,Pv,Pv/satprs(T1)*100,Pv/satprs(T2)*100
     write(70,'(10(G14.8,1x))') time,Ro1*Alpha_l,Ro2*Alpha,Ro_v*Alpha,Roa*Alpha, &
          Ro1*Alpha_l+Ro2*Alpha,Ro1*Alpha_l+Ro_v*Alpha,Ya
     write(75,'(10(G14.8,1x))') time,Gamma_Evap,Gamma_a,R1s,R2s
     write(80,'(10(G14.8,1x))') time,Ro1*Alpha_l*E1,Ro2*Alpha*E2,&
          Ro1*Alpha_l*E1+Ro2*Alpha*E2
     write(85,'(10(G14.8,1x))') time,Alpha_l,Alpha,Alpha_3,&
          & Alpha_l+Alpha+Alpha_3
     write(90,'(10(G14.8,1x))') time,YH2,XH2,W_a,RoH2*Alpha,RoH2

     write(*,'(10(G14.8,1x))') time,QVar(1:5),Ro1*Alpha_l,Ro2*Alpha,Ro_v*Alpha,Ya

     isteps_ts = isteps_ts + 1
     if(isteps_ts > nsteps_ts .and. tstep < tstepmax) then ! Try increasing time step
        tstep = min(tstepmax, tstep*ts_factor)
        isteps_ts = 0
     endif
  enddo
  close(55)
  close(60)
  close(65)
  close(70)
  close(75)
  close(80)
  close(85)
  close(90)
end program water_test
