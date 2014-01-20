!------------------------------------------------------------------!
!  Lagrangian model for dispersed phase                            !
!------------------------------------------------------------------!
!  $Id: MOD_Disp2D.f90 7 2014-01-14 10:07:24Z Sergey $
!------------------------------------------------------------------!
MODULE DISPERSED_PHASE_2D
  USE DISPERSED_PHASE
CONTAINS
  ! Main routine advancing dispersed particles
  SUBROUTINE SolveForDispersedPhase
    USE GLOBALS_2D
    USE GLOBALS
    USE VARIABLES_2D
    USE GRID_2D
    USE CORRELATIONS
    USE LOCAL_VARS
    USE PROBLEM_DATA
    IMPLICIT NONE
    !
    ! Local variables
    !
    REAL(KIND=8):: delz
    INTEGER:: i_new, n_new, i, j, k, id
    INTEGER:: iu, ju, iv, jv, ite, jte
    REAL(KIND=8):: xk,zk ! Particle coordinates
    REAL(KIND=8):: su1,su2,su3,su4
    REAL(KIND=8):: sv1,sv2,sv3,sv4
    REAL(KIND=8):: st1,st2,st3,st4
    REAL(KIND=8):: T1_D, T2_D, R1_D, R2_D, A1_D, A2_D, MUL_D, MUV_D
    REAL(KIND=8):: Vmod1, Vmod2, Denom, VMod13, VMod23
    REAL(KIND=8):: dD_k, dM_k, D_new, qmold, frac, zlead

    REAL(KIND=8):: surfwater,volwater,volvapor,voljet
    REAL(KIND=8):: qmwater,qmjeth
    REAL(KIND=8):: cdh,fi,RoRel,xt3old

    INTEGER:: iLocFr, jLocFr, iLocFr1,iLocFr2
    LOGICAL:: UnderWater
    REAL rand_xy(2) ! Used for randomly placing debris particle 
    REAL rand_num

    if(isThirdPhase == 0) return

    ! Indices for fragmentation rate assignment to 2D array
    iLocFr = 1
    DO i = 1,n-1
       IF(dxu(i+1) > 0.5*Diam .AND. dxu(i) >= 0.5*Diam) THEN
          iLocFr = i
          EXIT
       ENDIF
    ENDDO
    iLocFr1 = 1
    DO i = 1,n-1
       IF(dxu(i+1) > Diam .AND. dxu(i) >= Diam) THEN
          iLocFr1 = i
          EXIT
       ENDIF
    ENDDO
    iLocFr2 = 1
    DO i = 1,n-1
       IF(dxu(i+1) > 2*Diam .AND. dxu(i) >= 2*Diam) THEN
          iLocFr2 = i
          EXIT
       ENDIF
    ENDDO

    iLocFr  = MAX(2,MIN(iLocFr,n))
    iLocFr1 = MAX(2,MIN(iLocFr1,n))
    iLocFr2 = MAX(2,MIN(iLocFr2,n))

    IF(time_disp > time) THEN
       CALL DispersedToGrid
       RETURN
    ENDIF

    ! Integration loop for dispersed phase
1000 CONTINUE
    IF(time_disp > time) RETURN



    !----------------- RELEASE OF JET PARTICLES ----------------------
    IF (qmpart < qmpart0) THEN ! Release jet particles
       IF(ll == 0) THEN
          !-------Diameter of the jet particle
          xd3(1)=Djet
          !-------Radial location of the particle - 0.D0
          x(1)=0.D0
          !------Initial height of jet particles (release orifice length included)
          z(1)=hrelease
          xu3(1)=0.D0
          !-------Initial velocity of particle
          xv3(1)=-1.d-3
          !-------Temperature of the particle
          xt3(1)=Tin
          !--------Mass of the particle
          qmk(1)=qmpart_in
          !--------Type of the particle 1 - jet particle
          ind(1)=1
          !--------Mass of the melt released
          qmpart=qmpart+qmk(1)
          !------ll - total number of the particles
          ll=1
          !--------Index of leading edge particle
          lid=1
          !--------Index of the last particle which left the nozzle
          llj=1
       ELSE
          !------distance between initial height and jet particle
          !------which was the last to leave the nozzle
          delz = hrelease - z(llj)
          n_new = INT(delz/Z_JetPart)
          DO i_new = 1,n_new
             !-------If distance exceed Z_JetPart = qmpart_in/(rmelt*Spour)
             !-------next jet particle released from the nozzle
             ll=ll+1 ! count new particle
             xd3(ll)=Djet
             x(ll)=0.D0
             z(ll)=z(llj)+Z_JetPart
             qmk(ll)=qmpart_in
             xu3(ll)=0.D0
             !-------Velocity of the next particle equals to velocity
             !-------of the leading edge particle
             xv3(ll)= xv3(lid)
             xt3(ll)=Tin
             !-----Jet particle index
             ind(ll)=1
             !-----Index of the "last released" particle was changed
             llj=ll
             !------Total mass of the melt released
             qmpart=qmpart+qmk(ll)
          ENDDO
       ENDIF
    ENDIF

    !-----ll3 -auxiliary index - at this point it equals
    !-----to the number of jet particles
    ll3=ll

    !------ Mass, surface and volume of jet below water level
    qmwater= SUM(qmk,mask=(ind==1 .AND. z<WaterLevel))
    surfwater = pi*Z_JetPart*SUM(xd3,mask=(ind==1 .AND. z<WaterLevel))
    volwater=qmwater/rmelt
    !------ Total mass and volume of jet
    qmjeth = SUM(qmk,mask=ind==1)
    voljet=qmjeth/rmelt
    !------ Volume of the jet located in vapor
    volvapor=voljet-volwater
    !-------------------------- FRAGMENTATION ------------------------------
    ncount_frag=ncount_frag+1 !--------We count the number of time steps
    ! Fragmentation occurs after nskip_frag time steps
    FRAGMENT: IF(ncount_frag >  nskip_frag .AND. lid > 0) THEN 

       ncount_frag=1 ! Reset counter
       FRate = 0.D0      ! Reset total fragmentation rate
       FRateLocal = 0.D0 ! Reset fragmentation rate array
       !--------Leading edge is under water

       IF(z(lid) < WaterLevel) THEN
          FRAGMENTATION_LOOP: DO k=1,ll
             !------condition for fragmentation
             IF(ind(k) == 1 .AND.  We(k) > 1.D0) THEN

                ! Find out if the particle is under water

                UnderWater = z(k)<WaterLevel

                IF(.NOT.UnderWater) CYCLE

                ! (1/S)*(dM/dt)=v*rmelt/(2*(L/D)) ==> dDj/dt = tg(A)*v where tg(A)=1/(L/D)

                ! Reduce the diameter and mass of k-th jet macroparticle
                IF(isSaito) fragL_D = Sai(lid)

                dD_k = (tstep_disp*nskip_frag)*ABS(xv3(k))/fragL_D
                D_new = xd3(k) - dD_k

                IF(D_new <= 0.D0) THEN ! Total fragmentation, eliminate particle
                   dM_k = qmk(k)
                   xd3(k) = 0.D0; qmk(k) = 0.D0; ind(k) = 4
                ELSE
                   xd3(k) = D_new
                   qmold = qmk(k); qmk(k) = 0.25*pi*D_new**2*rmelt*Z_JetPart
                   dM_k = qmold-qmk(k)
                ENDIF

                FRate = FRate + dM_k/(tstep_disp*nskip_frag)

                ! Add the mass remaining from previous fragmentation
                dM_k = dM_k + dMfrag(k); dMfrag(k) = 0.D0

                IF(dM_k < DropM) THEN ! Insufficient droplet mass, postpone 
                   dMfrag(k) = dM_k
                ELSE
                   new_drops: DO id = 1,nfrag ! Distribute the mass 
                      ll3 = ll3 + 1
                      IF(ll3 > nump) THEN
                         PRINT *,'ERROR IN DISP: INSUFFICIENT 1D ARRAY SIZE ',nump
                         STOP
                      ENDIF
                      !------ll3 - the last index of jet particles
                      qmk(ll3)=dM_k/nfrag
                      xd3(ll3)=Ddrop
                      !-------initial radial location
                      x(ll3)=Diam/2.D0
                      !------intial height of the particles
                      z(ll3)= z(k)-0.5*Z_JetPart + Z_JetPart*float(id)/float(nfrag+1)
                      !-------droplet index
                      ind(ll3)=2
                      !-------intial droplet velocity and temperature
                      xu3(ll3)=ABS(xv3(k))*dropConeAngle
                      xv3(ll3)=xv3(k)
                      xt3(ll3)=xt3(k)

                      jLocFr=INT( (z(ll3)+0.5D0*hj)*hj_inv ) + 1

                      DO iLocFr=iLocFr1,iLocFr2
                         IF(z(ll3) >= dzt(jLocFr) .AND.&
                              z(ll3) < dzt(jLocFr+1)) THEN
                            frac = (dzt(jLocFr+1) - z(ll3))/(dzt(jLocFr+1)-dzt(jLocFr))

                            FRateLocal(iLocFr,jLocFr) = FRateLocal(iLocFr,jLocFr) +&
                                 frac*dM_k/(nfrag*tstep_disp*nskip_frag)*volc(iLocFr)/&
                                 (iLocFr2-iLocFr1+1)
                            FRateLocal(iLocFr,jLocFr+1) = FRateLocal(iLocFr,jLocFr+1) +&
                                 (1.D0-frac)*dM_k/(nfrag*tstep_disp*nskip_frag)*&
                                 volc(iLocFr)/(iLocFr2-iLocFr1+1)
                         ELSEIF(z(ll3) < dzt(jLocFr) .AND.&
                              z(ll3) >= dzt(jLocFr+1)) THEN
                            frac = (z(ll3) - dzt(jLocFr-1))/(dzt(jLocFr)-dzt(jLocFr-1))

                            FRateLocal(iLocFr,jLocFr) = FRateLocal(iLocFr,jLocFr) +&
                                 frac*dM_k/(nfrag*tstep_disp*nskip_frag)*volc(iLocFr)/&
                                 (iLocFr2-iLocFr1+1)
                            FRateLocal(iLocFr,jLocFr-1) = FRateLocal(iLocFr,jLocFr-1) +&
                                 (1.D0-frac)*dM_k/(nfrag*tstep_disp*nskip_frag)*&
                                 volc(iLocFr)/(iLocFr2-iLocFr1+1)
                         ELSE
                            PRINT *,'Bad index'
                            STOP
                         ENDIF
                         IF(frac < 0.D0) THEN
                            PRINT *,'Bad fraction'
                            STOP
                         ENDIF

                      ENDDO

                   ENDDO new_drops
                ENDIF
             ENDIF
          ENDDO FRAGMENTATION_LOOP
       ENDIF
       !------total number of the particles
       ll=ll3
       ! Update the index of the leading edge particle
       lid = -1 ! Initialize 
       zlead = 1.d10 ! Just make sure that search loop is properly started
       DO k=1,ll
          IF(ind(k)==1 .AND. z(k)<zlead) THEN
             zlead = z(k); lid = k
          ENDIF
       ENDDO
    ENDIF FRAGMENT
    !----------------------- DROPLET MOVEMENT ------------------------------
    !
    ! Perform one time step for each particle
    !
    LOOP_ALL_PARTICLES: DO k=1,ll
       IF(ind(k) == 4) CYCLE 
       zk=z(k); IF(zk <= 0.D0) CYCLE
       xk=x(k)

       iu=INT(xk*hi_inv) + 1
       ju=INT( (zk+0.5D0*hj)*hj_inv ) + 1

       su1=(dxu(iu+1)-xk)*(dzu(ju+1)-zk)
       su2=(xk-dxu(iu))*(dzu(ju+1)-zk)
       su3=(dxu(iu+1)-xk)*(zk-dzu(ju))
       su4=(xk-dxu(iu))*(zk-dzu(ju))

       VLiq(1)=(bu1(iu,ju)*su1 + bu1(iu+1,ju)*su2 + bu1(iu,ju+1)*su3 &
            + bu1(iu+1,ju+1)*su4)*s_inv

       VGas(1)=(bu2(iu,ju)*su1 + bu2(iu+1,ju)*su2 + bu2(iu,ju+1)*su3 &
            + bu2(iu+1,ju+1)*su4)*s_inv

       iv=INT( (xk+0.5D0*hi)*hi_inv ) + 1
       jv=INT(zk*hj_inv) + 1

       sv1=(dxv(iv+1)-xk)*(dzv(jv+1)-zk)
       sv2=(xk-dxv(iv))*(dzv(jv+1)-zk)
       sv3=(dxv(iv+1)-xk)*(zk-dzv(jv))
       sv4=(xk-dxv(iv))*(zk-dzv(jv))

       VLiq(2)=(bv1(iv,jv)*sv1 + bv1(iv+1,jv)*sv2 + bv1(iv,jv+1)*sv3 &
            + bv1(iv+1,jv+1)*sv4)*s_inv

       VGas(2)=(bv2(iv,jv)*sv1 + bv2(iv+1,jv)*sv2 + bv2(iv,jv+1)*sv3 &
            + bv2(iv+1,jv+1)*sv4)*s_inv


       ite=iv
       jte=ju

       st1=(dxt(ite+1)-xk)*(dzt(jte+1)-zk)
       st2=(xk-dxt(ite))*(dzt(jte+1)-zk)
       st3=(dxt(ite+1)-xk)*(zk-dzt(jte))
       st4=(xk-dxt(ite))*(zk-dzt(jte))

       T1_D=(at1(ite,jte)*st1 + at1(ite+1,jte)*st2 + at1(ite,jte+1)*st3 &
            + at1(ite+1,jte+1)*st4)*s_inv

       T2_D=(at2(ite,jte)*st1 + at2(ite+1,jte)*st2 + at2(ite,jte+1)*st3 &
            + at2(ite+1,jte+1)*st4)*s_inv

       R1_D=(ar1(ite,jte)*st1 + ar1(ite+1,jte)*st2 + ar1(ite,jte+1)*st3 &
            + ar1(ite+1,jte+1)*st4)*s_inv

       R2_D=(ar2(ite,jte)*st1 + ar2(ite+1,jte)*st2 + ar2(ite,jte+1)*st3 &
            + ar2(ite+1,jte+1)*st4)*s_inv

       A1_D=(aa1(ite,jte)*st1 + aa1(ite+1,jte)*st2 + aa1(ite,jte+1)*st3 &
            + aa1(ite+1,jte+1)*st4)*s_inv

       A2_D=(aa2(ite,jte)*st1 + aa2(ite+1,jte)*st2 + aa2(ite,jte+1)*st3 &
            + aa2(ite+1,jte+1)*st4)*s_inv

       MUL_D=(bmul(ite,jte)*st1 + bmul(ite+1,jte)*st2 + bmul(ite,jte+1)*st3 &
            + bmul(ite+1,jte+1)*st4)*s_inv

       MUV_D=(bmuv(ite,jte)*st1 + bmuv(ite+1,jte)*st2 + bmuv(ite,jte+1)*st3 &
            + bmuv(ite+1,jte+1)*st4)*s_inv

       VDisp(1) = xu3(k); VDisp(2) = xv3(k)



!!!DEBUG!!!
!!$       VLiq(1) = 0.D0
!!$       VLiq(2) = -0.1D0
!!$       VGas(1) = 0D0
!!$       VGas(2) = -0.2D0
!!$       
!!$       A1_d = 0.2D0
!!$       A2_d = 0.8D0
!!$       T1_d = 299.5D0
!!$       T2_d = 295D0
!!$       R1_d = 980D0
!!$       R2_d = 1.2D0
!!$       MUV_D = 2.D-5
!!$       MUL_D = 5.D-4
!!$


       fi=a2_d/(a1_d+a2_d)
       RoRel=((1-fi)*r1_d + fi*r2_d)/rmelt

       Vmod13=SQRT(DOT_PRODUCT(VLiq-VDisp,VLiq-VDisp))
       Vmod23=SQRT(DOT_PRODUCT(VGas-VDisp,VGas-VDisp))

       ! Set local parameters
       Visc_Liq = MUL_D
       Visc_Gas = MUV_D
       Ro1 = R1_D
       Ro2 = R2_D
       Alpha = A2_D
       Alpha_l = A1_D
       DDispPhase = xd3(k)
       T1 = T1_d
       T2 = T2_d
       T3 = xt3(k)
       CALL DragDispParticle            ! Calculate C13, C23
       CALL HeatExchangeDispParticle    ! Calculate 

       Vmod1=SQRT(DOT_PRODUCT(VLiq-VDisp,VLiq-VDisp))*c13*tstep_disp/rmelt
       Vmod2=SQRT(DOT_PRODUCT(VGas-VDisp,VGas-VDisp))*c23*tstep_disp/rmelt
       Denom=1.D0/(1.D0+Vmod1+Vmod2)
       xt3old=xt3(k)
       !
       !------Criteria used for fragmentation 
       !
       IF(z(k) > WaterLevel) THEN
          We(k) = rmelt*(Vmod23**2)*xd3(k)/sigme
          Sai(k) = 2.1D0*SQRT(1.D0/RoRel)*SQRT((Vmod23**2)/(Gravity*xd3(k)))
          Epf(k)= (SQRT(3.D0)/2.D0)*(1.D0 + R2_D/rmelt)*SQRT(rmelt/R2_D)
       ELSE
          We(k)= rmelt*(Vmod13**2)*xd3(k)/sigme
          Sai(k) = 2.1D0*SQRT(1.D0/RoRel)*SQRT((Vmod13**2)/(Gravity*xd3(k)))
          Epf(k)= (SQRT(3.D0)/2.D0)*(1.D0 + R1_D/rmelt)*SQRT(rmelt/R1_D)
       ENDIF

       IF(k == lid) THEN
          !------qpl=1.D0*(hr+hc)*prc6/dm *area
          !----- prc6=6.D0/(rmelt*cmelt) - qpl = ((hr+hc)/(rmelt*cmelt))*6/dm
          !------6/dm - surface concentration
          !------ correction  factor area=alvap/(alliq+alvap)
          !------correction for CYLINDER surface
          !------concentration  - for cylinder is 4/d
          !------for ball - 6/d
          qpg=4.D0*qpg/6.D0
          qpl=4.D0*qpl/6.D0
          ! Update particle velocity
          xu3(k)=(xu3(k) + VLiq(1)*Vmod1 + VGas(1)*Vmod2)*Denom

          cdh=0.06D0
          Denom=1.D0/(1.D0+tstep_disp*cdh*R1_D*surfwater*Vmod13/(rmelt*voljet))
          XV3(K)=(XV3(K) &
               -tstep_disp*Gravity*volwater/voljet &
               -tstep_disp*Gravity*(1.D0-RoRel)*volvapor/voljet+ &
               tstep_disp*cdh*R1_D*surfwater*vmod13*VLiq(2)/(voljet*rmelt))*Denom

          IF(isUseVin /= 0) THEN
             xu3(k) = 0.D0
             xv3(k) = MAX(xv3(k),-Vin)
          ENDIF

          xt3(K)=(xt3(K) + T1_D*tstep_disp*qpl + T2_D*tstep_disp*qpg) / &
               (1.D0+tstep_disp*qpl+tstep_disp*qpg)
       ENDIF

       !-----defining of the characteristics of the jet
       !-----particles
       IF(k > lid .AND. ind(k) == 1) THEN
          !------correction for CYLINDER
          qpg=4.D0*qpg/6.D0
          qpl=4.D0*qpl/6.D0
          xu3(k)=xu3(lid)
          xv3(k)=xv3(lid)
          xt3(K)=(xt3(K) + T1_D*tstep_disp*qpl + T2_D*tstep_disp*qpg) / &
               (1.D0+tstep_disp*qpl+tstep_disp*qpg)
       ENDIF

       !-----defining characteristics of the droplet particles
       IF(ind(k) == 2) THEN
          xu3(K)=(xu3(K)+ VLiq(1)*Vmod1 + VGas(1)*Vmod2)*Denom
          xv3(k)=(xv3(k) - tstep_disp*Gravity*(1.D0-RoRel) + &
               VLiq(2)*Vmod1 + VGas(2)*Vmod2)*Denom
          xt3(K)=(xt3(K) + T1_D*tstep_disp*qpl + T2_D*tstep_disp*qpg) / &
               (1.D0+tstep_disp*qpl+tstep_disp*qpg)

          ! Eliminate upward motion!
          xv3(k) = MIN(xv3(k),-VdropMin)
       ENDIF


       IF(ind(k) == 3) THEN
          !-----DEBRIS
          !-----if particles with ind=2 reach the bottom (fragmentation before
          !------the bottom reaching
          !------qpl = ((hr+hc)/(rmelt*cmelt))*6/dm
          !-----VesselRad - radius of the inner vessel
          !-----pi*(r**2)-
          !-----rmelt/qmdebris - total volume of debris
          !-----Sum of qpl = n * pi*(hx**2)/(n *vol of particle)
          !-----For  one particle it is pi*(hx**2)/(vol of particles)
          qpl=qpl*(xd3(k)/6.D0)*pi*(VesselRad**2)*rmelt/qmdebris
          qpg=qpg*(xd3(k)/6.D0)*pi*(VesselRad**2)*rmelt/qmdebris
          xv3(k)=1.d-10
          xu3(k)=1.d-10
          xt3(K)=(xt3(K) + T1_D*tstep_disp*qpl + T2_D*tstep_disp*qpg) / &
               (1.D0+tstep_disp*qpl+tstep_disp*qpg)
       ENDIF

       IF(ind(k) == 4) THEN
          xv3(k)=0.D0
          xu3(k)=0.D0
          xt3(k)=0.D0
          qmk(k)=0.D0
          z(k)=0.D0
          x(k)=0.D0
       ENDIF

       enerd=enerd+ &
            qpl*qmk(k)*cmelt*tstep_disp*(xt3(k)-T1_D) + &
            qpg*qmk(k)*cmelt*tstep_disp*(xt3(k)-T2_D)

       enerd1 = enerd1+ &
            qmk(k)*cmelt*(xt3(k)-xt3old)

       enerd2=enerd2+cmelt*(xt3(k)-tin)*qmk(k)

       r31k(k)=qpl*cmelt*rmelt
       r32k(k)=qpg*cmelt*rmelt


       c31k(k)=c13*Vmod1*qmk(k)*0.75D0/xd3(k)
       c32k(k)=c23*Vmod2*qmk(k)*0.75D0/xd3(k)
    ENDDO LOOP_ALL_PARTICLES

    MOVE_PARTICLES: DO k = 1,ll
       !
       ! Move particle
       !
       x(k) = x(k) + tstep_disp*xu3(k)
       z(k) = z(k) + tstep_disp*xv3(k)
       !
       ! Ensure all particles are inside the domain
       IF(x(k) < 0.D0) THEN
          x(k) = -x(k); xu3(k) = -xu3(k)
       ENDIF
       IF(x(k) > VesselRad) THEN
          x(k) = 2.D0*VesselRad-x(k); xu3(k) = -xu3(k)
       ENDIF
       IF(z(k) > VesselHeight) THEN
          z(k) = 2.D0*VesselHeight-z(k)
          xv3(k) = -xv3(k)
       ENDIF

       ! Avoid collection of particles in the stagnation zone near the axis!
       IF(ind(k) == 2 .AND. x(k) < 2*hi .AND. z(k)> dzv(m-1)) THEN
          xu3(k) = bu2(MIN(3,n),m)
       ENDIF
       !------For jet and droplet particles (ind =1 or 2)
       !------we determine particles
       !------reached the bottom - originating
       !------of the particles with index 3
       IF(ind(k) <3) THEN
          IF(z(k) < debrisHeight) THEN
             qmdebris = qmdebris + qmk(k)
             ind(k)=3
             xu3(k)=0.D0
             xv3(k)=0.D0
             CALL RANDOM_NUMBER(rand_xy)
             x(k)=VesselRad*SQRT(rand_xy(1)) !sqrt takes into account cylindric symmetry
             z(k)=rand_xy(2)*debrisHeight
          ENDIF
       ENDIF
    ENDDO MOVE_PARTICLES
    !------total number of the particles
    ll=ll3
    ! Update the index of the leading edge particle
    lid = -1 ! Initialize 
    zlead = 1.d10 ! Just make sure that search loop is properly started
    DO k=1,ll
       IF(ind(k)==1 .AND. z(k)<zlead) THEN
          zlead = z(k); lid = k
       ENDIF
    ENDDO
    !
    ! Check mass balance for the dispersed phase
    !
    qmjet=0.D0
    qmdrops=0.D0
    qmdebris=0.D0
    qmzero=0.D0
    lljet = 0; lldrops = 0; lldebris = 0; llzero = 0
    MASS_BALANCE: DO k=1,ll
       !-----jet mass
       IF(ind(k) == 1) THEN
          qmjet=qmjet+qmk(k); lljet = lljet+1
          !-----droplets mass
       ELSEIF(ind(k) == 2) THEN
          qmdrops=qmdrops+qmk(k); lldrops = lldrops +1
          !-----debris mass
       ELSEIF(ind(k) == 3) THEN
          qmdebris=qmdebris+qmk(k); lldebris = lldebris+1
          !----- "non-existing" particles - they are
          !------considered in order to avoid new
          !------numeration
       ELSEIF(ind(k) == 4) THEN
          qmzero=qmzero+qmk(k); llzero = llzero +1
       ENDIF
    END DO MASS_BALANCE
    !
    ! Transfer parameters of dispersed phase onto the grid
    !
    CALL DispersedToGrid
    !
    ! Check if we have to perform more time steps
    !
    time_disp = time_disp + tstep_disp
    GOTO 1000
  END SUBROUTINE SolveForDispersedPhase

  SUBROUTINE DispersedToGrid
    USE GLOBALS_2D
    USE GLOBALS
    USE VARIABLES_2D
    USE GRID_2D
    IMPLICIT NONE
    !
    ! Local variables
    !
    INTEGER:: i, j, k
    INTEGER:: iu, ju, iv, jv, ite, jte
    REAL(KIND=8):: xk,zk ! Particle coordinates
    REAL(KIND=8):: su1,su2,su3,su4
    REAL(KIND=8):: sv1,sv2,sv3,sv4
    REAL(KIND=8):: st1,st2,st3,st4
    REAL(KIND=8):: sst1,sst2,sst3,sst4,sumst

    REAL(KIND=8):: Ro3u(n,m9),Ro3v(n9,m)

    ! Zero out all 2D arrays
    al3 = 0.D0; at3 = 0.D0
    al3j = 0.D0; al3dr = 0.D0; al3de = 0.D0
    at3j = 0.D0; at3dr = 0.D0 
    al3Molten = 0.D0; al3drMolten = 0.D0; 
    at3Molten = 0.D0; at3drMolten = 0.D0; 
    r31c =0.D0; r32c = 0.D0

    Ro3u=0.D0; Ro3v=0.D0
    au3 = 0.D0; av3 = 0.D0
    dr3 = 0.D0
    !
    ! Save current values 
    !
    bl3 = al3 
    !
    PARTICLES_LOOP: DO k=1,ll

       xk=x(k)
       zk=z(k)

       IF(zk < 0.D0) CYCLE

       ! Find grid indices
       iu=INT(xk*hi_inv) + 1
       ju=INT( (zk+0.5D0*hj)*hj_inv ) + 1

       su1=(dxu(iu+1)-xk)*(dzu(ju+1)-zk)
       su2=(xk-dxu(iu))*(dzu(ju+1)-zk)
       su3=(dxu(iu+1)-xk)*(zk-dzu(ju))
       su4=(xk-dxu(iu))*(zk-dzu(ju))

       iv=INT( (xk+0.5D0*hi)*hi_inv ) + 1
       jv=INT(zk*hj_inv) + 1

       sv1=(dxv(iv+1)-xk)*(dzv(jv+1)-zk)
       sv2=(xk-dxv(iv))*(dzv(jv+1)-zk)
       sv3=(dxv(iv+1)-xk)*(zk-dzv(jv))
       sv4=(xk-dxv(iv))*(zk-dzv(jv))

       ite=iv
       jte=ju

       st1=(dxt(ite+1)-xk)*(dzt(jte+1)-zk)
       st2=(xk-dxt(ite))*(dzt(jte+1)-zk)
       st3=(dxt(ite+1)-xk)*(zk-dzt(jte))
       st4=(xk-dxt(ite))*(zk-dzt(jte))

       IF(ind(k) < 4) THEN ! All dispersed particles
          al3(ite,jte)=al3(ite,jte)+qmk(k)*st1*s_inv
          al3(ite+1,jte)=al3(ite+1,jte)+qmk(k)*st2*s_inv
          al3(ite,jte+1)=al3(ite,jte+1)+qmk(k)*st3*s_inv
          al3(ite+1,jte+1)=al3(ite+1,jte+1)+qmk(k)*st4*s_inv

          at3(ite,jte)=at3(ite,jte)+xt3(k)*qmk(k)*st1*s_inv
          at3(ite+1,jte)=at3(ite+1,jte)+xt3(k)*qmk(k)*st2*s_inv
          at3(ite,jte+1)=at3(ite,jte+1)+xt3(k)*qmk(k)*st3*s_inv
          at3(ite+1,jte+1)=at3(ite+1,jte+1)+xt3(k)*qmk(k)*st4*s_inv

          dr3(ite,jte)=dr3(ite,jte)+xd3(k)*qmk(k)*st1*s_inv
          dr3(ite+1,jte)=dr3(ite+1,jte)+xd3(k)*qmk(k)*st2*s_inv
          dr3(ite,jte+1)=dr3(ite,jte+1)+xd3(k)*qmk(k)*st3*s_inv
          dr3(ite+1,jte+1)=dr3(ite+1,jte+1)+xd3(k)*qmk(k)*st4*s_inv

          IF(xt3(K) > TMelting) THEN
             al3Molten(ite,jte)=al3Molten(ite,jte)+qmk(k)*st1*s_inv
             al3Molten(ite+1,jte)=al3Molten(ite+1,jte)+qmk(k)*st2*s_inv
             al3Molten(ite,jte+1)=al3Molten(ite,jte+1)+qmk(k)*st3*s_inv
             al3Molten(ite+1,jte+1)=al3Molten(ite+1,jte+1)+qmk(k)*st4*s_inv

             at3Molten(ite,jte)=at3Molten(ite,jte)+xt3(k)*qmk(k)*st1*s_inv
             at3Molten(ite+1,jte)=at3Molten(ite+1,jte)+xt3(k)*qmk(k)*st2*s_inv
             at3Molten(ite,jte+1)=at3Molten(ite,jte+1)+xt3(k)*qmk(k)*st3*s_inv
             at3Molten(ite+1,jte+1)=at3Molten(ite+1,jte+1)+xt3(k)*qmk(k)*st4*s_inv
          ENDIF
          ! Velocities
          Ro3u(iu,ju)=Ro3u(iu,ju)+qmk(k)*su1*s_inv
          Ro3u(iu+1,ju)=Ro3u(iu+1,ju)+qmk(k)*su2*s_inv
          Ro3u(iu,ju+1)=Ro3u(iu,ju+1)+qmk(k)*su3*s_inv
          Ro3u(iu+1,ju+1)=Ro3u(iu+1,ju+1)+qmk(k)*su4*s_inv

          Ro3v(iv,jv)=Ro3v(iv,jv)+qmk(k)*sv1*s_inv
          Ro3v(iv+1,jv)=Ro3v(iv+1,jv)+qmk(k)*sv2*s_inv
          Ro3v(iv,jv+1)=Ro3v(iv,jv+1)+qmk(k)*sv3*s_inv
          Ro3v(iv+1,jv+1)=Ro3v(iv+1,jv+1)+qmk(k)*sv4*s_inv

          au3(iu,ju)=au3(iu,ju) + XU3(k)*qmk(k)*su1*s_inv
          au3(iu+1,ju)=au3(iu+1,ju) + XU3(k)*qmk(k)*su2*s_inv
          au3(iu,ju+1)=au3(iu,ju+1) + XU3(k)*qmk(k)*su3*s_inv
          au3(iu+1,ju+1)=au3(iu+1,ju+1) + XU3(k)*qmk(k)*su4*s_inv

          av3(iv,jv)=av3(iv,jv) + xv3(k)*qmk(k)*sv1*s_inv
          av3(iv+1,jv)=av3(iv+1,jv) + xv3(k)*qmk(k)*sv2*s_inv
          av3(iv,jv+1)=av3(iv,jv+1) + xv3(k)*qmk(k)*sv3*s_inv
          av3(iv+1,jv+1)=av3(iv+1,jv+1) + xv3(k)*qmk(k)*sv4*s_inv
       ENDIF

       IF(ind(k) == 1) THEN ! Jet
          al3j(ite,jte)=al3j(ite,jte)+qmk(k)*st1*s_inv
          al3j(ite+1,jte)=al3j(ite+1,jte)+qmk(k)*st2*s_inv
          al3j(ite,jte+1)=al3j(ite,jte+1)+qmk(k)*st3*s_inv
          al3j(ite+1,jte+1)=al3j(ite+1,jte+1)+qmk(k)*st4*s_inv

          at3j(ite,jte)=at3j(ite,jte)+xt3(k)*qmk(k)*st1*s_inv
          at3j(ite+1,jte)=at3j(ite+1,jte)+xt3(k)*qmk(k)*st2*s_inv
          at3j(ite,jte+1)=at3j(ite,jte+1)+xt3(k)*qmk(k)*st3*s_inv
          at3j(ite+1,jte+1)=at3j(ite+1,jte+1)+xt3(k)*qmk(k)*st4*s_inv
       ENDIF

       IF(ind(k) == 2) THEN ! Droplets
          al3dr(ite,jte)=al3dr(ite,jte)+qmk(k)*st1*s_inv
          al3dr(ite+1,jte)=al3dr(ite+1,jte)+qmk(k)*st2*s_inv
          al3dr(ite,jte+1)=al3dr(ite,jte+1)+qmk(k)*st3*s_inv
          al3dr(ite+1,jte+1)=al3dr(ite+1,jte+1)+qmk(k)*st4*s_inv

          IF(xt3(k) > TMelting) THEN
             al3Molten(ite,jte)=al3Molten(ite,jte)+qmk(k)*st1*s_inv
             al3Molten(ite+1,jte)=al3Molten(ite+1,jte)+qmk(k)*st2*s_inv
             al3Molten(ite,jte+1)=al3Molten(ite,jte+1)+qmk(k)*st3*s_inv
             al3Molten(ite+1,jte+1)=al3Molten(ite+1,jte+1)+qmk(k)*st4*s_inv

             at3drMolten(ite,jte)=at3drMolten(ite,jte)+xt3(k)*qmk(k)*st1*s_inv
             at3drMolten(ite+1,jte)=at3drMolten(ite+1,jte)+xt3(k)*qmk(k)*st2*s_inv
             at3drMolten(ite,jte+1)=at3drMolten(ite,jte+1)+xt3(k)*qmk(k)*st3*s_inv
             at3drMolten(ite+1,jte+1)=at3drMolten(ite+1,jte+1)+xt3(k)*qmk(k)*st4*s_inv
          ENDIF
       ENDIF

       IF(ind(k) == 3) THEN ! Debris
          al3de(ite,jte)=al3de(ite,jte)+qmk(k)*st1*s_inv
          al3de(ite+1,jte)=al3de(ite+1,jte)+qmk(k)*st2*s_inv
          al3de(ite,jte+1)=al3de(ite,jte+1)+qmk(k)*st3*s_inv
          al3de(ite+1,jte+1)=al3de(ite+1,jte+1)+qmk(k)*st4*s_inv
       ENDIF

       IF(ind(k) < 4) THEN ! Heat transfer terms

          sst1 = st1*s_inv; sst2 = st2*s_inv; sst3 = st3*s_inv; sst4 = st4*s_inv
          IF(aa1(ite,jte)/(aa1(ite,jte)+aa2(ite,jte))<0.3D0) sst1 = 0.D0
          IF(aa1(ite+1,jte)/(aa1(ite+1,jte)+aa2(ite+1,jte))<0.3D0) sst2 = 0.D0
          IF(aa1(ite,jte+1)/(aa1(ite,jte+1)+aa2(ite,jte+1))<0.3D0) sst3 = 0.D0
          IF(aa1(ite+1,jte+1)/(aa1(ite+1,jte+1)+aa2(ite+1,jte+1))<0.3D0) sst4 = 0.D0
          sumst = sst1+sst2+sst3+sst4
          IF(sumst > 1.d-5) THEN
             sst1 = sst1/sumst;sst2 = sst2/sumst;sst3 = sst3/sumst;sst4 = sst4/sumst;

             r31c(ite,jte)=r31c(ite,jte) + r31k(k)*qmk(k)*sst1
             r31c(ite+1,jte)=r31c(ite+1,jte) + r31k(k)*qmk(k)*sst2
             r31c(ite,jte+1)=r31c(ite,jte+1) + r31k(k)*qmk(k)*sst3
             r31c(ite+1,jte+1)=r31c(ite+1,jte+1) + r31k(k)*qmk(k)*sst4

          ENDIF

          r32c(ite,jte)=r32c(ite,jte) + r32k(k)*qmk(k)*st1*s_inv
          r32c(ite+1,jte)=r32c(ite+1,jte) + r32k(k)*qmk(k)*st2*s_inv
          r32c(ite,jte+1)=r32c(ite,jte+1) + r32k(k)*qmk(k)*st3*s_inv
          r32c(ite+1,jte+1)=r32c(ite+1,jte+1) + r32k(k)*qmk(k)*st4*s_inv
       ENDIF

    ENDDO PARTICLES_LOOP

    ! Add the mass which was distributed across the axis to the 2-nd cell
    al3(2,:) = al3(1,:)+al3(2,:)
    al3j(2,:) = al3j(1,:)+al3j(2,:)
    al3dr(2,:) = al3dr(1,:)+al3dr(2,:)
    al3de(2,:) = al3de(1,:)+al3de(2,:)
    al3Molten(2,:) = al3Molten(1,:)+al3Molten(2,:)
    al3drMolten(2,:) = al3drMolten(1,:)+al3drMolten(2,:)
    at3Molten(2,:) = at3Molten(1,:)+at3Molten(2,:)
    at3drMolten(2,:) = at3drMolten(1,:)+at3drMolten(2,:)
    at3(2,:) = at3(1,:)+at3(2,:)
    dr3(2,:) = dr3(1,:)+dr3(2,:)
    at3j(2,:) = at3j(1,:)+at3j(2,:)
    r31c(2,:) = r31c(1,:)+r31c(2,:)
    r32c(2,:) = r32c(1,:)+r32c(2,:)

    ! Add the mass which was distributed across the right boundary
    al3(n9,:) = al3(n,:)+al3(n9,:)
    al3j(n9,:) = al3j(n,:)+al3j(n9,:)
    al3dr(n9,:) = al3dr(n,:)+al3dr(n9,:)
    al3de(n9,:) = al3de(n,:)+al3de(n9,:)
    al3Molten(n9,:) = al3Molten(n,:)+al3Molten(n9,:)
    al3drMolten(n9,:) = al3drMolten(n,:)+al3drMolten(n9,:)
    at3Molten(n9,:) = at3Molten(n,:)+at3Molten(n9,:)
    at3drMolten(n9,:) = at3drMolten(n,:)+at3drMolten(n9,:)
    at3(n9,:) = at3(n,:)+at3(n9,:)
    dr3(n9,:) = dr3(n,:)+dr3(n9,:)
    at3j(n9,:) = at3j(n,:)+at3j(n9,:)
    r31c(n9,:) = r31c(n,:)+r31c(n9,:)
    r32c(n9,:) = r32c(n,:)+r32c(n9,:)

    ! Add the mass which was distributed across the debris catcher
    al3(:,2) = al3(:,1)+al3(:,2)
    al3j(:,2) = al3j(:,1)+al3j(:,2)
    al3dr(:,2) = al3dr(:,1)+al3dr(:,2)
    al3de(:,2) = al3de(:,1)+al3de(:,2)
    al3Molten(:,2) = al3Molten(:,1)+al3Molten(:,2)
    al3drMolten(:,2) = al3drMolten(:,1)+al3drMolten(:,2)
    at3Molten(:,2) = at3Molten(:,1)+at3Molten(:,2)
    at3drMolten(:,2) = at3drMolten(:,1)+at3drMolten(:,2)
    at3(:,2) = at3(:,1)+at3(:,2)
    dr3(:,2) = dr3(:,1)+dr3(:,2)
    at3j(:,2) = at3j(:,1)+at3j(:,2)
    r31c(:,2) = r31c(:,1)+r31c(:,2)
    r32c(:,2) = r32c(:,1)+r32c(:,2)

    ! Add the mass which was distributed across the top boundary
    al3(:,m9) = al3(:,m)+al3(:,m9)
    al3j(:,m9) = al3j(:,m)+al3j(:,m9)
    al3dr(:,m9) = al3dr(:,m)+al3dr(:,m9)
    al3de(:,m9) = al3de(:,m)+al3de(:,m9)
    al3Molten(:,m9) = al3Molten(:,m)+al3Molten(:,m9)
    al3drMolten(:,m9) = al3drMolten(:,m)+al3drMolten(:,m9)
    at3Molten(:,m9) = at3Molten(:,m)+at3Molten(:,m9)
    at3drMolten(:,m9) = at3drMolten(:,m)+at3drMolten(:,m9)
    at3(:,m9) = at3(:,m)+at3(:,m9)
    dr3(:,m9) = dr3(:,m)+dr3(:,m9)
    at3j(:,m9) = at3j(:,m)+at3j(:,m9)
    r31c(:,m9) = r31c(:,m)+r31c(:,m9)
    r32c(:,m9) = r32c(:,m)+r32c(:,m9)

    DO j = 2, m
       DO i = 2, n
          ! Temperatures of dispersed particles
          at3(i,j)  = at3(i,j)/(al3(i,j)+1.D-10)
          at3j(i,j) = at3j(i,j)/(al3j(i,j)+1.D-10)
          at3Molten(i,j) = at3Molten(i,j)/(al3Molten(i,j)+1.D-10)
          at3drMolten(i,j) = at3drMolten(i,j)/(al3drMolten(i,j)+1.D-10)
          ! Diameter of dispersed phase
          IF(al3(i,j) > 1.d-10) THEN
             dr3(i,j)  = dr3(i,j)/(al3(i,j)+1.D-10)
          ELSE
             dr3(i,j) = 0.D0
          ENDIF
          !
          ! Convert masses to volume fractions
          al3(i,j) = al3(i,j)*volc(i)/rmelt
          al3j(i,j) = al3j(i,j)*volc(i)/rmelt
          al3dr(i,j) = al3dr(i,j)*volc(i)/rmelt
          al3de(i,j) = al3de(i,j)*volc(i)/rmelt
          al3Molten(i,j) = al3Molten(i,j)*volc(i)/rmelt
          al3drMolten(i,j) = al3drMolten(i,j)*volc(i)/rmelt
          !
          ! Heat transfer coefficients
          r31c(i,j) = r31c(i,j)*volc(i)/rmelt
          r32c(i,j) = r32c(i,j)*volc(i)/rmelt
       ENDDO
    ENDDO
    !
    ! Velocties
    !
    WHERE(Ro3u > 1.D-5)
       au3 = au3/Ro3u
    ELSEWHERE
       au3 = 0.D0
    endwhere

    WHERE(Ro3v > 1.D-5)
       av3 = av3/Ro3v
    ELSEWHERE
       av3 = 0.D0
    endwhere
    !
    ! Boundary conditions
    !
    al3(:,1) = al3(:,2); al3(:,m9) = al3(:,m)
    al3j(:,1) = al3j(:,2); al3j(:,m9) = al3j(:,m)
    al3dr(:,1) = al3dr(:,2); al3dr(:,m9) = al3dr(:,m)
    al3de(:,1) = al3de(:,2); al3de(:,m9) = al3de(:,m)
    at3(:,1) = at3(:,2); at3(:,m9) = at3(:,m)
    at3j(:,1) = at3j(:,2); at3j(:,m9) = at3j(:,m)
    al3Molten(:,1) = al3Molten(:,2); al3Molten(:,m9) = al3Molten(:,m)
    al3drMolten(:,1) = al3drMolten(:,2); al3drMolten(:,m9) = al3drMolten(:,m)
    at3Molten(:,1) = at3Molten(:,1); at3Molten(:,m9) = at3Molten(:,m)
    at3drMolten(:,1) = at3drMolten(:,2); at3drMolten(:,m9) = at3drMolten(:,m)
    r31c(:,1) = r31c(:,2); r31c(:,m9) = r31c(:,m)
    r32c(:,1) = r32c(:,2); r32c(:,m9) = r32c(:,m)

    al3(1,:) = al3(2,:); al3(n9,:) = al3(n,:)
    al3j(1,:) = al3j(2,:); al3j(n9,:) = al3j(n,:)
    al3dr(1,:) = al3dr(2,:); al3dr(n9,:) = al3dr(n,:)
    al3de(1,:) = al3de(2,:); al3de(n9,:) = al3de(n,:)
    at3(1,:) = at3(2,:); at3(n9,:) = at3(n,:)
    at3j(1,:) = at3j(2,:); at3j(n9,:) = at3j(n,:)
    al3Molten(1,:) = al3Molten(2,:); al3Molten(n9,:) = al3Molten(n,:)
    al3drMolten(1,:) = al3drMolten(2,:); al3drMolten(n9,:) = al3drMolten(n,:)
    at3Molten(1,:) = at3Molten(1,:); at3Molten(n9,:) = at3Molten(n,:)
    at3drMolten(1,:) = at3drMolten(2,:); at3drMolten(n9,:) = at3drMolten(n,:)
    r31c(1,:) = r31c(2,:); r31c(n9,:) = r31c(n,:)
    r32c(1,:) = r32c(2,:); r32c(n9,:) = r32c(n,:)

  END SUBROUTINE DispersedToGrid


END MODULE DISPERSED_PHASE_2D
