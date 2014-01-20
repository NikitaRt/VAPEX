!------------------------------------------------------------------!
!  Solver for 2D calculation: momentum equation                    !
!------------------------------------------------------------------!
!  $Id: MOD_Velocity2D.f90 8 2014-01-14 13:38:05Z Sergey $
!------------------------------------------------------------------!
MODULE VELOCITY_2D
  INTEGER,PRIVATE,PARAMETER::isVelDiagonalImpl=0!1
  REAL(KIND=8),PRIVATE, PARAMETER:: Cstab= 1.1d0

CONTAINS

  SUBROUTINE velocity_coefs( &
       u1n, v1n,   & ! u- and v-velocity of 1st phase on n-th time step
       u2n, v2n,   & ! u- and v-velocity of 2nd phase on n-th time step
       u3, v3,     & ! u- and v-velocity of 3rd phase
       bu1,bv1,    & ! u- and v-velocity of 1st phase on previous iteration 
       bu2,bv2,    & ! u- and v-velocity of 2nd phase on previous iteration
       ba1,ba2,    & ! phase volume fractions
       al3,dr3,    & ! 3rd phase volume fraction and diameter
       br1,br2,    & ! phase densities
       bp,         & ! pressure
       gam,        & ! phase change rate
       bmul,bmuv,  & ! viscosities of liquid and gas phases
       au1, av1,   & ! new (guessed) velocities on current iteration
       au2, av2,   & !
       fu1l, fu1r, & ! u1(i,j) = au1(i,j) - (fu1l(i,j)*dp(i,j)+fu1r(i,j)*dp(i+1,j))
       fv1l, fv1r, & ! v1(i,j) = av1(i,j) - (fv1l(i,j)*dp(i,j)+fv1r(i,j)*dp(i,j+1))
       fu2l, fu2r, & ! u2(i,j) = au2(i,j) - (fu2l(i,j)*dp(i,j)+fu2r(i,j)*dp(i+1,j))
       fv2l, fv2r, & ! v2(i,j) = av2(i,j) - (fv2l(i,j)*dp(i,j)+fv2r(i,j)*dp(i,j+1))
       tstep, hi, hj, Gravity, &
       n,m,n9,m9)

    USE GLOBALS
    USE CORRELATIONS
    USE LOCAL_VARS
    IMPLICIT NONE
    !
    ! Parmeters
    !
    REAL(KIND=8), INTENT(in):: u1n(n,m9), v1n(n9,m)
    REAL(KIND=8), INTENT(in):: u2n(n,m9), v2n(n9,m)
    REAL(KIND=8), INTENT(in):: u3(n,m9),  v3(n9,m)
    REAL(KIND=8), INTENT(in):: bu1(n,m9), bv1(n9,m)
    REAL(KIND=8), INTENT(in):: bu2(n,m9), bv2(n9,m)

    REAL(KIND=8), INTENT(in):: ba1(n9,m9),ba2(n9,m9)
    REAL(KIND=8), INTENT(in):: al3(n9,m9),dr3(n9,m9)
    REAL(KIND=8), INTENT(in):: br1(n9,m9),br2(n9,m9)

    REAL(KIND=8), INTENT(in):: bp(n9,m9), gam(n9,m9)
    REAL(KIND=8), INTENT(in):: bmuv(n9,m9), bmul(n9,m9)

    REAL(KIND=8), INTENT(out):: au1(n,m9), av1(n9,m)
    REAL(KIND=8), INTENT(out):: au2(n,m9), av2(n9,m)
    REAL(KIND=8), INTENT(out):: fu1l(n,m9), fu1r(n,m9)
    REAL(KIND=8), INTENT(out):: fu2l(n,m9), fu2r(n,m9)
    REAL(KIND=8), INTENT(out):: fv1l(n9,m), fv1r(n9,m)
    REAL(KIND=8), INTENT(out):: fv2l(n9,m), fv2r(n9,m)

    REAL(KIND=8),INTENT(in):: tstep, hi, hj, Gravity

    INTEGER, INTENT(in):: n,m,n9,m9
    !
    ! Local variables
    !
    REAL(KIND=8):: ur1, vz1
    REAL(KIND=8):: ur2, vz2
    REAL(KIND=8):: ur3, vz3
    REAL(KIND=8):: rr1, ra1, rr2, ra2, ra3, d3
    REAL(KIND=8):: gamma, gamp, gamm
    REAL(KIND=8):: vmod12, vmod13, vmod23
    REAL(KIND=8):: ur_diag, ur_off, vr_diag, vr_off, cr_diag
    REAL(KIND=8):: vz_diag, vz_off, uz_diag, uz_off, cz_diag 

    REAL(KIND=8):: AAU1, BBU1, CCU1
    REAL(KIND=8):: AAU2, BBU2, CCU2
    REAL(KIND=8):: AAV1, BBV1, CCV1
    REAL(KIND=8):: AAV2, BBV2, CCV2

    REAL(KIND=8):: x1, x2, Det

    INTEGER:: i,j

    REAL(KIND=8):: dF12dU, dF13dU, dF23dU
    REAL(KIND=8):: dF12dV, dF13dV, dF23dV
    REAL(KIND=8):: Stab
    REAL(KIND=8):: cpr, cpz, ctg

    cpr = tstep/hi; cpz = tstep/hj; ctg = tstep*Gravity

    !----------------------- U-velocity --------------------------------
    DO i=2,n-1
       DO j=2,m

          ur1= bu1(i,j)
          vz1=(bv1(i,j-1)+bv1(i,j)+bv1(i+1,j-1)+bv1(i+1,j))*0.25D0

          ur2= bu2(i,j)
          vz2=(bv2(i,j-1)+bv2(i,j)+bv2(i+1,j-1)+bv2(i+1,j))*0.25D0

          ur3= u3(i,j)
          vz3=(v3(i,j-1)+v3(i,j)+v3(i+1,j-1)+v3(i+1,j))*0.25D0

          rr1=(br1(i,j)+br1(i+1,j))*0.5D0
          ra1=(ba1(i,j)+ba1(i+1,j))*0.5D0
          ra1=MAX(ra1,0.D0)

          rr2=(br2(i,j)+br2(i+1,j))*0.5D0
          ra2=(ba2(i,j)+ba2(i+1,j))*0.5D0
          ra2=MAX(ra2,0.D0)

          ra3=(al3(i,j)+al3(i+1,j))*0.5D0
          ra3=MAX(ra3,0.D0)

          d3 = (al3(i,j)*dr3(i,j)+al3(i+1,j)*dr3(i+1,j))/&
               MAX(al3(i,j)+al3(i+1,j),allim)

          gamma = (gam(i,j)+gam(i+1,j))*0.5D0
          IF(gamma >= 0) THEN
             gamp = gamma; gamm = 0.D0
          ELSE
             gamp = 0.D0;  gamm = gamma
          ENDIF

          vmod12 = SQRT( (ur2-ur1)**2 + (vz2-vz1)**2 )
          vmod13 = SQRT( (ur3-ur1)**2 + (vz3-vz1)**2 )
          vmod23 = SQRT( (ur3-ur2)**2 + (vz3-vz2)**2 )
          !
          ! Liquid-gas drag coefficient
          !
          Alpha = ra2
          Alpha_l = ra1
          Alpha_3 = ra3
          Ro1 = rr1
          Ro2 = rr2
          VLiq(1) = ur1; VLiq(2) = vz1
          VGas(1) = ur2; VGas(2) = vz2
          VDisp(1) = ur3; VDisp(2) = vz3
          Visc_Liq = 0.5D0*(bmul(i+1,j)+bmul(i,j))
          Visc_Gas = 0.5D0*(bmuv(i+1,j)+bmuv(i,j))
          DDispPhase = d3
          !
          !
          CALL DragGasLiq
          CALL DragGasLiqDispersed
          !
          !
          !---------------------------- Water ---------------------------------
          ! AAU1*U1(n+1) + BBU1*U2(n+1) = CCU1 - cpr/Ro_1*(dP(i+1)-dP(i))
          ra1 = MAX(ra1,allim); ra2 = MAX(ra2,allim)
          ! Collecting coefficients AAU1, BBU1, CCU1       
          AAU1 = 0
          BBU1 = 0
          CCU1 = 0
          ! Adding mass transfer terms
          AAU1 = AAU1 - gamm/ra1
          BBU1 = BBU1 + gamm/ra1
          ! Adding drag between phases 1 and 2
          dF12dU = CDrag12_a1*(vmod12 + (ur2-ur1)**2/(vmod12+1.d-30))
          AAU1 = AAU1 + dF12dU
          BBU1 = BBU1 - dF12dU
          CCU1 = CCU1 - CDrag12_a1*(ur2-ur1)**3/(vmod12+1.d-30)
          ! Adding drag between phases 1 and 3
          dF13dU = C13_a1*(vmod13 + (ur3-ur1)**2/(vmod13+1.d-30))
          AAU1 = AAU1 + dF13dU
          CCU1 = CCU1 + C13_a1*(vmod13*ur3 + &
               ur1*(ur3-ur1)**2/(vmod13+1.d-30))
          ! Multiplying by t/(Alpha_1*Ro_1)
          AAU1 = AAU1*tstep/rr1
          BBU1 = BBU1*tstep/rr1
          CCU1 = CCU1*tstep/rr1
          ! Adding convective terms
          IF(ur1<=0.d0) THEN
             ur_off = bu1(i+1,j)
             ur_diag = -bu1(i,j)
             cr_diag = -1.D0
          ELSE
             ur_off = -bu1(i-1,j)
             ur_diag = bu1(i,j)
             cr_diag = 1.D0
          ENDIF

          IF(vz1<=0.d0) THEN
             uz_off = bu1(i,j+1)
             uz_diag = -bu1(i,j)
             cz_diag = -1.D0
          ELSE
             uz_off = -bu1(i,j-1)
             uz_diag = bu1(i,j)
             cz_diag = 1.D0
          ENDIF

          IF(isVelDiagonalImpl == 0) THEN ! Fully explicit scheme
             CCU1 = CCU1 - cpr*ur1*(ur_diag+ur_off) - &
                  cpz*vz1*(uz_diag+uz_off)
          ELSE                ! Implicit treatment of diagonal elements
             AAU1 = AAU1 + cpr*ur1*cr_diag + cpz*vz1*cz_diag
             CCU1 = CCU1 - cpr*ur1*ur_off - cpz*vz1*uz_off
          ENDIF
          ! Adding pressure gradient
          CCU1 = CCU1 - (bp(i+1,j)-bp(i,j))*cpr/rr1
          ! Adding time derivative
          CCU1 = CCU1 + u1n(i,j)

          ! Adding stabilizing terms:        
          IF(ra1+ra2 > allim) THEN
             Stab = ra2*rr2/(ra2*rr1+ra1*rr2)*vmod12**2*cpr*&
                  (ba2(i+1,j)/MAX(ba1(i+1,j)+ba2(i+1,j),allim)-&
                  ba2(i,j)/MAX(ba1(i,j)+ba2(i,j),allim))
             CCU1 = CCU1 - Stab*Cstab
          ENDIF

          !---------------------------- Vapor ---------------------------------
          ! AAU2*U1(n+1) + BBU2*U2(n+1) = CCU2 - cpr/Ro_2*(dP(i+1)-dP(i))

          ! Collecting coefficients AAU2, BBU2, CCU2
          AAU2 = 0
          BBU2 = 0
          CCU2 = 0
          ! Adding mass transfer terms
          AAU2 = AAU2 - gamp/ra2
          BBU2 = BBU2 + gamp/ra2
          ! Adding drag between phases 1 and 2
          ! Adding drag between phases 1 and 2
          dF12dU = CDrag12_a2*(vmod12 + (ur2-ur1)**2/(vmod12+1.d-30))
          AAU2 = AAU2 - dF12dU
          BBU2 = BBU2 + dF12dU
          CCU2 = CCU2 - CDrag12_a2*(ur1-ur2)**3/(vmod12+1.d-30)
          ! Adding drag between phases 2 and 3
          dF23dU = C23_a2*(vmod23 + (ur3-ur2)**2/(vmod23+1.d-30))
          BBU2 = BBU2 + dF23dU
          CCU2 = CCU2 + C23_a2*(vmod23*ur3 + &
               ur2*(ur3-ur2)**2/(vmod23+1.d-30))
          ! Multiplying by t/(Alpha_2*Ro_2)
          AAU2 = AAU2*tstep/rr2
          BBU2 = BBU2*tstep/rr2
          CCU2 = CCU2*tstep/rr2
          ! Adding convective terms
          IF(ur2<=0.d0) THEN
             ur_off = bu2(i+1,j)
             ur_diag = -bu2(i,j)
             cr_diag = -1.D0
          ELSE
             ur_off = -bu2(i-1,j)
             ur_diag = bu2(i,j)
             cr_diag = 1.D0
          ENDIF

          IF(vz2<=0.d0) THEN
             uz_off = bu2(i,j+1)
             uz_diag = -bu2(i,j)
             cz_diag = -1.D0
          ELSE
             uz_off = -bu2(i,j-1)
             uz_diag = bu2(i,j)
             cz_diag = 1.D0
          ENDIF
          IF(isVelDiagonalImpl == 0) THEN ! Fully explicit scheme
             CCU2 = CCU2 - cpr*ur2*(ur_diag+ur_off) - &
                  cpz*vz2*(uz_diag+uz_off)
          ELSE                ! Implicit treatment of diagonal elements
             AAU2 = AAU2 + cpr*ur2*cr_diag + cpz*vz2*cz_diag
             CCU2 = CCU2 - cpr*ur2*ur_off  - cpz*vz2*uz_off
          ENDIF
          ! Adding pressure gradient
          CCU2 = CCU2 - (bp(i+1,j)-bp(i,j))*cpr/rr2
          ! Adding time derivative
          CCU2 = CCU2 + u2n(i,j)

          ! Adding stabilizing terms:        
          IF(ra1+ra2 > allim) THEN
             Stab = ra1*rr1/(ra2*rr1+ra1*rr2)*vmod12**2*cpr*&
                  (ba2(i+1,j)/MAX(ba1(i+1,j)+ba2(i+1,j),allim)-&
                  ba2(i,j)/MAX(ba1(i,j)+ba2(i,j),allim))
             CCU2 = CCU2 + Stab*Cstab
          ENDIF
          !
          ! Calculation of coefficients 
          !
          x1=1.d0+AAU1
          x2=1.d0+BBU2
          Det=x1*x2 - BBU1*AAU2

          au1(i,j)=(CCU1*x2 - BBU1*CCU2) / Det
          fu1l(i,j)=(x2/rr1 - BBU1/rr2)*cpr/Det
          fu1r(i,j)=fu1l(i,j)

          au2(i,j)=(CCU2*x1 - AAU2*CCU1) / Det
          fu2l(i,j)=(x1/rr2 - AAU2/rr1)*cpr/Det
          fu2r(i,j)=fu2l(i,j)
       ENDDO
    ENDDO
    !
    ! Set boundary values for pressure coefficients
    !
    fu1l(1,:) = 2.*cpr/(br1(1,:)+br1(2,:))
    fu1r(1,:) = fu1l(1,:)
    fu1l(n,:) = 2.*cpr/(br1(n9,:)+br1(n,:))
    fu1r(n,:) = fu1l(n,:)
    !
    fu2l(1,:) = 2.*cpr/(br2(1,:)+br2(2,:))
    fu2r(1,:) = fu2l(1,:)
    fu2l(n,:) = 2.*cpr/(br2(n9,:)+br2(n,:))
    fu2r(n,:) = fu2l(n,:)

!!$    fu1l(1,:) = 0.D0
!!$    fu1r(1,:) = 0.D0
!!$    fu1l(n,:) = 0.D0
!!$    fu1r(n,:) = 0.D0
!!$    !
!!$    fu2l(1,:) = 0.D0
!!$    fu2r(1,:) = 0.D0
!!$    fu2l(n,:) = 0.D0
!!$    fu2r(n,:) = 0.D0
    !
    !----------------------- V-velocity --------------------------------
    !
    DO j=2,m-1
       DO i=2,n

          ur1=(bu1(i-1,j)+bu1(i,j)+bu1(i-1,j+1)+bu1(i,j+1))*0.25D0
          vz1=bv1(i,j)

          ur2=(bu2(i-1,j)+bu2(i,j)+bu2(i-1,j+1)+bu2(i,j+1))*0.25D0
          vz2=bv2(i,j)

          ur3=(u3(i-1,j)+u3(i,j)+u3(i-1,j+1)+u3(i,j+1))*0.25D0
          vz3=v3(i,j)

          rr1=(br1(i,j)+br1(i,j+1))*0.5D0
          ra1=(ba1(i,j)+ba1(i,j+1))*0.5D0
          ra1=MAX(ra1,0.D0)

          rr2=(br2(i,j)+br2(i,j+1))*0.5D0
          ra2=(ba2(i,j)+ba2(i,j+1))*0.5D0
          ra2=MAX(ra2,0.D0)

          ra3=(al3(i,j)+al3(i,j+1))*0.5D0
          ra3=MAX(ra3,0.D0)

          d3 = (al3(i,j)*dr3(i,j)+al3(i,j+1)*dr3(i,j+1))/&
               MAX(al3(i,j)+al3(i,j+1),allim)

          gamma=(gam(i,j)+gam(i,j+1))*0.5D0
          IF(gamma >= 0) THEN
             gamp = gamma; gamm = 0.D0
          ELSE
             gamp = 0.D0;  gamm = gamma
          ENDIF

          vmod12 = dsqrt( (ur2-ur1)**2 + (vz2-vz1)**2 )
          vmod13 = dsqrt( (ur3-ur1)**2 + (vz3-vz1)**2 )
          vmod23 = dsqrt( (ur3-ur2)**2 + (vz3-vz2)**2 )

          ! Liquid-gas drag coefficient

          Alpha = ra2
          Alpha_l = ra1
          Alpha_3 = ra3
          Ro1 = rr1
          Ro2 = rr2
          VLiq(1) = ur1; VLiq(2) = vz1
          VGas(1) = ur2; VGas(2) = vz2
          VDisp(1) = ur3; VDisp(2) = vz3
          Visc_Liq = 0.5D0*(bmul(i,j+1)+bmul(i,j))
          Visc_Gas = 0.5D0*(bmuv(i,j+1)+bmuv(i,j))
          DDispPhase = d3
          !
          !
          CALL DragGasLiq
          CALL DragGasLiqDispersed
          !
          !
          !----------------------------Water ---------------------------------
          ! AAV1*V1(n+1) + BBV1*V2(n+1) = CCV1 - cpz/Ro_1*(dP(i+1)-dP(i))
          ra1 = MAX(ra1,allim); ra2 = MAX(ra2,allim)
          ! Collecting coefficients AAV1, BBV1, CCV1
          AAV1 = 0
          BBV1 = 0
          CCV1 = 0
          ! Adding mass transfer terms
          AAV1 = AAV1 - gamm/ra1
          BBV1 = BBV1 + gamm/ra1
          ! Adding drag between phases 1 and 2
          dF12dV = CDrag12_a1*(vmod12 + (vz2-vz1)**2/(vmod12+1.d-30))
          AAV1 = AAV1 + dF12dV
          BBV1 = BBV1 - dF12dV
          CCV1 = CCV1 - CDrag12_a1*(vz2-vz1)**3/(vmod12+1.d-30)
          ! Adding drag between phases 1 and 3
          dF13dV = C13_a1*(vmod13 + (vz3-vz1)**2/(vmod13+1.d-30))
          AAV1 = AAV1 + dF13dV
          CCV1 = CCV1 + C13_a1*(vmod13*vz3 + &
               vz1*(vz3-vz1)**2/(vmod13+1.d-30))
          ! Multiplying by t/(Alpha_1*Ro_1)
          AAV1 = AAV1*tstep/rr1
          BBV1 = BBV1*tstep/rr1
          CCV1 = CCV1*tstep/rr1
          ! Adding convective terms
          IF(ur1<=0.d0) THEN
             vr_off  = bv1(i+1,j)
             vr_diag = -bv1(i,j)
             cr_diag = -1.D0
          ELSE
             vr_off  = -bv1(i-1,j)
             vr_diag = bv1(i,j)
             cr_diag = 1.D0
          ENDIF

          IF(vz1<=0.d0) THEN
             vz_off  = bv1(i,j+1)
             vz_diag = -bv1(i,j)
             cz_diag = -1.D0
          ELSE
             vz_off = -bv1(i,j-1)
             vz_diag = bv1(i,j)
             cz_diag = 1.D0
          ENDIF

          IF(isVelDiagonalImpl == 0) THEN ! Fully explicit scheme
             CCV1 = CCV1 - cpr*ur1*(vr_diag+vr_off) - &
                  cpz*vz1*(vz_diag+vz_off)
          ELSE                ! Implicit treatment of diagonal elements
             AAV1 = AAV1 + cpr*ur1*cr_diag + cpz*vz1*cz_diag
             CCV1 = CCV1 - cpr*ur1*vr_off - cpz*vz1*vz_off
          ENDIF
          ! Adding pressure gradient
          CCV1 = CCV1 - (bp(i,j+1)-bp(i,j))*cpz/rr1
          ! Adding gravity
          CCV1 = CCV1 - ctg
          ! Adding time derivative
          CCV1 = CCV1 + v1n(i,j)

          ! Adding stabilizing terms:        
          IF(ra1+ra2 > allim) THEN
             Stab = ra2*rr2/(ra2*rr1+ra1*rr2)*vmod12**2*cpz*&
                  (ba2(i,j+1)/MAX(ba1(i,j+1)+ba2(i,j+1),allim)-&
                  ba2(i,j)/MAX(ba1(i,j)+ba2(i,j),allim))
             CCV1 = CCV1 - Stab*Cstab
          ENDIF
          !---------------------------- Vapor ---------------------------------
          ! AAV2*V1(n+1) + BBV2*V2(n+1) = CCV2 - cpz/Ro_2*(dP(i+1)-dP(i))

          ! Collecting coefficients AAV2, BBV2, CCV2
          AAV2 = 0
          BBV2 = 0
          CCV2 = 0
          ! Adding mass transfer terms
          AAV2 = AAV2 - gamp/ra2
          BBV2 = BBV2 + gamp/ra2
          ! Adding drag between phases 1 and 2
          dF12dV = CDrag12_a2*(vmod12 + (vz2-vz1)**2/(vmod12+1.d-30))
          AAV2 = AAV2 - dF12dV
          BBV2 = BBV2 + dF12dV
          CCV2 = CCV2 - CDrag12_a2*(vz1-vz2)**3/(vmod12+1.d-30)
          ! Adding drag between phases 2 and 3
          dF23dV = C23_a2*(vmod23 + (vz3-vz2)**2/(vmod23+1.d-30))
          BBV2 = BBV2 + dF23dV
          CCV2 = CCV2 + C23_a2*(vmod23*vz3 + &
               vz2*(vz3-vz2)**2/(vmod23+1.d-30))
          ! Multiplying by t/(Alpha_2*Ro_2)
          AAV2 = AAV2*tstep/rr2
          BBV2 = BBV2*tstep/rr2
          CCV2 = CCV2*tstep/rr2
          ! Adding convective terms
          IF(ur2<=0.d0) THEN
             vr_off  = bv2(i+1,j)
             vr_diag = -bv2(i,j)
             cr_diag = -1.D0
          ELSE
             vr_off  = -bv2(i-1,j)
             vr_diag = bv2(i,j)
             cr_diag = 1.D0
          ENDIF

          IF(vz2<=0.d0) THEN
             vz_off  = bv2(i,j+1)
             vz_diag = -bv2(i,j)
             cz_diag = -1.D0
          ELSE
             vz_off = -bv2(i,j-1)
             vz_diag = bv2(i,j)
             cz_diag = 1.D0
          ENDIF
          IF(isVelDiagonalImpl == 0) THEN ! Fully explicit scheme
             CCV2 = CCV2 - cpr*ur2*(vr_diag+vr_off) - &
                  cpz*vz2*(vz_diag+vz_off)
          ELSE                ! Implicit treatment of diagonal elements
             AAV2 = AAV2 + cpr*ur2*cr_diag + cpz*vz2*cz_diag
             CCV2 = CCV2 - cpr*ur2*vr_off  - cpz*vz2*vz_off
          ENDIF
          ! Adding pressure gradient
          CCV2 = CCV2 - (bp(i,j+1)-bp(i,j))*cpz/rr2
          ! Adding gravity
          CCV2 = CCV2 - ctg
          ! Adding time derivative
          CCV2 = CCV2 + v2n(i,j)

          ! Adding stabilizing terms:        
          IF(ra1+ra2 > allim) THEN
             Stab = ra1*rr1/(ra2*rr1+ra1*rr2)*vmod12**2*cpz*&
                  (ba2(i,j+1)/MAX(ba1(i,j+1)+ba2(i,j+1),allim)-&
                  ba2(i,j)/MAX(ba1(i,j)+ba2(i,j),allim))
             CCV2 = CCV2 + Stab*Cstab
          ENDIF
          !------------------- Calculation of coefficients ----------------------
          x1=1.d0+AAV1
          x2=1.d0+BBV2
          Det=x1*x2 - BBV1*AAV2

          av1(i,j)=(CCV1*x2 - BBV1*CCV2) / Det
          fv1l(i,j)=( x2/rr1 - BBV1/rr2 )*cpz/Det
          fv1r(i,j)=fv1l(i,j)

          av2(i,j)=(CCV2*x1-AAV2*CCV1) / Det
          fv2l(i,j)=( x1/rr2 - AAV2/rr1 )*cpz/Det
          fv2r(i,j)=fv2l(i,j)
          !-----------------------------------------------------------------------
       ENDDO
    ENDDO
    !
    ! Set boundary values for pressure coefficients
    !
    fv1l(:,1) = 2.*cpz/(br1(:,1)+br1(:,2))
    fv1r(:,1) = fv1l(:,1)
    fv1l(:,m) = 2.*cpz/(br1(:,m9)+br1(:,m))
    fv1r(:,m) = fv1l(:,m)
    !
    fv2l(:,1) = 2.*cpz/(br2(:,1)+br2(:,2))
    fv2r(:,1) = fv2l(:,1)
    fv2l(:,m) = 2.*cpz/(br2(:,m9)+br2(:,m))
    fv2r(:,m) = fv2l(:,m)

  END SUBROUTINE velocity_coefs
  !-----------------------------------------------------------------------

END MODULE VELOCITY_2D
