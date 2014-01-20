!---------------------------------------------------------------------!
!  Calculation of thermodynamic properties and their                  !
!  derivatives for water/vapour/non-condensable mixture               !
!---------------------------------------------------------------------!
!  $Id: MOD_WaterProps.f90 5 2014-01-13 11:44:50Z Sergey $
!---------------------------------------------------------------------!
MODULE WATER_PROPS
  !
  ! Liquid water
  !
  REAL(kind=8):: Ro_1              ! Density(P,T1)
  REAL(kind=8):: dRo1_dT1, dRo1_dP ! and its derivatives
  REAL(kind=8):: E_1               ! Internal energy(T1)
  REAL(kind=8):: dE1_dT1, dE1_dP   ! and its derivatives
  !
  ! Combined gas (steam+non-condensable)
  !
  ! Ro2 = Roa + Rov; Roa = Roa(T2,Pa), Rov = Rov(T2,P,Pa)
  ! Ro2*E2 = Roa*Ea + Rov*Ev; Ea = Ea(T2,Pa), Ev = Ev(T2,P,Pa)
  !
  REAL(kind=8):: Ro_2              ! Density = Ro_v + Ro_a
  REAL(kind=8):: dRo2_dT2, dRo2_dP, dRo2_dPa ! and its derivatives
  REAL(kind=8):: E_2               ! Internal energy(T2)
  REAL(kind=8):: dE2_dT2, dE2_dP, dE2_dPa ! and its derivatives
  !
  ! Vapour (steam)
  !
  REAL(kind=8):: Ro_v
  REAL(kind=8):: dRov_dT2, dRov_dP, dRov_dPa
  REAL(kind=8):: E_v
  REAL(kind=8):: dEv_dT2, dEv_dP, dEv_dPa
  !
  ! Non-condensable gas
  !
  REAL(kind=8):: Ro_a              ! Density
  REAL(kind=8):: dRoa_dT2, dRoa_dPa! and its derivatives
  REAL(kind=8):: E_a               ! Internal energy
  REAL(kind=8):: dEa_dT2, dEa_dPa, dEa_dP ! and its derivatives
  !
  ! Mass fraction of non-condensible
  !
  REAL(kind=8):: Y_a ! Ro_a/Ro_2
  REAL(kind=8):: dYa_dT2, dYa_dP, dYa_dPa
  !
  ! Saturation temperature corresponding to total pressure P
  !
  REAL(kind=8):: T_sat, dTsat_dP
  !
  ! Saturation temperature corresponding to partial pressure Pv
  !
  REAL(kind=8):: Tsv,dTsv_dPv,dTsv_dP,dTsv_dPa
  REAL(kind=8), PARAMETER :: dpsdp = 1.D0
  REAL(kind=8), PARAMETER :: dpsdpa = -1.D0
  !
  ! Saturation enthalpies
  !
  REAL(kind=8):: hl_sat, dhl_sat_dP, dhl_sat_dPa
  REAL(kind=8):: hv_sat, dhv_sat_dP, dhv_sat_dPa
  REAL(kind=8):: dhv_sat_dPv
  !
  ! Work data
  !
  REAL(kind=8),PRIVATE:: dr(19)
  INTEGER, PARAMETER:: iop=0
  INTEGER, PARAMETER:: jstart=1
  INTEGER, PARAMETER:: ncells=1

  REAL(kind=8),PRIVATE:: ceoslp(40)
  REAL(kind=8),PRIVATE:: aeos14,ceos1 ,ceos2 , ceos3
  ! Correction to thermo: hl(P,Tsv) depends on P and Pa!
  REAL(kind=8),PRIVATE:: thermo_dhlsdpa

  INTEGER,PRIVATE:: igas, isGasH2mixture
  DATA igas/-1/,isGasH2mixture/-1/

  PRIVATE:: thermo,error,hev,star,seteos,e_H2,h_H2,Mixture_Gas_H2

  INTEGER,PARAMETER:: IGAS_AIR=1, IGAS_HYDROGEN=2, &
       & IGAS_HELIUM=3, IGAS_ARGON=4

  REAL(kind=8),PARAMETER,PRIVATE:: airmol_Air=28.96461d0
  REAL(kind=8),PARAMETER,PRIVATE:: airmol_H2=2.01594d0
  REAL(kind=8),PARAMETER,PRIVATE:: airmol_He=4.00260d0
  REAL(kind=8),PARAMETER,PRIVATE:: airmol_Ar=40.01594d0
  REAL(kind=8),PARAMETER,PRIVATE:: Cp_Air=1004.832d0
  REAL(kind=8),PARAMETER,PRIVATE:: Cp_H2=14533.2d0
  REAL(kind=8),PARAMETER,PRIVATE:: Cp_He=5234.0d0
  REAL(kind=8),PARAMETER,PRIVATE:: Cp_Ar=525.d0 

  REAL(kind=8),PARAMETER,PRIVATE:: Cp_NonCondensGas(4) = &
       (/Cp_Air,Cp_H2,Cp_He,Cp_Ar/)
  REAL(kind=8),PARAMETER,PRIVATE:: W_NonCondensGas(4)= &
       (/airmol_Air,airmol_H2,airmol_He,airmol_Ar/)

  REAL(kind=8),PARAMETER,PRIVATE:: gascon=8.314339d+03 ! Universal gas constant
  REAL(kind=8),PRIVATE:: airmol
CONTAINS
  !
  ! Public functions and subroutines
  !
  SUBROUTINE InitWaterProps
    IF(igas < 0) THEN
       PRINT *,'Wrong gas type in InitWaterProps: ',igas
       STOP
    ENDIF
    CALL seteos
  END SUBROUTINE InitWaterProps
  !
  SUBROUTINE SetNonCondensGasType(itype,isMixture)
    INTEGER, INTENT(in):: itype, isMixture
    IF(itype<IGAS_AIR .OR. itype>IGAS_ARGON) THEN
       PRINT *,'Wrong gas type in SetNonCondensGasType: ',itype
       STOP
    ENDIF
    igas = itype
    IF(isMixture /= 0) THEN
       isGasH2mixture = 1
    ELSE
       isGasH2mixture = 0
    ENDIF
  END SUBROUTINE SetNonCondensGasType
  !
  !
  SUBROUTINE Thermo_Props(Tliq,Tgas,Ptot,Pair)
    USE LOCAL_VARS
    USE HYDROGEN

    IMPLICIT NONE
    REAL(kind=8), OPTIONAL,INTENT(in):: Tliq
    REAL(kind=8), OPTIONAL,INTENT(in):: Tgas
    REAL(kind=8), OPTIONAL,INTENT(in):: Ptot
    REAL(kind=8), OPTIONAL,INTENT(in):: Pair
    !
    ! Dummy arrays for calling thermo
    !
    REAL(kind=8):: xP(1),xPa(1),xT1(1),xT2(1)
    REAL(kind=8):: xE_1(1),xE_2(1),xT_sat(1),&
         & xRo_1(1),xRo_2(1),xRo_a(1),&
         & xTsv(1),xE_a(1)
    !
    ! Use parameters if given
    !
    IF(PRESENT(Tliq)) T1 = Tliq
    IF(PRESENT(Tgas)) T2 = Tgas
    IF(PRESENT(Ptot)) P  = Ptot
    IF(PRESENT(Pair)) Pa = Pair
    !
    xP(1)  = P
    xPa(1) = Pa
    xT1(1) = T1
    xT2(1) = T2

    CALL Mixture_Gas_H2(YH2,XH2,W_a)
    CALL thermo(xP,xE_1,xE_2,xT1,xT2,xT_sat,xRo_1,xRo_2,xPa,xRo_a,&
         &       xTsv,xE_a)

    E_1 = xE_1(1); E_2 = xE_2(1); T_sat = xT_sat(1)
    Ro_1 = xRo_1(1); Ro_2 = xRo_2(1); Ro_a = xRo_a(1)
    Tsv = xTsv(1); E_a = xE_a(1)

    dTsat_dP = dr(1)

    Ro_v = Ro_2-Ro_a
    E_v = (Ro_2*E_2-Ro_a*E_a)/Ro_v

    dRo1_dT1 = dr(8)
    dRo1_dP  = dr(6)

    dE1_dT1 = dr(4)
    dE1_dP  = dr(2)

    dRoa_dT2 = dr(18)
    dEa_dT2  = dr(15)
    dRov_dT2 = dr(9)
    dEv_dT2  = dr(5)
    dRo2_dT2 = dRoa_dT2 + dRov_dT2

    dRoa_dPa = dr(17)
    dEa_dPa  = dr(16)
    dRov_dPa = dr(7)*dpsdpa
    dEv_dPa  = dr(3)*dpsdpa
    dRo2_dPa = dRoa_dPa + dRov_dPa

    dRov_dP  = dr(7)*dpsdp
    dEv_dP   = dr(3)*dpsdp
    dRo2_dP  = dRov_dP ! Roa does not depend on P !
    dEa_dP   = 0.D0

    dE2_dT2 = (Ro_v*dEv_dT2 + Ro_a*dEa_dT2 + &
         E_a*dRoa_dT2 + E_v*dRov_dT2 - E_2*dRo2_dT2)/Ro_2

    dE2_dP = (Ro_v*dEv_dP  + E_v*dRov_dP - E_2*dRo2_dP)/Ro_2

    dE2_dPa = (Ro_v*dEv_dPa + Ro_a*dEa_dPa + &
         E_a*dRoa_dPa + E_v*dRov_dPa - E_2*dRo2_dPa)/Ro_2

    Y_a = Ro_a/Ro_2
    dYa_dT2 = (dRoa_dT2 - dRo2_dT2*Ro_a/Ro_2)/Ro_2
    dYa_dP  = -dRo2_dP*Ro_a/Ro_2**2
    dYa_dPa = (dRoa_dPa - dRo2_dPa*Ro_a/Ro_2)/Ro_2

    hv_sat=dr(10)
    dhv_sat_dPv = dr(12)
    dhv_sat_dP  = dhv_sat_dPv*dpsdp
    dhv_sat_dPa = dhv_sat_dPv*dpsdpa

    hl_sat=dr(11)
    dhl_sat_dP = dr(13)
    dhl_sat_dPa = thermo_dhlsdpa

    dTsv_dPv = dr(14)
    dTsv_dP = dTsv_dPv*dpsdp
    dTsv_dPa = dTsv_dPv*dpsdpa

  END SUBROUTINE Thermo_Props

  REAL(kind=8) FUNCTION EOSPmin()
    EOSPmin = ceoslp(30)
  END FUNCTION EOSPmin
  REAL(kind=8) FUNCTION EOSPmax()
    EOSPmax = ceoslp(31)
  END FUNCTION EOSPmax
  REAL(kind=8) FUNCTION EOSTLiqMin()
    EOSTLiqMin = ceoslp(32)
  END FUNCTION EOSTLiqMin
  REAL(kind=8) FUNCTION EOSTLiqMax()
    EOSTLiqMax = ceoslp(33)
  END FUNCTION EOSTLiqMax
  REAL(kind=8) FUNCTION EOSTVapMin()
    EOSTVapMin = ceoslp(34)
  END FUNCTION EOSTVapMin
  REAL(kind=8) FUNCTION EOSTVapMax()
    EOSTVapMax = ceoslp(35)
  END FUNCTION EOSTVapMax
  REAL(kind=8) FUNCTION HEvaporation(T)
    IMPLICIT NONE
    REAL(kind=8), INTENT(in):: T
    HEvaporation = hev(T)
  END FUNCTION HEvaporation
  REAL(kind=8) FUNCTION SaturatedPressure(T)
    IMPLICIT NONE
    REAL(kind=8), INTENT(in):: T
    SaturatedPressure = satprs(T)
  END FUNCTION SaturatedPressure


  ! Private functions and subroutines  
  SUBROUTINE thermo(p,el,ev,tl,tv,tsat,rol,rov,pa,rova, &
       tssn,eva)
    !
    !Modified by S.Yakush: corrected calculation of dhlsp, added dhlspa
    !
    !     subroutine thermo evaluates the thermodynamic properties of h2o
    !
    !     Modified Jan. 1989 by S. Woodruff to include SNL extensions needed
    !     for P and T off the table.  *if,def,mel statements added to
    !     provide for alternate formulations for SNL.  *if,-def,mel for
    !     original LANL coding.  These changes should not affect results
    !     when "mel" is NOT defined.
    !     Corrections made to superheated vapor below critical point since
    !     original formulation was thermodynamically unstable.
    !     input variables
    !        1. p      pressure
    !        2. tl     liquid temperature
    !        3. tv     vapor temperature
    !        4. pa     partial pressure of the noncondensible
    !        5. iop    option selector - not in present version
    !
    !     output variables
    !        1. el     liquid internal energy
    !        2. ev     vapor (steam and noncondensable mixture) internal
    !                  energy
    !        3. tsat   saturation temperature corresponding to the total
    !                  pressure
    !        4. rol    liquid density
    !        5. rov    vapor (steam and noncondensable mixture) density
    !        6. rova   density of the noncondensable
    !        7. tssn   saturation temperature corresponding to the steam
    !                  partial pressure
    !        8. eva    internal energy of the noncondensable
    !        9. dtsdp  derivative of tsat wrt pressure: dr(1)
    !       10. deldp  derivative of el wrt pressure: dr(2)
    !       11. devdp  derivative of steam internal energy wrt steam
    !                  partial pressure: dr(3)
    !       12. deldt  derivative of el wrt tl: dr(4)
    !       13. devdt  derivative of steam internal energy wrt tv: dr(5)
    !       14. drolp  derivative of rol wrt pressure: dr(6)
    !       15. drovp  derivative of steam density wrt steam partial
    !                  pressure: dr(7)
    !       16. drolt  derivative of rol wrt tl : dr(8)
    !       17. drovt  derivative of steam density wrt tv: dr(9)
    !       18. hvst   saturated steam enthalpy (psteam,tssn): dr(10)
    !       19. hlst   saturated liquid enthalpy (p,tssn): dr(11)
    !       20. dhvsp  derivative of hvst wrt steam partial pressure: dr(12)
    !       21. dhlsp  derivative of hlst wrt pressure: dr(13)
    !       22. dtssp  derivative of tssn wrt steam partial pressure : dr(14)
    !       23. devat  derivative of eva wrt tv : dr(15)
    !       24. devap  derivative of eva wrt pa : dr(16)
    !       25. drvap  derivative of rova wrt pa: dr(17)
    !       26. drvat  derivative of rova wrt tv: dr(18)
    !
    IMPLICIT NONE

    REAL(kind=8), INTENT(inout):: p(*), pa(*)
    REAL(kind=8), INTENT(inout):: tl(*), tv(*)
    REAL(kind=8), INTENT(out):: el(*), ev(*), eva(*) 
    REAL(kind=8), INTENT(out):: rol(*), rov(*), rova(*)
    REAL(kind=8), INTENT(out):: tsat(*), tssn(*)

    REAL(kind=8):: ale(8),ble(8),cle(8),dle(8)
    REAL(kind=8):: ave(11),bve(11),cve(11),dve(11)
    REAL(kind=8):: avg(11),bvg(11),cvg(11),dvg(11)
    REAL(kind=8):: acp(10),bcp(10),ccp(10),dcp(10)

    REAL(kind=8), PARAMETER:: c26 = 0.30d0
    REAL(kind=8), PARAMETER:: ck0 = -8.329595d-04,ck2=-2.245825d-17,ck4=-1.450382d-06
    REAL(kind=8), PARAMETER:: a11=.10000887519691d-02,a12= .76916250454393d+03
    REAL(kind=8), PARAMETER:: a13=.13001153775598d-02
    !
    !     the following data-initialized variables not used in thermo:
    REAL(kind=8), PARAMETER:: a15 = .15448787270199d-02
    REAL(kind=8), PARAMETER:: delco= 3.7374522042975d+03
    REAL(kind=8), PARAMETER:: delc1= 7.8680714975650d+00
    REAL(kind=8), PARAMETER:: delc2=-4.5114568974853d-02
    REAL(kind=8), PARAMETER:: delc3= 1.2045547889272d-04
    REAL(kind=8), PARAMETER:: deld0=-.81454700000000d+05,deld1= .89319600000000d+03
    REAL(kind=8), PARAMETER:: deld2=-.31234800000000d+01,deld3= .37040880000000d-02
    REAL(kind=8), PARAMETER:: dele0=-.26221567700000d+08,dele1= .22589733400000d+06
    REAL(kind=8), PARAMETER:: dele2=-.64870195500000d+03,dele3= .62113375200000d+00
    !
    !     Pressure limits on transition regime between old and new
    !     formulations for superheated vapor.
    !     PLIM1    Pressure lower bound for transition region.
    !     PLIM2    Pressure upper bound for transition region. (Must be >
    !              lower bound.
    !    real(kind=8) plim1, plim2
    REAL(kind=8), PARAMETER:: plim1 = 190.d+05, plim2 = 200.d+05
    !     data for raw constants used in fits
    !
    DATA ale(1)/-1.1436668993222d+06/,ble(1)/ 4.1868000000000d+03/, &
         &   cle(1)/ 0.d0                 /,dle(1)/ 0.d0                 /
    DATA ale(2)/ 8.0957542810383d+06/,ble(2)/-5.7008855264640d+04/, &
         &   cle(2)/ 1.3443632119671d+02/,dle(2)/-9.7879669155946d-02/
    DATA ale(3)/-1.9373932457007d+06/,ble(3)/ 9.7492797103351d+03/, &
         &   cle(3)/-1.3299615999876d+01/,dle(3)/ 1.0879999999922d-02/
    DATA ale(4)/-5.3245827703670d+06/,ble(4)/ 2.9179372045334d+04/, &
         &   cle(4)/-5.0452192000967d+01/,dle(4)/ 3.4560000000583d-02/
    DATA ale(5)/-6.3583523639930d+07/,ble(5)/ 3.2873715263424d+05/, &
         &   cle(5)/-5.6371182000208d+02/,dle(5)/ 3.2760000000116d-01/
    DATA ale(6)/-6.6239163195929d+09/,ble(6)/ 3.1605562257270d+07/, &
         &   cle(6)/-5.0263730855532d+04/,dle(6)/ 2.6650075114186d+01/
    DATA ale(7)/-5.4759091078157d+09/,ble(7)/ 2.4635618770681d+07/, &
         &   cle(7)/-3.6931079506707d+04/,dle(7)/ 1.8454719393083d+01/
    DATA ale(8)/-7.1536399439453d+07/,ble(8)/ 3.0560801674842d+05/, &
         &   cle(8)/-4.2424553999630d+02/,dle(8)/ 1.9719999999823d-01/
    DATA ave(1)/ 2.4949771766385d+06/,bve(1)/ 2.0855856331827d-01/, &
         &   cve(1)/-1.3553894579716d-07/,dve(1)/ 2.8522684989198d-14/
    DATA ave(2)/ 2.5600870370371d+06/,bve(2)/ 3.1086111111026d-02/, &
         &   cve(2)/-6.8988888888580d-09/,dve(2)/ 4.3203703703379d-16/
    DATA ave(3)/ 2.5915500000006d+06/,bve(3)/ 8.7749999997567d-03/, &
         &   cve(3)/-1.7499999999663d-09/,dve(3)/ 4.2999999998503d-17/
    DATA ave(4)/ 2.6606000000024d+06/,bve(4)/-1.3545000000581d-02/, &
         &   cve(4)/ 6.4250000004682d-10/,dve(4)/-4.2100000001248d-17/
    DATA ave(5)/ 3.8201600000097d+06/,bve(5)/-2.3019900000170d-01/, &
         &   cve(5)/ 1.4068900000098d-08/,dve(5)/-3.1786000000187d-16/
    DATA ave(6)/-1.2103411633350d+08/,bve(6)/ 1.8018803375785d+01/, &
         &   cve(6)/-8.7442426507726d-07/,dve(6)/ 1.4091076856088d-14/
    DATA ave(7)/ 2.2000000000000d+06/,bve(7)/ 0.d0                 /, &
         &   cve(7)/ 0.d0                 /,dve(7)/ 0.d0               /
    DATA ave(8)/ 2.2000000000000d+06/,bve(8)/ 0.d0                 /, &
         &   cve(8)/ 0.d0                 /,dve(8)/ 0.d0               /
    DATA ave(9)/ 2.2000000000000d+06/,bve(9)/ 0.d0                 /, &
         &   cve(9)/ 0.d0                 /,dve(9)/ 0.d0               /
    DATA ave(10)/ 2.2000000000000d+06/,bve(10)/ 0.d0               /, &
         &   cve(10)/ 0.d0                 /,dve(10)/ 0.d0              /
    DATA ave(11)/ 2.2000000000000d+06/,bve(11)/ 0.d0                /, &
         &   cve(11)/ 0.d0                 /,dve(11)/ 0.d0              /
    DATA avg(1)/ 1.0666845123419d+00/,bvg(1)/ 2.8310838172462d-08/, &
         &   cvg(1)/-2.1151097428905d-14/,dvg(1)/ 4.7404001285964d-21/
    DATA avg(2)/ 1.0735412407407d+00/,bvg(2)/ 2.6518055555551d-09/, &
         &   cvg(2)/-6.3461111111128d-16/,dvg(2)/ 3.9824074074117d-23/
    DATA avg(3)/ 1.0777730000000d+00/,bvg(3)/-2.4300000008021d-11/, &
         &   cvg(3)/-7.1979999998378d-17/,dvg(3)/ 4.8799999990422d-25/
    DATA avg(4)/ 1.0851130000007d+00/,bvg(4)/-1.9307000001824d-09/, &
         &   cvg(4)/ 8.9100000014826d-17/,dvg(4)/-3.8960000003946d-24/
    DATA avg(5)/ 1.1639800000015d+00/,bvg(5)/-1.6338350000254d-08/, &
         &   cvg(5)/ 9.5856000001448d-16/,dvg(5)/-2.1194000000274d-23/
    DATA avg(6)/ 3.8898867259868d+00/,bvg(6)/-3.8595945559811d-07/, &
         &   cvg(6)/ 1.7476370114910d-14/,dvg(6)/-2.6377008249858d-22/
    DATA avg(7)/ 2.7168710524682d+00/,bvg(7)/-2.2832718294604d-07/, &
         &   cvg(7)/ 1.0417331983836d-14/,dvg(7)/-1.5842822199773d-22/
    DATA avg(8)/ 3.9749829999964d+00/,bvg(8)/-3.0657099999960d-07/, &
         &   cvg(8)/ 1.0637899999985d-14/,dvg(8)/-1.2257999999981d-22/
    DATA avg(9)/ 1.2946929999997d+00/,bvg(9)/-2.4834999999979d-08/, &
         &   cvg(9)/ 7.8979999999944d-16/,dvg(9)/-8.0799999999948d-24/
    DATA avg(10)/ 1.0590519999963d+00/,bvg(10)/-2.4615999996941d-09/, &
         &   cvg(10)/ 8.8399999991573d-17/,dvg(10)/-8.0799999992269d-25/
    DATA avg(11)/ 1.1430199999838d+00/,bvg(11)/-7.7095999988588d-09/, &
         &   cvg(11)/ 1.9335999997331d-16/,dvg(11)/-1.4639999997924d-24/
    DATA acp(1)/-7.9678485852270d+02/,bcp(1)/ 2.8187658437259d+01/, &
         &   ccp(1)/-1.0180624999920d-01/,dcp(1)/ 1.2499999999912d-04/
    DATA acp(2)/-9.7082632232795d+02/,bcp(2)/ 2.8324981030402d+01/, &
         &   ccp(2)/-9.7656200001157d-02/,dcp(2)/ 1.1600000000110d-04/
    DATA acp(3)/-1.6649701690752d+03/,bcp(3)/ 3.3159363169596d+01/, &
         &   ccp(3)/-1.0861179999898d-01/,dcp(3)/ 1.2399999999915d-04/
    DATA acp(4)/-6.1420486441088d+03/,bcp(4)/ 6.3630987079837d+01/, &
         &   ccp(4)/-1.7762319999965d-01/,dcp(4)/ 1.7599999999975d-04/
    DATA acp(5)/-8.2289951961933d+04/,bcp(5)/ 5.3773958896061d+02/, &
         &   ccp(5)/-1.1612491999609d+00/,dcp(5)/ 8.5599999997375d-04/
    DATA acp(6)/-6.5842104212475d+05/,bcp(6)/3.7934294783212d+03/, &
         &   ccp(6)/-7.2924928000022d+00/,dcp(6)/ 4.7040000000014d-03/
    DATA acp(7)/ 3.4561620732510d+05/,bcp(7)/-2.2129380791446d+02/, &
         &   ccp(7)/-2.4524285999925d+00/,dcp(7)/ 3.1479999999958d-03/
    DATA acp(8)/ 1.9798369474597d+06/,bcp(8)/-1.4782551342826d+04/, &
         &   ccp(8)/ 3.1656481897637d+01/,dcp(8)/-2.0843356864237d-02/
    DATA acp(9)/-9.6249385211359d+07/,bcp(9)/ 4.3633668884423d+05/, &
         &   ccp(9)/-6.5887615106930d+02/,dcp(9)/ 3.3146147264269d-01/
    DATA acp(10)/-1.1074934463333d+07/,bcp(10)/ 4.8073794630970d+04/, &
         &   ccp(10)/-6.9212173247881d+01/,dcp(10)/ 3.3091693999800d-02/
    !
    !
    !     definitions of combinations of constants
    !     (as initialized in data statements above)
    !
    !         a11=2.0*c26/(ceoslp(16)*ceoslp(12))=2.0/ceoslp(23)
    !         a12=1.0/a13=ceoslp(4)/2.0
    !         a13=a11*(1.0+c26)=2.0/ceoslp(4)
    !         aeos14=1.0/c28
    !         c26=ceoslp(16)-1.0
    !
    !    real(kind=8) thermo_dhlsdpa
    !    common/thermo_dhlsdpa_add/ thermo_dhlsdpa

    INTEGER j,jj,ieos,nthm,ii
    INTEGER ldtsdp,ldeldp,ldevdp,ldeldt,ldevdt,ldrolp,&
         & ldrovp,ldrolt,ldrovt,lhvst,lhlst,ldhvsp,ldhlsp,&
         & ldtssp,ldevat,ldevap,ldrvap,ldrvat
    REAL(kind=8):: pt,pg,ts,ps,tv1,tl2,dpes,cps,erp,elt,ert,elp
    REAL(kind=8):: psl, dps, expps
    REAL(kind=8):: t1,t2,t3,t4,t,tb
    REAL(kind=8):: rolst,drlsdp,drlsdt,rrolst
    REAL(kind=8):: drvde,drvde1,rov1,drov1dp,drov1dt
    REAL(kind=8):: rov2,drov2dt,drov2dp,drvde2
    REAL(kind=8):: es,dpcps,elsat,delsat
    REAL(kind=8):: ev1,de,dev1dt,dev1dp
    REAL(kind=8):: ev2,dev2dt,dev2dp
    REAL(kind=8):: gams, gamsm, dpgams
    REAL(kind=8):: beta, f,dfdp, weight,capk,dcapkp,dbetap


    ! FIL: Delete fragment
    !       i0=istrt3
    !       ix=iend3
    !       is=nvthm
    !     if(iop.ne.3) then
    !       i0=1
    !       ix=1
    !       is=1
    !     endif
    !     do 299 i=i0,ix,is
    ! FIL: End of fragment
    !
    !FIL: ADD 2 line
    ieos=0
    nthm=19

    DO 300 jj=jstart,ncells
       j=jj
       !FIL    if(iop.eq.3) then
       !FIL      j=i-i0+(jj-jstart)*ndimv1+1
       !FIL      ldtsdp=j
       !FIL    else
       ldtsdp=nthm*(j-1)+1
       !FIL    endif
       ldeldp=ldtsdp+1
       ldevdp=ldeldp+1
       ldeldt=ldevdp+1
       ldevdt=ldeldt+1
       ldrolp=ldevdt+1
       ldrovp=ldrolp+1
       ldrolt=ldrovp+1
       ldrovt=ldrolt+1
       lhvst=ldrovt+1
       lhlst=lhvst+1
       ldhvsp=lhlst+1
       ldhlsp=ldhvsp+1
       ldtssp=ldhlsp+1
       ldevat=ldtssp+1
       ldevap=ldevat+1
       ldrvap=ldevap+1
       ldrvat=ldrvap+1
       !
       pt = p(j)
       pg = pa(j)
       IF(p(j).LT.ceoslp(30).OR.p(j).GT.ceoslp(31)) THEN
          pt=MIN(ceoslp(31),MAX(ceoslp(30),p(j)))
          CALL error(2,'*thermo* pressure limit exceeded',4)
       ENDIF
       !
       !
       !     calculate saturation properties
       !
       ps=pt-pg
       ps = MAX(ps,ceoslp(30))
       tsat(j)=sattmp(pt)
       dr(ldtsdp)=satder(pt,tsat(j))
       tl2 = tl(j)
       tv1 = tv(j)
       IF(pg.GE.1.d-5) THEN
          tssn(j)=sattmp(ps)
          dr(ldtssp)=satder(ps,tssn(j))
       ELSE
          tssn(j)=tsat(j)
          dr(ldtssp)=dr(ldtsdp)
       ENDIF
       !
       IF (ieos .EQ. 0) THEN
          IF(ev(j).EQ.-1.0d0) THEN
             tv1=tsat(j)
             tv(j)=tv1
          ELSEIF(ev(j).EQ.-2.0d0) THEN
             tv1=tssn(j)
             tv(j)=tv1
          ENDIF
          IF(el(j).EQ.-1.0d0) THEN
             tl2=tsat(j)
             tl(j)=tl2
          ELSEIF(el(j).EQ.-2.0d0) THEN
             tl2=tssn(j)
             tl(j)=tl2
          ENDIF
          IF(tv1.LT.ceoslp(34).OR.tv1.GT.ceoslp(35)) THEN
             tv1=MIN(ceoslp(35),MAX(tv1,ceoslp(34)))
             PRINT *, 'steam temperature', tv(j)
             CALL error(2,'*thermo* vapor temp limit exceeded',4)
          ENDIF
       ENDIF
       !
       IF(tl2.LT.ceoslp(32).OR.tl2.GT.ceoslp(33)) THEN
          tl2=MIN(ceoslp(33),MAX(tl2,ceoslp(32)))
          PRINT *, 'liquid temperature', tl(j)
          CALL error(2,'*thermo* liquid temp limit exceeded',4)
       ENDIF
       !
       !     calculate liquid properties
       !
       !
       !        1. internal energy and its derivatives
       !
       ii=idint(MIN(MAX((tl2-373.15d0)/50.0d0,0.0d0)+1.0d0,7.0d0))
       IF(tl2.GT.645.15d0) ii=ii+1
       el(j)=ale(ii)+tl2*(ble(ii)+tl2*(cle(ii)+tl2*dle(ii)))
       dr(ldeldt)=ble(ii)+tl2*(2.0d0*cle(ii)+tl2*3.0d0*dle(ii))
       psl=satprs(tl2)
       dps=1.0d0/satder(psl,tl2)
       expps=EXP(ck4*psl)
       dr(ldeldp)=ck0*(1.0d0-expps)+ck2*psl*psl
       elp=(pt-psl)*dr(ldeldp)
       ert=dps*(ck0*(-1.0d0+expps*(1.0d0-ck4*(pt-psl))) &
            &      +ck2*(2.0d0*pt*psl-3.0d0*psl*psl))
       el(j)=el(j)+elp
       dr(ldeldt)=dr(ldeldt)+ert
       ii=idint(MIN(MAX((tssn(j)-373.15d0)/50.0d0,0.0d0)+1.0d0,7.0d0))
       IF(tssn(j).GT.645.15d0) ii=ii+1
       elsat=ale(ii)+tssn(j)*(ble(ii)+tssn(j)*(cle(ii)+tssn(j)*dle(ii)))
       delsat=ble(ii)+tssn(j)*(2.0d0*cle(ii)+tssn(j)*3.0d0*dle(ii))
       !
       CALL rholiq(pt,tssn(j),rolst,drlsdp,drlsdt)
       !
       rrolst=1.0d0/rolst
       dr(lhlst)=elsat+pt*rrolst
!!!      dr(ldhlsp)=delsat*dr(ldtsdp)+rrolst
!!!     1           -pt*rrolst*rrolst*(drlsdt*dr(ldtsdp)+drlsdp)
       dr(ldhlsp)=delsat*dr(ldtssp)+rrolst &
            &      -pt*rrolst*rrolst*(drlsdt*dr(ldtssp)+drlsdp)

       thermo_dhlsdpa = -delsat*dr(ldtssp)+ &
            &     pt*rrolst*rrolst*drlsdt*dr(ldtssp)

       !
       !        2. density and its derivatives
       !
       CALL rholiq(pt,tl2,rol(j),dr(ldrolp),dr(ldrolt))
       !
       !
       !
       !     calculate steam properties
       !
       !
       !
       !     properties at saturation
       !
       !     -----specific heat and its derivative
       ii=idint(MIN(MAX((tssn(j)-273.15d0)/50.0d0,0.0d0)+1.0d0,9.0d0))
       IF(tssn(j) > 647.3d0) ii=ii+1
       cps=acp(ii)+tssn(j)*(bcp(ii)+tssn(j)*(ccp(ii)+tssn(j)*dcp(ii)))
       dpcps=(bcp(ii)+tssn(j)*(2.0d0*ccp(ii)+tssn(j)*3.0d0*dcp(ii))) &
            &      *dr(ldtssp)
       !
       !     -----internal energy, enthalpy, and their derivatives
       IF(ps.LE.5.0d+05) THEN
          dr(lhvst)=ceoslp(26)+ceoslp(24)*(tssn(j)-ceoslp(5))+hev(tssn(j))
          dr(ldhvsp)=(ceoslp(24)-2470.2120d0)*dr(ldtssp)
          es=dr(lhvst)-ceoslp(12)*tssn(j)
          dpes=dr(ldhvsp)-ceoslp(12)*dr(ldtssp)
          gams=dr(lhvst)/es
          gamsm=gams-1.0d0
          dpgams=(dr(ldhvsp)-gams*dpes)/es
       ELSE
          ii=idint(MIN(ps/50.0d+05+1.0d0,9.0d0))
          IF(ps.GT.220.0d+05) THEN
             ii=ii+2
          ELSEIF(ps.GT.20.0d+05) THEN
             ii=ii+1
          ENDIF
          es=ave(ii)+ps*(bve(ii)+ps*(cve(ii)+ps*dve(ii)))
          dpes=(bve(ii)+ps*(2.0d0*cve(ii)+ps*3.0d0*dve(ii)))
          gams=avg(ii)+ps*(bvg(ii)+ps*(cvg(ii)+ps*dvg(ii)))
          dpgams=bvg(ii)+ps*(2.0d0*cvg(ii)+ps*3.0d0*dvg(ii))
          gamsm=gams-1.0d0
          dr(lhvst)=gams*es
          dr(ldhvsp)=gams*dpes+es*dpgams
       ENDIF
       !
       !     properties at actual steam temperature
       !
       !
       IF (ieos .EQ. 0) THEN
          !     Saturated and sub-cooled vapor.
          IF(tv1.LE.tssn(j)) THEN
             dr(ldevdt)=cps/ceoslp(16)
             de=dr(ldevdt)*(tv1-tssn(j))
             ev(j)=es+de
             dr(ldevdp)=dpes+de*dpcps/cps-dr(ldevdt)*dr(ldtssp)
             t4=1.0d0/(gamsm*es+c26*de)
             rov(j)=ps*t4
             drvde=-rov(j)*c26*t4
             dr(ldrovt)=drvde*dr(ldevdt)
             dr(ldrovp)=t4-rov(j)*(es*dpgams+(gamsm-c26)*dpes)*t4 &
                  & +drvde*dr(ldevdp)
             !     Superheated vapor.
          ELSE
             t1=1.0d0/(a11*cps-1.0d0)
             beta=tssn(j)*tssn(j)*(1.0d0-t1*t1 )
             t2=tssn(j)*t1
             weight = w (ps)
             !       weight = max(0., min(1., weight))
             f = fare (ps)
             !       write (59, '(''W = '', f10.3, ''F = '', f10.3)') weight, f
             !       If the pressure is in region 1 or in the transition region,
             !       compute the values using the original formulation.
             IF (weight .GE. 0.d0) THEN
                de=a12*(tv1-tssn(j)+SQRT(tv1*tv1-beta)-t2)
                ev1=es+de
                capk=a13*de+tssn(j)+t2
                dbetap=2.0d0*(beta*dr(ldtssp)+a11*dpcps*t2**3)/tssn(j)
                dcapkp=-a13*dpes+(1.0d0+t1)*dr(ldtssp)-a11*t1*t2*dpcps
                t3=1.0d0-beta/(capk*capk)
                dev1dt=ceoslp(4)/t3
                dev1dp=-0.5d0*(t3*dcapkp+dbetap/capk)*dev1dt
                t4=1.0d0/(gamsm*es+c26*de)
                rov1=ps*t4
                drvde1=-rov1*c26*t4
                drov1dt=drvde1*dev1dt
                drov1dp=t4-rov1*(es*dpgams+(gamsm-c26)*dpes)*t4 &
                     &  +drvde1*dev1dp
             ENDIF
             !     If the pressure is in region 2 or in the transition region,
             !     compute the values using the new formulation.
             IF ((1.d0-weight) .GE. 0.d0) THEN
                de = (1.d0/(a11*gams))*(tv1 - &
                     &  tssn(j) + &
                     &  SQRT (tv1**2-beta) - &
                     &  t2)
                ev2 = es + de
                dev2dt = (1.d0/(a11*gams))*(1.d0 + &
                     &  tv1/SQRT (tv1**2-beta))
                dbetap=2.0d0*(beta*dr(ldtssp)+a11*dpcps*t2**3)/tssn(j)
                dev2dp = dpes - &
                     &           (de/gams)*dpgams + &
                     &           (1.d0/(a11*gams))*(-dr(ldtssp) - &
                     &                            0.5d0*dbetap/SQRT (tv1**2-beta) - &
                     &                            t1*dr(ldtssp) + t1*t2*a11*dpcps)
                rov2   = ps/(gamsm*ev2)
                drvde2 = -(rov2/(gamsm*ev2))*gamsm
                drov2dt = drvde2*dev2dt
                drov2dp = (1.d0/(gamsm*ev2))*(1.d0 - rov2*ev2*dpgams) + &
                     &               drvde2*dev2dp
             ENDIF
             !
             !
             !       Store final values for superheated vapor properties.
             !       If the pressure is in region 1, use the original formulation.
             IF (weight .GE. 1.d0) THEN
                ev(j)       = ev1
                dr(ldevdp)  = dev1dp
                dr(ldevdt)  = dev1dt
                rov(j)      = rov1
                dr(ldrovp) = drov1dp
                dr(ldrovt) = drov1dt
                !       If the pressure is in the transition region, fare between
                !       the two formulations.
             ELSEIF (weight .GT. 0.d0)THEN
                dfdp = dfaredp (ps)
                ev(j)       = f*ev1     + (1.d0-f)*ev2
                dr(ldevdp)  = f*dev1dp +(1.d0-f)*dev2dp  + dfdp*(ev1 - ev2)
                dr(ldevdt)  = f*dev1dt  + (1.d0-f)*dev2dt
                rov(j)      = f*rov1    + (1.d0-f)*rov2
                dr(ldrovp) = f*drov1dp+(1.d0-f)*drov2dp + dfdp*(rov1 - rov2)
                dr(ldrovt) = f*drov1dt + (1.d0-f)*drov2dt
                !       If the pressure is in region 2, use the new formulation.
             ELSE
                ev(j)       = ev2
                dr(ldevdp)  = dev2dp
                dr(ldevdt)  = dev2dt
                rov(j)      = rov2
                dr(ldrovp) = drov2dp
                dr(ldrovt) = drov2dt
             ENDIF
          ENDIF
          IF(ps.GT.0.0d0.AND.rov(j).LE.0.0d0) THEN
             rov(j)=ps/(ceoslp(12)*tv1)
             dr(ldrovt)=-rov(j)/tv1
             dr(ldrovp)=rov(j)/ps
          ELSEIF(rov(j).GE.(0.999d0*rol(j))) THEN
             rov(j)=0.999d0*rol(j)
             dr(ldrovt)=0.999d0*dr(ldrolt)
             dr(ldrovp)=0.999d0*dr(ldrolp)
          ENDIF
       ENDIF
       !
       !
       !
       !
       !     calculate air properties
       !
       !
       !
       eva(j)=ceoslp(17)*tv1
       dr(ldevat)=ceoslp(17)
       dr(ldevap)=0.0d0
       dr(ldrvap)=1.0d0/(ceoslp(25)*tv1)
       rova(j)=dr(ldrvap)*pa(j)
       dr(ldrvat)=-ceoslp(25)*rova(j)*dr(ldrvap)
       !
       !
       !
       !     If "steam" is to be replaced by a non-condensible,
       !     copy non-condensible quantities into steam arrays.
       IF (ieos .NE. 0) THEN
          ev(j)=eva(j)
          rov(j)=dr(ldrvap)*ps
          dr(ldevdt)=dr(ldevat)
          dr(ldevdp)=dr(ldevap)
          dr(ldrovp)=dr(ldrvap)
          dr(ldrovt)=-ceoslp(25)*rov(j)*dr(ldrvap)
       ENDIF
       !     calculate air-steam mixture properties
       !
       !
       !
       ev(j)=(ev(j)*rov(j)+eva(j)*rova(j))
       rov(j)=rov(j)+rova(j)
       ev(j)=ev(j)/rov(j)
       p(j)=pt
       tl(j)=tl2
       tv(j)=tv1
300    CONTINUE
       ! FIL  299 continue
       RETURN

     CONTAINS
       !     ------------------------statement functions-----------------------
       !     Derivative of weighting function wrt p.
       REAL(kind=8) FUNCTION dwdp(p1)
         REAL(kind=8), INTENT(in):: p1
         dwdp = -1.d0/(plim2 - plim1)
       END FUNCTION dwdp
       !     weighting function.
       REAL(kind=8) FUNCTION w(p1)
         REAL(kind=8), INTENT(in):: p1
         w  = (plim2 - p1)*(-dwdp (p1))
       END FUNCTION w
       !     faring function (cubic).
       REAL(kind=8) FUNCTION fare (p1)
         REAL(kind=8), INTENT(in):: p1
         fare = 3.d0*(w (p1))**2 - 2.d0*(w (p1))**3
       END FUNCTION fare
       !     derivative of faring function.
       REAL(kind=8) FUNCTION dfaredp (p1)
         REAL(kind=8), INTENT(in):: p1
         dfaredp = 6.d0*w (p1)*dwdp (p1)*(1.d0 - w (p1))
       END FUNCTION dfaredp
     END SUBROUTINE thermo


     REAL(kind=8) FUNCTION cpll(h,p)
       !
       !     function cpll evaluates the specific heat of h2o liquid
       !     as a function of liquid enthalpy and total pressure
       !                (when b1, c1, and d1 .ne. 0.0)
       !
       !          liquid enthalpy             h      in (j/kg)
       !          total pressure              p      in (pa)
       !          liquid specific heat        cpll   in (j/kg/k)
       !
       IMPLICIT NONE
       SAVE

       REAL(kind=8),INTENT(in):: h, p

       REAL(kind=8),PARAMETER:: b0 = 2.394907d-04,b1 = -5.196250d-13
       REAL(kind=8),PARAMETER:: c0 = 1.193203d-11,c1 = 2.412704d-18
       REAL(kind=8),PARAMETER:: d0 = -3.944067d-17,d1=-1.680771d-24
       !
       cpll=1.0d0/MAX(2.5d-05,(b0+p*b1)+h*((c0+p*c1)+h*(d0+p*d1)))
     END FUNCTION cpll

     REAL(kind=8) FUNCTION cpvv1(t,p,pa)
       !
       !     function cpvv1 evaluates the specific heat of h2o vapor
       !     as a function of vapor temperature, total pressure, and
       !     noncondensable-gas pressure
       !
       !          vapor temperature           t      in (k)
       !          total pressure              p      in (pa)
       !          noncondensable gas pressure pa     in (pa)
       !          vapor specific heat         cpvv1  in (j/kg/k)
       !
       IMPLICIT NONE
       REAL(kind=8), INTENT(in):: t,p,pa
       SAVE

       REAL(kind=8),PARAMETER:: c1 = 1688.35968d0,c2 = 0.6029856d0, &
            & c3 = 482.0979623d0,c4 = 2.95317905d+7,c5 = 1.8d0,c6 = 460.d0
       REAL(kind=8) tb

       !      include 'function.i'
       !
       !     miscellaneous inline functions
       !
       !      stard(qanty1,expnt1,qanty2,expnt2)
       !     1     =exp(expnt1*dlog(qanty1)+expnt2*dlog(qanty2))
       !      starsr(qanty)=sqrt(sqrt(qanty))
       !


       !FIL  if(ieos.ne.0) go to 100
       !
       tb=t*c5-c6
       cpvv1=c1+c2*t+c3*p/star(tb,2.4d0)+c4*p*p*p/tb**9
       IF(pa.LT.1.0d-5) go to 110
       !
       !pvv1=((p-pa)*cpvv1+pa*ceoslp(22))/p
       !FIL      go to 110
       !
       !FIL  100 continue
       !FIL      cpvv1=ceoslp(22)
       !
110    CONTINUE
       RETURN
     END FUNCTION cpvv1

     SUBROUTINE error(code,mess,idum)
       INTEGER code
       CHARACTER *(*) mess

       PRINT *,mess
       !      if(code .ne. 2) stop 'ERROR'
       IF(code .EQ. 2) STOP 'ERROR'
       RETURN
     END SUBROUTINE error

     REAL(kind=8) FUNCTION hev(temp)
       !
       !     function hev calculates the heat of evaporation of h2o liquid
       !     as a function of liquid temperature for low pressures
       !
       !          liquid temperature          temp   in (k)
       !          heat of evaporation         hev    in (j/kg)
       !
       IMPLICIT NONE
       REAL(kind=8), INTENT(in):: temp
       SAVE

       REAL(kind=8),PARAMETER:: a = 3.18061959d+06,b = -2.4702120d+03
       REAL(kind=8):: t
       !
       t=MAX(temp,ceoslp(32))
       hev=a+b*t
       RETURN
     END FUNCTION hev


     REAL(kind=8) FUNCTION satder(pres,temp)
       !
       !     function satder evaluates the derivative of the h2o saturation
       !     temperature with respect to total pressure as a function
       !     of the saturation pressure and the saturation temperature
       !
       !          saturation pressure         pres   in (pa)
       !          saturation temperature      temp   in (k)
       !          dtsat/dp                    satder in (k/pa)
       !
       IMPLICIT REAL(kind=8) (a-h,o-z)
       SAVE
       !      include 'tsatcn.i'
       IF(temp.GE.ceoslp(21)) go to 10
       p=MAX(pres,ceoslp(36))
       t=MAX(temp,ceoslp(32))
       satder=ceoslp(12)*t*t/(p*hev(t))
       go to 40
       !
10     CONTINUE
       IF(temp.GE.ceoslp(40)) go to 20
       satder=ceos2*(temp-ceos3)/pres
       go to 40
       !
20     CONTINUE
       IF(temp.GE.ceoslp(38)) go to 30
       satder=-temp**2/(pres*(-8529.6481905883d0+2333338.6556656d0/temp))
       go to 40
       !
30     CONTINUE
       satder=2.0304886238506d-04*temp**2/pres
       !
40     CONTINUE
       RETURN
     END FUNCTION satder

     REAL(kind=8) FUNCTION satprs(temp)
       !
       !     function satprs evaluates the h2o saturation pressure
       !     as a function of the saturation temperature
       !
       !          saturation temperature      temp   in (k)
       !          saturation pressure         satprs in (pa)
       !
       IMPLICIT REAL(kind=8) (a-h,o-z)
       SAVE
       !      include 'tsatcn.i'
       !      include 'function.i'
       !      star(qanty,expnt)=exp(expnt*dlog(qanty))



       IF(temp.GE.ceoslp(21)) go to 10
       t=MAX(temp,ceoslp(32))
       satprs=24821.0d0*star(t/338.0d0,-5.3512d0)*EXP(20.387d0*&
            &(t-338.0d0)/t)
       satprs=MAX(satprs,ceoslp(36))
       go to 40
       !
10     CONTINUE
       IF(temp.GE.ceoslp(40)) go to 20
       satprs=(1.0d0/aeos14)*star((temp-ceos3)/ceos1,1.0d0/ceos2)
       go to 40
       !
20     CONTINUE
       IF(temp.GE.ceoslp(38)) go to 30
       satprs=7.2166948490268d+11 &
            &     *EXP((-8529.6481905883d0+1166669.3278328d0/temp)/temp)
       go to 40
       !
30     CONTINUE
       satprs=ceoslp(37) &
            &     *EXP(7.6084086799277d0-4924.9229385171d0/temp)
       !
40     CONTINUE
       RETURN
     END FUNCTION satprs

     REAL(kind=8) FUNCTION sattmp(pres)
       !
       !     function sattmp evaluates the h2o saturation temperature
       !     as a function of the saturation pressure
       !
       !          saturation pressure         pres   in (pa)
       !          saturation temperature      sattmp in (k)
       !
       IMPLICIT REAL(kind=8) (a-h,o-z)
       SAVE
       !      include 'tsatcn.i'
       !      include 'function.i'
       !      star(qanty,expnt)=exp(expnt*dlog(qanty))


       IF(pres.GE.ceoslp(20)) go to 10
       p=MAX(pres,ceoslp(36))
       sattmp=ceoslp(1)/(ceoslp(2)*LOG(p/ceoslp(11))+ceoslp(3))
       DO i=1,2
          hfgref=hev(sattmp)
          psref=satprs(sattmp)
          sattmp=sattmp/(1.0d0-ceoslp(12)*sattmp*LOG(p/psref)/hfgref)

       ENDDO
       go to 40
       !
10     CONTINUE
       IF(pres.GE.ceoslp(39)) go to 20
       sattmp=ceos1*star(aeos14*pres,ceos2)+ceos3
       go to 40
       !
20     CONTINUE
       IF(pres.GE.ceoslp(37)) go to 30
       plog=LOG(pres)
       sattmp=(4264.8240952941d0 &
            &        +SQRT(-13666986.708428d0+1166669.3278328d0*plog)) &
            &        /(27.304833093884d0-plog)
       go to 40
       !
30     CONTINUE
       sattmp=4924.9229385171d0/(24.520401414546d0-LOG(pres))
       !
40     CONTINUE
       RETURN
     END FUNCTION sattmp


     SUBROUTINE seteos
       !
       !     subroutine seteos initializes the h2o equation-of-state constants
       !
       !     airmol  = molecular weight of air
       !     ceoslp( 1) = 1st coeff of sat vap temp function
       !     ceoslp( 2) = 2nd coeff of sat vap temp function
       !     ceoslp( 3) = 3rd coeff of sat vap temp function
       !     ceoslp( 4) = vapor specific heat at constant volume
       !     ceoslp( 5) = reference temperature
       !     ceoslp( 6) = vapor specific internal energy at ref. temp.
       !     ceoslp( 7) = liquid specific heat at constant volume
       !     ceoslp( 8) = liquid specific internal energy at ref. temp.
       !     ceoslp( 9) = microscopic density of liquid
       !     ceoslp(10) = heat of evaporation at reference temperature
       !     ceoslp(11) = reference pressure
       !     ceoslp(12) = gas constant for vapor
       !     ceoslp(13) = noncondensable gas thermal conductivity coefficient
       !     ceoslp(14) = liquid thermal conductivity
       !     ceoslp(15) = square of reciprocal sound speed for liquid
       !     ceoslp(16) = vapor gamma, ratio of specific heats
       !     ceoslp(17) = air specific heat at constant volume
       !     ceoslp(18) = air gamma, ratio of specific heats
       !     ceoslp(19) = noncondensable gas thermal conductivity exponent
       !     ceoslp(20) = upper limit on pressure for low-pressure properties
       !     ceoslp(21) = upper limit on temp. for low-pressure properties
       !     ceoslp(22) = air specific heat at constant pressure
       !     ceoslp(23) = vapor specific heat at constant pressure
       !     ceoslp(24) = liquid specific heat at constant pressure
       !     ceoslp(25) = gas constant for air
       !     ceoslp(26) = liquid specific enthalpy at reference temperature
       !     ceoslp(27) = vapor specific enthalpy at reference temperature
       !     ceoslp(28) = (vap. gas const. - air gas const.) / vap. gas const.
       !     ceoslp(29) = datum temperature for the enthalpy of liquid
       !     ceoslp(30) = minimum allowable pressure
       !     ceoslp(31) = maximum allowable pressure
       !     ceoslp(32) = minimum allowable liquid temperature
       !     ceoslp(33) = maximum allowable liquid temperature
       !     ceoslp(34) = minimum allowable vapor temperature
       !     ceoslp(35) = maximum allowable vapor temperature
       !     ceoslp(36) = minimum saturation pressure
       !                  (corresponding to minimum liquid temperature)
       !     ceoslp(37) = critical pressure
       !     ceoslp(38) = critical temperature
       !     ceoslp(39) = lower limit on saturation pressure for the
       !                  high-pressure expression
       !     ceoslp(40) = lower limit on saturation temperature for the
       !                  high-pressure expression
       !     gascon     = universal gas constant
       !     igas       = 1 for air, 2 for hydrogen, 3 for helium
       !     vapmol     = molecular weight of vapor
       !
       IMPLICIT NONE

       REAL(kind=8):: vapmol

       aeos14=1.0d-05
       ceos1=117.8d0
       ceos2=0.223d0
       ceos3=255.2d0
       ceoslp(1)=-2263.0d0
       ceoslp(2)=0.434d0
       ceoslp(3)=-6.064d0
       !
       SELECT CASE(igas)
       CASE(IGAS_AIR) ! noncondensable gas is air
          airmol=28.96461d0
          !           = mole fraction * molecular weight sum
          !           =       0.78084 * 28.01352  for n2
          !                 + 0.20946 * 31.99874  for o2
          !                 + 0.00934 * 39.94766  for ar
          !                 + 0.00033 * 44.00989  for co2
          !                 + 0.00002 * 20.17135  for ne
          !                 + 0.00001 *  4.00260  for he
          ceoslp(13)=2.091d-04
          ceoslp(19)=0.846d0
          ceoslp(22)=1004.832d0
       CASE(IGAS_HYDROGEN) ! hydrogen
          airmol=2.01594d0
          ceoslp(13)=1.6355d-03
          ceoslp(19)=0.8213d0
          ceoslp(22)=14533.2d0
       CASE(IGAS_HELIUM) ! helium
          airmol=4.00260d0
          ceoslp(13)=3.366d-03
          ceoslp(19)=0.668d0
          ceoslp(22)=5234.0d0
       CASE(IGAS_ARGON) ! mixture of argon and hydrogen 
          ! here properties of pure argon are specified,
          ! mixture properties are defined by subsequent call to Mixture_Gas_H2()
          airmol=40.01594d0 ! Molecular mass
          ceoslp(13)=2.091d-04 ! Check values !!!
          ceoslp(19)=0.8d0
          ceoslp(22)=525.d0
       END SELECT

       !      gascon=8.314339d+03
       !           =6.022169e+26*1.380622e-23
       ceoslp(25)=gascon/airmol
       ceoslp(17)=ceoslp(22)-ceoslp(25)
       ceoslp(18)=ceoslp(22)/ceoslp(17)
       !    ceoslp(28)=(ceoslp(12)-ceoslp(25))/ceoslp(12)
       !
       ceoslp(5)=273.15d0
       ceoslp(29)=273.16d0
       ceoslp(30)=1.0d0
       ceoslp(31)=450.0d+5
       ceoslp(32)=ceoslp(5)
       ceoslp(33)=713.94025779311d0
       ceoslp(34)=ceoslp(5)
       ceoslp(35)=3000.0d0
       ceoslp(36)=610.8d0
       ceoslp(37)=221.2d+5
       ceoslp(38)=647.3d0
       ceoslp(39)=139.69971285053d+5
       ceoslp(40)=609.62462615967d0
       vapmol=18.016d0
       ceoslp(12)=gascon/vapmol
       ceoslp(16)=1.3d0
       ceoslp(4)=ceoslp(12)/(ceoslp(16)-1.0d0)
       ceoslp(23)=ceoslp(16)*ceoslp(4)
       !    ceoslp(28)=(ceoslp(12)-ceoslp(25))/ceoslp(12)
       ceoslp(24)=4186.800d0
       ceoslp(14)=0.65141d0
       ceoslp(7)=ceoslp(24)
       ceoslp(9)=990.0d0
       ceoslp(15)=0.0d0
       ceoslp(20)=9.056466d+4
       ceoslp(21)=370.4251d0
       ceoslp(11)=100000.0d0
       ceoslp(10)=hev(ceoslp(5))
       ceoslp(8)=-611.2d0*0.0010002d0+ceoslp(7)*(ceoslp(5)-ceoslp(29))
       ceoslp(26)=ceoslp(24)*(ceoslp(5)-ceoslp(29))
       ceoslp(27)=ceoslp(26)+ceoslp(10)
       ceoslp(6)=ceoslp(27)-ceoslp(12)*ceoslp(5)
       RETURN
     END SUBROUTINE seteos

     REAL(kind=8) FUNCTION viscl(h,p)
       !
       !     function viscl evaluates the h2o liquid dynamic viscosity
       !     as a function of liquid enthalpy and pressure
       !
       !          liquid enthalpy             h      in (j/kg)
       !          pressure                    p      in (pa)
       !          liquid viscosity            viscl  in (pa*s)
       !
       IMPLICIT REAL(kind=8) (a-h,o-z)
       SAVE
       DIMENSION a(5), b(4), d(5), e(4), f(4)
       DATA a/1.298102340d-3,-9.264032108d-4, 3.81047061d-4, &
            &     -8.219444458d-5, 7.022437984d-6/, ho/8.581289699d-6/, &
            &     con/4.265884d4/, econ/5.53588d4/, eho/6.484503981d-6/,&
            &     b/-6.5959d-12, 6.763d-12, -2.88825d-12, 4.4525d-13/,&
            &     d/3.026032306d-4, -1.836606896d-4, 7.567075775d-5,&
            &     -1.647878879d-5, 1.416457633d-6/, hoo/3.892077365d-6/,&
            &     e/1.4526052612d-3, -6.9880084985d-9, 1.5210230334d-14,&
            &     -1.2303194946d-20/, h1/.276d6/, h2/.394d6/, cn/4.014676d5/,&
            &     f/-3.8063507533d-11, 3.9285207677d-16, -1.2585799292d-21,&
            &     1.2860180788d-27/, pi/6.894575293d5/
       !
       poly3(a1,a2,a3,a4,x)=a1+x*(a2+x*(a3+x*a4))
       poly4(a1,a2,a3,a4,a5,x)=a1+x*(a2+x*(a3+x*(a4+x*a5)))
       xi      = (h-con)*ho
       eta     = (h-econ)*eho
       IF(h .GT. h1) go to 40
       !
       viscl=poly4(a(1),a(2),a(3),a(4),a(5),xi) &
            &     -poly3(b(1),b(2),b(3),b(4),eta)*(p-pi)
       go to 60
       !
40     CONTINUE
       IF(h .GT. h2) go to 50
       !
       viscl=poly3(e(1),e(2),e(3),e(4),h) &
            &     +poly3(f(1),f(2),f(3),f(4),h)*(p-pi)
       go to 60
       !
50     CONTINUE
       viscl=poly4(d(1),d(2),d(3),d(4),d(5),(h-cn)*hoo)
       !
60     CONTINUE
       RETURN
     END FUNCTION viscl

     SUBROUTINE rholiq(p,tl,rhol,drldp,drldt)
       !
       !     subroutine rholiq evaluates the density of h2o liquid and its
       !     derivatives with respect to total pressure and liquid temper-
       !     ature as a function of total pressure and liquid temperature
       !
       !          total pressure              p      in (pa)
       !          liquid temperature          tl     in (k)
       !          liquid density              rol    in (kg/m**3)
       !          drol/dp                     drldp  in (kg/m**3/pa)
       !          drol/dt                     drldt  in (kg/m**3/k)
       !
       IMPLICIT REAL(kind=8) (a-h,o-z)
       SAVE
       DIMENSION av0(12),bv0(12),cv0(12),dv0(12), &
            &     afn(12),bfn(12),cfn(12),dfn(12)
       DATA av0(1)/ 1.7057666777468d-03/,bv0(1)/-6.0320895569365d-06/,&
            &     cv0(1)/ 1.5944423965594d-08/,dv0(1)/-1.2149418561177d-11/
       DATA av0(2)/ 5.2145931517155d-04/,bv0(2)/ 3.5189228252915d-06/,&
            &     cv0(2)/-9.7304881862624d-09/,dv0(2)/ 1.0856688130631d-11/
       DATA av0(3)/-1.4931865836934d-02/,bv0(3)/ 9.7931556400429d-05/,&
            &     cv0(3)/-2.0172817692512d-07/,dv0(3)/ 1.4080475270259d-10/
       DATA av0(4)/-4.9334201381918d-01/,bv0(4)/ 2.5928571576499d-03/,&
            &     cv0(4)/-4.5387107397840d-06/,dv0(4)/ 2.6537936475365d-09/
       DATA av0(5)/-3.4558955902321d+00/,bv0(5)/ 1.7351793841884d-02/,&
            &     cv0(5)/-2.9047483637289d-05/,dv0(5)/ 1.6220227777320d-08/
       DATA av0(6)/-1.1952528427292d+01/,bv0(6)/ 5.8904962031842d-02/,&
            &     cv0(6)/-9.6786687447220d-05/,dv0(6)/ 5.3029284583415d-08/
       DATA av0(7)/-3.7446629978341d+01/,bv0(7)/ 1.8173474403006d-01/,&
            &     cv0(7)/-2.9404991620713d-04/,dv0(7)/ 1.5863005350824d-07/
       DATA av0(8)/-3.9713284923576d+02/,bv0(8)/ 1.8801824705202d+00/,&
            &     cv0(8)/-2.9673900150051d-03/,dv0(8)/ 1.5612171739106d-06/
       DATA av0(9)/-2.3142714272157d+03/,bv0(9)/ 1.0710216457395d+01/,&
            &     cv0(9)/-1.6521763202064d-02/,dv0(9)/ 8.4955209566212d-06/
       DATA av0(10)/ 2.0481569977849d+03/,bv0(10)/-9.3452783115489d+00/,&
            &     cv0(10)/ 1.4212077056589d-02/,dv0(10)/-7.2037202704367d-06/
       DATA av0(11)/-7.3864713248117d+01/,bv0(11)/ 3.3144939132191d-01/,&
            &     cv0(11)/-4.9608715522591d-04/,dv0(11)/ 2.4771793009809d-07/
       DATA av0(12)/-2.1891320674084d+01/,bv0(12)/ 9.6758467414310d-02/,&
            &     cv0(12)/-1.4289074953436d-04/,dv0(12)/ 7.0567217785700d-08/
       DATA afn(1)/-4.2486354144244d+09/,bfn(1)/ 3.7516769853867d+07/,&
            &     cfn(1)/-1.0064945851796d+05/,dfn(1)/ 8.7507285129715d+01/
       DATA afn(2)/-2.7936308563236d+08/,bfn(2)/ 5.5663179995300d+06/,&
            &     cfn(2)/-1.4921749894688d+04/,dfn(2)/ 1.0834095198280d+01/
       DATA afn(3)/-1.1761210016041d+08/,bfn(3)/ 4.3832221802974d+06/,&
            &     cfn(3)/-1.2088373365747d+04/,dfn(3)/ 8.6034520917150d+00/
       DATA afn(4)/-4.5415129389018d+09/,bfn(4)/ 2.7368608704680d+07/,&
            &     cfn(4)/-5.1894794477625d+04/,dfn(4)/ 3.1581281016141d+01/
       DATA afn(5)/-4.0104325667716d+10/,bfn(5)/ 2.0292575433752d+08/,&
            &     cfn(5)/-3.4075971373732d+05/,dfn(5)/ 1.9000660267975d+02/
       DATA afn(6)/-6.0173879922257d+10/,bfn(6)/ 2.9984925450490d+08/,&
            &     cfn(6)/-4.9675963282729d+05/,dfn(6)/ 2.7368658401451d+02/
       DATA afn(7)/ 2.0678826351719d+10/,bfn(7)/-8.9503807129603d+07/,&
            &     cfn(7)/ 1.2822787819385d+05/,dfn(7)/-6.0722291833340d+01/
       DATA afn(8)/ 8.3793557728900d+10/,bfn(8)/-3.8997180562867d+08/,&
            &     cfn(8)/ 6.0502628698976d+05/,dfn(8)/-3.1291965911464d+02/
       DATA afn(9)/ 9.2402374347985d+10/,bfn(9)/-4.2674923965292d+08/,&
            &     cfn(9)/ 6.5695613829284d+05/,dfn(9)/-3.3711122197289d+02/
       DATA afn(10)/-2.7547713637194d+10/,bfn(10)/ 1.2580004134443d+08/,&
            &     cfn(10)/-1.9147491048695d+05/,dfn(10)/ 9.7136148925404d+01/
       DATA afn(11)/ 6.8608195287374d+08/,bfn(11)/-3.0636028439513d+06/,&
            &     cfn(11)/ 4.5613625244005d+03/,dfn(11)/-2.2642074876391d+00/
       DATA afn(12)/ 4.3458430609231d+07/,bfn(12)/-1.8379937116289d+05/,&
            &     cfn(12)/ 2.5971646178490d+02/,dfn(12)/-1.2244044950391d-01/
       DATA an/7.146d0/
       !
       !     noac is the no-artificial-compressibility flag
       !          = 0 , no artificial compressibility is off
       !          = 1 , no artificial compressibility is on
       !
       DATA noac/1/
       !
       !     density and its derivatives
       !
       ii=idint(MIN((tl-273.15d0)*0.01d0+1.0d0,5.0d0)) &
            &     +idint(MIN(MAX((tl-593.15d0)*0.1d0,0.0d0),7.0d0))
       voll0=av0(ii)+tl*(bv0(ii)+tl*(cv0(ii)+tl*dv0(ii)))
       dvoll0=bv0(ii)+tl*(2.0d0*cv0(ii)+tl*3.0d0*dv0(ii))
       funct=afn(ii)+tl*(bfn(ii)+tl*(cfn(ii)+tl*dfn(ii)))
       dfunct=bfn(ii)+tl*(2.0d0*cfn(ii)+tl*3.0d0*dfn(ii))
       t1=1.0d0+p/funct
       rhol=1.0d0/(voll0*(1.0d0-LOG(t1)/an))
       drldp=rhol**2*voll0/((p+funct)*an)
       drldt=-rhol*dvoll0/voll0+(1.0d0-t1)*drldp*dfunct
       !
       !     artificial compressibility is applied when noac=0
       !
       IF(noac.EQ.1) go to 120
       !
       IF(p.GE.4.0d+05) go to 100
       !
       acf=5.0d-03-6.25d-09*p
       f=1.0d0-acf
       drldp=f*drldp+rhol*6.25d-09
       go to 110
       !
100    CONTINUE
       acf=1.0d+03/p
       f=1.0d0-acf
       drldp=f*drldp+rhol*1.0d-03*acf*acf
       !
110    CONTINUE
       drldt=f*drldt
       rhol=f*rhol
       !
120    CONTINUE
       RETURN
     END SUBROUTINE rholiq

     REAL(kind=8) FUNCTION viscv(h,p,rov,tv,pa)
       !
       !     function viscv evaluates the h2o vapor dynamic viscosity
       !     as a function of vapor temperature, vapor density, total
       !     pressure and noncondensable-gas pressure
       !
       !          vapor enthalpy              h      in (j/kg)
       !          total pressure              p      in (pa)
       !          vapor density               rov    in (kg/m**3)
       !          vapor temperature           tv     in (k)
       !          noncondensable gas pressure pa     in (pa)
       !          vapor viscosity             viscv  in (pa*s)
       !
       IMPLICIT REAL(kind=8) (a-h,o-z)
       SAVE
       ! $include contrllr.i


       !      include 'tsatcn.i'
       DATA vap1,vap2,vap3/-2.885d-06,2.427d-08,-6.789333333d-11/
       DATA vap4,vap5,vap6/6.317037d-14,1.76d+02,-1.6d0/
       DATA vap7,vap8,vap9/4.8d-03,-4.74074074d-06,3.53d-08/
       DATA vap10,vap11/6.765d-11,1.021d-14/
       DATA air1,air2,air3/1.707623d-05,5.927d-08,-8.14d-11/
       DATA air4,air5,air6/1.735d-05,4.193d-08,-1.09d-11/
       DATA   h1,  h2,  h3/4.175d-06,1.588d-08,7.6705d-13/
       DATA he1,he2,he3/5.9642d-06,5.2047d-08,-1.5345d-11/
       !
       !     requires 2nd-order and 3rd-order polynomials
       poly2(a1,a2,a3,x)=a1+x*(a2+x*a3)
       poly3(a1,a2,a3,a4,x)=a1+x*(a2+x*(a3+x*a4))
       !
       !FIL  if(ieos.ne.0) go to 130
       !
       !     viscosity of h2o vapor
       !
       t=tv-273.15d0
       v1=8.04d-06+4.07d-08*t
       IF(t.LE.300.0d0) go to 100
       IF(t.GE.375.0d0) go to 110
       viscv=v1+(poly3(vap1,vap2,vap3,vap4,t)+ &
            &     poly3(vap5,vap6,vap7,vap8,t)* &
            &     poly2(vap9,vap10,vap11,rov))*rov
       go to 120
       !
100    CONTINUE
       viscv=v1+(5.9d-10*t-1.858d-07)*rov
       IF(viscv.LT.1.0d-07) viscv=1.0d-07
       go to 120
       !
110    CONTINUE
       viscv=v1+poly2(vap9,vap10,vap11,rov)*rov
       !
120    CONTINUE
       IF(pa.LT.0.1d0) go to 180
       vsave=viscv
       !
130    CONTINUE
       IF(igas.EQ.3) go to 160
       IF(igas.EQ.2) go to 150
       !
       !     viscosity of noncondensable air
       !
       t=tv-273.15d0
       IF(t.GT.229.0d0) go to 140
       viscv=poly2(air1,air2,air3,t)
       go to 170
       !
140    CONTINUE
       viscv=poly2(air4,air5,air6,t)
       go to 170
       !
       !     viscosity of noncondensable hydrogen
       !       reference: handbook of thermodynamic tables & charts
       !                  kuzman, raznjevic, hemisphere publishing corp.
       !
150    CONTINUE
       viscv=poly2(h1,h2,h3,tv)
       go to 170
       !
       !     viscosity of noncondensable helium
       !
160    CONTINUE
       viscv=poly2(he1,he2,he3,tv)
       !
170    CONTINUE
       !FIL      if(ieos.ne.0) go to 180
       !
       f=MIN(1.0d0,pa/p)
       viscv=(1.0d0-f)*vsave+f*viscv
       !
180    CONTINUE
       RETURN
     END FUNCTION viscv

     REAL(kind=8) FUNCTION h_Ar_H2(T,Y_H2,Y_H2_restore)
       REAL(kind=8), INTENT(in):: T
       REAL(kind=8), INTENT(in):: Y_H2
       REAL(kind=8), INTENT(in),OPTIONAL:: Y_H2_restore

       CALL Mixture_Gas_H2(Y_H2)
       h_Ar_H2 = ceoslp(22)*T

       IF(PRESENT(Y_H2_restore)) THEN
          CALL Mixture_Gas_H2(Y_H2_restore)
       ENDIF

     END FUNCTION h_Ar_H2

     !------------------------------------Ar-H2 mixture -----------------------
     SUBROUTINE Mixture_Gas_H2(Y_H2,X_H2,W_a)
       !     
       !     subroutine Mixture_Gas_H2 calculates the properties of Ar-H2 mixture
       !     for the given mass fraction of hydrogen Y_H2
       !     it also returns molar (volume) fraction of hydrogen X_H2
       !
       IMPLICIT NONE

       REAL(kind=8),INTENT(in) :: Y_H2 ! Mass fraction of hydrogen
       REAL(kind=8),OPTIONAL,INTENT(out):: X_H2 ! Molar fraction of hydrogen
       REAL(kind=8),OPTIONAL,INTENT(out):: W_a  ! Moleculat mass of Ar-H2 mixture

       REAL(kind=8):: Y_A, X_A, Cp_A, airmol_A, X_hydrogen

       IF(isGasH2mixture == 0 .OR. igas == IGAS_HYDROGEN) THEN
          IF(PRESENT(W_a)) W_a = W_NonCondensGas(igas)
          IF(igas == IGAS_HYDROGEN) THEN
             IF(PRESENT(X_H2)) X_H2 = 1.D0
          ELSE
             IF(PRESENT(X_H2)) X_H2 = 0.D0
          ENDIF
          RETURN
       ENDIF
       !
       ! Mixture of gas with H2
       !
       Cp_A = Cp_NonCondensGas(igas)
       airmol_A = W_NonCondensGas(igas)
       Y_A = 1.d0 - Y_H2

       !     Calculate volume concentration
       X_hydrogen = (Y_H2/airmol_H2)/(Y_H2/airmol_H2+Y_A/airmol_A)
       X_A = 1.d0 - X_hydrogen
       !----------------------------------- Set equation of state ---------------

       airmol = X_A*airmol_A + X_hydrogen*airmol_H2

       ceoslp(25)=gascon/airmol  ! Gas constant for Ar-H2 mixture
       ceoslp(22) = Y_A*Cp_A + Y_H2*Cp_H2 ! Cp of Ar-H2 mixture

       ceoslp(17)=ceoslp(22)-ceoslp(25) ! Cv of Ar-H2 mixture
       ceoslp(18)=ceoslp(22)/ceoslp(17) ! Gamma of Ar-H2 mixture

       IF(PRESENT(W_a)) W_a = airmol
       IF(PRESENT(X_H2)) X_H2 = X_hydrogen
     END SUBROUTINE Mixture_Gas_H2

     REAL(kind=8) FUNCTION h_H2(T)
       REAL(kind=8), INTENT(in):: T ! Temperature
       REAL(kind=8), PARAMETER:: Cp_H2=14533.2d0 ! Specific heat at constant pressure

       h_H2 = Cp_H2*T
     END FUNCTION h_H2

     REAL(kind=8) FUNCTION e_H2(T)
       REAL(kind=8), INTENT(in):: T ! Temperature
       REAL(kind=8), PARAMETER:: Cv_H2=10163.d0 ! Specific heat at constant pressure

       e_H2 = Cv_H2*T
     END FUNCTION e_H2

     REAL(kind=8) FUNCTION star(qanty,expnt)
       REAL(kind=8), INTENT(in):: qanty,expnt
       star=EXP(expnt*dlog(qanty))
     END FUNCTION star


     SUBROUTINE SetLocalProps
       USE LOCAL_VARS
       IMPLICIT NONE
       !
       ! Set local variables to those from WATER_PROPS module
       !
       Ro1 = Ro_1; E1 = E_1 
       Ro2 = Ro_2; E2 = E_2 
       Roa = Ro_a
     END SUBROUTINE SetLocalProps


   END MODULE WATER_PROPS
