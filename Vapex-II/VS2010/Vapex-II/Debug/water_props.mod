	  ?!  |   k820309    �          12.0        h��R                                                                                                           
       D:\Projects\EREC\VAPEX-II\Vapex-II\src\MOD_WaterProps.f90 WATER_PROPS          THERMO ERROR HEV STAR SETEOS E_H2 H_H2 MIXTURE_GAS_H2 DR CEOSLP AEOS14 CEOS1 CEOS2 CEOS3 THERMO_DHLSDPA IGAS ISGASH2MIXTURE AIRMOL_AIR AIRMOL_H2 AIRMOL_HE AIRMOL_AR CP_AIR CP_H2 CP_HE CP_AR CP_NONCONDENSGAS W_NONCONDENSGAS GASCON AIRMOL        @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                	     
          @                                
     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                     
          @                                      
          @                                !     
          @                                "     
          @                                #     
          @                                $     
          @                                %     
          @                                &     
          @                                '     
                                           (     
         
               �?        1.D0                                        )     
           
                �               @                                *     
          @                                +     
          @                                ,     
          @                                -     
          @                                .     
          @                                /     
          @                                0     
                                            1                                       0                                         2                                      1                                         3                                      1                                         4                                      1                                         5                                      2                                         6                                      3                                         7                                      4#     @                                  8                    #     @                                  9                   #ITYPE :   #ISMIXTURE ;         
                                  :             
                                  ;       #     @                                  <                  #THERMO_PROPS%PRESENT =   #TLIQ >   #TGAS ?   #PTOT @   #PAIR A                                          =     PRESENT       
 @                              >     
        
 @                              ?     
        
 @                              @     
        
 @                              A     
  %     @                                B                     
   %     @                                C                     
   %     @                                D                     
   %     @                                E                     
   %     @                                F                     
   %     @                                G                     
   %     @                                H                    
   #T I         
  @                              I     
  %     @                                J                    
   #T K         
@ @                              K     
  %     @                               L                   
   #SATPRS%EXP M   #SATPRS%MAX N   #TEMP O                                          M     EXP                                        N     MAX         @                               O     
   %     @@ @                            P                   
   #SATTMP%LOG Q   #SATTMP%SQRT R   #SATTMP%MAX S   #PRES T                                          Q     LOG                                        R     SQRT                                        S     MAX         @                               T     
   %     @                               U                   
   #SATDER%MAX V   #PRES W   #TEMP X                                          V     MAX         @                               W     
           @                               X     
   #     @                                 Y                  #RHOLIQ%LOG Z   #RHOLIQ%IDINT [   #RHOLIQ%MAX \   #RHOLIQ%MIN ]   #P ^   #TL _   #RHOL `   #DRLDP a   #DRLDT b                                          Z     LOG                                        [     IDINT                                        \     MAX                                        ]     MIN                                         ^     
                                           _     
         D                                  `     
         D                                  a     
         D                                  b     
   %     @                                c                   
   #CPLL%MAX d   #H e   #P f                                          d     MAX       
                                 e     
        
                                 f     
  %     @                                g                    
   #T h   #P i   #PA j         
                                 h     
        
                                 i     
        
                                 j     
  %     @                                k                    
   #H l   #P m                                           l     
                                           m     
   %     @                                n                   
   #VISCV%MIN o   #H p   #P q   #ROV r   #TV s   #PA t                                          o     MIN                                         p     
                                           q     
                                           r     
                                           s     
                                           t     
   %     @                                u                   
   #H_AR_H2%PRESENT v   #T w   #Y_H2 x   #Y_H2_RESTORE y                                          v     PRESENT       
                                 w     
        
  @                              x     
        
 @                              y     
  #     @                                  z                       �   N      fn#fn !   �   �   b   uapp(WATER_PROPS    �  8       RO_1      8       DRO1_DT1    W  8       DRO1_DP    �  8       E_1    �  8       DE1_DT1    �  8       DE1_DP    7  8       RO_2    o  8       DRO2_DT2    �  8       DRO2_DP    �  8       DRO2_DPA      8       E_2    O  8       DE2_DT2    �  8       DE2_DP    �  8       DE2_DPA    �  8       RO_V    /  8       DROV_DT2    g  8       DROV_DP    �  8       DROV_DPA    �  8       E_V      8       DEV_DT2    G  8       DEV_DP      8       DEV_DPA    �  8       RO_A    �  8       DROA_DT2    '  8       DROA_DPA    _  8       E_A    �  8       DEA_DT2    �  8       DEA_DPA      8       DEA_DP    ?  8       Y_A    w  8       DYA_DT2    �  8       DYA_DP    �  8       DYA_DPA    	  8       T_SAT    W	  8       DTSAT_DP    �	  8       TSV    �	  8       DTSV_DPV    �	  8       DTSV_DP    7
  8       DTSV_DPA    o
  `       DPSDP    �
  \       DPSDPA    +  8       HL_SAT    c  8       DHL_SAT_DP    �  8       DHL_SAT_DPA    �  8       HV_SAT      8       DHV_SAT_DP    C  8       DHV_SAT_DPA    {  8       DHV_SAT_DPV    �  ]       IOP      ]       JSTART    m  ]       NCELLS    �  ]       IGAS_AIR    '  ]       IGAS_HYDROGEN    �  ]       IGAS_HELIUM    �  ]       IGAS_ARGON    >  D       INITWATERPROPS %   �  ^       SETNONCONDENSGASTYPE +   �  8   a   SETNONCONDENSGASTYPE%ITYPE /     8   a   SETNONCONDENSGASTYPE%ISMIXTURE    P  �       THERMO_PROPS %   �  <      THERMO_PROPS%PRESENT "     8   a   THERMO_PROPS%TLIQ "   J  8   a   THERMO_PROPS%TGAS "   �  8   a   THERMO_PROPS%PTOT "   �  8   a   THERMO_PROPS%PAIR    �  H       EOSPMIN    :  H       EOSPMAX    �  H       EOSTLIQMIN    �  H       EOSTLIQMAX      H       EOSTVAPMIN    Z  H       EOSTVAPMAX    �  O       HEVAPORATION    �  8   a   HEVAPORATION%T "   )  O       SATURATEDPRESSURE $   x  8   a   SATURATEDPRESSURE%T    �  r       SATPRS    "  8      SATPRS%EXP    Z  8      SATPRS%MAX    �  8   a   SATPRS%TEMP    �  �       SATTMP    M  8      SATTMP%LOG    �  9      SATTMP%SQRT    �  8      SATTMP%MAX    �  8   a   SATTMP%PRES    .  l       SATDER    �  8      SATDER%MAX    �  8   a   SATDER%PRES    
  8   a   SATDER%TEMP    B  �       RHOLIQ    �  8      RHOLIQ%LOG    /  :      RHOLIQ%IDINT    i  8      RHOLIQ%MAX    �  8      RHOLIQ%MIN    �  8   a   RHOLIQ%P      8   a   RHOLIQ%TL    I  8   a   RHOLIQ%RHOL    �  8   a   RHOLIQ%DRLDP    �  8   a   RHOLIQ%DRLDT    �  d       CPLL    U  8      CPLL%MAX    �  8   a   CPLL%H    �  8   a   CPLL%P    �  ^       CPVV1    [  8   a   CPVV1%T    �  8   a   CPVV1%P    �  8   a   CPVV1%PA      V       VISCL    Y  8   a   VISCL%H    �  8   a   VISCL%P    �  ~       VISCV    G  8      VISCV%MIN      8   a   VISCV%H    �  8   a   VISCV%P    �  8   a   VISCV%ROV    '  8   a   VISCV%TV    _  8   a   VISCV%PA    �  �       H_AR_H2        <      H_AR_H2%PRESENT    S   8   a   H_AR_H2%T    �   8   a   H_AR_H2%Y_H2 %   �   8   a   H_AR_H2%Y_H2_RESTORE    �   D       SETLOCALPROPS 