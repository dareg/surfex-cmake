!!   ########################
MODULE MODE_SLT_SURF
!!   ########################
!!

  USE MODD_SLT_SURF
  USE MODD_CSTS
  USE MODD_SLT_n

!
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
  USE PARKIND1  ,ONLY : JPRB
!
  IMPLICIT NONE
  PUBLIC

CONTAINS
  !!
  !!   ############################################################
  SUBROUTINE MASSFLUX2MOMENTFLUX_SLT(  &
            PFLUX,                        &![kg/m2/s] for M3, zero for other]
            PRHODREF,                    &![kg/m3] air density
            PEMISRADIUS,                 &![um] emitted radius for the different modes
            PEMISSIG                    &![-] emitted sigma for the different modes 
            )  
!!   ############################################################
!!
!!    PURPOSE
!!    -------
!!    Transform emissions in mass (kg/m2/sec) to emissions of moments which have
!!    a bit strange units
!!    MESONH carries the following units during transport:
!!    M0=#/molec_{air}
!!    M3=molec_{dst}/molec_{air}
!!    M6=um6/molec_{air}*1.d6
!!    The surface model should have (for sea salt)
!!    M0=#/m3*[kg_{dst}/mole_{dst}/XAVOGADRO]
!!    M3=kg/m3
!!    M6=um6/m3
!!
!!    REFERENCE
!!    ---------
!!    Tulet et al, ORILAM manuscript for transformation of modal parameters
!!    J. Geophys. Res., 110, D18201, doi:10.1029/2004JD005716
!!
!!    AUTHOR
!!    ------
!!    Alf Grini and Pierre TULET (CNRM/GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!    none
!!
!!    EXTERNAL
!!    --------
!!    None
!!
IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL,   DIMENSION(:,:),      INTENT(INOUT) :: PFLUX     !In; kg/m2/s (index #2, #5, #8 etc)
                                                        !Out: mole particles per mole air m/s *(MWdst/MWair*rhoair)(index #1)
                                                        !Out: kg/m2/s (index #2)
                                                        !Out: moles m6/moles air m/s *(MWdst/MWair*rhoair)(index #3)
REAL,   DIMENSION(:),        INTENT(IN)    :: PRHODREF  !I [kg/m3] density of air
REAL,   DIMENSION(:),        INTENT(IN)    :: PEMISRADIUS !I [um] emitted radius
REAL,   DIMENSION(:),        INTENT(IN)    :: PEMISSIG    !I [-] emitted sigma for the modes
!
!
!*      0.2    declarations local variables
!
REAL,DIMENSION(SIZE(PFLUX,1),3) :: ZFM               !Intermediate variable to get moments
REAL                            :: ZCONVERTFACM0     ![kg/mole*mole/molec] conversion factor for moment fluxes and used fluxes
REAL                            :: ZCONVERTFACM6

INTEGER   :: JMODE  ! Counter for sea salt modes
INTEGER   :: JSV_IDX ! Counter for sea salt scalar variables
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!MESONH carries the following units during transport:
!M0=#/molec_{air}
!M3=molec_{dst}/molec_{air}
!M6=um6/molec_{air}*1.d6

!The surface should get the following units
!M0=#/m3*MW_SLT/XAVOGADRO
!M3=kg/m3
!M6=um6/m3*1.d6 MW_SLT/XAVOGADRO

!Emissions of sea salt are in kg/m2/sec for mode 3 at this point

!Factor which is needed so that all gains normal units when leaving ground paramn
IF (LHOOK) CALL DR_HOOK('MODE_SLT_SURF:MASSFLUX2MOMENTFLUX_SLT',0,ZHOOK_HANDLE)
ZCONVERTFACM0 = XMOLARWEIGHT_SLT/XAVOGADRO !(kg_dst/mol_dst)/(molec/mol)

!Factor which is needed for moment 6, there is a factor 1.d6 transported around in M6 in MESONH
ZCONVERTFACM6 = XMOLARWEIGHT_SLT/XAVOGADRO*1.d6

!Initialize initermediate moments
ZFM(:,:)=0.

!
DO JMODE=1,NSLTMDE

   !Make index which is 0 for first mode, 3 for second, 6 for third etc
   IF (LVARSIG_SLT) THEN
       JSV_IDX = (JMODE-1)*3
   ELSE IF (LRGFIX_SLT) THEN
       JSV_IDX = JMODE-2
   ELSE
       JSV_IDX = (JMODE-1)*2
   END IF

   !IN THIS VERSION, MASS FLUX (kg/m2/sec) IS SENT IN INDEX #2, #5, #8 if 3 moments per mode
   !IF TWO MOMENTS PER MODE, MASS FLUX IS SENT AS INDEX #2, #4, #6

   !Get flux of number in #/m2/sec from flux of mass in kg/m2/sec
   ZFM(:,1)=PFLUX(:,JSV_IDX+2)                &!kg_{dst}/m2/sec
          /(4./3.*XPI*XDENSITY_SALT)            &!/(kg_{dst}/m^{3}_{dst)} ==> m^3_{dst}/m2/sec
          /PEMISRADIUS(JMODE)**3                &!*um^{-3} ==> #/m2/sec*(m3/um3)
          *1.d18                                &!um3/m3  ==> #/m2/sec
          *exp(-4.5*log(PEMISSIG(JMODE))*log(PEMISSIG(JMODE)))  !Take into account size distribution  
   
   ! Get flux of moment 6 consistent with the other moments
   ZFM(:,3) = ZFM(:,1)                          &![#/m3]
          * (PEMISRADIUS(JMODE)**6)               &!*um6 ==> um6/m2/sec 
          *EXP(18. *(LOG(PEMISSIG(JMODE)))**2)     !Take into account size distribution  

   !Get flux of Moment 0 in transport units
   IF (.NOT.(LRGFIX_SLT)) THEN
   PFLUX(:,JSV_IDX+1) = ZFM(:,1)            &!particles/m^2/sec
          *ZCONVERTFACM0                       !==> particles/m2/sec * kg_dst/m3_{air}  
   END IF
   
   ! Flux moment 6
   IF (LVARSIG_SLT) THEN
      PFLUX(:,JSV_IDX + 3) = ZFM(:,3)           &!um^6/m^2/sec
             *ZCONVERTFACM6                        !==>   
   ENDIF

   !Multiply with molecular weights so that you get back the units described above when
   !when multiply with the opposite variable in ground_paramn.f90
   !PFLUX(:,JSV_IDX+1) = PFLUX(:,JSV_IDX+1) * 100.E-3 * PRHODREF(:) / ZMD   !#_{aer}/molec_{air} m/s * kg_{aer}/m^3_{air}
   !IF (LVARSIG) PFLUX(:,JSV_IDX+3) = PFLUX(:,JSV_IDX+3) * 100.E-3 * PRHODREF(:) / ZMD   !um^6_{aer}/molec_{air}*cm^3/m^3 m/s kg_{aer}/m^3_{air}
   !
   !
ENDDO !Loop on modes
IF (LHOOK) CALL DR_HOOK('MODE_SLT_SURF:MASSFLUX2MOMENTFLUX_SLT',1,ZHOOK_HANDLE)
   !
END SUBROUTINE MASSFLUX2MOMENTFLUX_SLT

!**********************************************************************
!**********************************************************************
!**********************************************************************

SUBROUTINE SALTMOMENT2SIZE(       &
       PSVT,                          &!I [XX/m3] input scalar variables (moment of distribution)
        PRHODREF,                    &!I [kg/m3] density of air       
        PSIG1D,                      &!O [-] standard deviation of aerosol distribution
        PRG1D,                       &!O [um] number median diameter of aerosol distribution
        PN1D,                        &!O [#/m3] number concentration of aerosols
        PMASS1D,                     &!O [kg/m3] mass concentration of aerosol
        PM1D                        &!O aerosols moments 0, 3 and 6
       )  
  !!   ############################################################
  !!
  !!
  !!    PURPOSE
  !!    -------
  !!    Translate the three moments M0, M3 and M6 given in ppp into
  !!    Values which can be understood more easily (R, sigma, N, M)
  !!    At this point, M3 is in kg/m3, M0 in #/m3*(kg_{dst}/mole), M6 in um6/m3*1.d6*(kg_{dst}/mole)
  !!
  !!    All the moments have been transformed in MESONH (atmospheric model) so that the surface gets
  !!    M0 [#/m3] *XMOLARWEIGHT_SLT/XAVOGADRO
  !!    M3 [kg/m3]
  !!    M6 [um6/m3*1.d6] *XMOLARWEIGHT_SLT/XAVOGADRO
  !!   
  !!    REFERENCE
  !!    ---------
  !!    Tulet et al, ORILAM manuscript for transformation of modal parameters
  !!    J. Geophys. Res., 110, D18201, doi:10.1029/2004JD005716
  !!
  !!    AUTHOR
  !!    ------
  !!    Pierre TULET (LA)
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!    Alf Grini (CNRM)
  !!
  !!    EXTERNAL
  !!    --------
  !!    None
  !!
  IMPLICIT NONE
  !!
  !-------------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  !*      0.1    declarations of arguments
  !
  !INPUT
  REAL,       DIMENSION(:,:),  INTENT(IN)     :: PSVT      !I [ppp] moments in surface units
  REAL,       DIMENSION(:),    INTENT(IN)     :: PRHODREF  !I [kg/m3] density of air
  
  !OUTPUT
  REAL,       DIMENSION(:,:),  OPTIONAL, INTENT(OUT)     :: PSIG1D   !O [-] standard deviation
  REAL,       DIMENSION(:,:),  OPTIONAL, INTENT(OUT)     :: PRG1D    !O [um] number median diameter
  REAL,       DIMENSION(:,:),  OPTIONAL, INTENT(OUT)     :: PN1D     !O [#/m3] number concentration
  REAL,       DIMENSION(:,:),  OPTIONAL, INTENT(OUT)     :: PMASS1D  !O [kg_{aer}/m3] mass concentration
  REAL,       DIMENSION(:,:),  OPTIONAL, INTENT(OUT)     :: PM1D     !O aerosols moments (MESONH units)
  !
  !*      0.2    declarations local variables
  !
  REAL,DIMENSION(SIZE(PSVT,1), SIZE(PSVT,2)) :: ZSV       ! [sea salts moment concentration]
  REAL,DIMENSION(SIZE(PSVT,1))               :: ZSIGMA    ! [-] standard deviation
  REAL,DIMENSION(SIZE(PSVT,1))               :: ZRG       ! [um] number median diameter
  REAL,DIMENSION(SIZE(PSVT,1), JPMODE_SLT*3) :: ZM        ! [moments] local array for moments
  INTEGER,DIMENSION(JPMODE_SLT)              :: NM0       ! [idx] index for Mode 0 in passed variables
  INTEGER,DIMENSION(JPMODE_SLT)              :: NM3       ! [idx] indexes for Mode 3 in passed variables
  INTEGER,DIMENSION(JPMODE_SLT)              :: NM6       ! [idx] indexes for Mode 6 in passed variables
  INTEGER                               :: JN,JMODEIDX,JJ ! [idx] loop counters
  !-------------------------------------------------------------------------------
  !
  !Conversion factor for M0
  REAL                :: ZCONVERTFACM0         !==> to obtain #/m3
  !Conversion factor for M6 (there is always a factor 1.d6 floating around in MESONH for M6)
  REAL                :: ZCONVERTFACM6         !==> to get um6/m3
  integer             :: I !counter
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  !
  !        1.1    initialisation 

  
  !Get the conversion factors
  IF (LHOOK) CALL DR_HOOK('MODE_SLT_SURF:SALTMOMENT2SIZE',0,ZHOOK_HANDLE)
  ZCONVERTFACM0 = XAVOGADRO/XMOLARWEIGHT_SLT         !==> to obtain #/m3
  ZCONVERTFACM6 = XAVOGADRO/XMOLARWEIGHT_SLT *1.D-6  !==> to get um6/m3
  
  !Get scalar variable indexes
  
  !Save the moments in a local array
  ZSV(:,:) = MAX(PSVT(:,:), 1E-80)
  
  DO JN=1,JPMODE_SLT
     !Set counter for number, M3 and M6
     NM0(JN) = 1+(JN-1)*3
     NM3(JN) = 2+(JN-1)*3
     NM6(JN) = 3+(JN-1)*3
  END DO
  
  DO JN=1,JPMODE_SLT
     
     IF (LVARSIG_SLT) THEN ! give M6 (case of variable standard deviation)
        
        !Get number concentration (#/molec_{air}==>#/m3)
        ZM(:,NM0(JN))=                         &
               ZSV(:,1+(JN-1)*3)                  &!#/m3air*M_{dst}/avogadro 
               *ZCONVERTFACM0                      ! ==> #/m3  
        
        !calculate moment 3 from total aerosol mass in kg/m3 ==> um3/m3
        ZM(:,NM3(JN)) =                        &
               ZSV(:,2+(JN-1)*3)                  &!kg_{aer}/m3_{air}
               /(XPI*4./3.*XDENSITY_SLT)          &!==> m3_{dst}/m3_{air}
               *1.d18                              !==> um3_{dst}/m3_{air} (volume ==> 3rd moment)  
             
        !Calculate moment 6 from the sent value
        ZM(:,NM6(JN)) = ZSV(:,3+(JN-1)*3)      &!um6/m3_{air}*(cm3/m3)*M_{dst}/Avogadro
               *ZCONVERTFACM6                     !==> um6/m3  


        !Get sigma (only if sigma is allowed to vary)
        !Get intermediate values for sigma M3^2/(M0*M6) (ORILAM paper, eqn 8)
        ZSIGMA(:)=ZM(:,NM3(JN))**2/(ZM(:,NM0(JN))*ZM(:,NM6(JN)))
        !Limit the intermediate value, can not be larger than 1
        ZSIGMA(:)=MIN(1-1E-10,ZSIGMA(:))
        !Limit the value for intermediate, can not be smaller than 0
        ZSIGMA(:)=MAX(1E-10,ZSIGMA(:))
        !Calculate log(sigma)
        ZSIGMA(:)= LOG(ZSIGMA(:))
        !Finally get the real sigma the negative sign is because of 
        !The way the equation is written (M3^2/(M0*M6)) instead of (M0*M6)/M3^3
        ZSIGMA(:)= EXP(1./3.*SQRT(-ZSIGMA(:)))
        
     ELSE IF (LRGFIX_SLT) THEN ! compute M6 from M3, Rg and SIGMA
        !calculate moment 3 from total aerosol mass in kg/m3 ==> um3/m3
        ZM(:,NM3(JN)) =                        &
               ZSV(:,JN)                          &!kg_{aer}/m3_{air}
               /(XPI*4./3.*XDENSITY_SLT)          &!==> m3_{dst}/m3_{air}
               *1.d18                              !==> um3_{dst}/m3_{air} (volume ==> 3rd moment)  
        !Get the emitted sigma for this mode
        ZSIGMA(:) = XEMISSIG_SLT(JN)

        ZM(:,NM0(JN))= ZM(:,NM3(JN)) /&
         ((XEMISRADIUS_SLT(JN)**3)*EXP(4.5 * LOG(ZSIGMA(:))**2))  


     ELSE ! compute M6 from M0, M3 and SIGMA
       
        !Get number concentration (#/molec_{air}==>#/m3)
        ZM(:,NM0(JN))=                         &
               ZSV(:,1+(JN-1)*2)                  &!#/m3air*M_{dst}/avogadro 
               *ZCONVERTFACM0                      ! ==> #/m3  

        !calculate moment 3 from total aerosol mass in kg/m3 ==> um3/m3
        ZM(:,NM3(JN)) =                        &
               ZSV(:,2+(JN-1)*2)                  &!kg_{aer}/m3_{air}
               /(XPI*4./3.*XDENSITY_SLT)          &!==> m3_{dst}/m3_{air}
               *1.d18                              !==> um3_{dst}/m3_{air} (volume ==> 3rd moment)  


     END IF
        !Get the emitted sigma for this mode
        ZSIGMA(:) = XEMISSIG_SLT(JN)

        !Calculate moment 6 from this emitted sigma
        ZM(:,NM6(JN)) = ZM(:,NM0(JN)) &
               * ( (ZM(:,NM3(JN))/ZM(:,NM0(JN)))**(1./3.) &
               * exp(-(3./2.)*log(ZSIGMA(:))**2))**6 &
               * exp(18.*log(ZSIGMA(:))**2)  

     
     
     !Get number median radius (eqn. 7 in Orilam manuscript)
     ZRG(:)=     &
            (      &
            ZM(:,NM3(JN))*ZM(:,NM3(JN))*ZM(:,NM3(JN))*ZM(:,NM3(JN)) &
            /(ZM(:,NM6(JN))*ZM(:,NM0(JN))*ZM(:,NM0(JN))*ZM(:,NM0(JN))) &
            )                                                                        &
            ** XSIXTH   

     !
     !Give the sigma-values to the passed array
     IF(PRESENT(PSIG1D))THEN
        PSIG1D(:,JN) = ZSIGMA(:)
     ENDIF

     !Set the number concentrations in the passed array
     IF(PRESENT(PN1D))THEN
        PN1D(:,JN) = ZM(:,NM0(JN))
     ENDIF

    !Get the number median radius
    IF(PRESENT(PRG1D))THEN
       PRG1D(:,JN)= ZRG(:)
    ENDIF
    !
    
    IF(PRESENT(PMASS1D))THEN
       PMASS1D(:,JN)=     &
              ZM(:,NM0(JN))         &!#/m^3_{air}
              * XPI*4./3.           &
              * XDENSITY_SALT       &!==>kg/m^3_{aeros}/m^3_{air}
              * ZRG(:) * ZRG(:) * ZRG(:) &
              * XUM3TOM3             &!==>kg/m^3_{air}
              * exp(4.5*log(ZSIGMA(:))*log(ZSIGMA(:)))  
    ENDIF
!
 END DO  !Loop on modes

IF(PRESENT(PM1D)) PM1D(:,:) = ZM(:,:)

IF (LHOOK) CALL DR_HOOK('MODE_SLT_SURF:SALTMOMENT2SIZE',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE SALTMOMENT2SIZE


END MODULE MODE_SLT_SURF
