!!   ############################################################
     SUBROUTINE CH_AER_EMISSION(PFLUX, PRHODREF, HSV, KSV_CHSBEG,  PFCO)
!!   ############################################################
!!
!!    PURPOSE
!!    -------
!!    Transforme les emissions  d'aérosol en masse kg.kg-1.m.s-1 en molecules.m-2.s-1 : flux du moment m3
!!    Calcule les flux des moments m0 et m6 à partir de sigma  et Rg (um)
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (CNRM/GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!    none
!!
!!    EXTERNAL
!!    --------
!!    None
!!
USE MODD_CHS_AEROSOL
USE MODD_DST_SURF, ONLY : XDENSITY_DST
USE MODI_ABOR1_SFX
!!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PFLUX
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF
INTEGER,              INTENT(IN)    :: KSV_CHSBEG
 CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSV      ! name of chemical species
REAL, DIMENSION(:), OPTIONAL,    INTENT(IN) :: PFCO   ! CO flux 
!
!*      0.2    declarations local variables
!
REAL,DIMENSION(SIZE(PFLUX,1),NSP+NCARB+NSOA,JPMODE) :: ZFCTOTA
REAL,DIMENSION(SIZE(PFLUX,1),JPIN) :: ZFM
REAL,DIMENSION(SIZE(PFLUX,1)) :: ZFCO
REAL,DIMENSION(SIZE(PFLUX,1)) :: ZCONVERSION
REAL,DIMENSION(NSP+NCARB+NSOA) :: ZFAC, ZRHOI
REAL :: ZDEN2MOL
        !  ZDEN2MOL = 6.0221367E+23 * 1E-6 / 28.9644E-3
        !  conversion factor density to mol/cm3
        !  n_molec (moelc./cm3):  M = 1E-6*RHO(kg/m3) * XAVOGADRO / XMD
!
!
REAL   :: ZEMISRADIUSI, ZEMISRADIUSJ
REAL   :: ZVALBC, ZVALOC
!
INTEGER :: I_CH_M0i, I_CH_M0j, I_CH_M6i, I_CH_M6j, I_CH_H2Oi, I_CH_H2Oj,&
                  I_CH_SO4i,I_CH_SO4j, I_CH_NO3i, I_CH_NO3j, I_CH_NH3i, I_CH_NH3j,&
                  I_CH_OCi, I_CH_OCj, I_CH_BCi, I_CH_BCj  , I_CH_DSTi, I_CH_DSTj   
INTEGER :: JJ, JSV  ! loop counter
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
!-------------------------------------------------------------------------------
!
!*       1.     TRANSFER FROM GAS TO AEROSOL MODULE
!               ------------------------------------
!        1.1    initialisation 
!
IF (LHOOK) CALL DR_HOOK('CH_AER_EMISSION',0,ZHOOK_HANDLE)
!
! Initial aerosols fluxes have been transformed into molecu.m-2.s-1, 
! conversion into are in kg.kg-1.m.s-1 
!  conversion in kg.kg-1.m.s-1
ZCONVERSION(:) =  XAVOGADRO * PRHODREF(:)
!
!*       0      conversion into  kg.kg-1.m.s-1 (due to PCONVERSION)
!
ZVALBC  = 0.
ZVALOC  = 0.
ZFCO(:) = 0.
IF (LCO2PM .AND. PRESENT(PFCO)) THEN
  !ZVALBC=2.748549E-02  ! CO / BC conversion factor
  !ZVALOC=4.947248E-02  ! CO / POM conversion factor
  ZVALBC = 5.* 0.6E-10 / 0.4E-8  ! CO / BC conversion factor
  ZVALOC = 5.* 0.3E-10 / 0.4E-8  ! CO / POM conversion factor
  ZFCO(:) = PFCO(:)
END IF
!
!
PFLUX(:,:) = MAX(PFLUX(:,:),0.)
ZFCTOTA(:,:,:) = 0.
!
CALL GET_FCTOTA("SO4I  ",98.,I_CH_SO4i,ZFCTOTA(:,JP_AER_SO4,1))
CALL GET_FCTOTA("SO4J  ",98.,I_CH_SO4j,ZFCTOTA(:,JP_AER_SO4,2))
!
CALL GET_FCTOTA("NH3I  ",17.,I_CH_NH3i,ZFCTOTA(:,JP_AER_NH3,1))
CALL GET_FCTOTA("NH3J  ",17.,I_CH_NH3j,ZFCTOTA(:,JP_AER_NH3,2))
!
CALL GET_FCTOTA("NO3I  ",63.,I_CH_NO3i,ZFCTOTA(:,JP_AER_NO3,1))
CALL GET_FCTOTA("NO3J  ",63.,I_CH_NO3j,ZFCTOTA(:,JP_AER_NO3,2))
!
CALL GET_FCTOTA("H2OI  ",18.,I_CH_H2Oi,ZFCTOTA(:,JP_AER_H2O,1))
CALL GET_FCTOTA("H2OJ  ",18.,I_CH_H2Oj,ZFCTOTA(:,JP_AER_H2O,2))
!
CALL GET_FCTOTA("OCI   ",250.,I_CH_OCi,ZFCTOTA(:,JP_AER_OC,1))
CALL GET_FCTOTA("OCJ   ",250.,I_CH_OCj,ZFCTOTA(:,JP_AER_OC,2))
!
CALL GET_FCTOTA("BCI   ",250.,I_CH_BCi,ZFCTOTA(:,JP_AER_BC,1))
CALL GET_FCTOTA("BCJ   ",250.,I_CH_BCj,ZFCTOTA(:,JP_AER_BC,2))
!
CALL GET_FCTOTA("DSTI  ",100.,I_CH_DSTi,ZFCTOTA(:,JP_AER_DST,1))
CALL GET_FCTOTA("DSTJ  ",100.,I_CH_DSTj,ZFCTOTA(:,JP_AER_DST,2))
!
!
!*       1.1    calculate moment 3 flux from total aerosol mass
!
! Aerosol Density
! Cf Ackermann (all to black carbon except water)
ZRHOI(:)          = 1.8e3
ZRHOI(JP_AER_H2O) = 1.0e3   ! water
ZRHOI(JP_AER_DST) = XDENSITY_DST
DO JJ = 1,NSP+NCARB+NSOA
  ZFAC(JJ)=(4./3.) * 3.14292654 * ZRHOI(JJ) * 1.e-9
ENDDO
!
!
ZFM(:,2) = 0.
ZFM(:,5) = 0.
DO JJ = 1,NSP+NCARB+NSOA
  ZFM(:,2) = ZFM(:,2) + ZFCTOTA(:,JJ,1) / ZFAC(JJ)
  ZFM(:,5) = ZFM(:,5) + ZFCTOTA(:,JJ,2) / ZFAC(JJ)
ENDDO
!
!
IF (CRGUNIT=="MASS") THEN
  ZEMISRADIUSI = XEMISRADIUSI * EXP(-3.*(LOG(XEMISSIGI))**2)
  ZEMISRADIUSJ = XEMISRADIUSJ * EXP(-3.*(LOG(XEMISSIGJ))**2)
ELSE
  ZEMISRADIUSI = XEMISRADIUSI
  ZEMISRADIUSJ = XEMISRADIUSJ
END IF
!
!*       1.2    calculate moment 0 flux from dispersion and mean radius Rg
!
ZFM(:,1)= ZFM(:,2) / ((ZEMISRADIUSI**3)*EXP(4.5 * (LOG(XEMISSIGI))**2)) 
ZFM(:,4)= ZFM(:,5) / ((ZEMISRADIUSJ**3)*EXP(4.5 * (LOG(XEMISSIGJ))**2)) 
!
!*       1.3    calculate moment 6 flux from dispersion and mean diameter
!
ZFM(:,3) = ZFM(:,1) * (ZEMISRADIUSI**6) *EXP(18 *(LOG(XEMISSIGI))**2)
ZFM(:,6) = ZFM(:,4) * (ZEMISRADIUSJ**6) *EXP(18 *(LOG(XEMISSIGJ))**2)
!
!
I_CH_M0i=-999
I_CH_M0j=-999
I_CH_M6i=-999
I_CH_M6j=-999
DO JSV=1, size(HSV)
 IF (TRIM(HSV(JSV)) == "M0I") I_CH_M0i = JSV-KSV_CHSBEG+1
 IF (TRIM(HSV(JSV)) == "M0J") I_CH_M0j = JSV-KSV_CHSBEG+1
 IF (TRIM(HSV(JSV)) == "M6I") I_CH_M6i = JSV-KSV_CHSBEG+1
 IF (TRIM(HSV(JSV)) == "M6J") I_CH_M6j = JSV-KSV_CHSBEG+1
ENDDO
IF (I_CH_M0i ==-999) CALL ABOR1_SFX ('WRONG VALUE FOR I_CH_M0i ')
IF (I_CH_M0j ==-999) CALL ABOR1_SFX ('WRONG VALUE FOR I_CH_M0j ')
IF (I_CH_M6i ==-999) CALL ABOR1_SFX ('WRONG VALUE FOR I_CH_M6i ')
IF (I_CH_M6j ==-999) CALL ABOR1_SFX ('WRONG VALUE FOR I_CH_M6j ')
!
!*       1.4    conversion en ppp.m.s-1
!
ZDEN2MOL = 1E-6 * XAVOGADRO / XMD
!
! conversion in atmospheric unit only for moments 0 and 6 
! then conversion in  molecules.m-2.s-1
PFLUX(:,I_CH_M0i) = ZFM(:,1) * ZCONVERSION(:) * 1E-6 / (ZDEN2MOL * PRHODREF(:))
PFLUX(:,I_CH_M0j) = ZFM(:,4) * ZCONVERSION(:) * 1E-6 / (ZDEN2MOL * PRHODREF(:))
!
IF (LVARSIGI) PFLUX(:,I_CH_M6i) = ZFM(:,3) * ZCONVERSION(:) / (ZDEN2MOL * PRHODREF(:))
IF (LVARSIGJ) PFLUX(:,I_CH_M6j) = ZFM(:,6) * ZCONVERSION(:) / (ZDEN2MOL * PRHODREF(:))
!
IF (LHOOK) CALL DR_HOOK('CH_AER_EMISSION',1,ZHOOK_HANDLE)
!
CONTAINS 
!
SUBROUTINE GET_FCTOTA(H_CH,PMI,K_CH,PFCTOTA)
!
 CHARACTER(LEN=6), INTENT(IN) :: H_CH
REAL, INTENT(IN) :: PMI
INTEGER, INTENT(OUT) :: K_CH
REAL, DIMENSION(:), INTENT(INOUT) :: PFCTOTA
!
K_CH = -999
DO JSV=1,SIZE(HSV)
  IF ( TRIM(HSV(JSV)) == TRIM(H_CH) ) THEN
    K_CH = JSV - KSV_CHSBEG + 1
    EXIT
  ENDIF
ENDDO
IF ( K_CH == -999 ) CALL ABOR1_SFX ("WRONG VALUE FOR I_CH_"//H_CH)
!
IF ( TRIM(H_CH)=="OCI" ) THEN
  PFLUX(:,K_CH) = PFLUX(:,K_CH) + ( ZFCO(:) * ZVALOC / 2. )
ELSEIF ( TRIM(H_CH)=="OCJ" ) THEN
  PFLUX(:,K_CH) = PFLUX(:,K_CH) + ( ZFCO(:) * ZVALOC )
ELSEIF ( TRIM(H_CH)=="BCI" ) THEN
  PFLUX(:,K_CH) = PFLUX(:,K_CH) + ( ZFCO(:) * ZVALBC / 2. )
ELSEIF ( TRIM(H_CH)=="BCJ" ) THEN
  PFLUX(:,K_CH) = PFLUX(:,K_CH) + ( ZFCO(:) * ZVALBC )
ENDIF
!
!*       1.0    transfer aerosol mass from gas to aerosol variables
!               (and conversion of kg.kg-1.m.s-1 --> microgram.m-2.s-1)
!
! Initial aerosols fluxes have been transformed into molecu.m-2.s-1, 
! conversion into are in kg.kg-1.m.s-1
PFCTOTA(:) = PFLUX(:,K_CH) / ZCONVERSION(:) * PMI * 1E-3 * 1E+9 * PRHODREF(:)
!
! aerosol phase conversion molecu.m-2.s-1 into ppp.m.s-1
!
PFLUX(:,K_CH) = PFLUX(:,K_CH) * XMD / ( PMI * 1E-3 )
!
END SUBROUTINE GET_FCTOTA
!
END SUBROUTINE CH_AER_EMISSION
