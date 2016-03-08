!#############################################################
SUBROUTINE INIT_VEG_PGD_n (MSS, DTCO, U, IO, IP, P, MX, MT, AG, &
                           HPROGRAM, HSURF, KLUOUT, KI, KMONTH,               &
                           ODEEPSOIL, OPHYSDOMC, PTDEEP_CLI, PGAMMAT_CLI,     &
                           OAGRIP, PTHRESHOLD, HINIT, PCO2, PRHOA     )  
!#############################################################
!
!!****  *INIT_VEG_PGD_n_n* - routine to initialize ISBA
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      23/07/13     (Decharme) Surface / Water table depth coupling
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t
USE MODD_AGRI_n, ONLY : AGRI_t
!
USE MODD_SURF_ATM,       ONLY : LCPL_ARP
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_CSTS,           ONLY : XCPD, XLVTT, XLSTT
USE MODD_SNOW_PAR,       ONLY : XEMISSN
USE MODD_ISBA_PAR,       ONLY : XTAU_ICE
!
USE MODD_SGH_PAR,        ONLY : XICE_DEPH_MAX
!
USE MODE_COTWO,          ONLY : GAULEG
!
USE MODI_SURF_PATCH
USE MODI_GET_1D_MASK
USE MODI_CO2_INIT_n
USE MODI_SUBSCALE_Z0EFF
!
USE MODE_SOIL
!
USE MODI_HEATCAPZ
USE MODI_THRMCONDZ
USE MODI_ABOR1_SFX
USE MODI_DIF_LAYER
USE MODI_DRY_WET_SOIL_ALBEDOS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              --²-----------------------
!
!
TYPE(SSO_t), INTENT(INOUT) :: MSS
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PGD_t), INTENT(INOUT) :: P
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: MX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: MT
TYPE(AGRI_t), INTENT(INOUT) :: AG
!
 CHARACTER(LEN=6), INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=6), INTENT(IN)  :: HSURF     ! Type of surface
INTEGER, INTENT(IN)  :: KLUOUT
!
INTEGER, INTENT(IN)  :: KI
!
INTEGER, INTENT(IN)  :: KMONTH
!
LOGICAL, INTENT(IN) :: ODEEPSOIL
LOGICAL, INTENT(IN) :: OPHYSDOMC
REAL, DIMENSION(:), INTENT(IN) :: PTDEEP_CLI
REAL, DIMENSION(:), INTENT(IN) :: PGAMMAT_CLI
!
LOGICAL, INTENT(IN) :: OAGRIP
REAL, DIMENSION(:), INTENT(IN) :: PTHRESHOLD
!
 CHARACTER(LEN=3), INTENT(IN) :: HINIT
 !
REAL, DIMENSION(:), INTENT(IN) :: PCO2
REAL, DIMENSION(:), INTENT(IN) :: PRHOA
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JPATCH  ! loop counter on tiles
INTEGER :: JILU,JP, JMAXLOC    ! loop increment
INTEGER :: JLAYER  ! loop counter on layers
!
INTEGER :: ISIZE
!
REAL, DIMENSION(SIZE(PCO2))       :: ZCO2  ! CO2 concentration  (kg/kg)
!
INTEGER, DIMENSION(:), ALLOCATABLE :: IR_NATURE_P
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('INIT_VEG_PGD_n',0,ZHOOK_HANDLE)
!
!*       2.4    Fraction of each tile
!               ---------------------
!
ALLOCATE(IP%XPATCH         (KI,IO%NPATCH))
ALLOCATE(IP%XVEGTYPE_PATCH (KI,NVEGTYPE,IO%NPATCH))
ALLOCATE(IP%NSIZE_NATURE_P (IO%NPATCH))
ALLOCATE(IP%NR_NATURE_P    (KI,IO%NPATCH))
!
 CALL SURF_PATCH(IO%NPATCH,MX%XVEGTYPE,IP%XPATCH,IP%XVEGTYPE_PATCH)
!
!*       2.5    Masks for tiles
!               ---------------
!
IF (IO%XRM_PATCH/=0.) THEN
  !
  WRITE(KLUOUT,*) " REMOVE PATCH below 5 % add to dominant patch " 
  ! remove small fraction of PATCHES and add to MAIN PATCH
  DO JP = 1,KI
    !1) find most present patch maximum value 
    JMAXLOC = MAXVAL(MAXLOC(IP%XPATCH(JP,:)))
    !2) FIND small value of cover 
    DO JPATCH = 1,IO%NPATCH
      IF ( IP%XPATCH(JP,JPATCH)<IO%XRM_PATCH ) THEN
        IP%XPATCH(JP,JMAXLOC) = IP%XPATCH(JP,JMAXLOC) + IP%XPATCH(JP,JPATCH)
        IP%XPATCH(JP,JPATCH) = 0.0
       ENDIF
    ENDDO
  ENDDO
  !
ENDIF
!
DO JPATCH=1,IO%NPATCH
  IP%NSIZE_NATURE_P(JPATCH) = COUNT(IP%XPATCH(:,JPATCH) > 0.0)
ENDDO
!
IP%NR_NATURE_P(:,:) = 0
DO JPATCH=1,IO%NPATCH
  ALLOCATE(IR_NATURE_P(IP%NSIZE_NATURE_P(JPATCH)))
  CALL GET_1D_MASK(IP%NSIZE_NATURE_P(JPATCH),KI,IP%XPATCH(:,JPATCH),IR_NATURE_P)
  IP%NR_NATURE_P(:IP%NSIZE_NATURE_P(JPATCH),JPATCH) = IR_NATURE_P(:)
  DEALLOCATE(IR_NATURE_P)
ENDDO
!
!
!*       2.6    Miscellaneous fields for ISBA:
!               -----------------------------
!
!* default value for:
! lateral water flux, deep soil temperature climatology and its relaxation time-scale
!
ALLOCATE(IP%XTDEEP (KI))
ALLOCATE(IP%XGAMMAT(KI))
IP%XTDEEP (:) = XUNDEF
IP%XGAMMAT(:) = XUNDEF
!
IF (ODEEPSOIL) THEN
   DO JILU = 1, KI
      IP%XTDEEP (JILU) = PTDEEP_CLI (KMONTH)
      IP%XGAMMAT(JILU) = 1. / PGAMMAT_CLI(KMONTH)
   END DO
   !
   WRITE(KLUOUT,*)' LDEEPSOIL = ',ODEEPSOIL,' LPHYSDOMC = ',OPHYSDOMC
   WRITE(KLUOUT,*)' XTDEEP    = ',MINVAL(IP%XTDEEP(:)),MAXVAL(IP%XTDEEP(:))
   WRITE(KLUOUT,*)' XGAMMAT   = ',MINVAL(IP%XGAMMAT(:)),MAXVAL(IP%XGAMMAT(:))
ENDIF
!
!
!*       2.7    Irrigation
!               ----------
!
IF (OAGRIP) THEN
   ALLOCATE(AG%NIRRINUM(KI,IO%NPATCH))
   ALLOCATE(AG%LIRRIDAY(KI,IO%NPATCH))
   ALLOCATE(AG%LIRRIGATE(KI,IO%NPATCH))
   ALLOCATE(AG%XTHRESHOLDSPT(KI,IO%NPATCH))
   !
   AG%NIRRINUM (:,:) = 1
   AG%LIRRIDAY (:,:) = .FALSE.                          
   AG%LIRRIGATE(:,:) = .FALSE.                          
   !
   DO JILU = 1, KI
      DO JPATCH = 1, IO%NPATCH
         AG%XTHRESHOLDSPT(JILU,JPATCH) = PTHRESHOLD(AG%NIRRINUM(JILU,JPATCH))
      END DO
   END DO
ELSE
   ALLOCATE(AG%NIRRINUM(0,0))
   ALLOCATE(AG%LIRRIDAY(0,0))
   ALLOCATE(AG%LIRRIGATE(0,0))
   ALLOCATE(AG%XTHRESHOLDSPT(0,0))
ENDIF
!
!
!*       2.8    Additional fields for ISBA-AGS:
!               ------------------------------                        
!
IF(IO%CPHOTO /= 'NON' .AND. HINIT == 'ALL') THEN
  IF (IO%LTR_ML) THEN
    ISIZE = 10
  ELSE
    ISIZE = 3
  ENDIF
  ALLOCATE(IP%XABC(ISIZE))
  ALLOCATE(IP%XPOI(ISIZE))
  IP%XABC(:) = 0.
  IP%XPOI(:) = 0.          
  ZCO2(:) = PCO2(:) / PRHOA(:)
  ALLOCATE(IP%XANMAX        (KI,IO%NPATCH))
  ALLOCATE(IP%XFZERO        (KI,IO%NPATCH))
  ALLOCATE(IP%XEPSO         (KI,IO%NPATCH))
  ALLOCATE(IP%XGAMM         (KI,IO%NPATCH))
  ALLOCATE(IP%XQDGAMM       (KI,IO%NPATCH))
  ALLOCATE(IP%XQDGMES       (KI,IO%NPATCH))
  ALLOCATE(IP%XT1GMES       (KI,IO%NPATCH))
  ALLOCATE(IP%XT2GMES       (KI,IO%NPATCH))
  ALLOCATE(IP%XAMAX         (KI,IO%NPATCH))
  ALLOCATE(IP%XQDAMAX       (KI,IO%NPATCH))
  ALLOCATE(IP%XT1AMAX       (KI,IO%NPATCH))
  ALLOCATE(IP%XT2AMAX       (KI,IO%NPATCH))
  ALLOCATE(IP%XAH           (KI,IO%NPATCH))
  ALLOCATE(IP%XBH           (KI,IO%NPATCH))
  ALLOCATE(IP%XTAU_WOOD     (KI,IO%NPATCH))
  ALLOCATE(IP%XINCREASE     (KI,IO%NNBIOMASS,IO%NPATCH))
  ALLOCATE(IP%XTURNOVER     (KI,IO%NNBIOMASS,IO%NPATCH))
  CALL CO2_INIT_n(IO, IP, ZCO2, MT%XGMES, MT%XGC, MX%XDMAX    )

ELSEIF(IO%CPHOTO == 'NON' .AND. IO%LTR_ML)THEN ! Case for MEB
   ISIZE = 10
   ALLOCATE (IP%XABC(ISIZE))
   ALLOCATE (IP%XPOI(ISIZE)) ! Working
   IP%XABC(:) = 0.
   IP%XPOI(:) = 0.
   CALL GAULEG(0.0,1.0,IP%XABC,IP%XPOI,SIZE(IP%XABC))
   DEALLOCATE (IP%XPOI)
   ALLOCATE   (IP%XPOI(0))
ELSE
  ALLOCATE(IP%XABC(0))
  ALLOCATE(IP%XPOI(0))
  ALLOCATE(IP%XANMAX        (0,0))
  ALLOCATE(IP%XFZERO        (0,0))
  ALLOCATE(IP%XEPSO         (0,0))
  ALLOCATE(IP%XGAMM         (0,0))
  ALLOCATE(IP%XQDGAMM       (0,0))
  ALLOCATE(IP%XQDGMES       (0,0))
  ALLOCATE(IP%XT1GMES       (0,0))
  ALLOCATE(IP%XT2GMES       (0,0))
  ALLOCATE(IP%XAMAX         (0,0))
  ALLOCATE(IP%XQDAMAX       (0,0))
  ALLOCATE(IP%XT1AMAX       (0,0))
  ALLOCATE(IP%XT2AMAX       (0,0))
  ALLOCATE(IP%XAH           (0,0))
  ALLOCATE(IP%XBH           (0,0))
  ALLOCATE(IP%XTAU_WOOD     (0,0))
  ALLOCATE(IP%XINCREASE     (0,0,0))
  ALLOCATE(IP%XTURNOVER     (0,0,0))  
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       4.     Orographic roughness length
!               ---------------------------
!
ALLOCATE(MSS%XZ0EFFIP(KI,IO%NPATCH))
ALLOCATE(MSS%XZ0EFFIM(KI,IO%NPATCH))
ALLOCATE(MSS%XZ0EFFJP(KI,IO%NPATCH))
ALLOCATE(MSS%XZ0EFFJM(KI,IO%NPATCH))
!
MSS%XZ0EFFIP(:,:) = XUNDEF
MSS%XZ0EFFIM(:,:) = XUNDEF
MSS%XZ0EFFJP(:,:) = XUNDEF
MSS%XZ0EFFJM(:,:) = XUNDEF
!
ALLOCATE(MSS%XZ0REL  (KI))
!
IF (SIZE(MSS%XAOSIP)>0) CALL SUBSCALE_Z0EFF(MSS,MT%XZ0,.TRUE.) 
!
!-------------------------------------------------------------------------------
!
!*       5.1     Soil hydraulic characteristics:
!                -------------------------------
!
ALLOCATE(IP%XCONDSAT (KI,IO%NGROUND_LAYER,IO%NPATCH))
ALLOCATE(IP%XMPOTSAT (KI,IO%NGROUND_LAYER))
ALLOCATE(IP%XBCOEF   (KI,IO%NGROUND_LAYER))
ALLOCATE(IP%XWWILT   (KI,IO%NGROUND_LAYER)) ! wilting point
ALLOCATE(IP%XWFC     (KI,IO%NGROUND_LAYER)) ! field capacity
ALLOCATE(IP%XWSAT    (KI,IO%NGROUND_LAYER)) ! saturation
ALLOCATE(IP%XTAUICE  (KI))
!        
DO JLAYER=1,IO%NGROUND_LAYER
   IP%XBCOEF  (:,JLAYER) = BCOEF_FUNC     (P%XCLAY(:,JLAYER),P%XSAND(:,JLAYER),IO%CPEDOTF)
   IP%XMPOTSAT(:,JLAYER) = MATPOTSAT_FUNC (P%XCLAY(:,JLAYER),P%XSAND(:,JLAYER),IO%CPEDOTF)
   DO JPATCH=1,IO%NPATCH
      IP%XCONDSAT(:,JLAYER,JPATCH) = HYDCONDSAT_FUNC(P%XCLAY(:,JLAYER),P%XSAND(:,JLAYER),IO%CPEDOTF)
   ENDDO   
   IP%XWSAT (:,JLAYER) = WSAT_FUNC (P%XCLAY(:,JLAYER),P%XSAND(:,JLAYER),IO%CPEDOTF)
   IP%XWWILT(:,JLAYER) = WWILT_FUNC(P%XCLAY(:,JLAYER),P%XSAND(:,JLAYER),IO%CPEDOTF)
END DO
!
IF (IO%CISBA=='2-L' .OR. IO%CISBA=='3-L') THEN
  !  field capacity at hydraulic conductivity = 0.1mm/day
  IP%XWFC(:,:) = WFC_FUNC(P%XCLAY(:,:),P%XSAND(:,:),IO%CPEDOTF)
ELSE IF (IO%CISBA=='DIF') THEN
  !  field capacity at water potential = 0.33bar        
  IP%XWFC(:,:) = W33_FUNC(P%XCLAY(:,:),P%XSAND(:,:),IO%CPEDOTF)
END IF
!
IP%XTAUICE(:) = XTAU_ICE
!
IF (IO%CISBA=='2-L' .OR. IO%CISBA=='3-L') THEN
  ALLOCATE(IP%XCGSAT (KI))
  ALLOCATE(IP%XC1SAT (KI,IO%NPATCH))
  ALLOCATE(IP%XC2REF (KI,IO%NPATCH))
  ALLOCATE(IP%XC3    (KI,2,IO%NPATCH))
  ALLOCATE(IP%XC4B   (KI))
  ALLOCATE(IP%XACOEF (KI))
  ALLOCATE(IP%XPCOEF (KI))
  ALLOCATE(IP%XC4REF (KI,IO%NPATCH))
  IP%XCGSAT(:)  = CGSAT_FUNC(P%XCLAY(:,1),P%XSAND(:,1))
  IP%XC4B(:)    = C4B_FUNC(P%XCLAY(:,1))
  !
  IP%XACOEF(:)  = ACOEF_FUNC(P%XCLAY(:,1))
  IP%XPCOEF(:)  = PCOEF_FUNC(P%XCLAY(:,1))
  !
  DO JPATCH=1,IO%NPATCH
    IP%XC1SAT(:,JPATCH) = C1SAT_FUNC(P%XCLAY(:,1))
    IP%XC2REF(:,JPATCH) = C2REF_FUNC(P%XCLAY(:,1))         
    IP%XC4REF(:,JPATCH) = C4REF_FUNC(P%XCLAY(:,1),P%XSAND(:,1),       &
                                  MX%XDG(:,2,            JPATCH), &
                                  MX%XDG(:,IO%NGROUND_LAYER,JPATCH)  )
    IP%XC3     (:,1,JPATCH) = C3_FUNC(P%XCLAY(:,1))
    IP%XC3     (:,2,JPATCH) = C3_FUNC(P%XCLAY(:,2))

  END DO
  !
ELSE IF (IO%CISBA=='DIF') THEN
  !
  ALLOCATE(IP%XCGSAT (0))
  ALLOCATE(IP%XC1SAT (0,0))
  ALLOCATE(IP%XC2REF (0,0))
  ALLOCATE(IP%XC3    (0,0,0))
  ALLOCATE(IP%XC4B   (0))
  ALLOCATE(IP%XC4REF (0,0))
  ALLOCATE(IP%XACOEF (0))
  ALLOCATE(IP%XPCOEF (0))
  !
END IF
!
IF(IO%CRUNOFF=='SGH')THEN
!
  ALLOCATE(IP%XWD0   (KI,IO%NGROUND_LAYER))
  ALLOCATE(IP%XKANISO(KI,IO%NGROUND_LAYER))
!
  IF(IO%CISBA=='DIF')THEN
     IP%XWD0(:,:) = WFC_FUNC(P%XCLAY(:,:),P%XSAND(:,:),IO%CPEDOTF)
  ELSE
     IP%XWD0(:,:) = IP%XWWILT(:,:)
  ENDIF
  IP%XKANISO(:,:) = ANISO_FUNC(P%XCLAY(:,:))
!
ELSE
!
  ALLOCATE(IP%XWD0   (0,0))
  ALLOCATE(IP%XKANISO(0,0))
!
ENDIF
!
!*       5.2     Soil thermal characteristics:
!               --------------------------------
!
ALLOCATE(IP%XPCPS (KI,IO%NPATCH))
ALLOCATE(IP%XPLVTT(KI,IO%NPATCH))
ALLOCATE(IP%XPLSTT(KI,IO%NPATCH))
IP%XPCPS (:,:) = XCPD
IP%XPLVTT(:,:) = XLVTT
IP%XPLSTT(:,:) = XLSTT
!
!CSCOND used in soil.F90 and soildif.F90
!
IF (IO%CSCOND=='NP89'.AND.IO%CISBA=='DIF') THEN
   WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   WRITE(KLUOUT,*)'IF CISBA=DIF, CSCOND=NP89 is not available'
   WRITE(KLUOUT,*)'because not physic. CSCOND is put to PL98 '
   WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   WRITE(KLUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
ENDIF
!
IF (IO%CSCOND=='PL98'.OR.IO%CISBA=='DIF') THEN
  ALLOCATE(IP%XHCAPSOIL(KI,IO%NGROUND_LAYER))
  ALLOCATE(IP%XCONDDRY (KI,IO%NGROUND_LAYER))
  ALLOCATE(IP%XCONDSLD (KI,IO%NGROUND_LAYER))
  ! 
  CALL HEATCAPZ(P%XSAND,IP%XHCAPSOIL)
  CALL THRMCONDZ(P%XSAND,IP%XWSAT,IP%XCONDDRY,IP%XCONDSLD)
  !
ELSE
  ALLOCATE(IP%XHCAPSOIL(0,0))
  ALLOCATE(IP%XCONDDRY (0,0))
  ALLOCATE(IP%XCONDSLD (0,0))
END IF
!
!-------------------------------------------------------------------------------
!CPSURF used in drag.F90
!CPL_ARP used in drag.F90 and e_budget.F90
IF(IO%CCPSURF=='DRY'.AND.LCPL_ARP) THEN
  CALL ABOR1_SFX('CCPSURF=DRY must not be used with LCPL_ARP')
ENDIF
!
!*       6.1    Initialize hydrology
!               --------------------
!
ALLOCATE(IP%XRUNOFFD (KI,IO%NPATCH))
IP%XRUNOFFD(:,:)=XUNDEF
!
IF (IO%CISBA == 'DIF') THEN
!  
  ALLOCATE(IP%XDZG       (KI,IO%NGROUND_LAYER,IO%NPATCH))
  ALLOCATE(IP%XDZDIF     (KI,IO%NGROUND_LAYER,IO%NPATCH))
  ALLOCATE(IP%XSOILWGHT  (KI,IO%NGROUND_LAYER,IO%NPATCH))
  CALL DIF_LAYER(KI, IO, IP, MX )
!
   ALLOCATE(IP%XFWTD(KI))
   ALLOCATE(IP%XWTD (KI))
   IP%XFWTD(:) = 0.0
   IP%XWTD (:) = XUNDEF
!
ELSE
!    
  ALLOCATE(IP%XDZG       (0,0,0))
  ALLOCATE(IP%XDZDIF     (0,0,0))
  ALLOCATE(IP%XSOILWGHT  (0,0,0))
  DO JPATCH=1,IO%NPATCH
    WHERE(IP%XPATCH(:,JPATCH)>0.0)
      IP%XRUNOFFD(:,JPATCH) = MX%XDG(:,2,JPATCH)
    ENDWHERE
  END DO
!  
  IO%NLAYER_DUN=2
  IO%NLAYER_HORT=2
!
  ALLOCATE(IP%XFWTD(0))
  ALLOCATE(IP%XWTD (0))
!   
ENDIF
!
!Horton (also used by the flooding sheme)
! 
ALLOCATE(IP%XKSAT_ICE(KI,IO%NPATCH))
!
IF(IO%CISBA/='DIF')THEN
  MX%XD_ICE   (:,:)=MIN(MX%XDG(:,2,:),MX%XD_ICE(:,:))
  MX%XD_ICE   (:,:)=MAX(XICE_DEPH_MAX,MX%XD_ICE(:,:))
  IP%XKSAT_ICE(:,:)=IP%XCONDSAT(:,1,:)
ELSE
  MX%XD_ICE   (:,:)=0.0
  IP%XKSAT_ICE(:,:)=0.0
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       8.     Physiographic Radiative fields:  
!               ------------------------------
!
!
!* dry and wet bare soil albedos
!
ALLOCATE(IP%XALBNIR_DRY  (KI))
ALLOCATE(IP%XALBVIS_DRY  (KI))
ALLOCATE(IP%XALBUV_DRY   (KI))
ALLOCATE(IP%XALBNIR_WET  (KI))
ALLOCATE(IP%XALBVIS_WET  (KI))
ALLOCATE(IP%XALBUV_WET   (KI))
!
 CALL DRY_WET_SOIL_ALBEDOS(P%XSAND(:,1),P%XCLAY(:,1), MX%XVEGTYPE, IP )  
!
!
!
!*       2.9    Nitrogen version for isbaAgs
!               ------------------------------                        
!
IF (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='NCB') THEN
  ALLOCATE(IP%XBSLAI_NITRO            (KI,IO%NPATCH              ))
  WHERE ((MT%XCE_NITRO (:,:)*MT%XCNA_NITRO(:,:)+MT%XCF_NITRO (:,:)) /= 0. )
      IP%XBSLAI_NITRO(:,:) = 1. / (MT%XCE_NITRO (:,:)*MT%XCNA_NITRO(:,:)+MT%XCF_NITRO (:,:))
  ELSEWHERE
      IP%XBSLAI_NITRO(:,:) = XUNDEF
  ENDWHERE
ELSE
  ALLOCATE(IP%XBSLAI_NITRO (0,0))
ENDIF
!
IF (LHOOK) CALL DR_HOOK('INIT_VEG_PGD_n',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_VEG_PGD_n
