!     #########
SUBROUTINE UNPACK_ISBA_PATCH_n (AG, IO, IP, IMT, IMA, IR, ISS, PK, KMASK,KSIZE,KPATCH)
!##############################################
!
!!****  *UNPACK_ISBA_PATCH_n* - unpacks ISBA prognostic variables
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     A. Boone
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      B. Decharme    2008 Floodplains
!!      A.L. Gibelin 04/2009 : BIOMASS and RESP_BIOMASS arrays 
!!      A.L. Gibelin 04/2009 : TAU_WOOD for NCB option 
!!      A.L. Gibelin 05/2009 : Add carbon spinup
!!      A.L. Gibelin 06/2009 : Soil carbon variables for CNT option
!!      A.L. Gibelin 07/2009 : Suppress RDK and transform GPP as a diagnostic
!!      A.L. Gibelin 07/2009 : Suppress PPST and PPSTF as outputs
!!      B. Decharme  06/2013 : add lateral drainage flux diag for DIF
!!                             water table / surface coupling
!!      P. Samuelsson 02/2012 : MEB
!!
!!------------------------------------------------------------------
!

!
USE MODD_AGRI_n, ONLY : AGRI_t, AGRI_INIT
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_TIME_t, ISBA_PARAM_ALB_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t, ISBA_INIT
USE MODD_GRID_n, ONLY : GRID_INIT
USE MODD_SSO_n, ONLY : SSO_INIT, SSO_t
USE MODD_PACK_ISBA, ONLY : PACK_ISBA_t
!
USE MODD_AGRI,     ONLY :  LAGRIP

!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
TYPE(AGRI_t), INTENT(INOUT) :: AG
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: IMT
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: IMA
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(SSO_t), INTENT(INOUT) :: ISS
TYPE(PACK_ISBA_t), INTENT(INOUT) :: PK
!
INTEGER, INTENT(IN)               :: KSIZE, KPATCH
!
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
!
INTEGER JJ, JI, JK, JL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('UNPACK_ISBA_PATCH_N',0,ZHOOK_HANDLE)
IF (IO%NPATCH==1) THEN
!  I = PK%I
  !
  IF(LAGRIP .AND. (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='LAI' .OR. IO%CPHOTO=='LST' .OR. IO%CPHOTO=='NCB') ) THEN
    AG%LIRRIDAY (:,1)  =    PK%AG%LIRRIDAY (:,1)
  END IF

  IR%TSNOW%WSNOW     (:, :, 1) = PK%I%R%TSNOW%WSNOW    (:, :,1)
  IR%TSNOW%RHO       (:, :, 1) = PK%I%R%TSNOW%RHO    (:, :,1)
  IR%TSNOW%ALB       (:, 1)    = PK%I%R%TSNOW%ALB    (:,1)
  IR%XWR             (:, 1)    = PK%I%R%XWR         (:,1)
  IR%XTG             (:, :, 1) = PK%I%R%XTG         (:, :,1)
  IR%XWG             (:, :, 1) = PK%I%R%XWG         (:, :,1)
  IR%XWGI            (:, :, 1) = PK%I%R%XWGI        (:, :,1)
  IR%XRESA           (:, 1)    = PK%I%R%XRESA       (:,1) 
  IP%XPCPS           (:, 1)    = PK%I%IP%XPCPS        (:,1) 
  IP%XPLVTT          (:, 1)    = PK%I%IP%XPLVTT       (:,1) 
  IP%XPLSTT          (:, 1)    = PK%I%IP%XPLSTT       (:,1) 
  IMT%XALBNIR         (:, 1)    = PK%I%M%T%XALBNIR     (:,1) 
  IMT%XALBVIS         (:, 1)    = PK%I%M%T%XALBVIS     (:,1) 
  IMT%XALBUV          (:, 1)    = PK%I%M%T%XALBUV      (:,1) 
  IMT%XALBNIR_VEG     (:, 1)    = PK%I%M%T%XALBNIR_VEG (:,1) 
  IMT%XALBVIS_VEG     (:, 1)    = PK%I%M%T%XALBVIS_VEG (:,1) 
  IMT%XALBUV_VEG      (:, 1)    = PK%I%M%T%XALBUV_VEG  (:,1) 
  IMA%XALBNIR_SOIL    (:, 1)    = PK%I%M%A%XALBNIR_SOIL(:,1) 
  IMA%XALBVIS_SOIL    (:, 1)    = PK%I%M%A%XALBVIS_SOIL(:,1) 
  IMA%XALBUV_SOIL     (:, 1)    = PK%I%M%A%XALBUV_SOIL (:,1) 
  IMT%XEMIS           (:, 1)    = PK%I%M%T%XEMIS       (:,1) 
  ISS%XZ0EFFIP        (:, 1)    = PK%ISS%XZ0EFFIP    (:,1) 
  ISS%XZ0EFFIM        (:, 1)    = PK%ISS%XZ0EFFIM    (:,1) 
  ISS%XZ0EFFJP        (:, 1)    = PK%ISS%XZ0EFFJP    (:,1) 
  ISS%XZ0EFFJM        (:, 1)    = PK%ISS%XZ0EFFJM    (:,1) 
  IR%XLE             (:, 1)    = PK%I%R%XLE         (:,1)
  !
   IF(IO%LMEB_PATCH(KPATCH))THEN
     IR%XWRL            (:, 1)    = PK%I%R%XWRL        (:,1)
     IR%XWRLI           (:, 1)    = PK%I%R%XWRLI       (:,1)
     IR%XWRVN           (:, 1)    = PK%I%R%XWRVN       (:,1)
     IR%XTV             (:, 1)    = PK%I%R%XTV         (:,1)
     IR%XTL             (:, 1)    = PK%I%R%XTL         (:,1)
     IR%XTC             (:, 1)    = PK%I%R%XTC         (:,1)
     IR%XQC             (:, 1)    = PK%I%R%XQC         (:,1)
   ELSE
! Please note that XLAI, XVEG, and XZ0 are not unpacked
! in the case of MEB.
     IMT%XLAI            (:, 1)    = PK%I%M%T%XLAI        (:,1) 
     IMT%XVEG            (:, 1)    = PK%I%M%T%XVEG        (:,1) 
     IMT%XZ0             (:, 1)    = PK%I%M%T%XZ0         (:,1) 
   ENDIF
  !
  IF (IO%LTR_ML) THEN
    IR%XFAPARC         (:, 1)    = PK%I%R%XFAPARC     (:,1)
    IR%XFAPIRC         (:, 1)    = PK%I%R%XFAPIRC     (:,1)
    IR%XLAI_EFFC       (:, 1)    = PK%I%R%XLAI_EFFC   (:,1)
    IR%XMUS            (:, 1)    = PK%I%R%XMUS        (:,1)
  ENDIF   
  !
  IF (IO%CPHOTO/='NON') THEN
     IR%XAN             (:, 1)    = PK%I%R%XAN         (:,1)
     IR%XANDAY          (:, 1)    = PK%I%R%XANDAY      (:,1)
     IR%XANFM           (:, 1)    = PK%I%R%XANFM       (:,1)
     IR%XBIOMASS        (:,:,1)   = PK%I%R%XBIOMASS        (:,:,1)
     IR%XRESP_BIOMASS   (:,:,1)   = PK%I%R%XRESP_BIOMASS   (:,:,1)
  END IF
  !
  IF(IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='NCB') THEN
     IP%XBSLAI_NITRO    (:,1)    =    PK%I%IP%XBSLAI_NITRO    (:,1)          
  END IF
  !
    IF(IO%CPHOTO=='NCB') THEN
     IP%XINCREASE       (:,:,1)   =    PK%I%IP%XINCREASE       (:,:,1)
  END IF
  !
  IF(IO%CRESPSL=='CNT') THEN
     IR%XLITTER         (:,:,:,1) =    PK%I%R%XLITTER         (:,:,:,1)
     IR%XSOILCARB       (:,:,1)   =    PK%I%R%XSOILCARB       (:,:,1)
     IR%XLIGNIN_STRUC   (:,:,1)   =    PK%I%R%XLIGNIN_STRUC   (:,:,1)
     IP%XTURNOVER       (:,:,1)   =    PK%I%IP%XTURNOVER       (:,:,1)
  END IF
  !
  IF(LAGRIP .AND. (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='LAI' .OR. IO%CPHOTO=='LST' .OR. IO%CPHOTO=='NCB') ) THEN
    AG%LIRRIDAY (:,1)  =    PK%AG%LIRRIDAY (:,1)
  END IF
  !
  IF (IR%TSNOW%SCHEME=='3-L' .OR. IR%TSNOW%SCHEME=='CRO') THEN
     IR%TSNOW%HEAT      (:, :, 1) = PK%I%R%TSNOW%HEAT   (:, :,1)
     IR%TSNOW%EMIS      (:, 1)    = PK%I%R%TSNOW%EMIS   (:,1)
     IR%TSNOW%AGE       (:, :, 1) = PK%I%R%TSNOW%AGE    (:, :,1)
     IR%TSNOW%ALBVIS    (:, 1)    = PK%I%R%TSNOW%ALBVIS (:,1)
     IR%TSNOW%ALBNIR    (:, 1)    = PK%I%R%TSNOW%ALBNIR (:,1)
     IR%TSNOW%ALBFIR    (:, 1)    = PK%I%R%TSNOW%ALBFIR (:,1)     
  END IF

  IF (IR%TSNOW%SCHEME=='CRO') THEN
     IR%TSNOW%GRAN1     (:, :, 1) = PK%I%R%TSNOW%GRAN1   (:, :,1)
     IR%TSNOW%GRAN2     (:, :, 1) = PK%I%R%TSNOW%GRAN2   (:, :,1)
     IR%TSNOW%HIST      (:, :, 1) = PK%I%R%TSNOW%HIST    (:, :,1)
  END IF
  !
  IF(IO%LGLACIER)THEN
     IR%XICE_STO        (:,1)     = PK%I%R%XICE_STO    (:,1)
  ENDIF

ELSE
!
! Only save values for patches which are in use:
!
  DO JJ=1,KSIZE
    JI                              = KMASK         (JJ)
    IR%TSNOW%ALB       (JI, KPATCH)    = PK%I%R%TSNOW%ALB    (JJ,1)
    IR%XWR             (JI, KPATCH)    = PK%I%R%XWR         (JJ,1)
    IR%XRESA           (JI, KPATCH)    = PK%I%R%XRESA       (JJ,1) 
    IP%XPCPS           (JI, KPATCH)    = PK%I%IP%XPCPS        (JJ,1) 
    IP%XPLVTT          (JI, KPATCH)    = PK%I%IP%XPLVTT       (JJ,1) 
    IP%XPLSTT          (JI, KPATCH)    = PK%I%IP%XPLSTT       (JJ,1) 
    IMT%XALBNIR         (JI, KPATCH)    = PK%I%M%T%XALBNIR     (JJ,1) 
    IMT%XALBVIS         (JI, KPATCH)    = PK%I%M%T%XALBVIS     (JJ,1) 
    IMT%XALBUV          (JI, KPATCH)    = PK%I%M%T%XALBUV      (JJ,1) 
    IMT%XALBNIR_VEG     (JI, KPATCH)    = PK%I%M%T%XALBNIR_VEG (JJ,1) 
    IMT%XALBVIS_VEG     (JI, KPATCH)    = PK%I%M%T%XALBVIS_VEG (JJ,1) 
    IMT%XALBUV_VEG      (JI, KPATCH)    = PK%I%M%T%XALBUV_VEG  (JJ,1) 
    IMA%XALBNIR_SOIL    (JI, KPATCH)    = PK%I%M%A%XALBNIR_SOIL(JJ,1) 
    IMA%XALBVIS_SOIL    (JI, KPATCH)    = PK%I%M%A%XALBVIS_SOIL(JJ,1) 
    IMA%XALBUV_SOIL     (JI, KPATCH)    = PK%I%M%A%XALBUV_SOIL (JJ,1) 
    IMT%XEMIS           (JI, KPATCH)    = PK%I%M%T%XEMIS       (JJ,1) 
    ISS%XZ0EFFIP        (JI, KPATCH)    = PK%ISS%XZ0EFFIP    (JJ,1) 
    ISS%XZ0EFFIM        (JI, KPATCH)    = PK%ISS%XZ0EFFIM    (JJ,1) 
    ISS%XZ0EFFJP        (JI, KPATCH)    = PK%ISS%XZ0EFFJP    (JJ,1) 
    ISS%XZ0EFFJM        (JI, KPATCH)    = PK%ISS%XZ0EFFJM    (JJ,1) 
    IR%XLE             (JI, KPATCH)    = PK%I%R%XLE         (JJ,1)
  !
  END DO
  !
  IF(IO%LMEB_PATCH(KPATCH))THEN
    DO JJ=1,KSIZE
      JI                              = KMASK         (JJ)
      IR%XWRL            (JI, KPATCH)    = PK%I%R%XWRL        (JJ,1)
      IR%XWRLI           (JI, KPATCH)    = PK%I%R%XWRLI       (JJ,1)
      IR%XWRVN           (JI, KPATCH)    = PK%I%R%XWRVN       (JJ,1)
      IR%XTV             (JI, KPATCH)    = PK%I%R%XTV         (JJ,1)
      IR%XTL             (JI, KPATCH)    = PK%I%R%XTL         (JJ,1)
      IR%XTC             (JI, KPATCH)    = PK%I%R%XTC         (JJ,1)
      IR%XQC             (JI, KPATCH)    = PK%I%R%XQC         (JJ,1)
    END DO
  ELSE
! Please note that XLAI, XVEG, and XZ0 are not unpacked
! in the case of MEB yet. This must be done when interactive/carbon
! vegetation is activated for MEB.
    DO JJ=1,KSIZE
      JI                              = KMASK         (JJ)
      IMT%XLAI            (JI, KPATCH)    = PK%I%M%T%XLAI        (JJ,1) 
      IMT%XVEG            (JI, KPATCH)    = PK%I%M%T%XVEG        (JJ,1) 
      IMT%XZ0             (JI, KPATCH)    = PK%I%M%T%XZ0         (JJ,1) 
    END DO
  ENDIF
  !
  DO JK=1,SIZE(IR%XTG,2)
    DO JJ=1,KSIZE
      JI                      =    KMASK(JJ)
      IR%XTG             (JI, JK, KPATCH) = PK%I%R%XTG         (JJ, JK,1)
    ENDDO
  ENDDO
!  
  DO JK=1,SIZE(IR%XWG,2)
    DO JJ=1,KSIZE
      JI                      =    KMASK(JJ)
      IR%XWG             (JI, JK, KPATCH) = PK%I%R%XWG         (JJ, JK,1)
      IR%XWGI            (JI, JK, KPATCH) = PK%I%R%XWGI        (JJ, JK,1)
    ENDDO
  ENDDO
!  
  DO JK=1,SIZE(PK%I%R%TSNOW%WSNOW,2)
    DO JJ=1,KSIZE
      JI                      =    KMASK(JJ)
      IR%TSNOW%WSNOW     (JI, JK, KPATCH) = PK%I%R%TSNOW%WSNOW    (JJ, JK,1)
      IR%TSNOW%RHO       (JI, JK, KPATCH) = PK%I%R%TSNOW%RHO    (JJ, JK,1)
    ENDDO
  ENDDO
  !
  IF (IO%LTR_ML) THEN
    DO JJ=1,KSIZE
      JI                      =    KMASK(JJ)          
      IR%XFAPARC         (JI, KPATCH)    = PK%I%R%XFAPARC     (JJ,1)
      IR%XFAPIRC         (JI, KPATCH)    = PK%I%R%XFAPIRC     (JJ,1)
      IR%XLAI_EFFC       (JI, KPATCH)    = PK%I%R%XLAI_EFFC   (JJ,1)
      IR%XMUS            (JI, KPATCH)    = PK%I%R%XMUS        (JJ,1)
    ENDDO
  ENDIF  
  !
  IF (IO%CPHOTO/='NON') THEN
    DO JJ=1,KSIZE
      JI                              = KMASK         (JJ)
      IR%XAN             (JI, KPATCH)    = PK%I%R%XAN         (JJ,1)
      IR%XANDAY          (JI, KPATCH)    = PK%I%R%XANDAY      (JJ,1)
      IR%XANFM           (JI, KPATCH)    = PK%I%R%XANFM       (JJ,1)
    ENDDO
    DO JK=1,SIZE(IR%XBIOMASS,2)
      DO JJ=1,KSIZE
        JI                              = KMASK         (JJ)       
        IR%XBIOMASS        (JI, JK, KPATCH) = PK%I%R%XBIOMASS        (JJ, JK,1)
        IR%XRESP_BIOMASS   (JI, JK, KPATCH) = PK%I%R%XRESP_BIOMASS   (JJ, JK,1)
      ENDDO
    END DO
  END IF
  !
  IF (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='NCB') THEN
    DO JJ=1,KSIZE
      JI                                 = KMASK             (JJ)
      IP%XBSLAI_NITRO    (JI, KPATCH)       = PK%I%IP%XBSLAI_NITRO    (JJ,1)
    END DO
  END IF
  !
  IF (IO%CPHOTO=='NCB') THEN
    DO JK=1,SIZE(IP%XINCREASE,2)
      DO JJ=1,KSIZE
        JI                                 = KMASK             (JJ)
        IP%XINCREASE       (JI, JK, KPATCH)   = PK%I%IP%XINCREASE       (JJ, JK,1)
      ENDDO
    END DO
  END IF
  !
  IF (IO%CRESPSL=='CNT') THEN
    DO JL=1,SIZE(PK%I%R%XLITTER,3)
      DO JK=1,SIZE(PK%I%R%XLITTER,2)
        DO JJ=1,KSIZE
          JI                                 = KMASK             (JJ)
          IR%XLITTER       (JI, JK, JL, KPATCH) = PK%I%R%XLITTER         (JJ, JK, JL,1)
        ENDDO
      ENDDO
    ENDDO
    DO JK=1,SIZE(PK%I%R%XSOILCARB,2)
      DO JJ=1,KSIZE
        JI                                 = KMASK             (JJ)
        IR%XSOILCARB       (JI, JK, KPATCH)   = PK%I%R%XSOILCARB       (JJ, JK,1)
      ENDDO
    ENDDO
    DO JK=1,SIZE(PK%I%R%XLIGNIN_STRUC,2)
      DO JJ=1,KSIZE
        JI                                  = KMASK             (JJ)
        IR%XLIGNIN_STRUC   (JI, JK, KPATCH)    = PK%I%R%XLIGNIN_STRUC   (JJ, JK,1)
      ENDDO
    ENDDO
    DO JK=1,SIZE(PK%I%IP%XTURNOVER,2)
      DO JJ=1,KSIZE
        JI                      =    KMASK(JJ)
        IP%XTURNOVER       (JI, JK, KPATCH)    = PK%I%IP%XTURNOVER       (JJ, JK,1)
      ENDDO
    END DO
  END IF
  !
  IF(LAGRIP .AND. (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='LAI' .OR. IO%CPHOTO=='LST' .OR. IO%CPHOTO=='NCB') ) THEN
     DO JJ=1,KSIZE
       JI                    =  KMASK             (JJ)
       AG%LIRRIDAY (JI,KPATCH)  =  PK%AG%LIRRIDAY       (JJ,1)
     END DO
  END IF
  !
  IF (IR%TSNOW%SCHEME=='3-L' .OR. IR%TSNOW%SCHEME=='CRO') THEN
    DO JK=1,SIZE(PK%I%R%TSNOW%HEAT,2)
      DO JJ=1,KSIZE
        JI                              = KMASK         (JJ)
        IR%TSNOW%HEAT      (JI, JK, KPATCH) = PK%I%R%TSNOW%HEAT  (JJ, JK,1)
        IR%TSNOW%AGE       (JI, JK, KPATCH) = PK%I%R%TSNOW%AGE   (JJ, JK,1)
      ENDDO
    ENDDO
    DO JJ=1,KSIZE
      JI                              = KMASK         (JJ)
      IR%TSNOW%EMIS      (JI, KPATCH)    = PK%I%R%TSNOW%EMIS   (JJ,1)
      IR%TSNOW%ALBVIS    (JI, KPATCH)    = PK%I%R%TSNOW%ALBVIS (JJ,1)
      IR%TSNOW%ALBNIR    (JI, KPATCH)    = PK%I%R%TSNOW%ALBNIR (JJ,1)
      IR%TSNOW%ALBFIR    (JI, KPATCH)    = PK%I%R%TSNOW%ALBFIR (JJ,1)     
    END DO
  END IF

  IF (IR%TSNOW%SCHEME=='CRO') THEN
    DO JK=1,SIZE(PK%I%R%TSNOW%GRAN1,2)
      DO JJ=1,KSIZE
        JI                              = KMASK         (JJ)
        IR%TSNOW%GRAN1     (JI, JK, KPATCH) = PK%I%R%TSNOW%GRAN1   (JJ, JK,1)
        IR%TSNOW%GRAN2     (JI, JK, KPATCH) = PK%I%R%TSNOW%GRAN2   (JJ, JK,1)
        IR%TSNOW%HIST      (JI, JK, KPATCH) = PK%I%R%TSNOW%HIST    (JJ, JK,1)
      ENDDO
    END DO
  END IF
  !
  IF(IO%LGLACIER)THEN
    DO JJ=1,KSIZE
       JI                   = KMASK     (JJ)
       IR%XICE_STO(JI, KPATCH) = PK%I%R%XICE_STO(JJ,1)
    ENDDO
  ENDIF
!
END IF
!
!------------------------------------------------------------------
!
CALL ISBA_INIT(PK%I)
CALL GRID_INIT(PK%G)
CALL AGRI_INIT(PK%AG)
CALL SSO_INIT(PK%ISS)
!
!PK%I%M%X%XZ0_O_Z0H     => NULL()
!PK%I%M%T%XEMIS         => NULL()
!PK%I%M%T%XALBNIR       => NULL()
!PK%I%M%T%XALBVIS       => NULL()
!PK%I%M%T%XALBUV        => NULL()
!PK%I%M%T%XALBNIR_VEG   => NULL()
!PK%I%M%T%XALBVIS_VEG   => NULL()
!PK%I%M%T%XALBUV_VEG    => NULL()
!PK%I%M%A%XALBNIR_SOIL  => NULL()
!PK%I%M%A%XALBVIS_SOIL  => NULL()
!PK%I%M%A%XALBUV_SOIL   => NULL()
!PK%I%M%T%XZ0           => NULL()
!PK%I%M%T%XWRMAX_CF     => NULL()
!PK%I%M%T%XGAMMA        => NULL()
!PK%I%M%T%XCV           => NULL()
!PK%I%M%T%XRGL          => NULL()
!PK%I%IP%XRUNOFFD      => NULL()
!PK%I%IP%XZ0EFFIP      => NULL()
!PK%I%IP%XZ0EFFIM      => NULL()
!PK%I%IP%XZ0EFFJP      => NULL()
!PK%I%IP%XZ0EFFJM      => NULL()
!PK%I%R%XWR           => NULL() 
!PK%I%M%T%XLAI          => NULL() 
!PK%I%R%XRESA         => NULL()
!PK%I%IP%XPCPS          => NULL()
!PK%I%IP%XPLVTT         => NULL()
!PK%I%IP%XPLSTT         => NULL()
!PK%I%M%T%XVEG          => NULL()
!PK%I%R%TSNOW%ALB      => NULL()
!PK%I%R%TSNOW%ALBVIS   => NULL()
!PK%I%R%TSNOW%ALBNIR   => NULL()
!PK%I%R%TSNOW%ALBFIR   => NULL()
!PK%I%R%XLE           => NULL() 
!PK%I%R%XPSN          => NULL()
!PK%I%R%XPSNG         => NULL()
!PK%I%R%XPSNV         => NULL()
!PK%I%IP%XALBNIR_DRY   => NULL()
!PK%I%IP%XALBVIS_DRY   => NULL()
!PK%I%IP%XALBUV_DRY    => NULL()
!PK%I%IP%XALBNIR_WET   => NULL()
!PK%I%IP%XALBVIS_WET   => NULL()
!PK%I%IP%XALBUV_WET    => NULL()
!PK%I%P%XRUNOFFB      => NULL()
!PK%I%P%XWDRAIN       => NULL()
!PK%I%IP%XTAUICE       => NULL()
!PK%I%IP%XZ0REL        => NULL()
!PK%I%P%XAOSIP        => NULL()
!PK%I%P%XAOSIM        => NULL()
!PK%I%P%XAOSJP        => NULL()
!PK%I%P%XAOSJM        => NULL()
!PK%I%P%XHO2IP        => NULL()
!PK%I%P%XHO2IM        => NULL()
!PK%I%P%XHO2JP        => NULL()
!PK%I%P%XHO2JM        => NULL()
!PK%I%P%XSSO_SLOPE    => NULL()
!PK%I%IP%XGAMMAT       => NULL()
!PK%I%IP%XTDEEP        => NULL() 
!!
!PK%I%P%XCLAY         => NULL() 
!PK%I%P%XSAND         => NULL() 
!PK%I%IP%XWFC          => NULL()
!PK%I%IP%XWWILT        => NULL()
!PK%I%IP%XWSAT         => NULL()
!PK%I%IP%XCONDSAT      => NULL()
!PK%I%M%X%XDG           => NULL()
!PK%I%R%XWG           => NULL()
!PK%I%R%XWGI          => NULL()
!!
!PK%I%IP%XKSAT_ICE     => NULL()
!PK%I%M%X%XD_ICE        => NULL()
!!
!PK%I%IP%XVEGTYPE_PATCH=> NULL()
!!
!PK%I%R%XTG           => NULL()
!!
!PK%I%R%TSNOW%WSNOW      => NULL()
!PK%I%R%TSNOW%RHO      => NULL()
!!
!PK%I%I%XDIR_ALB_WITH_SNOW=> NULL()
!PK%I%I%XSCA_ALB_WITH_SNOW=> NULL()
!!
!PK%I%I%XFFLOOD       => NULL()
!PK%I%I%XPIFLOOD      => NULL()
!PK%I%I%XFF           => NULL()
!PK%I%I%XFFG          => NULL()
!PK%I%I%XFFV          => NULL()
!PK%I%I%XFFROZEN      => NULL()
!PK%I%I%XALBF         => NULL()
!PK%I%I%XEMISF        => NULL()
!!
!PK%I%R%XPSNV_A       => NULL()
!!
!PK%I%R%TSNOW%HEAT     => NULL()
!PK%I%R%TSNOW%EMIS     => NULL() 
!!
!PK%I%R%TSNOW%GRAN1    => NULL()
!PK%I%R%TSNOW%GRAN2    => NULL()
!PK%I%R%TSNOW%HIST     => NULL()
!PK%I%R%TSNOW%AGE      => NULL()
!!
!PK%I%R%XICE_STO      => NULL()
!!
!PK%I%IP%XFWTD         => NULL()
!PK%I%IP%XWTD          => NULL()
!!
!PK%I%IP%XHCAPSOIL     => NULL()
!!
!PK%I%IP%XCONDDRY      => NULL()
!PK%I%IP%XCONDSLD      => NULL()
!!
!PK%I%IP%XC4B          => NULL() 
!PK%I%IP%XACOEF        => NULL() 
!PK%I%IP%XPCOEF        => NULL()
!PK%I%IP%XCGSAT        => NULL() 
!PK%I%IP%XC1SAT        => NULL() 
!PK%I%IP%XC2REF        => NULL() 
!PK%I%IP%XC4REF        => NULL()
!PK%I%IP%XC3           => NULL() 
!!
!PK%I%IP%XMPOTSAT      => NULL()
!PK%I%IP%XBCOEF        => NULL()
!!
!PK%I%M%X%XROOTFRAC     => NULL()
!PK%I%IP%XDZG          => NULL()
!PK%I%IP%XDZDIF        => NULL()
!PK%I%M%X%NWG_LAYER     => NULL()
!PK%I%IP%XSOILWGHT     => NULL()
!!
!PK%I%M%T%XRSMIN        => NULL()
!!
!PK%I%M%T%XBSLAI        => NULL()
!PK%I%M%T%XLAIMIN       => NULL()
!PK%I%M%T%XSEFOLD       => NULL()
!PK%I%M%X%XH_TREE       => NULL()
!PK%I%IP%XANMAX        => NULL()
!PK%I%IP%XFZERO        => NULL()
!PK%I%IP%XEPSO         => NULL()
!PK%I%IP%XGAMM         => NULL()
!PK%I%IP%XQDGAMM       => NULL()
!PK%I%M%T%XGMES         => NULL()
!PK%I%M%X%XRE25         => NULL()
!PK%I%IP%XQDGMES       => NULL()
!PK%I%IP%XT1GMES       => NULL()
!PK%I%IP%XT2GMES       => NULL()
!PK%I%IP%XAMAX         => NULL()
!PK%I%IP%XQDAMAX       => NULL()
!PK%I%IP%XT1AMAX       => NULL()
!PK%I%IP%XT2AMAX       => NULL()
!PK%I%R%XFAPARC       => NULL()
!PK%I%R%XFAPIRC       => NULL()
!PK%I%R%XLAI_EFFC     => NULL()
!PK%I%R%XMUS          => NULL()
!PK%I%R%XAN           => NULL() 
!PK%I%R%XANDAY        => NULL() 
!PK%I%R%XANFM         => NULL() 
!PK%I%M%T%XGC           => NULL()
!PK%G%XLAT          => NULL()
!PK%G%XLON          => NULL()
!PK%I%R%XBIOMASS      => NULL()
!PK%I%R%XRESP_BIOMASS => NULL()
!!
!PK%I%M%T%LSTRESS       => NULL()
!PK%I%M%T%XF2I          => NULL()
!PK%I%IP%XAH           => NULL()
!PK%I%IP%XBH           => NULL()
!PK%I%M%X%XDMAX         => NULL()
!!
!PK%I%M%I%TSEED         => NULL()
!PK%I%M%I%TREAP         => NULL()
!PK%I%M%I%XIRRIG        => NULL()
!PK%I%M%I%XWATSUP       => NULL()
!!
!PK%AG%LIRRIDAY     => NULL()
!PK%AG%XTHRESHOLDSPT    => NULL()
!PK%AG%LIRRIGATE    => NULL()
!!
!PK%I%M%T%XCE_NITRO     => NULL()
!PK%I%M%T%XCF_NITRO     => NULL()
!PK%I%M%T%XCNA_NITRO    => NULL()
!PK%I%IP%XBSLAI_NITRO  => NULL()
!!
!PK%I%IP%XINCREASE     => NULL()
!PK%I%IP%XTAU_WOOD     => NULL()
!!
!PK%I%R%XLITTER       => NULL()
!PK%I%R%XSOILCARB     => NULL()
!PK%I%R%XLIGNIN_STRUC => NULL()
!PK%I%IP%XTURNOVER     => NULL()
!!
!PK%I%I%XFSAT=> NULL()
!PK%I%I%XTOPQS=> NULL()
!!
!PK%I%I%XMUF=> NULL()
!!
!PK%I%R%XWRL          => NULL()
!PK%I%R%XWRLI         => NULL()
!PK%I%R%XWRVN         => NULL() 
!PK%I%R%XTV           => NULL() 
!PK%I%R%XTL           => NULL() 
!PK%I%R%XTC           => NULL() 
!PK%I%R%XQC           => NULL() 
!
!PK%I%M%M%XH_VEG        => NULL()
!PK%I%M%M%XRGLGV         => NULL()
!PK%I%M%M%XGAMMAGV       => NULL()
!PK%I%M%M%XWRMAX_CFGV    => NULL()
!PK%I%M%M%XLAIGV         => NULL()
!PK%I%M%M%XZ0GV          => NULL()
!PK%I%M%M%XRSMINGV       => NULL()
!PK%I%M%M%XROOTFRACGV    => NULL()
!PK%I%M%M%XGNDLITTER    => NULL()
!PK%I%M%M%XZ0LITTER     => NULL()
!
!
DEALLOCATE(PK%LBLOCK_SIMPLE)
DEALLOCATE(PK%LBLOCK_0)
DEALLOCATE(PK%NBLOCK_SIMPLE)
DEALLOCATE(PK%NBLOCK_0)
DEALLOCATE(PK%TBLOCK_SIMPLE)
DEALLOCATE(PK%TBLOCK_0)
DEALLOCATE(PK%XBLOCK_SIMPLE)
DEALLOCATE(PK%XBLOCK_GROUND)
DEALLOCATE(PK%XBLOCK_VEGTYPE)
DEALLOCATE(PK%XBLOCK_TG)
DEALLOCATE(PK%XBLOCK_SNOW)
DEALLOCATE(PK%XBLOCK_ALB)
DEALLOCATE(PK%XBLOCK_2)
DEALLOCATE(PK%XBLOCK_BIOMASS)
DEALLOCATE(PK%XBLOCK_SOILCARB)
DEALLOCATE(PK%XBLOCK_LITTLEVS)
DEALLOCATE(PK%XBLOCK_LITTER)
DEALLOCATE(PK%XBLOCK_0)
DEALLOCATE(PK%XBLOCK_00)
DEALLOCATE(PK%XBLOCK_000)
DEALLOCATE(PK%XBLOCK_01)
!
IF (LHOOK) CALL DR_HOOK('UNPACK_ISBA_PATCH_N',1,ZHOOK_HANDLE)
!------------------------------------------------------------------
!
END SUBROUTINE UNPACK_ISBA_PATCH_n
