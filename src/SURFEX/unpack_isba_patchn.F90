!     #########
SUBROUTINE UNPACK_ISBA_PATCH_n (AG, I, PKI, &
                                KMASK,KSIZE,KPATCH)
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
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_ISBA_n, ONLY : ISBA_t
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
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(PACK_ISBA_t), INTENT(INOUT) :: PKI
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
IF (I%O%NPATCH==1) THEN
  I%R%TSNOW%WSNOW     (:, :, 1) = PKI%XP_SNOWSWE    (:, :)
  I%R%TSNOW%RHO       (:, :, 1) = PKI%XP_SNOWRHO    (:, :)
  I%R%TSNOW%ALB       (:, 1)    = PKI%XP_SNOWALB    (:)
  I%R%XWR             (:, 1)    = PKI%XP_WR         (:)
  I%R%XTG             (:, :, 1) = PKI%XP_TG         (:, :)
  I%R%XWG             (:, :, 1) = PKI%XP_WG         (:, :)
  I%R%XWGI            (:, :, 1) = PKI%XP_WGI        (:, :)
  I%R%XRESA           (:, 1)    = PKI%XP_RESA       (:) 
  I%IP%XPCPS           (:, 1)    = PKI%XP_CPS        (:) 
  I%IP%XPLVTT          (:, 1)    = PKI%XP_LVTT       (:) 
  I%IP%XPLSTT          (:, 1)    = PKI%XP_LSTT       (:) 
  I%M%T%XALBNIR         (:, 1)    = PKI%XP_ALBNIR     (:) 
  I%M%T%XALBVIS         (:, 1)    = PKI%XP_ALBVIS     (:) 
  I%M%T%XALBUV          (:, 1)    = PKI%XP_ALBUV      (:) 
  I%M%T%XALBNIR_VEG     (:, 1)    = PKI%XP_ALBNIR_VEG (:) 
  I%M%T%XALBVIS_VEG     (:, 1)    = PKI%XP_ALBVIS_VEG (:) 
  I%M%T%XALBUV_VEG      (:, 1)    = PKI%XP_ALBUV_VEG  (:) 
  I%M%A%XALBNIR_SOIL    (:, 1)    = PKI%XP_ALBNIR_SOIL(:) 
  I%M%A%XALBVIS_SOIL    (:, 1)    = PKI%XP_ALBVIS_SOIL(:) 
  I%M%A%XALBUV_SOIL     (:, 1)    = PKI%XP_ALBUV_SOIL (:) 
  I%M%T%XEMIS           (:, 1)    = PKI%XP_EMIS       (:) 
  I%IP%XZ0EFFIP        (:, 1)    = PKI%XP_Z0EFFIP    (:) 
  I%IP%XZ0EFFIM        (:, 1)    = PKI%XP_Z0EFFIM    (:) 
  I%IP%XZ0EFFJP        (:, 1)    = PKI%XP_Z0EFFJP    (:) 
  I%IP%XZ0EFFJM        (:, 1)    = PKI%XP_Z0EFFJM    (:) 
  I%R%XLE             (:, 1)    = PKI%XP_LE         (:)
  !
   IF(I%O%LMEB_PATCH(KPATCH))THEN
     I%R%XWRL            (:, 1)    = PKI%XP_WRL        (:)
     I%R%XWRLI           (:, 1)    = PKI%XP_WRLI       (:)
     I%R%XWRVN           (:, 1)    = PKI%XP_WRVN       (:)
     I%R%XTV             (:, 1)    = PKI%XP_TV         (:)
     I%R%XTL             (:, 1)    = PKI%XP_TL         (:)
     I%R%XTC             (:, 1)    = PKI%XP_TC         (:)
     I%R%XQC             (:, 1)    = PKI%XP_QC         (:)
   ELSE
! Please note that XLAI, XVEG, and XZ0 are not unpacked
! in the case of MEB.
     I%M%T%XLAI            (:, 1)    = PKI%XP_LAI        (:) 
     I%M%T%XVEG            (:, 1)    = PKI%XP_VEG        (:) 
     I%M%T%XZ0             (:, 1)    = PKI%XP_Z0         (:) 
   ENDIF
  !
  IF (I%O%LTR_ML) THEN
    I%R%XFAPARC         (:, 1)    = PKI%XP_FAPARC     (:)
    I%R%XFAPIRC         (:, 1)    = PKI%XP_FAPIRC     (:)
    I%R%XLAI_EFFC       (:, 1)    = PKI%XP_LAI_EFFC   (:)
    I%R%XMUS            (:, 1)    = PKI%XP_MUS        (:)
  ENDIF   
  !
  IF (I%O%CPHOTO/='NON') THEN
     I%R%XAN             (:, 1)    = PKI%XP_AN         (:)
     I%R%XANDAY          (:, 1)    = PKI%XP_ANDAY      (:)
     I%R%XANFM           (:, 1)    = PKI%XP_ANFM       (:)
     I%R%XBIOMASS        (:,:,1)   = PKI%XP_BIOMASS        (:,:)
     I%R%XRESP_BIOMASS   (:,:,1)   = PKI%XP_RESP_BIOMASS   (:,:)
  END IF
  !
  IF(I%O%CPHOTO=='NIT' .OR. I%O%CPHOTO=='NCB') THEN
     I%IP%XBSLAI_NITRO    (:,1)    =    PKI%XP_BSLAI_NITRO    (:)          
  END IF
  !
    IF(I%O%CPHOTO=='NCB') THEN
     I%IP%XINCREASE       (:,:,1)   =    PKI%XP_INCREASE       (:,:)
  END IF
  !
  IF(I%O%CRESPSL=='CNT') THEN
     I%R%XLITTER         (:,:,:,1) =    PKI%XP_LITTER         (:,:,:)
     I%R%XSOILCARB       (:,:,1)   =    PKI%XP_SOILCARB       (:,:)
     I%R%XLIGNIN_STRUC   (:,:,1)   =    PKI%XP_LIGNIN_STRUC   (:,:)
     I%IP%XTURNOVER       (:,:,1)   =    PKI%XP_TURNOVER       (:,:)
  END IF
  !
  IF(LAGRIP .AND. (I%O%CPHOTO=='NIT' .OR. I%O%CPHOTO=='LAI' .OR. I%O%CPHOTO=='LST' .OR. I%O%CPHOTO=='NCB') ) THEN
    AG%LIRRIDAY (:,1)  =    PKI%XP_LIRRIDAY (:)
  END IF
  !
  IF (I%R%TSNOW%SCHEME=='3-L' .OR. I%R%TSNOW%SCHEME=='CRO') THEN
     I%R%TSNOW%HEAT      (:, :, 1) = PKI%XP_SNOWHEAT   (:, :)
     I%R%TSNOW%EMIS      (:, 1)    = PKI%XP_SNOWEMIS   (:)
     I%R%TSNOW%AGE       (:, :, 1) = PKI%XP_SNOWAGE    (:, :)
     I%R%TSNOW%ALBVIS    (:, 1)    = PKI%XP_SNOWALBVIS (:)
     I%R%TSNOW%ALBNIR    (:, 1)    = PKI%XP_SNOWALBNIR (:)
     I%R%TSNOW%ALBFIR    (:, 1)    = PKI%XP_SNOWALBFIR (:)     
  END IF

  IF (I%R%TSNOW%SCHEME=='CRO') THEN
     I%R%TSNOW%GRAN1     (:, :, 1) = PKI%XP_SNOWGRAN1   (:, :)
     I%R%TSNOW%GRAN2     (:, :, 1) = PKI%XP_SNOWGRAN2   (:, :)
     I%R%TSNOW%HIST      (:, :, 1) = PKI%XP_SNOWHIST    (:, :)
  END IF
  !
  IF(I%O%LGLACIER)THEN
     I%R%XICE_STO        (:,1)     = PKI%XP_ICE_STO    (:)
  ENDIF
!
ELSE
!
! Only save values for patches which are in use:
!
  DO JJ=1,KSIZE
    JI                              = KMASK         (JJ)
    I%R%TSNOW%ALB       (JI, KPATCH)    = PKI%XP_SNOWALB    (JJ)
    I%R%XWR             (JI, KPATCH)    = PKI%XP_WR         (JJ)
    I%R%XRESA           (JI, KPATCH)    = PKI%XP_RESA       (JJ) 
    I%IP%XPCPS           (JI, KPATCH)    = PKI%XP_CPS        (JJ) 
    I%IP%XPLVTT          (JI, KPATCH)    = PKI%XP_LVTT       (JJ) 
    I%IP%XPLSTT          (JI, KPATCH)    = PKI%XP_LSTT       (JJ) 
    I%M%T%XALBNIR         (JI, KPATCH)    = PKI%XP_ALBNIR     (JJ) 
    I%M%T%XALBVIS         (JI, KPATCH)    = PKI%XP_ALBVIS     (JJ) 
    I%M%T%XALBUV          (JI, KPATCH)    = PKI%XP_ALBUV      (JJ) 
    I%M%T%XALBNIR_VEG     (JI, KPATCH)    = PKI%XP_ALBNIR_VEG (JJ) 
    I%M%T%XALBVIS_VEG     (JI, KPATCH)    = PKI%XP_ALBVIS_VEG (JJ) 
    I%M%T%XALBUV_VEG      (JI, KPATCH)    = PKI%XP_ALBUV_VEG  (JJ) 
    I%M%A%XALBNIR_SOIL    (JI, KPATCH)    = PKI%XP_ALBNIR_SOIL(JJ) 
    I%M%A%XALBVIS_SOIL    (JI, KPATCH)    = PKI%XP_ALBVIS_SOIL(JJ) 
    I%M%A%XALBUV_SOIL     (JI, KPATCH)    = PKI%XP_ALBUV_SOIL (JJ) 
    I%M%T%XEMIS           (JI, KPATCH)    = PKI%XP_EMIS       (JJ) 
    I%IP%XZ0EFFIP        (JI, KPATCH)    = PKI%XP_Z0EFFIP    (JJ) 
    I%IP%XZ0EFFIM        (JI, KPATCH)    = PKI%XP_Z0EFFIM    (JJ) 
    I%IP%XZ0EFFJP        (JI, KPATCH)    = PKI%XP_Z0EFFJP    (JJ) 
    I%IP%XZ0EFFJM        (JI, KPATCH)    = PKI%XP_Z0EFFJM    (JJ) 
    I%R%XLE             (JI, KPATCH)    = PKI%XP_LE         (JJ)
  !
  END DO
  !
  IF(I%O%LMEB_PATCH(KPATCH))THEN
    DO JJ=1,KSIZE
      JI                              = KMASK         (JJ)
      I%R%XWRL            (JI, KPATCH)    = PKI%XP_WRL        (JJ)
      I%R%XWRLI           (JI, KPATCH)    = PKI%XP_WRLI       (JJ)
      I%R%XWRVN           (JI, KPATCH)    = PKI%XP_WRVN       (JJ)
      I%R%XTV             (JI, KPATCH)    = PKI%XP_TV         (JJ)
      I%R%XTL             (JI, KPATCH)    = PKI%XP_TL         (JJ)
      I%R%XTC             (JI, KPATCH)    = PKI%XP_TC         (JJ)
      I%R%XQC             (JI, KPATCH)    = PKI%XP_QC         (JJ)
    END DO
  ELSE
! Please note that XLAI, XVEG, and XZ0 are not unpacked
! in the case of MEB yet. This must be done when interactive/carbon
! vegetation is activated for MEB.
    DO JJ=1,KSIZE
      JI                              = KMASK         (JJ)
      I%M%T%XLAI            (JI, KPATCH)    = PKI%XP_LAI        (JJ) 
      I%M%T%XVEG            (JI, KPATCH)    = PKI%XP_VEG        (JJ) 
      I%M%T%XZ0             (JI, KPATCH)    = PKI%XP_Z0         (JJ) 
    END DO
  ENDIF
  !
  DO JK=1,SIZE(I%R%XTG,2)
    DO JJ=1,KSIZE
      JI                      =    KMASK(JJ)
      I%R%XTG             (JI, JK, KPATCH) = PKI%XP_TG         (JJ, JK)
    ENDDO
  ENDDO
!  
  DO JK=1,SIZE(I%R%XWG,2)
    DO JJ=1,KSIZE
      JI                      =    KMASK(JJ)
      I%R%XWG             (JI, JK, KPATCH) = PKI%XP_WG         (JJ, JK)
      I%R%XWGI            (JI, JK, KPATCH) = PKI%XP_WGI        (JJ, JK)
    ENDDO
  ENDDO
!  
  DO JK=1,SIZE(PKI%XP_SNOWSWE,2)
    DO JJ=1,KSIZE
      JI                      =    KMASK(JJ)
      I%R%TSNOW%WSNOW     (JI, JK, KPATCH) = PKI%XP_SNOWSWE    (JJ, JK)
      I%R%TSNOW%RHO       (JI, JK, KPATCH) = PKI%XP_SNOWRHO    (JJ, JK)
    ENDDO
  ENDDO
  !
  IF (I%O%LTR_ML) THEN
    DO JJ=1,KSIZE
      JI                      =    KMASK(JJ)          
      I%R%XFAPARC         (JI, KPATCH)    = PKI%XP_FAPARC     (JJ)
      I%R%XFAPIRC         (JI, KPATCH)    = PKI%XP_FAPIRC     (JJ)
      I%R%XLAI_EFFC       (JI, KPATCH)    = PKI%XP_LAI_EFFC   (JJ)
      I%R%XMUS            (JI, KPATCH)    = PKI%XP_MUS        (JJ)
    ENDDO
  ENDIF  
  !
  IF (I%O%CPHOTO/='NON') THEN
    DO JJ=1,KSIZE
      JI                              = KMASK         (JJ)
      I%R%XAN             (JI, KPATCH)    = PKI%XP_AN         (JJ)
      I%R%XANDAY          (JI, KPATCH)    = PKI%XP_ANDAY      (JJ)
      I%R%XANFM           (JI, KPATCH)    = PKI%XP_ANFM       (JJ)
    ENDDO
    DO JK=1,SIZE(I%R%XBIOMASS,2)
      DO JJ=1,KSIZE
        JI                              = KMASK         (JJ)       
        I%R%XBIOMASS        (JI, JK, KPATCH) = PKI%XP_BIOMASS        (JJ, JK)
        I%R%XRESP_BIOMASS   (JI, JK, KPATCH) = PKI%XP_RESP_BIOMASS   (JJ, JK)
      ENDDO
    END DO
  END IF
  !
  IF (I%O%CPHOTO=='NIT' .OR. I%O%CPHOTO=='NCB') THEN
    DO JJ=1,KSIZE
      JI                                 = KMASK             (JJ)
      I%IP%XBSLAI_NITRO    (JI, KPATCH)       = PKI%XP_BSLAI_NITRO    (JJ)
    END DO
  END IF
  !
  IF (I%O%CPHOTO=='NCB') THEN
    DO JK=1,SIZE(I%IP%XINCREASE,2)
      DO JJ=1,KSIZE
        JI                                 = KMASK             (JJ)
        I%IP%XINCREASE       (JI, JK, KPATCH)   = PKI%XP_INCREASE       (JJ, JK)
      ENDDO
    END DO
  END IF
  !
  IF (I%O%CRESPSL=='CNT') THEN
    DO JL=1,SIZE(PKI%XP_LITTER,3)
      DO JK=1,SIZE(PKI%XP_LITTER,2)
        DO JJ=1,KSIZE
          JI                                 = KMASK             (JJ)
          I%R%XLITTER       (JI, JK, JL, KPATCH) = PKI%XP_LITTER         (JJ, JK, JL)
        ENDDO
      ENDDO
    ENDDO
    DO JK=1,SIZE(PKI%XP_SOILCARB,2)
      DO JJ=1,KSIZE
        JI                                 = KMASK             (JJ)
        I%R%XSOILCARB       (JI, JK, KPATCH)   = PKI%XP_SOILCARB       (JJ, JK)
      ENDDO
    ENDDO
    DO JK=1,SIZE(PKI%XP_LIGNIN_STRUC,2)
      DO JJ=1,KSIZE
        JI                                  = KMASK             (JJ)
        I%R%XLIGNIN_STRUC   (JI, JK, KPATCH)    = PKI%XP_LIGNIN_STRUC   (JJ, JK)
      ENDDO
    ENDDO
    DO JK=1,SIZE(PKI%XP_TURNOVER,2)
      DO JJ=1,KSIZE
        JI                      =    KMASK(JJ)
        I%IP%XTURNOVER       (JI, JK, KPATCH)    = PKI%XP_TURNOVER       (JJ, JK)
      ENDDO
    END DO
  END IF
  !
  IF(LAGRIP .AND. (I%O%CPHOTO=='NIT' .OR. I%O%CPHOTO=='LAI' .OR. I%O%CPHOTO=='LST' .OR. I%O%CPHOTO=='NCB') ) THEN
     DO JJ=1,KSIZE
       JI                    =  KMASK             (JJ)
       AG%LIRRIDAY (JI,KPATCH)  =  PKI%XP_LIRRIDAY       (JJ)
     END DO
  END IF
  !
  IF (I%R%TSNOW%SCHEME=='3-L' .OR. I%R%TSNOW%SCHEME=='CRO') THEN
    DO JK=1,SIZE(PKI%XP_SNOWHEAT,2)
      DO JJ=1,KSIZE
        JI                              = KMASK         (JJ)
        I%R%TSNOW%HEAT      (JI, JK, KPATCH) = PKI%XP_SNOWHEAT  (JJ, JK)
        I%R%TSNOW%AGE       (JI, JK, KPATCH) = PKI%XP_SNOWAGE   (JJ, JK)
      ENDDO
    ENDDO
    DO JJ=1,KSIZE
      JI                              = KMASK         (JJ)
      I%R%TSNOW%EMIS      (JI, KPATCH)    = PKI%XP_SNOWEMIS   (JJ)
      I%R%TSNOW%ALBVIS    (JI, KPATCH)    = PKI%XP_SNOWALBVIS (JJ)
      I%R%TSNOW%ALBNIR    (JI, KPATCH)    = PKI%XP_SNOWALBNIR (JJ)
      I%R%TSNOW%ALBFIR    (JI, KPATCH)    = PKI%XP_SNOWALBFIR (JJ)     
    END DO
  END IF

  IF (I%R%TSNOW%SCHEME=='CRO') THEN
    DO JK=1,SIZE(PKI%XP_SNOWGRAN1,2)
      DO JJ=1,KSIZE
        JI                              = KMASK         (JJ)
        I%R%TSNOW%GRAN1     (JI, JK, KPATCH) = PKI%XP_SNOWGRAN1   (JJ, JK)
        I%R%TSNOW%GRAN2     (JI, JK, KPATCH) = PKI%XP_SNOWGRAN2   (JJ, JK)
        I%R%TSNOW%HIST      (JI, JK, KPATCH) = PKI%XP_SNOWHIST    (JJ, JK)
      ENDDO
    END DO
  END IF
  !
  IF(I%O%LGLACIER)THEN
    DO JJ=1,KSIZE
       JI                   = KMASK     (JJ)
       I%R%XICE_STO(JI, KPATCH) = PKI%XP_ICE_STO(JJ)
    ENDDO
  ENDIF
!
END IF
!
!------------------------------------------------------------------
!
PKI%XP_Z0_O_Z0H     => NULL()
PKI%XP_EMIS         => NULL()
PKI%XP_ALBNIR       => NULL()
PKI%XP_ALBVIS       => NULL()
PKI%XP_ALBUV        => NULL()
PKI%XP_ALBNIR_VEG   => NULL()
PKI%XP_ALBVIS_VEG   => NULL()
PKI%XP_ALBUV_VEG    => NULL()
PKI%XP_ALBNIR_SOIL  => NULL()
PKI%XP_ALBVIS_SOIL  => NULL()
PKI%XP_ALBUV_SOIL   => NULL()
PKI%XP_Z0           => NULL()
PKI%XP_WRMAX_CF     => NULL()
PKI%XP_GAMMA        => NULL()
PKI%XP_CV           => NULL()
PKI%XP_RGL          => NULL()
PKI%XP_RUNOFFD      => NULL()
PKI%XP_Z0EFFIP      => NULL()
PKI%XP_Z0EFFIM      => NULL()
PKI%XP_Z0EFFJP      => NULL()
PKI%XP_Z0EFFJM      => NULL()
PKI%XP_WR           => NULL() 
PKI%XP_LAI          => NULL() 
PKI%XP_RESA         => NULL()
PKI%XP_CPS          => NULL()
PKI%XP_LVTT         => NULL()
PKI%XP_LSTT         => NULL()
PKI%XP_VEG          => NULL()
PKI%XP_SNOWALB      => NULL()
PKI%XP_SNOWALBVIS   => NULL()
PKI%XP_SNOWALBNIR   => NULL()
PKI%XP_SNOWALBFIR   => NULL()
PKI%XP_LE           => NULL() 
PKI%XP_PSN          => NULL()
PKI%XP_PSNG         => NULL()
PKI%XP_PSNV         => NULL()
PKI%XP_ALBNIR_DRY   => NULL()
PKI%XP_ALBVIS_DRY   => NULL()
PKI%XP_ALBUV_DRY    => NULL()
PKI%XP_ALBNIR_WET   => NULL()
PKI%XP_ALBVIS_WET   => NULL()
PKI%XP_ALBUV_WET    => NULL()
PKI%XP_RUNOFFB      => NULL()
PKI%XP_WDRAIN       => NULL()
PKI%XP_TAUICE       => NULL()
PKI%XP_Z0REL        => NULL()
PKI%XP_AOSIP        => NULL()
PKI%XP_AOSIM        => NULL()
PKI%XP_AOSJP        => NULL()
PKI%XP_AOSJM        => NULL()
PKI%XP_HO2IP        => NULL()
PKI%XP_HO2IM        => NULL()
PKI%XP_HO2JP        => NULL()
PKI%XP_HO2JM        => NULL()
PKI%XP_SSO_SLOPE    => NULL()
PKI%XP_GAMMAT       => NULL()
PKI%XP_TDEEP        => NULL() 
!
PKI%XP_CLAY         => NULL() 
PKI%XP_SAND         => NULL() 
PKI%XP_WFC          => NULL()
PKI%XP_WWILT        => NULL()
PKI%XP_WSAT         => NULL()
PKI%XP_CONDSAT      => NULL()
PKI%XP_DG           => NULL()
PKI%XP_WG           => NULL()
PKI%XP_WGI          => NULL()
!
PKI%XP_KSAT_ICE     => NULL()
PKI%XP_D_ICE        => NULL()
!
PKI%XP_VEGTYPE_PATCH=> NULL()
!
PKI%XP_TG           => NULL()
!
PKI%XP_SNOWSWE      => NULL()
PKI%XP_SNOWRHO      => NULL()
!
PKI%XP_DIR_ALB_WITH_SNOW=> NULL()
PKI%XP_SCA_ALB_WITH_SNOW=> NULL()
!
PKI%XP_FFLOOD       => NULL()
PKI%XP_PIFLOOD      => NULL()
PKI%XP_FF           => NULL()
PKI%XP_FFG          => NULL()
PKI%XP_FFV          => NULL()
PKI%XP_FFROZEN      => NULL()
PKI%XP_ALBF         => NULL()
PKI%XP_EMISF        => NULL()
!
PKI%XP_PSNV_A       => NULL()
!
PKI%XP_SNOWHEAT     => NULL()
PKI%XP_SNOWEMIS     => NULL() 
!
PKI%XP_SNOWGRAN1    => NULL()
PKI%XP_SNOWGRAN2    => NULL()
PKI%XP_SNOWHIST     => NULL()
PKI%XP_SNOWAGE      => NULL()
!
PKI%XP_ICE_STO      => NULL()
!
PKI%XP_FWTD         => NULL()
PKI%XP_WTD          => NULL()
!
PKI%XP_HCAPSOIL     => NULL()
!
PKI%XP_CONDDRY      => NULL()
PKI%XP_CONDSLD      => NULL()
!
PKI%XP_C4B          => NULL() 
PKI%XP_ACOEF        => NULL() 
PKI%XP_PCOEF        => NULL()
PKI%XP_CGSAT        => NULL() 
PKI%XP_C1SAT        => NULL() 
PKI%XP_C2REF        => NULL() 
PKI%XP_C4REF        => NULL()
PKI%XP_C3           => NULL() 
!
PKI%XP_MPOTSAT      => NULL()
PKI%XP_BCOEF        => NULL()
!
PKI%XP_ROOTFRAC     => NULL()
PKI%XP_DZG          => NULL()
PKI%XP_DZDIF        => NULL()
PKI%NK_WG_LAYER     => NULL()
PKI%XP_SOILWGHT     => NULL()
!
PKI%XP_RSMIN        => NULL()
!
PKI%XP_BSLAI        => NULL()
PKI%XP_LAIMIN       => NULL()
PKI%XP_SEFOLD       => NULL()
PKI%XP_H_TREE       => NULL()
PKI%XP_ANF          => NULL()
PKI%XP_ANMAX        => NULL()
PKI%XP_FZERO        => NULL()
PKI%XP_EPSO         => NULL()
PKI%XP_GAMM         => NULL()
PKI%XP_QDGAMM       => NULL()
PKI%XP_GMES         => NULL()
PKI%XP_RE25         => NULL()
PKI%XP_QDGMES       => NULL()
PKI%XP_T1GMES       => NULL()
PKI%XP_T2GMES       => NULL()
PKI%XP_AMAX         => NULL()
PKI%XP_QDAMAX       => NULL()
PKI%XP_T1AMAX       => NULL()
PKI%XP_T2AMAX       => NULL()
PKI%XP_FAPARC       => NULL()
PKI%XP_FAPIRC       => NULL()
PKI%XP_LAI_EFFC     => NULL()
PKI%XP_MUS          => NULL()
PKI%XP_AN           => NULL() 
PKI%XP_ANDAY        => NULL() 
PKI%XP_ANFM         => NULL() 
PKI%XP_GC           => NULL()
PKI%XP_LAT          => NULL()
PKI%XP_LON          => NULL()
PKI%XP_BIOMASS      => NULL()
PKI%XP_RESP_BIOMASS => NULL()
!
PKI%LP_STRESS       => NULL()
PKI%XP_F2I          => NULL()
PKI%XP_AH           => NULL()
PKI%XP_BH           => NULL()
PKI%XP_DMAX         => NULL()
!
PKI%TP_SEED         => NULL()
PKI%TP_REAP         => NULL()
PKI%XP_IRRIG        => NULL()
PKI%XP_WATSUP       => NULL()
!
PKI%XP_LIRRIDAY     => NULL()
PKI%XP_THRESHOLD    => NULL()
PKI%XP_LIRRIGATE    => NULL()
!
PKI%XP_CE_NITRO     => NULL()
PKI%XP_CF_NITRO     => NULL()
PKI%XP_CNA_NITRO    => NULL()
PKI%XP_BSLAI_NITRO  => NULL()
!
PKI%XP_INCREASE     => NULL()
PKI%XP_TAU_WOOD     => NULL()
!
PKI%XP_LITTER       => NULL()
PKI%XP_SOILCARB     => NULL()
PKI%XP_LIGNIN_STRUC => NULL()
PKI%XP_TURNOVER     => NULL()
!
PKI%XP_FSAT=> NULL()
PKI%XP_TOPQS=> NULL()
!
PKI%XP_MUF=> NULL()
!
PKI%XP_WRL          => NULL()
PKI%XP_WRLI         => NULL()
PKI%XP_WRVN         => NULL() 
PKI%XP_TV           => NULL() 
PKI%XP_TL           => NULL() 
PKI%XP_TC           => NULL() 
PKI%XP_QC           => NULL() 

PKI%XP_H_VEG        => NULL()
PKI%XP_RGLV         => NULL()
PKI%XP_GAMMAV       => NULL()
PKI%XP_WRMAX_CFV    => NULL()
PKI%XP_LAIV         => NULL()
PKI%XP_Z0V          => NULL()
PKI%XP_RSMINV       => NULL()
PKI%XP_ROOTFRACV    => NULL()
PKI%XP_GNDLITTER    => NULL()
PKI%XP_Z0LITTER     => NULL()
!
!
DEALLOCATE(PKI%LBLOCK_SIMPLE)
DEALLOCATE(PKI%LBLOCK_0)
DEALLOCATE(PKI%NBLOCK_SIMPLE)
DEALLOCATE(PKI%NBLOCK_0)
DEALLOCATE(PKI%TBLOCK_SIMPLE)
DEALLOCATE(PKI%TBLOCK_0)
DEALLOCATE(PKI%XBLOCK_SIMPLE)
DEALLOCATE(PKI%XBLOCK_GROUND)
DEALLOCATE(PKI%XBLOCK_VEGTYPE)
DEALLOCATE(PKI%XBLOCK_TG)
DEALLOCATE(PKI%XBLOCK_SNOW)
DEALLOCATE(PKI%XBLOCK_ALB)
DEALLOCATE(PKI%XBLOCK_2)
DEALLOCATE(PKI%XBLOCK_BIOMASS)
DEALLOCATE(PKI%XBLOCK_SOILCARB)
DEALLOCATE(PKI%XBLOCK_LITTLEVS)
DEALLOCATE(PKI%XBLOCK_LITTER)
DEALLOCATE(PKI%XBLOCK_0)
DEALLOCATE(PKI%XBLOCK_00)
DEALLOCATE(PKI%XBLOCK_000)
DEALLOCATE(PKI%XBLOCK_01)
!
IF (LHOOK) CALL DR_HOOK('UNPACK_ISBA_PATCH_N',1,ZHOOK_HANDLE)
!------------------------------------------------------------------
!
END SUBROUTINE UNPACK_ISBA_PATCH_n
