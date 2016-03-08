!     #########
SUBROUTINE UNPACK_DIAG_PATCH_n(DIO, DGIP, GB, IO, TPSNOW, PKD, PK,  &
                               KMASK, KSIZE, KNPATCH, KPATCH,       &
                               PCPL_DRAIN, PCPL_RUNOFF, PCPL_EFLOOD,&
                               PCPL_PFLOOD, PCPL_IFLOOD, PCPL_ICEFLUX )  
!##############################################
!
!!****  *UNPACK_DIAG_PATCH_n* - unpacks ISBA diagnostics
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
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified    10/2004 by P. Le Moigne: Halstead coefficient
!!      Modified    10/2005 by P. Le Moigne: Deallocation (EBA)
!!      Modified    05/2008 by B. Decharme : Flooding scheme
!!      Modified    01/2010 by B. Decharme : new diag
!!      Modified      04-09 by A.L. Gibelin : Add carbon diagnostics
!!      Modified      05-09 by A.L. Gibelin : Add carbon spinup
!!      Modified    08/2012 by B. Decharme : optimization
!!      Modified    06/2013 by B. Decharme : add lateral drainage flux diag for DIF
!!                                           add tiotale sublimation flux
!!      Modified    10/2014 by P. Samuelsson: MEB
!!
!!------------------------------------------------------------------
!
USE MODD_TYPE_SNOW
!
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_INIT
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_INIT
USE MODD_DIAG_n, ONLY : DIAG_OPTIONS_t, DIAG_PATCH_t, DIAG_INIT
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_PACK_DIAG_ISBA, ONLY : PACK_DIAG_ISBA_t
USE MODD_PACK_ISBA, ONLY : PACK_ISBA_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DIO
TYPE(DIAG_PATCH_t), INTENT(INOUT) :: DGIP
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(SURF_SNOW), INTENT(INOUT) :: TPSNOW
TYPE(PACK_DIAG_ISBA_t), INTENT(INOUT) :: PKD
TYPE(PACK_ISBA_t), INTENT(INOUT) :: PK
!
INTEGER, INTENT(IN)                :: KSIZE, KPATCH, KNPATCH
INTEGER, DIMENSION(:), INTENT(IN)  :: KMASK
!
!Coupling variable
REAL, DIMENSION(:,:),  INTENT(OUT) :: PCPL_DRAIN
REAL, DIMENSION(:,:),  INTENT(OUT) :: PCPL_RUNOFF
REAL, DIMENSION(:,:),  INTENT(OUT) :: PCPL_EFLOOD
REAL, DIMENSION(:,:),  INTENT(OUT) :: PCPL_PFLOOD
REAL, DIMENSION(:,:),  INTENT(OUT) :: PCPL_IFLOOD
REAL, DIMENSION(:,:),  INTENT(OUT) :: PCPL_ICEFLUX
!
INTEGER :: JJ, JI, JSW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('UNPACK_DIAG_PATCH_N',0,ZHOOK_HANDLE)
!
IF (KNPATCH==1) THEN
  !
  DGIP%AL(KPATCH)%XTS   (:) = PKD%DGIP%XTS   (:)
  DGIP%AL(KPATCH)%XTSRAD(:) = PKD%DGIP%XTSRAD(:)
  DGIP%AL(KPATCH)%XALBT (:) = PKD%DGIP%XALBT (:)
  !
  IF (DIO%N2M>=1) THEN
    DGIP%AL(KPATCH)%XT2M    (:)    = PKD%DGIP%XT2M    (:)
    DGIP%AL(KPATCH)%XQ2M    (:)    = PKD%DGIP%XQ2M    (:)
    DGIP%AL(KPATCH)%XHU2M   (:)    = PKD%DGIP%XHU2M   (:)
    DGIP%AL(KPATCH)%XZON10M (:)    = PKD%DGIP%XZON10M (:)
    DGIP%AL(KPATCH)%XMER10M (:)    = PKD%DGIP%XMER10M (:)
    DGIP%AL(KPATCH)%XRI     (:)    = PKD%DGIP%XRI     (:)
!    
    DGIP%AL(KPATCH)%XWIND10M(:) = SQRT(PKD%DGIP%XZON10M(:)**2+PKD%DGIP%XMER10M(:)**2)
!    
  END IF
  !
  IF (DIO%LSURF_BUDGET) THEN
    DGIP%AL(KPATCH)%XRN    (:)    = PKD%DGIP%XRN         (:)
    DGIP%AL(KPATCH)%XH     (:)    = PKD%DGIP%XH          (:)
    DGIP%AL(KPATCH)%XGFLUX (:)    = PKD%DGIP%XGFLUX      (:)
    DGIP%AL(KPATCH)%XLEI   (:)    = PKD%DGIP%XLEI        (:)
    DGIP%AL(KPATCH)%XEVAP  (:)    = PKD%DGIP%XEVAP       (:)
    DGIP%AL(KPATCH)%XSUBL  (:)    = PKD%DGIP%XSUBL       (:)    
    DGIP%AL(KPATCH)%XSWD   (:)    = PKD%DGIP%XSWD        (:)
    DGIP%AL(KPATCH)%XSWU   (:)    = PKD%DGIP%XSWU        (:)
    DGIP%AL(KPATCH)%XLWD   (:)    = PKD%DGIP%XLWD        (:)
    DGIP%AL(KPATCH)%XLWU   (:)    = PKD%DGIP%XLWU        (:)
    DGIP%AL(KPATCH)%XFMU   (:)    = PKD%DGIP%XFMU        (:)
    DGIP%AL(KPATCH)%XFMV   (:)    = PKD%DGIP%XFMV        (:)
    !
    DGIP%AL(KPATCH)%XSWBD   (:, :) = PKD%DGIP%XSWBD  (:,:)
    DGIP%AL(KPATCH)%XSWBU   (:, :) = PKD%DGIP%XSWBU  (:,:)
    !
  END IF
  !
  IF (DIO%LCOEF) THEN
    DGIP%AL(KPATCH)%XCD  (:) = PKD%DGIP%XCD (:)
    DGIP%AL(KPATCH)%XCH  (:) = PKD%DGIP%XCH (:)
    DGIP%AL(KPATCH)%XCE  (:) = PKD%DGIP%XCE (:)
    DGIP%AL(KPATCH)%XZ0  (:) = PKD%DGIP%XZ0 (:)
    DGIP%AL(KPATCH)%XZ0H (:) = PKD%DGIP%XZ0H(:)
    DGIP%AL(KPATCH)%XZ0EFF(:)= PKD%DGIP%XZ0EFF(:)
  END IF
  !
  IF (DIO%LSURF_VARS) THEN
    DGIP%AL(KPATCH)%XQS(:)    = PKD%DGIP%XQS(:)
  END IF
  !
  IF (IO%LCPL_RRM) THEN
    PCPL_DRAIN (:,KPATCH) = PKD%DGEIP%XDRAIN (:)
    PCPL_RUNOFF(:,KPATCH) = PKD%DGEIP%XRUNOFF(:)
  END IF
  !
  IF (IO%LFLOOD) THEN
    PCPL_EFLOOD(:,KPATCH) = PKD%DGEIP%XLE_FLOOD (:) / PK%I%IP%XPLVTT(:,1) &
                          + PKD%DGEIP%XLEI_FLOOD(:) / PK%I%IP%XPLSTT(:,1)
    PCPL_PFLOOD(:,KPATCH) = PKD%DGEIP%XPFLOOD(:)
    PCPL_IFLOOD(:,KPATCH) = PKD%DGEIP%XIFLOOD(:)
  END IF    
  !
  IF(IO%LCPL_RRM.AND.IO%LGLACIER)THEN
    PCPL_ICEFLUX   (:,KPATCH)    = PKD%DGEIP%XICEFLUX       (:)
  ENDIF
  !
  IF(TPSNOW%SCHEME=='3-L' .OR.TPSNOW%SCHEME=='CRO')THEN
   TPSNOW%TEMP(:,:,KPATCH) = PKD%DGMI%XSNOWTEMP(:,:)
   TPSNOW%TS  (:,KPATCH)   = PKD%DGMI%XSNOWTEMP(:,1)
  ENDIF
  !
  IF (IO%CPHOTO/='NON') THEN
    GB%XIACAN(:,:,KPATCH) = PKD%GB%XIACAN(:,:,1)
  ENDIF
  !
ELSE
  !
  DO JJ=1,KSIZE
     JI                      = KMASK     (JJ)
     DGIP%AL(KPATCH)%XALBT    (JI)  =  PKD%DGIP%XALBT     (JJ)
     DGIP%AL(KPATCH)%XTS    (JI)     = PKD%DGIP%XTS    (JJ)  
     DGIP%AL(KPATCH)%XTSRAD (JI)     = PKD%DGIP%XTSRAD    (JJ)  
  END DO
  IF (DIO%N2M>=1) THEN
    DO JJ=1,KSIZE
      JI                      = KMASK     (JJ)
      DGIP%AL(KPATCH)%XT2M    (JI)    = PKD%DGIP%XT2M    (JJ)
      DGIP%AL(KPATCH)%XQ2M    (JI)    = PKD%DGIP%XQ2M    (JJ)
      DGIP%AL(KPATCH)%XHU2M   (JI)    = PKD%DGIP%XHU2M   (JJ)
      DGIP%AL(KPATCH)%XZON10M (JI)    = PKD%DGIP%XZON10M (JJ)
      DGIP%AL(KPATCH)%XMER10M (JI)    = PKD%DGIP%XMER10M (JJ)
      DGIP%AL(KPATCH)%XRI     (JI)    = PKD%DGIP%XRI     (JJ)
      !     
      DGIP%AL(KPATCH)%XWIND10M(JI)  = SQRT(PKD%DGIP%XZON10M(JJ)**2+PKD%DGIP%XMER10M(JJ)**2)
      !      
    END DO
  END IF
  !
  IF (DIO%LSURF_BUDGET) THEN
    DO JJ=1,KSIZE
      JI                     = KMASK         (JJ)
      DGIP%AL(KPATCH)%XRN    (JI)    = PKD%DGIP%XRN         (JJ)
      DGIP%AL(KPATCH)%XH     (JI)    = PKD%DGIP%XH          (JJ)
      DGIP%AL(KPATCH)%XGFLUX (JI)    = PKD%DGIP%XGFLUX      (JJ)
      DGIP%AL(KPATCH)%XLEI   (JI)    = PKD%DGIP%XLEI        (JJ)
      DGIP%AL(KPATCH)%XEVAP      (JI)  =  PKD%DGIP%XEVAP       (JJ)
      DGIP%AL(KPATCH)%XSUBL      (JI)  =  PKD%DGIP%XSUBL       (JJ)      
      DGIP%AL(KPATCH)%XSWD   (JI)    = PKD%DGIP%XSWD        (JJ)
      DGIP%AL(KPATCH)%XSWU   (JI)    = PKD%DGIP%XSWU        (JJ)
      DGIP%AL(KPATCH)%XLWD   (JI)    = PKD%DGIP%XLWD        (JJ)
      DGIP%AL(KPATCH)%XLWU   (JI)    = PKD%DGIP%XLWU        (JJ)
      DGIP%AL(KPATCH)%XFMU   (JI)    = PKD%DGIP%XFMU        (JJ)
      DGIP%AL(KPATCH)%XFMV   (JI)    = PKD%DGIP%XFMV        (JJ)
      !
      DO JSW=1,SIZE(DGIP%AL(KPATCH)%XSWBD,2)
        DGIP%AL(KPATCH)%XSWBD   (JI, JSW) = PKD%DGIP%XSWBD  (JJ,JSW)
        DGIP%AL(KPATCH)%XSWBU   (JI, JSW) = PKD%DGIP%XSWBU  (JJ,JSW)
      END DO
      !
    END DO
  END IF
  !
  IF (DIO%LCOEF) THEN
    DO JJ=1,KSIZE
      JI                               = KMASK             (JJ)
      DGIP%AL(KPATCH)%XCD              (JI)    = PKD%DGIP%XCD             (JJ)
      DGIP%AL(KPATCH)%XCH              (JI)    = PKD%DGIP%XCH             (JJ)
      DGIP%AL(KPATCH)%XCE              (JI)    = PKD%DGIP%XCE             (JJ)
      DGIP%AL(KPATCH)%XZ0    (JI)    = PKD%DGIP%XZ0   (JJ)
      DGIP%AL(KPATCH)%XZ0H   (JI)    = PKD%DGIP%XZ0H  (JJ)
      DGIP%AL(KPATCH)%XZ0EFF           (JI)    = PKD%DGIP%XZ0EFF          (JJ)
    END DO
  END IF
  !
  IF (DIO%LSURF_VARS) THEN
    DO JJ=1,KSIZE
      JI                               = KMASK             (JJ)
      DGIP%AL(KPATCH)%XQS              (JI)    = PKD%DGIP%XQS             (JJ)
    END DO
  END IF
  !
  IF (IO%LCPL_RRM) THEN
    DO JJ=1,KSIZE
      JI                               = KMASK             (JJ)
      PCPL_DRAIN       (JI,KPATCH)    = PKD%DGEIP%XDRAIN          (JJ)
      PCPL_RUNOFF      (JI,KPATCH)    = PKD%DGEIP%XRUNOFF         (JJ)
    END DO
  END IF
  !
  IF (IO%LFLOOD) THEN
    DO JJ=1,KSIZE
      JI                               = KMASK                     (JJ)
      PCPL_EFLOOD      (JI,KPATCH)    = PKD%DGEIP%XLE_FLOOD (JJ) / PK%I%IP%XPLVTT(JJ,1) &
                                       + PKD%DGEIP%XLEI_FLOOD(JJ) / PK%I%IP%XPLSTT(JJ,1)
      PCPL_PFLOOD      (JI,KPATCH)    = PKD%DGEIP%XPFLOOD                 (JJ)
      PCPL_IFLOOD      (JI,KPATCH)    = PKD%DGEIP%XIFLOOD                 (JJ)
    END DO
  END IF
  !
  IF(IO%LGLACIER)THEN
    DO JJ=1,KSIZE
      JI                              = KMASK             (JJ)
      PCPL_ICEFLUX    (JI,KPATCH)    = PKD%DGEIP%XICEFLUX        (JJ)
    END DO          
  ENDIF
  !
  IF(TPSNOW%SCHEME=='3-L' .OR.TPSNOW%SCHEME=='CRO')THEN
    DO JJ=1,KSIZE
      JI                       = KMASK             (JJ)
     TPSNOW%TS    (JI,KPATCH)  = PKD%DGMI%XSNOWTEMP(JJ,1)
      DO JSW=1,SIZE(TPSNOW%TEMP,2)
       TPSNOW%TEMP(JI,JSW,KPATCH)  = PKD%DGMI%XSNOWTEMP(JJ,JSW)
      ENDDO
    ENDDO          
  ENDIF
  !  
  IF (IO%CPHOTO/='NON') THEN
    DO JJ=1,KSIZE
      JI                  = KMASK   (JJ)
      DO JSW=1,SIZE(GB%XIACAN,2)
         GB%XIACAN(JI,JSW,KPATCH) = PKD%GB%XIACAN(JJ,JSW,1)
      ENDDO
    ENDDO
  ENDIF
  !
ENDIF
!
!------------------------------------------------------------------------
!
CALL DIAG_INIT(PKD%DGI)
CALL DIAG_INIT(PKD%DGIP)
CALL DIAG_MISC_ISBA_INIT(PKD%DGMI)
CALL DIAG_EVAP_ISBA_INIT(PKD%DGEI)
CALL DIAG_EVAP_ISBA_INIT(PKD%DGEIP)
!
!PKD%DGIP%XCH           => NULL()
!PKD%DGIP%XCE           => NULL()
!PKD%DGIP%XCD           => NULL()
!PKD%DGIP%XCDN          => NULL()
!PKD%DGIP%XRI           => NULL()
!PKD%DGIP%XHU           => NULL()
!PKD%DGIP%XHUG          => NULL()
!PKD%DGIP%XALBT         => NULL()
!!
!PKD%DGIP%XRN           => NULL()
!PKD%DGIP%XH            => NULL()
!PKD%DGIP%XLEI          => NULL()
!PKD%DGIP%XGFLUX        => NULL()
!!
!PKD%DGEIP%XLEG          => NULL()
!PKD%DGEIP%XLEGI         => NULL()
!PKD%DGEIP%XLEV          => NULL()
!PKD%DGEIP%XLES          => NULL()
!PKD%DGEIP%XLER          => NULL()
!PKD%DGEIP%XLETR         => NULL()
!PKD%DGIP%XEVAP         => NULL()
!PKD%DGIP%XSUBL         => NULL()
!PKD%DGEIP%XRESTORE      => NULL()
!PKD%DGEIP%XDRAIN        => NULL()
!PKD%DGEIP%XQSB          => NULL()
!PKD%DGEIP%XRUNOFF       => NULL()
!PKD%DGEIP%XMELT         => NULL()
!PKD%DGEIP%XMELTADV      => NULL()
!PKD%DGMI%XSRSFC        => NULL()
!PKD%DGMI%XRRSFC        => NULL()
!PK%I%R%XSNOWFREE_ALB => NULL()
!!
!PKD%DGEIP%XHORT         => NULL()
!PKD%DGEIP%XDRIP         => NULL()
!PKD%DGEIP%XRRVEG        => NULL()
!PKD%DGEIP%XIRRIG_FLUX   => NULL()
!!
!PKD%DGIP%XSWBD         => NULL()
!PKD%DGIP%XSWBU         => NULL()
!!
!PKD%DGIP%XSWD          => NULL()
!PKD%DGIP%XSWU          => NULL()
!PKD%DGIP%XLWD          => NULL()
!PKD%DGIP%XLWU          => NULL()
!PKD%DGIP%XFMU          => NULL()
!PKD%DGIP%XFMV          => NULL()
!!
!PKD%DGIP%XZ0 => NULL()
!PKD%DGIP%XZ0H => NULL()
!PKD%DGIP%XZ0EFF        => NULL()
!!
!PKD%DGMI%XCG           => NULL()
!PKD%DGMI%XC1           => NULL()
!PKD%DGMI%XC2           => NULL()
!PKD%DGMI%XWGEQ         => NULL()
!PKD%DGMI%XCT           => NULL()
!PKD%DGMI%XRS           => NULL()
!PKD%DGIP%XHV           => NULL()
!PKD%DGIP%XQS           => NULL()
!!
!PKD%DGIP%XTS           => NULL()
!PKD%DGIP%XTSRAD        => NULL()
!!
!PKD%DGEIP%XRESP_AUTO    => NULL()
!PKD%DGEIP%XRESP_ECO     => NULL()
!PKD%DGEIP%XGPP          => NULL()
!PKD%DGMI%XFAPAR        => NULL()
!PKD%DGMI%XFAPIR        => NULL()
!PKD%DGMI%XFAPAR_BS     => NULL()
!PKD%DGMI%XFAPIR_BS     => NULL()
!!
!PKD%DGEIP%XIFLOOD       => NULL()
!PKD%DGEIP%XPFLOOD       => NULL()
!PKD%DGEIP%XLE_FLOOD     => NULL()
!PKD%DGEIP%XLEI_FLOOD    => NULL()
!!
!PKD%DGMI%XRNSNOW       => NULL()
!PKD%DGMI%XHSNOW        => NULL()
!PKD%DGMI%XHPSNOW       => NULL()
!PKD%DGMI%XGFLUXSNOW    => NULL()
!PKD%DGMI%XUSTARSNOW    => NULL()
!PKD%DGMI%XGRNDFLUX     => NULL()
!PKD%DGEIP%XLESL         => NULL()
!PKD%DGEIP%XSNDRIFT      => NULL()
!PKD%DGMI%XCDSNOW       => NULL()
!PKD%DGMI%XCHSNOW       => NULL()
!PKD%DGMI%XSNOWHMASS    => NULL()
!
!PKD%DIO%XRN      => NULL()
!PKD%DIO%XH       => NULL()
!PKD%DIO%XLE      => NULL()
!PKD%DIO%XLEI     => NULL()
!PKD%DIO%XGFLUX   => NULL()
!!
!PKD%DGEI%XLEG     => NULL()
!PKD%DGEI%XLEGI    => NULL()
!PKD%DGEI%XLEV     => NULL()
!PKD%DGEI%XLETR    => NULL()
!PKD%DGEI%XUSTAR   => NULL()
!PKD%DGEI%XLER     => NULL()
!
!PKD%DGMI%XSNOWLIQ      => NULL()
!PKD%DGMI%XSNOWDZ       => NULL()
!!
!PKD%DGMI%XSNOWTEMP     => NULL()
!!
!PK%I%R%XSNOWFREE_ALB_VEG=> NULL()
!PK%I%R%XSNOWFREE_ALB_SOIL=> NULL()
!!
!PKD%GB%XIACAN        => NULL()
!!
!PKD%DGIP%XT2M          => NULL()
!PKD%DGIP%XQ2M          => NULL()
!PKD%DGIP%XHU2M         => NULL()
!PKD%DGIP%XZON10M       => NULL()
!PKD%DGIP%XMER10M       => NULL()
!!
!PKD%DGMI%XSWI          => NULL()
!PKD%DGMI%XTSWI         => NULL()
!PKD%DGMI%XTWSNOW       => NULL()
!PKD%DGMI%XTDSNOW       => NULL()
!!
!PKD%DGEIP%XICEFLUX      => NULL()
!!
!PKD%DGEIP%XDWG          => NULL()
!PKD%DGEIP%XDWGI         => NULL()
!PKD%DGEIP%XDSWE         => NULL()
!PKD%DGEIP%XWATBUD       => NULL()
!!
!PKD%DGEIP%XSWUP       => NULL()
!! MEB stuff
!PKD%DGEIP%XSWNET_V       => NULL()
!PKD%DGEIP%XSWNET_G       => NULL()
!PKD%DGEIP%XSWNET_N       => NULL()
!PKD%DGEIP%XSWNET_NS       => NULL()
!PKD%DGEIP%XLWUP       => NULL()
!PKD%DGEIP%XLWNET_V       => NULL()
!PKD%DGEIP%XLWNET_G       => NULL()
!PKD%DGEIP%XLWNET_N       => NULL()
!PKD%DGEIP%XLEVCV       => NULL()
!PKD%DGEIP%XLESC       => NULL()
!PKD%DGEIP%XH_V_C       => NULL()
!PKD%DGEIP%XH_G_C       => NULL()
!PKD%DGEIP%XLETRGV       => NULL()
!PKD%DGEIP%XLETRCV       => NULL()
!PKD%DGEIP%XLERGV       => NULL()
!PKD%DGEIP%XLELITTER     => NULL()
!PKD%DGEIP%XLELITTERI    => NULL()
!PKD%DGEIP%XDRIPLIT      => NULL()
!PKD%DGEIP%XRRLIT       => NULL()
!PKD%DGEIP%XLERCV       => NULL()
!PKD%DGEIP%XH_C_A       => NULL()
!PKD%DGEIP%XH_N_C       => NULL()
!PKD%DGEIP%XLE_C_A       => NULL()
!PKD%DGEIP%XLE_V_C       => NULL()
!PKD%DGEIP%XLE_G_C       => NULL()
!PKD%DGEIP%XLE_N_C       => NULL()
!PKD%DGEIP%XEVAP_N_C       => NULL()
!PKD%DGEIP%XEVAP_G_C       => NULL()
!PKD%DGEIP%XSR_GN       => NULL()
!PKD%DGEIP%XMELTCV       => NULL()
!PKD%DGEIP%XFRZCV       => NULL()
!PKD%DGEIP%XSWDOWN_GN       => NULL()
!PKD%DGEIP%XLWDOWN_GN       => NULL()
!!
DEALLOCATE(PKD%XBLOCK_SIMPLE)
DEALLOCATE(PKD%XBLOCK_GROUND)
DEALLOCATE(PKD%XBLOCK_SNOW)
DEALLOCATE(PKD%XBLOCK_KSW)
DEALLOCATE(PKD%XBLOCK_ABC)
DEALLOCATE(PKD%XBLOCK_0)
DEALLOCATE(PKD%XBLOCK_00)
!
IF (LHOOK) CALL DR_HOOK('UNPACK_DIAG_PATCH_N',1,ZHOOK_HANDLE)
!------------------------------------------------------------------------
!
END SUBROUTINE UNPACK_DIAG_PATCH_n
