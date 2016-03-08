!     #########
SUBROUTINE PREP_TEB (DTCO, UG, U, USS, TOP, BOP, B, TG, TCP, T, GD, GR, &
                     HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!     #################################################################################
!
!!****  *PREP_TEB* - prepares TEB fields
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
!!      S. Riette   06/2009 PREP_TEB_CANOPY has no more argument
!!------------------------------------------------------------------
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SSO_n, ONLY : SSO_t
!
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_CANOPY_n, ONLY : CANOPY_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_t
!
USE MODI_PREP_HOR_TEB_FIELD
USE MODI_PREP_VER_TEB
USE MODI_PREP_OUTPUT_GRID
USE MODI_GET_LUOUT
USE MODI_PREP_TEB_CANOPY
USE MODI_PREP_TEB_GARDEN
USE MODI_PREP_TEB_GREENROOF
USE MODI_GOTO_WRAPPER_TEB_PATCH
!
USE MODN_PREP_TEB
!
USE MODD_READ_NAMELIST, ONLY : LNAM_READ
!
USE MODD_PREP,       ONLY : XZS_LS
!
USE MODD_PREP_TEB_GARDEN, ONLY : XWSNOW_GD, XRSNOW_GD, XTSNOW_GD, XLWCSNOW_GD, &
                                 XAGESNOW_GD
!
USE MODD_PREP_TEB_GREENROOF, ONLY : XWSNOW_GR, XRSNOW_GR, XTSNOW_GR, XLWCSNOW_GR, &
                                    XAGESNOW_GR
!
USE MODD_SURF_ATM,   ONLY : LVERTSHIFT
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_CLEAN_PREP_OUTPUT_GRID
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SSO_t), INTENT(INOUT) :: USS
!
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(GRID_t), INTENT(INOUT) :: TG
TYPE(CANOPY_t), INTENT(INOUT) :: TCP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_VEG_t), INTENT(INOUT) :: GD
TYPE(TEB_VEG_t), INTENT(INOUT) :: GR
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
!
!*      0.2    declarations of local variables
!
INTEGER :: ILUOUT
INTEGER :: JPATCH         ! TEB patch number
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
!*      1.     Default of configuration
!
!
IF (LHOOK) CALL DR_HOOK('PREP_TEB',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
 CALL PREP_OUTPUT_GRID(UG%G, TG, U%NSIZE_FULL, ILUOUT)
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading and horizontal interpolations
!
!
!
!* option for roads
!
TOP%CROAD_DIR = CROAD_DIR
TOP%CWALL_OPT = CWALL_OPT
!
DO JPATCH=1,TOP%NTEB_PATCH
  !
  CALL GOTO_WRAPPER_TEB_PATCH(JPATCH, B=B, T=T, TGDR=GD%R, TGDMT=GD%M%T, TGRR=GR%R, TGRMT=GR%M%T)    
  !*      2.0    Large scale orography
  !
  IF (JPATCH==1) &
    CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                            HPROGRAM,'ZS     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,1)


  !*      2.1    Water reservoirs
  !
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'WS_ROOF',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'WS_ROAD',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  !
  !*      2.2    Building temperature
  !
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'TI_BLD ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  !
  !*      2.3    Road deep temperature
  !
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'TI_ROAD',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  !
  !*      2.4    Temperature profiles
  !
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'T_ROAD ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'T_WALLA',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'T_WALLB',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'T_ROOF ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  !
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'T_WIN1 ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  IF (TOP%CBEM == 'BEM') THEN
    CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'QI_BLD ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
    CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'T_WIN2 ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
    CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'T_FLOOR',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
    CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'T_MASS ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  ENDIF  
  !*      2.5    Snow variables
  !
  T%CUR%TSNOW_ROOF%SCHEME='1-L'
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'SN_ROOF',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  T%CUR%TSNOW_ROAD%SCHEME='1-L'
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'SN_ROAD',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  !
  !*      2.6    Canyon air variables
  !
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'T_CAN  ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  CALL PREP_HOR_TEB_FIELD(B%CUR, BOP, DTCO, U, TG, T%CUR, TOP, &
                         HPROGRAM,'Q_CAN  ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  !
  !-------------------------------------------------------------------------------------
  !
  !*      3.     Vertical interpolations of all variables
  !
  IF(LVERTSHIFT)THEN
    CALL PREP_VER_TEB(B%CUR, T%CUR, TOP%XZS, TOP%CBEM)
  ENDIF
  !
  !-------------------------------------------------------------------------------------
  !
  !*      4.     Urban green areas
  !
  
  IF (TOP%LGARDEN)    CALL PREP_TEB_GARDEN(DTCO, UG, U, USS, TG, TOP, GD%R%CUR, GD%O, &
                                              GD%M%X, GD%M%T%CUR, GD%IP,  &
                                              HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
  IF (TOP%LGREENROOF) THEN
     CALL PREP_TEB_GREENROOF(DTCO, UG, U, USS, TG, TOP, GR%R%CUR, GR%O, GR%M%X, GR%M%T%CUR, GR%IP, &
                             HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE, JPATCH)
    !
    ! Initializing deep GR temp. with that of the outer layer of the structural roof 
    !
    GR%IP%XTDEEP(:) = T%CUR%XT_ROOF(:,1)  
    !
  ENDIF
  !
ENDDO
!
DEALLOCATE(XWSNOW_GD,XRSNOW_GD,XTSNOW_GD,XLWCSNOW_GD,XAGESNOW_GD)
DEALLOCATE(XWSNOW_GR,XRSNOW_GR,XTSNOW_GR,XLWCSNOW_GR,XAGESNOW_GR)
!
!-------------------------------------------------------------------------------------
!
!*      5.     Preparation of canopy air variables
!
TOP%LCANOPY = LTEB_CANOPY
IF (TOP%LCANOPY) CALL PREP_TEB_CANOPY(TCP, TG%NDIM)
!
DEALLOCATE(XZS_LS)
!
!-------------------------------------------------------------------------------------
 CALL CLEAN_PREP_OUTPUT_GRID
IF (LHOOK) CALL DR_HOOK('PREP_TEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_TEB
