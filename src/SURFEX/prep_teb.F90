!     #########
SUBROUTINE PREP_TEB(HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
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
USE MODD_TEB_GREENROOF_n, ONLY : TGR => TEB_GREENROOF
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TGRPE => TEB_GREENROOF_PGD_EVOL
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TGRP => TEB_GREENROOF_PGD
!
USE MODD_ISBA_n, ONLY : I => ISBA
USE MODD_SURF_ATM_SSO_n, ONLY : USS => SURF_ATM_SSO
!
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
!
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
!
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_BEM_n, ONLY : B => BEM
!
USE MODD_TEB_CANOPY_n, ONLY : TCP => TEB_CANOPY
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
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_GRID_n, ONLY : TG => TEB_GRID
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
 CALL PREP_OUTPUT_GRID(UG, U, &
                       ILUOUT,TG%CGRID,TG%XGRID_PAR,TG%XLAT,TG%XLON)
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading and horizontal interpolations
!
!
!*      2.0    Large scale orography
!
 CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'ZS     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!* option for roads
!
TOP%CROAD_DIR = CROAD_DIR
TOP%CWALL_OPT = CWALL_OPT
!
DO JPATCH=1,TOP%NTEB_PATCH
  !
  CALL GOTO_WRAPPER_TEB_PATCH(JPATCH)
  !*      2.1    Water reservoirs
  !
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'WS_ROOF',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'WS_ROAD',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  !
  !*      2.2    Building temperature
  !
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'TI_BLD ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  !
  !*      2.3    Road deep temperature
  !
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'TI_ROAD',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  !
  !*      2.4    Temperature profiles
  !
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'T_ROAD ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'T_WALLA',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'T_WALLB',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'T_ROOF ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  !
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'T_WIN1 ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  IF (TOP%CBEM == 'BEM') THEN
    CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'QI_BLD ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
    CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'T_WIN2 ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
    CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'T_FLOOR',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
    CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'T_MASS ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  ENDIF  
  !*      2.5    Snow variables
  !
  T%CUR%TSNOW_ROOF%SCHEME='1-L'
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'SN_ROOF',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  T%CUR%TSNOW_ROAD%SCHEME='1-L'
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'SN_ROAD',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  !
  !*      2.6    Canyon air variables
  !
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'T_CAN  ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_TEB_FIELD(B, BOP, DTCO, IOB, IG, U, TG, T, TOP, &
                         HPROGRAM,'Q_CAN  ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  !
  !-------------------------------------------------------------------------------------
  !
  !*      3.     Vertical interpolations of all variables
  !
  IF(LVERTSHIFT)THEN
    CALL PREP_VER_TEB(B, T, TOP)
  ENDIF
  !
  !-------------------------------------------------------------------------------------
  !
  !*      4.     Urban green areas
  !
  
  IF (TOP%LGARDEN)    CALL PREP_TEB_GARDEN(DTCO, IOB, IG, I, UG, U, USS, &
                                           TGD, TGDO, TGDPE, TGDP, TG, TOP, TVG, &
                                              HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  IF (TOP%LGREENROOF) CALL PREP_TEB_GREENROOF(DTCO, IOB, IG, I, UG, U, USS, TGR, TGRO, TGRPE, TGRP, &
                               TG, T, TOP, TVG, &
                                              HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
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
IF (TOP%LCANOPY) CALL PREP_TEB_CANOPY(TCP, TG)
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
