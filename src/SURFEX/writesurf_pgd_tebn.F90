!     #########
      SUBROUTINE WRITESURF_PGD_TEB_n(HPROGRAM)
!     ###############################################
!
!!****  *WRITE_PGD_TEB_n* - writes TEB fields
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
!!      Original    01/2003 
!!      B. Decharme 07/2011 : delete argument HWRITE
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_BLD_DESCRIPTION_n, ONLY : BDD => BLD_DESC
USE MODD_DATA_BEM_n, ONLY : DTB => DATA_BEM
USE MODD_DATA_TEB_GREENROOF_n, ONLY : DTGR => DATA_TEB_GREENROOF
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
USE MODD_TEB_IRRIG_n, ONLY : TIR => TEB_IRRIG
!
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_DATA_TEB_GARDEN_n, ONLY : DTGD => DATA_TEB_GARDEN
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
!
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TGRP => TEB_GREENROOF_PGD
!
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_TEB_GRID_n, ONLY : TG => TEB_GRID
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
!
USE MODE_WRITE_SURF_COV, ONLY : WRITE_SURF_COV
!
USE MODI_WRITE_SURF
USE MODI_WRITE_GRID
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_WRITESURF_PGD_TEB_PAR_n
USE MODI_WRITESURF_PGD_TEB_VEG_n
USE MODI_WRITESURF_PGD_TEB_GREENROOF_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_TEB_N',0,ZHOOK_HANDLE)
!
!*       1.     Dimension initializations:
!               -------------------------
!
!
!* number of TEB patches
!
YRECFM='TEB_PATCH'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TOP%NTEB_PATCH,IRESP,HCOMMENT=YCOMMENT)
!
!
!* number of roof layers
!
YRECFM='ROOF_LAYER'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TOP%NROOF_LAYER,IRESP,HCOMMENT=YCOMMENT)
!
!* number of road layers
!
YRECFM='ROAD_LAYER'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TOP%NROAD_LAYER,IRESP,HCOMMENT=YCOMMENT)
!
!* number of wall layers
!
YRECFM='WALL_LAYER'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TOP%NWALL_LAYER,IRESP,HCOMMENT=YCOMMENT)
!
!* flag indicating if fields are computed from ecoclimap or not
!
YRECFM='ECOCLIMAP'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TOP%LECOCLIMAP,IRESP,HCOMMENT=YCOMMENT)
!
!
!* Type of Building Energy Model
!
YRECFM='BEM'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TOP%CBEM,IRESP,HCOMMENT=YCOMMENT) 
!
IF (TOP%CBEM=='BEM') THEN
  YRECFM='COOL_COIL'
  YCOMMENT=YRECFM
  CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,BOP%CCOOL_COIL,IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='HEAT_COIL'
  YCOMMENT=YRECFM
  CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,BOP%CHEAT_COIL,IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='AUTOSIZE'
  YCOMMENT=YRECFM
  CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,BOP%LAUTOSIZE,IRESP,HCOMMENT=YCOMMENT)
END IF
!
!* Type of averaging of buildings characteristics
!
YRECFM='BLD_ATYPE'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TOP%CBLD_ATYPE,IRESP,HCOMMENT=YCOMMENT)
!
!
!
!* number of floor layers
!
IF (TOP%CBEM=="BEM") THEN
  YRECFM='FLOOR_LAYER'
  YCOMMENT=YRECFM
  CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,BOP%NFLOOR_LAYER,IRESP,HCOMMENT=YCOMMENT)
ENDIF
!
!
!* Use of solar panels
!
YRECFM='SOLAR_PANEL'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TOP%LSOLAR_PANEL,IRESP,HCOMMENT=YCOMMENT)
!
!------------------------------------------------------------------------------
!
! * ISBA fields for urban green areas
! 
IF (TOP%LGARDEN) THEN
!
! * Greenroofs and hydrology (only activated if LGARDEN)
!
YRECFM='LGREENROOF'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TOP%LGREENROOF,IRESP,HCOMMENT=YCOMMENT) 
!
YRECFM='LURBAN_HYDRO'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TOP%LHYDRO,IRESP,HCOMMENT=YCOMMENT) 
!
! * General ISBA options for urban vegetation
!
! * Pedo-transfert function
!
YRECFM='GD_PEDOTF'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TVG%CPEDOTF,IRESP,HCOMMENT=YCOMMENT)
!
! * type of photosynthesis
!
YRECFM='GD_PHOTO'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TVG%CPHOTO,IRESP,HCOMMENT=YCOMMENT)
!
!* new radiative transfert
!
YRECFM='GD_TR_ML'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TVG%LTR_ML,IRESP,HCOMMENT=YCOMMENT)
!
! * ISBA fields specific to urban gardens
!
 CALL WRITESURF_PGD_TEB_VEG_n(DGU, IOB, U, &
                              DTGD, TGDO, TGDP, TVG, &
                              HPROGRAM)
!
! * ISBA fields specific to urban greenroofs
!
IF (TOP%LGREENROOF) CALL WRITESURF_PGD_TEB_GREENROOF_n(DGU, IOB, U, &
                                                       TGRO, TGRP, &
                                                       HPROGRAM)
!
ENDIF
!
!------------------------------------------------------------------------------
!
!*       2.     Physiographic data fields:
!               -------------------------
!
!* cover classes
!
YRECFM='COVER_LIST'
YCOMMENT='(LOGICAL LIST)'
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TOP%LCOVER(:),IRESP,HCOMMENT=YCOMMENT,HDIR='-')
!
YCOMMENT='COVER FIELDS'
 CALL WRITE_SURF_COV(DGU, IOB, U, &
                     HPROGRAM,'COVER',TOP%XCOVER(:,:),TOP%LCOVER,IRESP,HCOMMENT=YCOMMENT)
!
!* orography
!
YRECFM='ZS'
YCOMMENT='ZS'
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,TOP%XZS(:),IRESP,HCOMMENT=YCOMMENT)
!
!* latitude, longitude
!
 CALL WRITE_GRID(DGU, IOB, U, &
                 HPROGRAM,TG%CGRID,TG%XGRID_PAR,TG%XLAT,TG%XLON,TG%XMESH_SIZE,IRESP)
!
!-------------------------------------------------------------------------------
 CALL WRITESURF_PGD_TEB_PAR_n(BDD, DTB, DTGD, DTGR, DTT, DGU, IOB, U, TGDO, TGDP, TGRO, TIR, TOP, &
                              HPROGRAM)
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_TEB_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_PGD_TEB_n
