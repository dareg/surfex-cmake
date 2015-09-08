!     #################################################################################
SUBROUTINE PREP_SURF_ATM (YSC, &
                          HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!     #################################################################################
!
!!****  *PREP_SURF_ATM* - driver for surface fields preparation
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
!!      P. Le Moigne 10/2005, Phasage Arome
!!------------------------------------------------------------------
!

!
USE MODD_SURFEX_n, ONLY : SURFEX_t
!
USE MODI_PREP_NATURE
USE MODI_PREP_SEA
USE MODI_PREP_INLAND_WATER
USE MODI_PREP_TOWN
!
USE MODE_READ_GRIB
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_SURF_VERSION
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(SURFEX_t), INTENT(INOUT) :: YSC
!
 CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM    ! program calling surf. schemes
 CHARACTER(LEN=28), INTENT(IN) :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),  INTENT(IN) :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28), INTENT(IN) :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),  INTENT(IN) :: HPGDFILETYPE! type of the Atmospheric file
!
!*      0.2    declarations of local variables
 CHARACTER(LEN=28)               :: YATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6)                :: YATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28)               :: YPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6)                :: YPGDFILETYPE! type of the Atmospheric file
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
!-------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PREP_SURF_ATM',0,ZHOOK_HANDLE)
 CALL SURF_VERSION
!-------------------------------------------------------------------------------------
!
IF ( LEN_TRIM(HATMFILE)>0 ) THEN
  YATMFILE=HATMFILE
ELSE
  YATMFILE='                            '
ENDIF
!
IF ( LEN_TRIM(HPGDFILE)>0 ) THEN
  YPGDFILE=HPGDFILE
ELSE
  YPGDFILE='                            '
ENDIF
!
IF (  LEN_TRIM(HATMFILETYPE)>0 ) THEN
  YATMFILETYPE=HATMFILETYPE
ELSE
  YATMFILETYPE='      '
ENDIF
!
IF (  LEN_TRIM(HPGDFILETYPE)>0 ) THEN
  YPGDFILETYPE=HPGDFILETYPE
ELSE
  YPGDFILETYPE='      '
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! SEA Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
IF(YSC%U%NDIM_SEA>0) CALL PREP_SEA(YSC%IOB, &
   YSC%DTCO, YSC%DTS, YSC%O, YSC%OR, YSC%SG, YSC%S, YSC%SSB, YSC%UG, YSC%U, &
   HPROGRAM,YATMFILE,YATMFILETYPE,YPGDFILE,YPGDFILETYPE)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! INLAND WATER Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
IF(YSC%U%NDIM_WATER>0) CALL PREP_INLAND_WATER(YSC%DTCO, YSC%IOB, YSC%USS, &
   YSC%FG, YSC%F, YSC%FSB, YSC%UG, YSC%U, YSC%WG, YSC%W, YSC%WSB, &
   HPROGRAM,YATMFILE,YATMFILETYPE,YPGDFILE,YPGDFILETYPE)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! NATURAL SURFACE Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
IF(YSC%U%NDIM_NATURE>0) CALL PREP_NATURE(YSC%DTCO, YSC%IOB, YSC%ICP, YSC%IG, &
               YSC%I, YSC%UG, YSC%U, YSC%USS, &
   HPROGRAM,YATMFILE,YATMFILETYPE,YPGDFILE,YPGDFILETYPE)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! URBAN Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
IF(YSC%U%NDIM_TOWN>0) CALL PREP_TOWN(YSC%DGCT, YSC%DGMT, YSC%B, YSC%BOP, YSC%DTCO, YSC%IOB, YSC%IG, &
                      YSC%I, YSC%UG, YSC%U, YSC%USS, YSC%TCP, YSC%TGD, &
                      YSC%TGDO, YSC%TGDPE, YSC%TGDP, YSC%TGR, YSC%TGRO, YSC%TGRPE, &
                      YSC%TGRP, YSC%TG, YSC%T, YSC%TOP, YSC%TVG, &
                      HPROGRAM,YATMFILE,YATMFILETYPE,YPGDFILE,YPGDFILETYPE)
!
 CALL CLEAR_GRIB_INDEX
!
IF (LHOOK) CALL DR_HOOK('PREP_SURF_ATM',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_SURF_ATM
