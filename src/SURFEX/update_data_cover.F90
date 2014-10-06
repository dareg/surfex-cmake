!     #########
      SUBROUTINE UPDATE_DATA_COVER(KYEAR)
!     #########################
!
!!**** *INI_DATA_COVER* initializes cover-field correspondance arrays
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    09/2008
!!    P. Samuelsson 10/2014 MEB
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_ISBA_n,         ONLY : NPATCH, LMEB_PATCH
USE MODD_DATA_COVER_n,   ONLY : XDATA_VEGTYPE, NYEAR, XDATA_NATURE, XDATA_GARDEN
USE MODD_DATA_COVER,     ONLY :   XDATA_LAI, XDATA_H_TREE, &
                                  XDATA_VEG, XDATA_GREEN, XDATA_Z0, XDATA_EMIS_ECO, &
                                  XDATA_LAIGV, XDATA_Z0GV, XDATA_H_VEG, XDATA_LAIMIN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ECOCLIMAP2_LAI
!
USE MODI_INI_DATA_PARAM
USE MODI_FIX_MEB_VEG
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
INTEGER,             INTENT(IN)    :: KYEAR        ! new year
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER           :: ISIZE_LMEB_PATCH  ! Number of patches with MEB=true
!
!*    0.3    Declaration of namelists
!            ------------------------
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('UPDATE_DATA_COVER',0,ZHOOK_HANDLE)
IF (KYEAR /= NYEAR) THEN
  NYEAR = KYEAR
  CALL ECOCLIMAP2_LAI
  CALL INI_DATA_PARAM(XDATA_VEGTYPE, PSURF=XDATA_NATURE, PSURF2=XDATA_GARDEN, &
             PLAI=XDATA_LAI, PH_TREE=XDATA_H_TREE, PVEG_OUT=XDATA_VEG,        &
             PGREEN=XDATA_GREEN, PZ0=XDATA_Z0, PEMIS_ECO=XDATA_EMIS_ECO,      &
             PLAIMIN_OUT=XDATA_LAIMIN,                                        &
             PLAIGV_OUT=XDATA_LAIGV, PZ0GV=XDATA_Z0GV, PH_VEG=XDATA_H_VEG     )
!
  ISIZE_LMEB_PATCH=COUNT(LMEB_PATCH(:))
  IF (ISIZE_LMEB_PATCH>0)  THEN
    CALL FIX_MEB_VEG(NPATCH)
  ENDIF
!
END IF
IF (LHOOK) CALL DR_HOOK('UPDATE_DATA_COVER',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------

END SUBROUTINE UPDATE_DATA_COVER
