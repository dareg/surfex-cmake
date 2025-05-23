!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE DEFAULT_DIAG_TEB (K2M,OSURF_BUDGET,O2M_MIN_ZS,ORAD_BUDGET,OCOEF,OSURF_VARS, &
                                   OSURF_MISC_BUDGET,OSURF_DIAG_ALBEDO,OUTCI,OPGD,PDIAG_TSTEP )  
!     #################################################################################################################
!
!!****  *DEFAULT_DIAG_TEB * - routine to set default values for the choice of diagnostics
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
!!      Original    01/2004 
!!      Modified by P. Le Moigne, 11/2004: add budget switch 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
INTEGER,  INTENT(OUT) :: K2M                ! flag for operational 2m quantities
LOGICAL,  INTENT(OUT) :: OSURF_BUDGET       ! flag for surface budget
LOGICAL,  INTENT(OUT) :: O2M_MIN_ZS
LOGICAL,  INTENT(OUT) :: ORAD_BUDGET        ! flag for radiative budget
LOGICAL,  INTENT(OUT) :: OCOEF
LOGICAL,  INTENT(OUT) :: OSURF_VARS
LOGICAL,  INTENT(OUT) :: OSURF_MISC_BUDGET  ! flag for surface miscellaneous budget
LOGICAL,  INTENT(OUT) :: OSURF_DIAG_ALBEDO  ! flag for albedo
LOGICAL,  INTENT(OUT) :: OUTCI              ! flag for UTCI fields
LOGICAL,  INTENT(OUT) :: OPGD               ! flag for PGD fields
REAL,     INTENT(OUT) :: PDIAG_TSTEP        ! time-step for writing
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DEFAULT_DIAG_TEB',0,ZHOOK_HANDLE)
K2M               = 0
OSURF_BUDGET      = .FALSE.
!
O2M_MIN_ZS        = .FALSE.
!
ORAD_BUDGET       = .FALSE.
!
OCOEF             = .FALSE.
OSURF_VARS        = .FALSE.
!
OSURF_MISC_BUDGET = .FALSE.
OSURF_DIAG_ALBEDO = .FALSE.
!
OUTCI             = .FALSE.
!
OPGD              = .FALSE.
!
PDIAG_TSTEP       = XUNDEF
IF (LHOOK) CALL DR_HOOK('DEFAULT_DIAG_TEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_DIAG_TEB
