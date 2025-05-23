!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###################################################################
      SUBROUTINE UPDATE_ESM_SEAFLUX_n (S, KI,KSW,PZENITH,PDIR_ALB,     &
                                      PSCA_ALB,PEMIS,PTSRAD,PTSURF )
!     ##############################################################
!
!!****  *UPDATE_ESM_SEAFLUX_n* - routine to update SEAFLUX radiative and physical properties in
!!                               Earth System Model after the call to OASIS coupler in order 
!!                               to close the energy budget between radiative scheme and surfex
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
!!     B. Decharme 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2009
!!      B. Decharme 06/2013 new coupling variables
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
!
USE MODD_CSTS,           ONLY : XTTS
!
USE MODI_UPDATE_RAD_SEA
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
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
!
INTEGER,                            INTENT(IN)  :: KI        ! number of points
INTEGER,                            INTENT(IN)  :: KSW       ! number of short-wave spectral bands
!
REAL,             DIMENSION(KI),    INTENT(IN)  :: PZENITH   ! solar zenithal angle
!
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(KI),    INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(KI),    INTENT(OUT) :: PTSRAD    ! radiative temperature
REAL,             DIMENSION(KI),    INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!*            Albedo and emissivity on open sea and sea ice
!             ---------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('UPDATE_ESM_SEAFLUX_N',0,ZHOOK_HANDLE)
!
CALL UPDATE_RAD_SEA(S,PZENITH,XTTS,PDIR_ALB,PSCA_ALB,PEMIS,PTSRAD   ) 
!
PTSURF(:) = S%XSST(:)
!
IF (LHOOK) CALL DR_HOOK('UPDATE_ESM_SEAFLUX_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE UPDATE_ESM_SEAFLUX_n
