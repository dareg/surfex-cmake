!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###########################################################
      SUBROUTINE PGD_GRID_SURF_ATM (UG, KDIM_FULL, &
                                    HPROGRAM,HFILE,HFILETYPE,OGRID,HDIR)
!     ###########################################################
!!
!!    PURPOSE
!!    -------
!!   This program prepares the physiographic data fields.
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson                   Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     13/10/03
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURFEX_MPI, ONLY : NRANK
!
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
!
USE MODD_SURF_PAR,        ONLY : NVERSION, NBUGFIX
USE MODD_SURF_CONF,       ONLY : CPROGNAME
USE MODD_PGD_GRID,        ONLY : LLATLONMASK, NL, NGRID_PAR
!
USE MODI_PGD_GRID
USE MODI_INI_CSTS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_PGD_GRID_IO_INIT
USE MODI_SURF_VERSION
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
!
!
!
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
!
INTEGER, INTENT(INOUT) :: KDIM_FULL
 CHARACTER(LEN=6),     INTENT(IN)  :: HPROGRAM ! program calling
 CHARACTER(LEN=28),    INTENT(IN)  :: HFILE    ! atmospheric file name
 CHARACTER(LEN=6),     INTENT(IN)  :: HFILETYPE! atmospheric file type
LOGICAL,              INTENT(IN)  :: OGRID    ! .true. if grid is imposed by atm. model
 CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: HDIR
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
 CHARACTER(LEN=1) :: YDIR
 CHARACTER(LEN=100) :: YCOMMENT
INTEGER :: IRESP ! error return code
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PGD_GRID_SURF_ATM',0,ZHOOK_HANDLE)
!
YDIR = 'A'
IF (PRESENT(HDIR)) YDIR=HDIR
!
CPROGNAME=HPROGRAM
!
!*    1.      Set default constant values 
!             ---------------------------
!
 CALL SURF_VERSION
!
 CALL INI_CSTS
!
!-------------------------------------------------------------------------------
!
!*    2.      Initialisation of output grid
!             -----------------------------
!
 CALL PGD_GRID(UG, KDIM_FULL, HPROGRAM,HFILE,HFILETYPE,OGRID,YDIR)
! 
!
!-------------------------------------------------------------------------------
!
!
!*    3.      Additional actions for I/O
!
 CALL PGD_GRID_IO_INIT(HPROGRAM)
!
IF (LHOOK) CALL DR_HOOK('PGD_GRID_SURF_ATM',1,ZHOOK_HANDLE)
!_______________________________________________________________________________
!
END SUBROUTINE PGD_GRID_SURF_ATM
