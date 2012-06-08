!     ##################
      MODULE MODN_PGD_SCHEMES
!     ##################
!
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!       
!!    AUTHOR
!!    ------
!!	V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2003                    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
CHARACTER(LEN=6) :: CNATURE  ! type of scheme for natural covers: NONE, ISBA
CHARACTER(LEN=6) :: CSEA     ! type of scheme for oceans        : NONE, CSEAFLX
CHARACTER(LEN=6) :: CTOWN    ! type of scheme for towns         : NONE, CTEB
CHARACTER(LEN=6) :: CWATER   ! type of scheme for inland waters : NONE, WATFLX, FLAKE
LOGICAL          :: LGARDEN  ! T : green urban areas treated inside TOWN by CTOWN scheme
!                            ! F : green urban areas treated inside NATURE aggregated
!                            !     with other vegetated surfaces
!
!
NAMELIST/NAM_PGD_SCHEMES/CNATURE, CSEA, CTOWN, CWATER, LGARDEN
!
END MODULE MODN_PGD_SCHEMES
