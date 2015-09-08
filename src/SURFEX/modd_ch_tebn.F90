!     #####################
      MODULE MODD_CH_TEB_n
!     ######################
!
!!
!!    PURPOSE
!!    -------
!     
!   
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!
!!    AUTHOR
!!    ------
!!  P. Tulet   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!  16/07/03 (P. Tulet)  restructured for externalization
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE CH_TEB_t
!
  CHARACTER(LEN=28)  :: CCHEM_SURF_FILE  ! name of general (chemical) purpose ASCII input file
  CHARACTER(LEN=6)                :: CCH_DRY_DEP            !  deposition scheme
  REAL, DIMENSION(:,:),   POINTER :: XDEP                   ! final dry deposition  
                                                            ! velocity  for nature
  REAL, DIMENSION(:,:),   POINTER :: XSOILRC_SO2            ! for SO2
  REAL, DIMENSION(:,:),   POINTER :: XSOILRC_O3             ! for O3                                                            
  LOGICAL                         :: LCH_BIO_FLUX           ! flag for the calculation of
                                                            ! biogenic fluxes
  LOGICAL                         :: LCH_NO_FLUX            ! flag for the calculation of
                                                            ! biogenic NO fluxes
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSV            ! name of the scalar var.

  INTEGER    :: NSV_CHSBEG, NSV_CHSEND                      ! chemical begin and ending
                                                            ! index of the HSV/CSV array
  INTEGER    :: NBEQ                                        ! number of chemical species
                                                            ! in the surface scheme
  INTEGER    :: NSV_DSTBEG, NSV_DSTEND                      ! index of first and last dust
!                                                           ! related scalar variable
  INTEGER    :: NSV_SLTBEG, NSV_SLTEND                      ! index of first and last sea salt
!                                                           ! related scalar variable
  INTEGER    :: NSV_AERBEG, NSV_AEREND                      ! index of first and last dust
!                                                           ! related scalar variable
  INTEGER    :: NDSTEQ                                      ! number of dust related species
!                                                           ! in scalar variables list
  INTEGER    :: NAEREQ                                      ! number of aerosol species
  INTEGER    :: NSLTEQ                                      ! number of sea salt related species
!                                                           ! in scalar variables list
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES      ! NAME OF CHEMICAL SPECIES
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES     ! NAME OF CHEMICAL SPECIES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CDSTNAMES      ! NAME OF CHEMICAL SPECIES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSLTNAMES      ! NAME OF CHEMICAL SPECIES
!
END TYPE CH_TEB_t



CONTAINS

!





SUBROUTINE CH_TEB_INIT(YCH_TEB)
TYPE(CH_TEB_t), INTENT(INOUT) :: YCH_TEB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_TEB_N:CH_TEB_INIT",0,ZHOOK_HANDLE)
NULLIFY(YCH_TEB%XDEP)
NULLIFY(YCH_TEB%XSOILRC_SO2)
NULLIFY(YCH_TEB%XSOILRC_O3)
NULLIFY(YCH_TEB%CSV)
NULLIFY(YCH_TEB%CCH_NAMES)
NULLIFY(YCH_TEB%CAER_NAMES)
NULLIFY(YCH_TEB%CDSTNAMES)
NULLIFY(YCH_TEB%CSLTNAMES)
YCH_TEB%CCHEM_SURF_FILE=' '
YCH_TEB%CCH_DRY_DEP=' '
YCH_TEB%LCH_BIO_FLUX=.FALSE.
YCH_TEB%LCH_NO_FLUX=.FALSE.
YCH_TEB%NBEQ=0
YCH_TEB%NDSTEQ=0
YCH_TEB%NAEREQ=0
YCH_TEB%NSLTEQ=0
YCH_TEB%NSV_CHSBEG=0
YCH_TEB%NSV_CHSEND=0
YCH_TEB%NSV_DSTBEG=0
YCH_TEB%NSV_DSTEND=0
YCH_TEB%NSV_SLTBEG=0
YCH_TEB%NSV_SLTEND=0
YCH_TEB%NSV_AERBEG=0
YCH_TEB%NSV_AEREND=0
IF (LHOOK) CALL DR_HOOK("MODD_CH_TEB_N:CH_TEB_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE CH_TEB_INIT


END MODULE MODD_CH_TEB_n
