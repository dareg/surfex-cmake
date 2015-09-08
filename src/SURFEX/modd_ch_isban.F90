!     #####################
      MODULE MODD_CH_ISBA_n
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

TYPE CH_ISBA_t
!
  CHARACTER(LEN=28)  :: CCHEM_SURF_FILE  ! name of general (chemical) purpose ASCII input file
  CHARACTER(LEN=6)                :: CCH_DRY_DEP            !  deposition scheme
  REAL, DIMENSION(:,:,:), POINTER :: XDEP                   ! final dry deposition  
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
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CCH_NAMES      ! NAME OF CHEMICAL SPECIES
                                                            ! (FOR DIAG ONLY)
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CAER_NAMES     ! NAME OF CHEMICAL SPECIES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CDSTNAMES      ! NAME OF CHEMICAL SPECIES
  CHARACTER(LEN=6), DIMENSION(:), POINTER :: CSLTNAMES      ! NAME OF CHEMICAL SPECIES                                                            
!
END TYPE CH_ISBA_t



CONTAINS

!





SUBROUTINE CH_ISBA_INIT(YCH_ISBA)
TYPE(CH_ISBA_t), INTENT(INOUT) :: YCH_ISBA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_CH_ISBA_N:CH_ISBA_INIT",0,ZHOOK_HANDLE)
NULLIFY(YCH_ISBA%XDEP)
NULLIFY(YCH_ISBA%XSOILRC_SO2)
NULLIFY(YCH_ISBA%XSOILRC_O3)
NULLIFY(YCH_ISBA%CSV)
NULLIFY(YCH_ISBA%CCH_NAMES)
NULLIFY(YCH_ISBA%CAER_NAMES)
NULLIFY(YCH_ISBA%CDSTNAMES)
NULLIFY(YCH_ISBA%CSLTNAMES)
YCH_ISBA%CCHEM_SURF_FILE=' '
YCH_ISBA%CCH_DRY_DEP=' '
YCH_ISBA%LCH_BIO_FLUX=.FALSE.
YCH_ISBA%LCH_NO_FLUX=.FALSE.
YCH_ISBA%NBEQ=0
YCH_ISBA%NDSTEQ=0
YCH_ISBA%NAEREQ=0
YCH_ISBA%NSLTEQ=0
YCH_ISBA%NSV_CHSBEG=0
YCH_ISBA%NSV_CHSEND=0
YCH_ISBA%NSV_DSTBEG=0
YCH_ISBA%NSV_DSTEND=0
YCH_ISBA%NSV_SLTBEG=0
YCH_ISBA%NSV_SLTEND=0
YCH_ISBA%NSV_AERBEG=0
YCH_ISBA%NSV_AEREND=0
IF (LHOOK) CALL DR_HOOK("MODD_CH_ISBA_N:CH_ISBA_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE CH_ISBA_INIT


END MODULE MODD_CH_ISBA_n
