!     ########################################
      SUBROUTINE PUT_ZS_n(HPROGRAM,KI,PZS)
!     ########################################
!
!!****  *PUT_ZS_n* - routine to modify surface oropgraphy of each tile using atmospheric
!                    model orography
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
!!      Original     01/2004
!!      P. Le Moigne 05/2007: write model orography over each tile
!-------------------------------------------------------------------------------
!      
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_PUT_ZS_INLAND_WATER_n
!
USE MODI_PUT_ZS_NATURE_n
!
USE MODI_PUT_ZS_SEA_n
!
USE MODI_PUT_ZS_SURF_ATM_n
!
USE MODI_PUT_ZS_TOWN_n
USE MODI_GET_SIZE_FULL_n
USE MODI_GET_1D_MASK
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM
INTEGER,             INTENT(IN)  :: KI      ! horizontal dim. of cover
REAL, DIMENSION(KI), INTENT(IN)  :: PZS     ! orography
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PUT_ZS_N',0,ZHOOK_HANDLE)
!
!*       1. Full surface
!           ------------
!
 CALL PUT_ZS_SURF_ATM_n(HPROGRAM,KI,PZS)
!
!*       2. inland water
!           ------------
!
IF (U%NSIZE_WATER > 0 .AND. U%CWATER/='NONE' .AND. U%CWATER/='FLUX') CALL PACK_ZS(U%NSIZE_WATER,U%NR_WATER,'W')
!
!*       3. nature
!           ------
!
IF (U%NSIZE_NATURE > 0 .AND. U%CNATURE/='NONE' .AND. U%CNATURE/='FLUX') CALL PACK_ZS(U%NSIZE_NATURE,U%NR_NATURE,'N')
!
!*       4. town
!           ----
!
IF (U%NSIZE_TOWN > 0 .AND. U%CTOWN/='NONE' .AND. U%CTOWN/='FLUX') CALL PACK_ZS(U%NSIZE_TOWN,U%NR_TOWN,'T')
!
!        5.sea
!           ----
!
IF (U%NSIZE_SEA > 0 .AND. U%CSEA/='NONE' .AND. U%CSEA/='FLUX') CALL PACK_ZS(U%NSIZE_SEA,U%NR_SEA,'S')
!
IF (LHOOK) CALL DR_HOOK('PUT_ZS_N',1,ZHOOK_HANDLE)
!
CONTAINS
!=======================================================================================
SUBROUTINE PACK_ZS(KSIZE,KMASK,YTYPE)
!
INTEGER, INTENT(IN)               :: KSIZE
INTEGER, POINTER, DIMENSION(:)    :: KMASK
 CHARACTER(LEN=1), INTENT(IN)      :: YTYPE
!
REAL, DIMENSION(KSIZE) :: ZP_ZS
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! input arguments:
!
IF (LHOOK) CALL DR_HOOK('PUT_ZS_N:PACK_ZS',0,ZHOOK_HANDLE)
!
IF (.NOT.ASSOCIATED(KMASK)) THEN
  ALLOCATE(KMASK (KSIZE))
  IF (KSIZE>0) THEN
    CALL GET_SIZE_FULL_n(HPROGRAM,U%NDIM_FULL,U%NSIZE_FULL)
    IF (YTYPE=='W') THEN
      CALL GET_1D_MASK( KSIZE, U%NSIZE_FULL, U%XWATER, KMASK)
    ELSEIF (YTYPE=='N') THEN
      CALL GET_1D_MASK( KSIZE, U%NSIZE_FULL, U%XNATURE, KMASK)
    ELSEIF (YTYPE=='T') THEN
      CALL GET_1D_MASK( KSIZE, U%NSIZE_FULL, U%XTOWN, KMASK)
    ELSEIF (YTYPE=='S') THEN
      CALL GET_1D_MASK( KSIZE, U%NSIZE_FULL, U%XSEA, KMASK)
    ENDIF     
  ENDIF
ENDIF
!
DO JJ=1,KSIZE
  ZP_ZS(JJ)         = PZS         (KMASK(JJ))
ENDDO
!
IF (YTYPE=='W') THEN
  CALL PUT_ZS_INLAND_WATER_n(HPROGRAM,KSIZE,ZP_ZS,U%CWATER)
ELSEIF (YTYPE=='N') THEN
  CALL PUT_ZS_NATURE_n(HPROGRAM,KSIZE,ZP_ZS)
ELSEIF (YTYPE=='T') THEN
  CALL PUT_ZS_TOWN_n(HPROGRAM,KSIZE,ZP_ZS)
ELSEIF (YTYPE=='S') THEN
  CALL PUT_ZS_SEA_n(HPROGRAM,KSIZE,ZP_ZS)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('PUT_ZS_N:PACK_ZS',1,ZHOOK_HANDLE)
!
END SUBROUTINE PACK_ZS
!=======================================================================================
!
END SUBROUTINE PUT_ZS_n
