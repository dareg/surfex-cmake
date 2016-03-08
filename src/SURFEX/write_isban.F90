!     #########
      SUBROUTINE WRITE_ISBA_n (DTCO, HSELECT, U, CHI, DGI, ICP, I, DST, HPROGRAM,HWRITE,OLAND_USE)
!     ####################################
!
!!****  *WRITE_ISBA_n* - routine to write surface variables in their respective files
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
!!      B. Decharme 07/2011 : Suppress pgd output
!       B. Decharme 07/2011 : land_use key for writing semi-prognostic variables
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_DST_n, ONLY : DST_t
!
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_SURFEX_n, ONLY : ISBA_DIAG_t
USE MODD_CANOPY_n, ONLY : CANOPY_t
USE MODD_ISBA_n, ONLY : ISBA_t
!
USE MODD_WRITE_SURF_ATM, ONLY : LNOWRITE_CANOPY
USE MODI_INIT_IO_SURF_n
USE MODI_WRITESURF_ISBA_n
USE MODI_WRITESURF_ISBA_CONF_n
USE MODI_END_IO_SURF_n
USE MODI_WRITESURF_ISBA_CANOPY_n
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
 CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSELECT 
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(ISBA_DIAG_t), INTENT(INOUT) :: DGI
TYPE(CANOPY_t), INTENT(INOUT) :: ICP
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(DST_t), INTENT(INOUT) :: DST
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),    INTENT(IN)  :: HWRITE    ! 'PREP' : does not write SBL XUNDEF fields
!                                             ! 'ALL' : all fields are written
LOGICAL,             INTENT(IN)  :: OLAND_USE !
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_ISBA_N',0,ZHOOK_HANDLE)
CALL INIT_IO_SURF_n(DTCO, U, HPROGRAM,'NATURE','ISBA  ','WRITE')
!
!*       1.     Selection of surface scheme
!               ---------------------------
!        
 CALL WRITESURF_ISBA_CONF_n(CHI, DGI%DE, DGI%O, DGI%DM, I%O, HPROGRAM)
 CALL WRITESURF_ISBA_n(HSELECT, CHI, DST, I%O, I%IP%XPATCH, I%R, I%M%T%XLAI, &
                       I%M%X%XDG, I%I%TTIME, HPROGRAM,OLAND_USE)
!
IF ((.NOT.LNOWRITE_CANOPY).OR.SIZE(HSELECT)>0) THEN
  CALL WRITESURF_ISBA_CANOPY_n(HSELECT, ICP, I%O%LCANOPY, HPROGRAM,HWRITE)
ENDIF
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_ISBA_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_ISBA_n
