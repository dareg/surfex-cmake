!     #########
      SUBROUTINE WRITESURF_PGD_TEB_GREENROOF_n (DGU, U, &
                                                 DTGR, TGRO, TGRP, &
                                                HPROGRAM)
!     ###############################################
!
!!****  *WRITESURF_PGD_TEB_GREENROOF_n* - writes ISBA fields describing urban greenroofs
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
!!     A. Lemonsu & C. de Munck   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2011 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
!
USE MODI_WRITE_SURF
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
!
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTGR
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(ISBA_PGD_t), INTENT(INOUT) :: TGRP
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_TEB_GREENROOF_N',0,ZHOOK_HANDLE)
!
!* soil scheme option
!
YRECFM='GR_ISBA'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,TGRO%CISBA,IRESP,HCOMMENT=YCOMMENT)
!
!* thermal conductivity option
!
YRECFM='GR_SCOND'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,TGRO%CSCOND,IRESP,HCOMMENT=YCOMMENT)
!
!* number of soil layers
!
YRECFM='GR_LAYER'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,TGRO%NGROUND_LAYER,IRESP,HCOMMENT=YCOMMENT)
!
!* number of time data for green roof chacteristics (VEG, LAI, EMIS, Z0) 
!
YRECFM='GR_NTIME'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,DTGR%NTIME,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GR_RUNOFFB'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,TGRP%XRUNOFFB,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GR_WDRAIN'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,TGRP%XWDRAIN,IRESP,HCOMMENT=YCOMMENT)
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_TEB_GREENROOF_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_PGD_TEB_GREENROOF_n
