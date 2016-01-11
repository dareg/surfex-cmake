!     #########
      SUBROUTINE WRITESURF_PGD_SEAF_PAR_n (DGU, U, &
                                            DTS, &
                                           HPROGRAM)
!     ################################################
!
!!****  *WRITESURF_PGD_SEAF_PAR_n* - writes SEAFLUX sst
!!                        
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
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2007 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_DATA_SEAFLUX_n, ONLY : DATA_SEAFLUX_t
!
USE MODD_TYPE_DATE_SURF
!
USE MODI_WRITE_SURF
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
!
TYPE(DIAG_t), INTENT(INOUT) :: DGU
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
TYPE(DATA_SEAFLUX_t), INTENT(INOUT) :: DTS
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
INTEGER           :: JTIME          ! loop index
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_SEAF_PAR_N',0,ZHOOK_HANDLE)
DTS%NTIME = SIZE(DTS%XDATA_SST,2)
YRECFM='ND_SEA_TIME'
YCOMMENT='(-)'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,DTS%NTIME,IRESP,HCOMMENT=YCOMMENT)
!
DO JTIME=1,DTS%NTIME
  WRITE(YRECFM,FMT='(A7,I3.3)') 'D_SST_T',JTIME
  YCOMMENT='X_Y_DATA_SST'
  CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,DTS%XDATA_SST(:,JTIME),IRESP,HCOMMENT=YCOMMENT)
END DO
!
YRECFM='TD_SST'
YCOMMENT='(-)'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,DTS%TDATA_SST,IRESP,HCOMMENT=YCOMMENT)
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_SEAF_PAR_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_PGD_SEAF_PAR_n
