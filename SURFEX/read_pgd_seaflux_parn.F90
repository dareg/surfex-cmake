!     #########
      SUBROUTINE READ_PGD_SEAFLUX_PAR_n(HPROGRAM)
!     ################################################
!
!!****  *READ_PGD_SEAFLUX_PAR_n* - reads SEAFLUX sst
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
!!	P. Le Moigne   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2007 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SEAFLUX_GRID_n,    ONLY : NDIM
USE MODD_TYPE_DATE_SURF
USE MODD_DATA_SEAFLUX_n,    ONLY : NTIME, XDATA_SST, TDATA_SST
!
USE MODI_READ_SURF
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
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
INTEGER           :: JTIME          ! loop index
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('READ_PGD_SEAFLUX_PAR_N',0,ZHOOK_HANDLE)
YRECFM='ND_SEA_TIME'
CALL READ_SURF(HPROGRAM,YRECFM,NTIME,IRESP,HCOMMENT=YCOMMENT)
!
ALLOCATE(XDATA_SST       (NDIM,NTIME))
!
DO JTIME=1,NTIME
  WRITE(YRECFM,FMT='(A7,I3.3)') 'D_SST_T',JTIME
  CALL READ_SURF(HPROGRAM,YRECFM,XDATA_SST(:,JTIME),IRESP,HCOMMENT=YCOMMENT)
END DO
!
ALLOCATE(TDATA_SST       (NTIME))
!
YRECFM='TDATA_SST'
YCOMMENT='(-)'
CALL READ_SURF(HPROGRAM,YRECFM,TDATA_SST,IRESP,HCOMMENT=YCOMMENT)
IF (LHOOK) CALL DR_HOOK('READ_PGD_SEAFLUX_PAR_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_PGD_SEAFLUX_PAR_n
