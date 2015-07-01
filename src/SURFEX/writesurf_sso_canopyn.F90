!     ####################################
      SUBROUTINE WRITESURF_SSO_CANOPY_n (DGU, IOB, U, &
                                          SSCP, &
                                         HPROGRAM,OWRITE)
!     ####################################
!
!!****  *WRITE_SSO_n* - writes SSO fields
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
!!      E. Martin   01/2012 avoid write of XUNDEF fields
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
!
!
!
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_SSO_CANOPY_n, ONLY : SSO_CANOPY_t
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
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
TYPE(SSO_CANOPY_t), INTENT(INOUT) :: SSCP
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
LOGICAL,           INTENT(IN)  :: OWRITE   ! flag to write canopy terms
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
!
INTEGER :: JLAYER  ! loop counter on layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       1.     Prognostic fields:
!               -----------------
!
!
!* flag to define if canopy is computed
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_SSO_CANOPY_N',0,ZHOOK_HANDLE)
YRECFM='SSO_CANOPY'
YCOMMENT='flag to use canopy levels'
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,OWRITE,IRESP,HCOMMENT=YCOMMENT)
!
IF (.NOT. OWRITE .AND. LHOOK) CALL DR_HOOK('WRITESURF_SSO_CANOPY_N',1,ZHOOK_HANDLE)
IF (.NOT. OWRITE) RETURN
!
!* number of levels
!
YRECFM='SSO_CAN_LVL'
YCOMMENT='number of canopy levels'
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,SSCP%NLVL,IRESP,HCOMMENT=YCOMMENT)
!
!* altitudes
!
DO JLAYER=1,SSCP%NLVL
  WRITE(YRECFM,'(A9,I2.2,A1)') 'SSO_CAN_Z',JLAYER,' '
  YCOMMENT='altitudes of canopy levels (m)'
  CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,SSCP%XZ(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!* wind in canopy
!
DO JLAYER=1,SSCP%NLVL
  WRITE(YRECFM,'(A9,I2.2,A1)') 'SSO_CAN_U',JLAYER,' '
  YCOMMENT='wind at canopy levels (m/s)'
  CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,SSCP%XU(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!* Tke in canopy
!
DO JLAYER=1,SSCP%NLVL
  WRITE(YRECFM,'(A9,I2.2,A1)') 'SSO_CAN_E',JLAYER,' '
  YCOMMENT='Tke at canopy levels (m2/s2)'
  CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,SSCP%XTKE(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_SSO_CANOPY_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_SSO_CANOPY_n
