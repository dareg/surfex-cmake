!     #########################################
      SUBROUTINE READ_SSO_CANOPY_n(HPROGRAM)
!     #########################################
!
!!****  *READ_SSO_CANOPY_n* - reads SSO fields
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2010 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SSO_CANOPY_n,   ONLY : NLVL, XZ, XU, XTKE, XDZ, XZF, XDZF
!
USE MODI_READ_SURF
USE MODI_PREP_SSO_CANOPY
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_CANOPY_GRID
USE MODI_GET_TYPE_DIM_n
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! calling program
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: ILU          ! 1D physical dimension
!
INTEGER           :: IRESP          ! Error code after redding
!
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
!
INTEGER :: JLAYER  ! loop counter on layers
!
INTEGER           :: IVERSION   ! surface version
!
LOGICAL           :: GCANOPY    ! flag to test if SSO canopy fields are in the file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_SSO_CANOPY_N',0,ZHOOK_HANDLE)
CALL GET_TYPE_DIM_n('FULL  ',ILU)
!
!* flag to use or not canopy levels
!
YRECFM='VERSION'
CALL READ_SURF(HPROGRAM,YRECFM,IVERSION,IRESP)
!
IF (IVERSION<6) THEN
  GCANOPY = .FALSE.
ELSE
  YRECFM='SSO_CANOPY'
  CALL READ_SURF(HPROGRAM,YRECFM,GCANOPY,IRESP)
END IF
!
!*       2.     Allocation of Prognostic fields:
!               --------------------------------
!
!
!* number of vertical levels
!
IF (.NOT. GCANOPY) THEN
  CALL PREP_SSO_CANOPY(ILU)
ELSE

YRECFM='SSO_CAN_LVL'
CALL READ_SURF(HPROGRAM,YRECFM,NLVL,IRESP)
ALLOCATE(XZ(ILU,NLVL))
ALLOCATE(XU(ILU,NLVL))
ALLOCATE(XTKE(ILU,NLVL))
!
!
!*       3.     Reading of Prognostic fields:
!               -----------------------------
!
!* altitudes
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A4)') 'SSO_CAN_Z',JLAYER,'    '
  CALL READ_SURF(HPROGRAM,YRECFM,XZ(:,JLAYER),IRESP)
END DO
!
!* wind in canopy
!
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A4)') 'SSO_CAN_U',JLAYER,'    '
  CALL READ_SURF(HPROGRAM,YRECFM,XU(:,JLAYER),IRESP)
END DO
!
!* Tke in canopy
!
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A4)') 'SSO_CAN_E',JLAYER,'    '
  CALL READ_SURF(HPROGRAM,YRECFM,XTKE(:,JLAYER),IRESP)
END DO
!
END IF
!
!* Grid characteristics
!
!
!  --------------------------------- XZ(k+1)                     XDZ(k+1)
!                                                                           ^
!                                                                           |
!                                                                           |
!  - - - - - - - - - - - - - - - - - XZf(k+1)                               | XDZf(k+1)
!                                                              ^            |
!                                                              |            |
!  --------------------------------- XZ(k), XU, XT, XQ, XTKE   | XDZ(k)     V
!                                                              |            ^
!  - - - - - - - - - - - - - - - - - XZf(k)                    V            | XDZf(k)
!  --------------------------------- XZ(k-1)                     XDZ(k-1)   V
!  - - - - - - - - - - - - - - - - - XZf(k-1)
!
ALLOCATE(XDZ (ILU,NLVL))
ALLOCATE(XZF (ILU,NLVL))
ALLOCATE(XDZF(ILU,NLVL))
CALL CANOPY_GRID(ILU,NLVL,XZ,XZF,XDZ,XDZF)
IF (LHOOK) CALL DR_HOOK('READ_SSO_CANOPY_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_SSO_CANOPY_n
