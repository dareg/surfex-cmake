!     #########
      SUBROUTINE READ_TEB_CANOPY_n(HPROGRAM)
!     #########################################
!
!!****  *READ_TEB_CANOPY_n* - reads TEB fields
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
!!      Original    01/2003 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TEB_n,          ONLY : LCANOPY
USE MODD_TEB_CANOPY_n,   ONLY : NLVL, XZ, XU, XT, XQ, XTKE, XLMO, XLM,XLEPS, XDZ, XZF, XDZF, XP, XQ0
!
USE MODI_READ_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_CANOPY_GRID
!
USE MODI_GET_TYPE_DIM_n
!
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_TEB_CANOPY_N',0,ZHOOK_HANDLE)
YRECFM='SIZE_TOWN'
CALL GET_TYPE_DIM_n('TOWN  ',ILU)
!
!* flag to use or not canopy levels
!
YRECFM='VERSION'
CALL READ_SURF(HPROGRAM,YRECFM,IVERSION,IRESP)
!
IF (IVERSION<3) THEN
  LCANOPY = .FALSE.
ELSE
  YRECFM='TEB_CANOPY'
  CALL READ_SURF(HPROGRAM,YRECFM,LCANOPY,IRESP)
END IF
!
!
IF (.NOT. LCANOPY .AND. LHOOK) CALL DR_HOOK('READ_TEB_CANOPY_N',1,ZHOOK_HANDLE)
IF (.NOT. LCANOPY) RETURN
!
!* number of vertical levels
!
YRECFM='TEB_CAN_LVL'
CALL READ_SURF(HPROGRAM,YRECFM,NLVL,IRESP)
!
!*       2.     Prognostic fields:
!               -----------------
!
!* altitudes
!
ALLOCATE(XZ(ILU,NLVL))
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A5)') 'TEB_CAN_Z',JLAYER,'     '
  CALL READ_SURF(HPROGRAM,YRECFM,XZ(:,JLAYER),IRESP)
END DO
!
!* wind in canopy
!
ALLOCATE(XU(ILU,NLVL))
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A5)') 'TEB_CAN_U',JLAYER,'     '
  CALL READ_SURF(HPROGRAM,YRECFM,XU(:,JLAYER),IRESP)
END DO
!
!* theta in canopy
!
ALLOCATE(XT(ILU,NLVL))
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A5)') 'TEB_CAN_T',JLAYER,'     '
  CALL READ_SURF(HPROGRAM,YRECFM,XT(:,JLAYER),IRESP)
END DO
!
!* humidity in canopy
!
ALLOCATE(XQ(ILU,NLVL))
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A5)') 'TEB_CAN_Q',JLAYER,'     '
  CALL READ_SURF(HPROGRAM,YRECFM,XQ(:,JLAYER),IRESP)
END DO
!
!* Tke in canopy
!
ALLOCATE(XTKE(ILU,NLVL))
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A5)') 'TEB_CAN_E',JLAYER,'     '
  CALL READ_SURF(HPROGRAM,YRECFM,XTKE(:,JLAYER),IRESP)
END DO
!
!* Monin-Obhukov length
!
ALLOCATE(XLMO(ILU,NLVL))
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A10,I2.2,A4)') 'TEB_CAN_MO',JLAYER,'    '
  CALL READ_SURF(HPROGRAM,YRECFM,XLMO(:,JLAYER),IRESP)
ENDDO
!
!* mixing length
!
ALLOCATE(XLM(ILU,NLVL))
!
!* dissipative length
!
ALLOCATE(XLEPS(ILU,NLVL))
!
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
!
!
!* Pressure
!
ALLOCATE(XP(ILU,NLVL))
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A5)') 'TEB_CAN_P',JLAYER,'     '
  CALL READ_SURF(HPROGRAM,YRECFM,XP(:,JLAYER),IRESP)
END DO
!
!
!* Canyon heat flux
!
ALLOCATE(XQ0(ILU))
!
CALL READ_SURF(HPROGRAM,'TEB_CAN_Q0',XQ0(:),IRESP)
IF (LHOOK) CALL DR_HOOK('READ_TEB_CANOPY_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_TEB_CANOPY_n
