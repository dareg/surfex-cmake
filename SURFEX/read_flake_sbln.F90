!     #########
      SUBROUTINE READ_FLAKE_SBL_n(HPROGRAM)
!     #########################################
!
!!****  *READ_FLAKE_SBL_n* - reads FLAKE fields
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
USE MODD_FLAKE_n,       ONLY : LSBL
USE MODD_FLAKE_SBL_n,   ONLY : NLVL, XZ, XU, XT, XQ, XTKE, XLMO, XDZ, XZF, XDZF, XP
!
USE MODI_READ_SURF
!
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_FLAKE_SBL_N',0,ZHOOK_HANDLE)
CALL GET_TYPE_DIM_n('WATER ',ILU)
!
!
!* flag to use or not SBL levels
!
YRECFM='WAT_SBL'
CALL READ_SURF(HPROGRAM,YRECFM,LSBL,IRESP)
!
IF (.NOT. LSBL .AND. LHOOK) CALL DR_HOOK('READ_FLAKE_SBL_N',1,ZHOOK_HANDLE)
IF (.NOT. LSBL) RETURN
!
!* number of vertical levels
!
YRECFM='WAT_SBL_LVL'
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
  WRITE(YRECFM,'(A9,I2.2,A4)') 'WAT_SBL_Z',JLAYER,'    '
  CALL READ_SURF(HPROGRAM,YRECFM,XZ(:,JLAYER),IRESP)
END DO
!
!* wind in SBL
!
ALLOCATE(XU(ILU,NLVL))
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A4)') 'WAT_SBL_U',JLAYER,'    '
  CALL READ_SURF(HPROGRAM,YRECFM,XU(:,JLAYER),IRESP)
END DO
!
!* theta in SBL
!
ALLOCATE(XT(ILU,NLVL))
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A4)') 'WAT_SBL_T',JLAYER,'    '
  CALL READ_SURF(HPROGRAM,YRECFM,XT(:,JLAYER),IRESP)
END DO
!
!* humidity in SBL
!
ALLOCATE(XQ(ILU,NLVL))
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A4)') 'WAT_SBL_Q',JLAYER,'    '
  CALL READ_SURF(HPROGRAM,YRECFM,XQ(:,JLAYER),IRESP)
END DO
!
!* Tke in SBL
!
ALLOCATE(XTKE(ILU,NLVL))
!
DO JLAYER=1,NLVL
  WRITE(YRECFM,'(A9,I2.2,A4)') 'WAT_SBL_E',JLAYER,'    '
  CALL READ_SURF(HPROGRAM,YRECFM,XTKE(:,JLAYER),IRESP)
END DO
!
!* Monin-Obhukov length
!
ALLOCATE(XLMO(ILU))
!
YRECFM='WAT_SBL_LMO     '
CALL READ_SURF(HPROGRAM,YRECFM,XLMO(:),IRESP)
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
  WRITE(YRECFM,'(A9,I2.2,A4)') 'WAT_SBL_P',JLAYER,'    '
  CALL READ_SURF(HPROGRAM,YRECFM,XP(:,JLAYER),IRESP)
END DO
IF (LHOOK) CALL DR_HOOK('READ_FLAKE_SBL_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_FLAKE_SBL_n
