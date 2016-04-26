!     #########
      SUBROUTINE WRITESURF_SEAFLUX_SBL_n (HSELECT, OSBL, SB, HPROGRAM,HWRITE)
!     ####################################
!
!!****  *WRITE_SEAFLUX_n* - writes SEAFLUX fields
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
USE MODD_CANOPY_n, ONLY : CANOPY_t
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
 CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSELECT
!
LOGICAL, INTENT(IN) :: OSBL
TYPE(CANOPY_t), INTENT(INOUT) :: SB
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
 CHARACTER(LEN=3),  INTENT(IN)  :: HWRITE   ! 'PREP' : does not write SBL XUNDEF fields
!                                          ! 'ALL' : all fields are written
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
!* flag to define if SBL is computed
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_SEAFLUX_SBL_N',0,ZHOOK_HANDLE)
YRECFM='SEA_SBL'
YCOMMENT='flag to use SBL levels'
 CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,OSBL,IRESP,HCOMMENT=YCOMMENT)
!
IF (.NOT. OSBL .AND. LHOOK) CALL DR_HOOK('WRITESURF_SEAFLUX_SBL_N',1,ZHOOK_HANDLE)
IF (.NOT. OSBL) RETURN
!
!* number of levels
!
YRECFM='SEA_SBL_LVL'
YCOMMENT='number of SBL levels'
 CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,SB%NLVL,IRESP,HCOMMENT=YCOMMENT)
!
!* altitudes
!
DO JLAYER=1,SB%NLVL
  WRITE(YRECFM,'(A9,I2.2,A1)') 'SEA_SBL_Z',JLAYER,' '
  YCOMMENT='altitudes of SBL levels (m)'
  CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,SB%XZ(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
IF (HWRITE/='PRE') THEN
  !
  !* wind in SBL
  !
  DO JLAYER=1,SB%NLVL
    WRITE(YRECFM,'(A9,I2.2,A1)') 'SEA_SBL_U',JLAYER,' '
    YCOMMENT='wind at SBL levels (m/s)'
    CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,SB%XU(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
  END DO
  !
  !* temperature in SBL
  !
  DO JLAYER=1,SB%NLVL
    WRITE(YRECFM,'(A9,I2.2,A1)') 'SEA_SBL_T',JLAYER,' '
    YCOMMENT='temperature at SBL levels (K)'
    CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,SB%XT(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
  END DO
  !
  !* humidity in SBL
  !
  DO JLAYER=1,SB%NLVL
    WRITE(YRECFM,'(A9,I2.2,A1)') 'SEA_SBL_Q',JLAYER,' '
    YCOMMENT='humidity at SBL levels (kg/m3)'
    CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,SB%XQ(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
  END DO
  !
  !* Tke in SBL
  !
  DO JLAYER=1,SB%NLVL
    WRITE(YRECFM,'(A9,I2.2,A1)') 'SEA_SBL_E',JLAYER,' '
    YCOMMENT='Tke at SBL levels (m2/s2)'
    CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,SB%XTKE(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
  END DO
  !
  !* Monin-Obhukov length
  !
  YRECFM='SEA_SBL_LMO '
  CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,SB%XLMO(:,SB%NLVL),IRESP,HCOMMENT=YCOMMENT)
  !
  !* Air pressure in SBL
  !
  DO JLAYER=1,SB%NLVL
    WRITE(YRECFM,'(A9,I2.2,A1)') 'SEA_SBL_P',JLAYER,' '
    YCOMMENT='Pressure at SBL levels (Pa)'
    CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,SB%XP(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
  END DO
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_SEAFLUX_SBL_N',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_SEAFLUX_SBL_n
