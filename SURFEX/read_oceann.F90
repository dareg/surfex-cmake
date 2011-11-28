!     #########
      SUBROUTINE READ_OCEAN_n(HPROGRAM)
!     #########################################
!
!!****  *READ_OCEAN_n* - read oceanic variables
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
!!	C. Lebeaupin Brossier   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2007 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_OCEAN_n
USE MODD_DIAG_OCEAN_n
USE MODD_OCEAN_GRID_n, ONLY : NOCKMIN,NOCKMAX
!
!
USE MODI_READ_SURF
USE MODI_OCEAN_MERCATORVERGRID
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
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
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: ILU          ! 1D physical dimension
!
INTEGER           :: IRESP          ! Error code after redding
!
CHARACTER(LEN=4)  :: YLVL
!
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=14) :: YFORM          ! Writing format
REAL, DIMENSION(:),ALLOCATABLE  :: ZWORK      ! 1D array to write data in file
!
INTEGER :: JLEVEL ! loop counter on oceanic levels
INTEGER :: J      ! loop counter on sea grid points
!
INTEGER           :: IVERSION       ! surface version
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_OCEAN_N',0,ZHOOK_HANDLE)
NOCTCOUNT=0
!
YRECFM='VERSION'
CALL READ_SURF(HPROGRAM,YRECFM,IVERSION,IRESP)
!
!* flag to use or not Ocean model
!
IF (IVERSION<=3) THEN
   LMERCATOR=.FALSE.
ELSE
   YRECFM='SEA_OCEAN'
   CALL READ_SURF(HPROGRAM,YRECFM,LMERCATOR,IRESP)
ENDIF
!
IF (.NOT. LMERCATOR .AND. LHOOK) CALL DR_HOOK('READ_OCEAN_N',1,ZHOOK_HANDLE)
IF (.NOT. LMERCATOR) RETURN
!
!-------------------------------------------------------------------------------
CALL OCEAN_MERCATORVERGRID
!
!* 1D physical dimension
!
YRECFM='SIZE_SEA'
CALL GET_TYPE_DIM_n('SEA   ',ILU)
!
!*       2.     Prognostic fields:
!               -----------------
!
ALLOCATE(ZWORK(ILU))
!* oceanic temperature
!
ALLOCATE(XSEAT(ILU,NOCKMIN:NOCKMAX))
!
DO JLEVEL=NOCKMIN+1,NOCKMAX
  WRITE(YLVL,'(I4)') JLEVEL
  YRECFM='TEMP_OC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:),IRESP)
  XSEAT(:,JLEVEL)=ZWORK(:)
END DO
XSEAT(:,NOCKMIN)=XSEAT(:,NOCKMIN+1)
!
!* oceanic salinity
!
ALLOCATE(XSEAS(ILU,NOCKMIN:NOCKMAX))
!
DO JLEVEL=NOCKMIN+1,NOCKMAX
  WRITE(YLVL,'(I4)') JLEVEL
  YRECFM='SALT_OC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:),IRESP)
  XSEAS(:,JLEVEL)=ZWORK(:)
END DO
XSEAS(:,NOCKMIN)=XSEAS(:,NOCKMIN+1)
!
!* oceanic current
!
ALLOCATE(XSEAU(ILU,NOCKMIN:NOCKMAX))
ALLOCATE(XSEAV(ILU,NOCKMIN:NOCKMAX))
!
DO JLEVEL=NOCKMIN+1,NOCKMAX
  WRITE(YLVL,'(I4)') JLEVEL
  YRECFM='UCUR_OC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:),IRESP)
  XSEAU(:,JLEVEL)=ZWORK(:)
END DO
DO JLEVEL=NOCKMIN+1,NOCKMAX
  WRITE(YLVL,'(I4)') JLEVEL
  YRECFM='VCUR_OC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:),IRESP)
  XSEAV(:,JLEVEL)=ZWORK(:)
END DO
XSEAU(:,NOCKMIN)=XSEAU(:,NOCKMIN+1)
XSEAV(:,NOCKMIN)=XSEAV(:,NOCKMIN+1)
!
!* oceanic turbulent kinetic energy
!
ALLOCATE(XSEAE(ILU,NOCKMIN:NOCKMAX))
!
DO JLEVEL=NOCKMIN+1,NOCKMAX
  WRITE(YLVL,'(I4)') JLEVEL
  YRECFM='TKE_OC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:),IRESP)
  XSEAE(:,JLEVEL)=ZWORK(:)
END DO
XSEAE(:,NOCKMIN)=XSEAE(:,NOCKMIN+1)
!
!
!-------------------------------------------------------------------------------
!
!*       4.     Semi-prognostic fields:
!               ----------------------
!
!* bathymetry indice
!
ALLOCATE(XSEABATH(ILU,NOCKMIN:NOCKMAX))
!
DO JLEVEL=NOCKMIN+1,NOCKMAX
  WRITE(YLVL,'(I4)') JLEVEL
  YRECFM='SEAINDBATH'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:),IRESP)
  XSEABATH(:,JLEVEL)=ZWORK(:)
END DO
XSEABATH(:,NOCKMIN)=1.
!
!-------------------------------------------------------------------------------
!Complete undefined values for the oceanic 1D model convergence
DO J=1,ILU
  DO JLEVEL=NOCKMIN+2,NOCKMAX
    IF (XSEABATH(J,JLEVEL)==0.) THEN
      XSEAT(J,JLEVEL)=XSEAT(J,JLEVEL-1)
      XSEAS(J,JLEVEL)=XSEAS(J,JLEVEL-1)
      XSEAU(J,JLEVEL)=XSEAU(J,JLEVEL-1)
      XSEAV(J,JLEVEL)=XSEAV(J,JLEVEL-1)
      XSEAE(J,JLEVEL)=XSEAE(J,JLEVEL-1)
    ENDIF
  ENDDO
ENDDO
!
DEALLOCATE(ZWORK)
!-------------------------------------------------------------------------------
ALLOCATE(XSEAHMO(ILU))
YRECFM='SEA_HMO'
CALL READ_SURF(HPROGRAM,YRECFM,XSEAHMO(:),IRESP)
!
ALLOCATE(XHMOKMEL(ILU))
YRECFM='SEA_KMEL'
CALL READ_SURF(HPROGRAM,YRECFM,XHMOKMEL(:),IRESP)
!
!-------------------------------------------------------------------------------
ALLOCATE(XLE        (SIZE(XSEAT,1),NOCKMIN:NOCKMAX))
ALLOCATE(XLK        (SIZE(XSEAT,1),NOCKMIN:NOCKMAX))
ALLOCATE(XKMEL      (SIZE(XSEAT,1),NOCKMIN:NOCKMAX))
ALLOCATE(XKMELM     (SIZE(XSEAT,1),NOCKMIN:NOCKMAX))
XLE(:,:)    =XUNDEF
XLK(:,:)    =XUNDEF
XKMEL(:,:)  =XUNDEF
XKMELM(:,:) =XUNDEF
!
ALLOCATE(XSEATEND   (SIZE(XSEAT,1)))
XSEATEND(:) =XUNDEF
IF (LHOOK) CALL DR_HOOK('READ_OCEAN_N',1,ZHOOK_HANDLE)
!
!------------------------------------------------------------------------------
END SUBROUTINE READ_OCEAN_n
