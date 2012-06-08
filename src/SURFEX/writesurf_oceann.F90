!     #########
      SUBROUTINE WRITESURF_OCEAN_n(HPROGRAM)
!     ########################################
!
!!****  *WRITE_OCEAN_n* - writes OCEAN fields
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
!!      Original    03/2007
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_OCEAN_n, ONLY : XSEAT, XSEAS, XSEAU, XSEAV, XSEAE, XSEABATH,&
                         XSEAHMO, LMERCATOR  
USE MODD_OCEAN_GRID_n
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
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=4 ) :: YLVL
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
CHARACTER(LEN=14) :: YFORM          ! Writing format
!
INTEGER :: JLEVEL ! loop counter on oceanic levels
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WRITESURF_OCEAN_N',0,ZHOOK_HANDLE)
!
!* flag to define if OCEAN is used
!
YRECFM='SEA_OCEAN'
YCOMMENT='flag to use OCEAN model'
CALL WRITE_SURF(HPROGRAM,YRECFM,LMERCATOR,IRESP,HCOMMENT=YCOMMENT)
!
IF (.NOT. LMERCATOR .AND. LHOOK) CALL DR_HOOK('WRITESURF_OCEAN_N',1,ZHOOK_HANDLE)
IF (.NOT. LMERCATOR) RETURN
!
!-------------------------------------------------------------------------------
!
!*       3.     Prognostic fields:
!               -----------------
!* oceanic temperature
!
DO JLEVEL=NOCKMIN+1,NOCKMAX
  WRITE(YLVL,'(I4)') JLEVEL
  YRECFM='TEMP_OC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A11,I1.1,A5)'
  IF (JLEVEL >= 10)  YFORM='(A11,I2.2,A5)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_TEMP_OC',JLEVEL,' (°C)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XSEAT(:,JLEVEL),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!* oceanic salinity
!
DO JLEVEL=NOCKMIN+1,NOCKMAX
  WRITE(YLVL,'(I4)') JLEVEL
  YRECFM='SALT_OC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A11,I1.1,A5)'
  IF (JLEVEL >= 10)  YFORM='(A11,I2.2,A5)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_SALT_OC',JLEVEL,'(psu)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XSEAS(:,JLEVEL),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!* oceanic zonal current
!
DO JLEVEL=NOCKMIN+1,NOCKMAX
  WRITE(YLVL,'(I4)') JLEVEL
  YRECFM='UCUR_OC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A11,I1.1,A5)'
  IF (JLEVEL >= 10)  YFORM='(A11,I2.2,A5)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_UCUR_OC',JLEVEL,' (m/s)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XSEAU(:,JLEVEL),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!* oceanic meridian current
!
DO JLEVEL=NOCKMIN+1,NOCKMAX
  WRITE(YLVL,'(I4)') JLEVEL
  YRECFM='VCUR_OC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A11,I1.1,A5)'
  IF (JLEVEL >= 10)  YFORM='(A11,I2.2,A5)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_VCUR_OC',JLEVEL,'(m/s)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XSEAV(:,JLEVEL),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!* oceanic turbulent kinetic energy
!
DO JLEVEL=NOCKMIN+1,NOCKMAX
  WRITE(YLVL,'(I4)') JLEVEL
  YRECFM='TKE_OC'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A11,I1.1,A5)'
  IF (JLEVEL >= 10)  YFORM='(A11,I2.2,A5)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_TKE_OC',JLEVEL,' (J)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XSEAE(:,JLEVEL),IRESP,HCOMMENT=YCOMMENT)
END DO
!-------------------------------------------------------------------------------
!
!*       4.     Semi-prognostic fields:
!               ----------------------
!* bathymetry indice
!
DO JLEVEL=NOCKMIN+1,NOCKMAX
  WRITE(YLVL,'(I4)') JLEVEL
  YRECFM='SEAINDBATH'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A11,I1.1,A5)'
  IF (JLEVEL >= 10)  YFORM='(A11,I2.2,A5)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_SEAINDBATH',JLEVEL,' (J)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XSEABATH(:,JLEVEL),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!-------------------------------------------------------------------------------
!Sea Surface Salinity
!
YRECFM='SSS'
YCOMMENT='SSS'
CALL WRITE_SURF(HPROGRAM,YRECFM,XSEAS(:,NOCKMIN),IRESP,HCOMMENT=YCOMMENT)
!-------------------------------------------------------------------------------
!Sea Surface Salinity
!
YRECFM='SEA_HMO'
YCOMMENT='X_Y_'
CALL WRITE_SURF(HPROGRAM,YRECFM,XSEAHMO(:),IRESP,HCOMMENT=YCOMMENT)
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WRITESURF_OCEAN_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITESURF_OCEAN_n
