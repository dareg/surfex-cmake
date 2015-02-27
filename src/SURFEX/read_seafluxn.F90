!     #########
      SUBROUTINE READ_SEAFLUX_n(HPROGRAM,KLUOUT)
!     #########################################
!
!!****  *READ_SEAFLUX_n* - read SEAFLUX varaibles
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003 
!!      Modified    02/2008 Add oceanic variables initialisation
!!      S. Belamari 04/2014 Suppress LMERCATOR
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_SEAFLUX_n,      ONLY : XSST, XZ0, LINTERPOL_SST, CINTERPOL_SST, &
                                XPERTFLUX, LPERTFLUX,                    &
                                XSST_MTH, TTIME, LHANDLE_SIC,            &
                                XSSS, XSSS_MTH, LINTERPOL_SSS,           &
                                CINTERPOL_SSS
!
USE MODI_READ_SURF
USE MODI_INTERPOL_SST_MTH
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_TYPE_DIM_n
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! calling program
INTEGER,           INTENT(IN)  :: KLUOUT
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: JMTH, INMTH
CHARACTER(LEN=2 ) :: YMTH
!
INTEGER           :: ILU          ! 1D physical dimension
!
INTEGER           :: IRESP          ! Error code after redding
!
CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
!
INTEGER           :: IVERSION       ! surface version
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_SEAFLUX_N',0,ZHOOK_HANDLE)
!
YRECFM='SIZE_SEA'
CALL GET_TYPE_DIM_n('SEA   ',ILU)
!
!*       2.     Prognostic fields:
!               -----------------
!
!* water temperature
!
ALLOCATE(XSST(ILU))
!
IF(LINTERPOL_SST)THEN
!
! Precedent, Current and Next Monthly/annual SST
  INMTH=3
  IF(TRIM(CINTERPOL_SST)=='ANNUAL')INMTH=14
!
  ALLOCATE(XSST_MTH(SIZE(XSST),INMTH))
  DO JMTH=1,INMTH
     WRITE(YMTH,'(I2)') (JMTH-1)
     YRECFM='SST_MTH'//ADJUSTL(YMTH(:LEN_TRIM(YMTH)))
     CALL READ_SURF(HPROGRAM,YRECFM,XSST_MTH(:,JMTH),IRESP)
  ENDDO
!
  CALL INTERPOL_SST_MTH(TTIME%TDATE%YEAR,TTIME%TDATE%MONTH,TTIME%TDATE%DAY,'T',XSST)
!
ELSE
! 
  ALLOCATE(XSST_MTH(0,0))
!
  YRECFM='SST'
  CALL READ_SURF(HPROGRAM,YRECFM,XSST(:),IRESP)
!
ENDIF
!
!* stochastic flux perturbation pattern
!
ALLOCATE(XPERTFLUX(ILU))
IF( LPERTFLUX ) THEN
   CALL READ_SURF(HPROGRAM,'PERTSEAFLUX',XPERTFLUX(:),IRESP)
ELSE
  XPERTFLUX(:) = 0.
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.     Semi-prognostic fields:
!               ----------------------
!
!* roughness length
!
ALLOCATE(XZ0(ILU))
YRECFM='Z0SEA'
XZ0(:) = 0.001
CALL READ_SURF(HPROGRAM,YRECFM,XZ0(:),IRESP)
!
!* flag to use or not the SeaIce model 
!
CALL READ_SURF(HPROGRAM,'VERSION',IVERSION,IRESP)
IF (IVERSION <8) THEN
   LHANDLE_SIC=.FALSE.
ELSE
   CALL READ_SURF(HPROGRAM,'HANDLE_SIC',LHANDLE_SIC,IRESP)
ENDIF
!
!
! * sea surface salinity
!
ALLOCATE(XSSS(ILU))
XSSS(:)=0.0
!
!* Sea surface salinity nudging data
!
IF(LINTERPOL_SSS)THEN
   !
   ! Precedent, Current and Next Monthly/Annual SSS
   INMTH=3
   ! Precedent, Current and Next Annual Monthly SSS
   IF(TRIM(CINTERPOL_SSS)=='ANNUAL')INMTH=14
   !
   ALLOCATE(XSSS_MTH(ILU,INMTH))
   DO JMTH=1,INMTH
      WRITE(YMTH,'(I2)') (JMTH-1)
      YRECFM='SSS_MTH'//ADJUSTL(YMTH(:LEN_TRIM(YMTH)))
      CALL READ_SURF(HPROGRAM,YRECFM,XSSS_MTH(:,JMTH),IRESP)
      CALL CHECK_SEA(YRECFM,XSSS_MTH(:,JMTH))
   ENDDO
   !
   CALL INTERPOL_SST_MTH(TTIME%TDATE%YEAR,TTIME%TDATE%MONTH,TTIME%TDATE%DAY,'S',XSSS)
   !
ELSEIF (IVERSION>=8) THEN
   ! 
   ALLOCATE(XSSS_MTH(0,0))
   !
   YRECFM='SSS'
   CALL READ_SURF(HPROGRAM,YRECFM,XSSS,IRESP)
   IF(LHANDLE_SIC)THEN
     CALL CHECK_SEA(YRECFM,XSSS(:))
   ENDIF
   !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('READ_SEAFLUX_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE CHECK_SEA(HFIELD,PFIELD)
!
USE MODD_SEAFLUX_GRID_n, ONLY : XLON, XLAT
!
IMPLICIT NONE
!
CHARACTER(LEN=12),  INTENT(IN) :: HFIELD
REAL, DIMENSION(:), INTENT(IN) :: PFIELD
!
REAL            :: ZMAX,ZMIN
INTEGER         :: JI, IERRC
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('READ_SEAFLUX_N:CHECK_SEA',0,ZHOOK_HANDLE)
!
ZMIN=0.0
ZMAX=1.0E10
!
IERRC=0
!
DO JI=1,ILU
   IF(PFIELD(JI)>ZMAX.OR.PFIELD(JI)<ZMIN)THEN
      IERRC=IERRC+1
      WRITE(KLUOUT,*)'PROBLEM FIELD '//TRIM(HFIELD)//' =',PFIELD(JI),&
                     'NOT REALISTIC AT LOCATION (LAT/LON)',XLAT(JI),XLON(JI)
   ENDIF
ENDDO
!         
IF(IERRC>0) CALL ABOR1_SFX('READ_SEAFLUX_N: FIELD '//TRIM(HFIELD)//' NOT REALISTIC')
!
IF (LHOOK) CALL DR_HOOK('READ_SEAFLUX_N:CHECK_SEA',1,ZHOOK_HANDLE)

END SUBROUTINE CHECK_SEA
!
!------------------------------------------------------------------------------
END SUBROUTINE READ_SEAFLUX_n
