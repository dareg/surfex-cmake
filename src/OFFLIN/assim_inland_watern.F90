!     ###############################################################################
SUBROUTINE ASSIM_INLAND_WATER_n(YPROGRAM,KI,PTS_O,PITM,HTEST)

!     ###############################################################################
!
!!****  *ASSIM_INLAND_WATER_n * - Chooses the surface assimilation schemes for INLAND_WATER parts  
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     T. Aspelien
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2012
!!--------------------------------------------------------------------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_ASSIM,          ONLY : LPRINT,LEXTRAP_WATER,LWATERTG2
USE MODN_IO_OFFLINE,     ONLY : CPREPFILE,CPGDFILE
!
USE MODD_SURF_ATM_n,     ONLY : CWATER,NR_WATER
USE MODD_SURF_ATM_GRID_n,ONLY : XLAT, XLON
!
#ifdef LFI
USE MODD_IO_SURF_LFI,    ONLY : CFILEIN_LFI, CFILE_LFI,CFILEOUT_LFI
#endif
!
USE YOMHOOK,             ONLY : LHOOK,DR_HOOK
USE PARKIND1,            ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_INIT_IO_SURF_n
USE MODI_READ_SURF
USE MODI_END_IO_SURF_n
USE MODI_IO_BUFF_CLEAN_n
USE MODI_OI_HOR_EXTRAPOL_SURF
USE MODI_FLAG_UPDATE
USE MODI_WRITE_SURF
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN) :: YPROGRAM  ! program calling surf. schemes
INTEGER,            INTENT(IN) :: KI
REAL,DIMENSION(:), INTENT(IN) :: PTS_O
REAL,DIMENSION(:), INTENT(IN) :: PITM
CHARACTER(LEN=2),   INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
REAL, DIMENSION(SIZE(PTS_O))              :: ZLSTINC
REAL, DIMENSION(SIZE(PTS_O))              :: ZLST
REAL, DIMENSION(SIZE(PTS_O))              :: ZTP
!
REAL,ALLOCATABLE,DIMENSION(:)    :: ZLON
REAL,ALLOCATABLE,DIMENSION(:)    :: ZLAT
REAL,ALLOCATABLE,DIMENSION(:)    :: ZALT
REAL,ALLOCATABLE,DIMENSION(:)    :: ZLST_IN
LOGICAL,ALLOCATABLE,DIMENSION(:) :: OINTERP_LST
!
CHARACTER(LEN=10)                :: YVAR    ! Name of the prognostic variable (in LFI file)
CHARACTER(LEN=100)               :: YPREFIX ! Prefix of the prognostic variable  (in LFI file)
INTEGER                          :: IRESP,I
!
REAL(KIND=JPRB)                  :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ASSIM_INLAND_WATER_N',0,ZHOOK_HANDLE)

IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_INLAND_WATER_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF

WRITE(*,*) 'UPDATING LST FOR INLAND_WATER: ',TRIM(CWATER)

IF ( LEXTRAP_WATER ) THEN
  ALLOCATE(ZALT(KI))

#ifdef LFI
  CFILEIN_LFI = CPGDFILE        ! PGD file orography
  CFILE_LFI=CFILEIN_LFI
#endif

  CALL INIT_IO_SURF_n(YPROGRAM,'WATER ','SURF  ','READ ')
  !
  !  Read orography
  !
  CALL READ_SURF(YPROGRAM,'ZS',        ZALT,  IRESP)
  CALL END_IO_SURF_n(YPROGRAM)
  CALL IO_BUFF_CLEAN_n

  ALLOCATE(OINTERP_LST(KI))
  ALLOCATE(ZLON(KI))
  ALLOCATE(ZLAT(KI))
  ALLOCATE(ZLST_IN(KI))

  ! Set longitudes/latitudes for water point
  DO I=1,KI
    ZLON(I)=XLON(NR_WATER(I))
    ZLAT(I)=XLAT(NR_WATER(I))
  ENDDO
  OINTERP_LST(:) = .FALSE.
ENDIF

!
!------------------------------------------------------------
! READ PREP FILE
!------------------------------------------------------------
!
!   File handling definition
!
#ifdef LFI
CFILEIN_LFI = CPREPFILE        ! input PREP file (surface fields)
CFILE_LFI=CFILEIN_LFI
#endif
!
!   Read grid dimension for allocation
!
CALL INIT_IO_SURF_n(YPROGRAM,'WATER ','SURF  ','READ ')
!
!  Read prognostic variables
!
CALL READ_SURF(YPROGRAM,'TS_WATER',  ZLST,  IRESP)
IF (LWATERTG2) THEN
  CALL READ_SURF(YPROGRAM,'TG2',       ZTP,   IRESP)
ENDIF
CALL END_IO_SURF_n(YPROGRAM)
CALL IO_BUFF_CLEAN_n

ZLSTINC(:) = 0.0
ZLSTINC(:) = ZLST(:)

!
!*     ZLST updated
!
DO I=1,KI
  IF ( LWATERTG2 ) THEN
    !
    !*     ZLST updated from LAND values of climatological TS
    !
    IF (ZTP(I)/=XUNDEF .AND. PITM(I) > 0.5 ) THEN
      ZLST(I) = ZTP(I)
    ELSE
      ! Keep ZLST or do extrapolation from neighbour points
      IF ( LEXTRAP_WATER ) THEN
        OINTERP_LST(I) = .TRUE.
        ZLST(I) = XUNDEF
      ENDIF
    ENDIF

  ELSE
    !
    !*     ZLST updated from from CANARI analysis
    !
    IF ( PITM(I) < 0.5 ) THEN
      ZLST(I) = PTS_O(I)
    ELSE
      ! Keep ZLST or do extrapolation from neighbour points
      IF ( LEXTRAP_WATER ) THEN
        OINTERP_LST(I) = .TRUE.
        ZLST(I) = XUNDEF
      ENDIF
    ENDIF
  ENDIF
ENDDO

IF ( LEXTRAP_WATER ) THEN
  !
  !*     Extrapolation
  !
  ZLST_IN = ZLST
  CALL OI_HOR_EXTRAPOL_SURF(KI,ZLAT,ZLON,ZLST_IN,ZLAT,ZLON,ZLST,OINTERP_LST,ZALT)

  !
  !*     Print values produced by OI_HO_EXTRAPOL_SURF
  !
  IF (LPRINT) THEN
    DO I=1,KI
      IF (OINTERP_LST(I)) THEN
        PRINT *,'Lake surface temperature set to ',ZLST(I),'from nearest neighbour at I=',NR_WATER(I)
      ENDIF
    ENDDO
  ENDIF

  DEALLOCATE(OINTERP_LST)
  DEALLOCATE(ZLON)
  DEALLOCATE(ZLAT)
  DEALLOCATE(ZALT)
  DEALLOCATE(ZLST_IN)
ENDIF

! Sum the increments
ZLSTINC(:) = ZLST(:) - ZLSTINC(:)

WRITE(*,*) 'Mean LST increments over inland water   ',SUM(ZLSTINC)/KI

! Write updated LST field
WRITE(*,*) 'WRITING UPDATED LST'
!
#ifdef LFI
CFILEOUT_LFI=CPREPFILE
#endif
!
CALL FLAG_UPDATE(.FALSE.,.TRUE.,.FALSE.,.FALSE.)
CALL INIT_IO_SURF_n(YPROGRAM,'WATER ','SURF  ','WRITE')

YVAR='TS_WATER'
YPREFIX='X_Y_TS_WATER (K)                                  '
CALL WRITE_SURF(YPROGRAM,YVAR,ZLST,IRESP,HCOMMENT=YPREFIX)

CALL END_IO_SURF_n(YPROGRAM)
CALL IO_BUFF_CLEAN_n

IF (LHOOK) CALL DR_HOOK('ASSIM_INLAND_WATER_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ASSIM_INLAND_WATER_n
