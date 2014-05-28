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
!!      Trygve Aspelien, Separating IO  06/2013
!!--------------------------------------------------------------------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_ASSIM,          ONLY : NPRINTLEV,LEXTRAP_WATER,LWATERTG2
USE MODN_IO_OFFLINE,     ONLY : CPREPFILE,CPGDFILE
USE MODD_ISBA_n,         ONLY : XTG
USE MODD_WATFLUX_n,      ONLY : XTS,XZS
USE MODD_SURF_ATM_n,     ONLY : CWATER,NR_WATER,NSIZE_WATER,NR_NATURE,NSIZE_NATURE
USE MODD_SURF_ATM_GRID_n,ONLY : XLAT, XLON
USE YOMHOOK,             ONLY : LHOOK,DR_HOOK
USE PARKIND1,            ONLY : JPRB

USE MODI_ABOR1_SFX
USE MODI_OI_HOR_EXTRAPOL_SURF
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN) :: YPROGRAM  ! program calling surf. schemes
INTEGER,            INTENT(IN) :: KI
REAL,DIMENSION(KI), INTENT(IN) :: PTS_O
REAL,DIMENSION(KI), INTENT(IN) :: PITM
CHARACTER(LEN=2),   INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
INTEGER                          :: IRESP,I,J
REAL(KIND=JPRB)                  :: ZHOOK_HANDLE
REAL, DIMENSION(KI)              :: ZLSTINC
REAL, DIMENSION(KI)              :: ZLST
REAL,ALLOCATABLE,DIMENSION(:)    :: ZTP
REAL,ALLOCATABLE,DIMENSION(:)    :: ZLON
REAL,ALLOCATABLE,DIMENSION(:)    :: ZLAT
REAL,ALLOCATABLE,DIMENSION(:)    :: ZALT
REAL,ALLOCATABLE,DIMENSION(:)    :: ZLST_IN
LOGICAL,ALLOCATABLE,DIMENSION(:) :: OINTERP_LST
INTEGER                          :: NPATCH=1

IF (LHOOK) CALL DR_HOOK('ASSIM_INLAND_WATER_N',0,ZHOOK_HANDLE)

IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_INLAND_WATER_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF

WRITE(*,*) 'UPDATING LST FOR INLAND_WATER: ',TRIM(CWATER)

! Set local array from global
ZLST=XTS

IF ( LEXTRAP_WATER ) THEN
  ALLOCATE(ZALT(KI))
  ! Set local array from global
  ZALT=XZS

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

ZLSTINC(:) = 0.0
ZLSTINC(:) = ZLST(:)

!
!*     ZLST updated
!

IF ( LWATERTG2 ) ALLOCATE(ZTP(KI))

DO I=1,KI
  !
  IF ( LWATERTG2 ) THEN

    ! Set TG2 from global array
    ZTP(I)=XUNDEF
    loop: DO J=1,NSIZE_NATURE
      IF ( NR_WATER(I) == NR_NATURE(J) ) THEN
        ZTP(I)=XTG(J,2,NPATCH)
        CYCLE loop
      ENDIF
    ENDDO loop
    !
    !*     ZLST updated from LAND values of climatological TS
    !
    IF (ZTP(I)/=XUNDEF .AND. PITM(I) > 0.5 ) THEN
      ZLST(I) = ZTP(I)
    ELSEIF ( LEXTRAP_WATER ) THEN
      ! Keep ZLST or do extrapolation from neighbour points
      OINTERP_LST(I) = .TRUE.
      ZLST(I) = XUNDEF
    ENDIF
    !
  ELSE
    !
    !*     ZLST updated from from CANARI analysis
    !
    IF ( PITM(I) < 0.5 ) THEN
      ZLST(I) = PTS_O(I)
    ELSEIF ( LEXTRAP_WATER ) THEN
      ! Keep ZLST or do extrapolation from neighbour points
      OINTERP_LST(I) = .TRUE.
      ZLST(I) = XUNDEF
    ENDIF
    !
  ENDIF
  !
ENDDO
!
IF ( LWATERTG2 ) DEALLOCATE(ZTP)


IF ( LEXTRAP_WATER ) THEN
  !
  !*     Extrapolation
  !
  ZLST_IN = ZLST
  CALL OI_HOR_EXTRAPOL_SURF(KI,ZLAT,ZLON,ZLST_IN,ZLAT,ZLON,ZLST,OINTERP_LST,ZALT)

  !
  !*     Print values produced by OI_HO_EXTRAPOL_SURF
  !
  IF ( NPRINTLEV > 2 ) THEN
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

! Setting modified variables
XTS=ZLST

IF (LHOOK) CALL DR_HOOK('ASSIM_INLAND_WATER_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ASSIM_INLAND_WATER_n
