!     ###############################################################################
SUBROUTINE ASSIM_SEA_n(YPROGRAM,KI,PTS_IN,PSST_IN,PSIC_IN,PITM,HTEST)

!     ###############################################################################
!
!!****  *ASSIM_SEA_n * - Chooses the surface assimilation schemes for SEA tile
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
!!      Original       04/2012
!!      Trygve Aspelien, Separating IO  06/2013 
!!--------------------------------------------------------------------
!
USE MODD_ASSIM,          ONLY : NPRINTLEV,LAESST,LEXTRAP_SEA
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_SURF_ATM_n,     ONLY : CSEA,NR_SEA,XZS,XNATURE
USE MODD_SEAFLUX_n,      ONLY : XSST
USE MODD_SURF_ATM_GRID_n,ONLY : XLAT, XLON
USE MODN_IO_OFFLINE,     ONLY : CPGDFILE,CPREPFILE
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
REAL,DIMENSION(KI), INTENT(IN) :: PTS_IN
REAL,DIMENSION(KI), INTENT(IN) :: PSST_IN
REAL,DIMENSION(KI), INTENT(IN) :: PSIC_IN
REAL,DIMENSION(KI), INTENT(IN) :: PITM
CHARACTER(LEN=2),   INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
REAL(KIND=JPRB)             :: ZHOOK_HANDLE
INTEGER                     :: IRESP,I
REAL                        :: ZFMAX,ZFMIN,ZFMEAN
REAL, DIMENSION (KI)        :: ZT2INC
REAL, DIMENSION (KI)        :: ZTCLS
REAL, DIMENSION (KI)        :: ZSST
REAL, DIMENSION (KI)        :: PSST
REAL, DIMENSION (KI)        :: PTS
REAL, DIMENSION (KI)        :: ZSSTINC
REAL,ALLOCATABLE,DIMENSION(:)    :: PLON
REAL,ALLOCATABLE,DIMENSION(:)    :: PLAT
REAL,ALLOCATABLE,DIMENSION(:)    :: ZALT
LOGICAL,ALLOCATABLE,DIMENSION(:) :: OINTERP_SST

IF (LHOOK) CALL DR_HOOK('ASSIM_SEA_N',0,ZHOOK_HANDLE)

IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_SEA_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF

WRITE(*,*) 'UPDATING SST FOR SCHEME: ',TRIM(CSEA)

IF ( LEXTRAP_SEA ) THEN
  ALLOCATE(ZALT(KI))

  ! Set local array from global
  DO I=1,KI
    ZALT(I)=XZS(NR_SEA(I))
  ENDDO

  ALLOCATE(OINTERP_SST(KI))
  ALLOCATE(PLON(KI))
  ALLOCATE(PLAT(KI))

  ! Set longitudes/latitudes for sea point
  DO I=1,KI
    PLON(I)=XLON(NR_SEA(I))
    PLAT(I)=XLAT(NR_SEA(I))
  ENDDO

  OINTERP_SST(:) = .FALSE.
ENDIF

! Set SST from watfluxn
PSST = XSST

! Read SST from file or set it to input SST
IF ( .NOT. LAESST ) THEN

  ! Set SST to input
  ZSST(:) = PSST_IN(:)
ELSE

  ! SST analysed in CANARI 
  ZSST(:)    = XUNDEF
  DO I=1,KI
    IF (PITM(I)< 0.5 .AND. XNATURE(NR_SEA(I)) == 0. ) THEN
     ZSST(I) = PTS_IN(I)   ! set SST analysis from CANARI
    ENDIF
  END DO
  !
  ZFMIN = MINVAL(ZSST)
  ZFMAX = MAXVAL(ZSST)
  ZFMEAN = SUM(ZSST)/FLOAT(KI)
  WRITE(*,*) '  SST analysis from CANARI '
  WRITE(*,'("  ZSST            - min, mean, max: ",3E13.4)') ZFMIN, ZFMEAN, ZFMAX
ENDIF

ZSSTINC(:) = PSST(:)

!*     PSST updated at all sea points with ZSST where ZSST is available

DO I=1,KI
  !
  IF (ZSST(I)/=XUNDEF) THEN
    PSST(I) = ZSST(I)
  ELSEIF ( LEXTRAP_SEA ) THEN
    OINTERP_SST(I) = .TRUE.
    PSST(I) = XUNDEF
  ENDIF
  !
ENDDO

IF ( LEXTRAP_SEA ) THEN
  !
  !*     Extrapolation
  !
  ZSST(:) = PSST(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,PLAT,PLON,ZSST,PLAT,PLON,PSST,OINTERP_SST,ZALT)

  !
  !*     Print values produced by OI_HO_EXTRAPOL_SURF
  !
  IF ( NPRINTLEV > 2 ) THEN
    DO I=1,KI
      IF (OINTERP_SST(I)) THEN
        PRINT *,'Sea surface temperature set to ',PSST(I),'from nearest neighbour at I=',NR_SEA(I)
      ENDIF
    ENDDO
  ENDIF

  DEALLOCATE(OINTERP_SST)
  DEALLOCATE(PLON)
  DEALLOCATE(PLAT)
  DEALLOCATE(ZALT)
ENDIF

! Sum the increments
ZSSTINC(:) = PSST(:) - ZSSTINC(:)

WRITE(*,*) 'Mean SST increments over SEA   ',SUM(ZSSTINC)/KI

! Setting modified variables
XSST=PSST

IF (LHOOK) CALL DR_HOOK('ASSIM_SEA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ASSIM_SEA_n
