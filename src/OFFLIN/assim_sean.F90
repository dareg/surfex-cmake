!     ###############################################################################
SUBROUTINE ASSIM_SEA_n(YPROGRAM,KI,PTS_IN,PITM,HTEST)

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
!!      Original    04/2012
!!--------------------------------------------------------------------
!
USE MODD_ASSIM,          ONLY : LPRINT,LAESST,LEXTRAP_SEA
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODD_SURF_ATM_n,     ONLY : CSEA,NR_SEA
USE MODD_SURF_ATM_GRID_n,ONLY : XLAT, XLON
USE MODN_IO_OFFLINE,     ONLY : CPGDFILE,CPREPFILE
!
#ifdef LFI
USE MODD_IO_SURF_LFI,    ONLY : CFILEIN_LFI,CFILE_LFI,CFILEOUT_LFI
#endif
!
USE YOMHOOK,             ONLY : LHOOK,DR_HOOK
USE PARKIND1,            ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_INIT_IO_SURF_n
USE MODI_READ_SURF
USE MODI_IO_BUFF_CLEAN_n
USE MODI_OI_HOR_EXTRAPOL_SURF
USE MODI_FLAG_UPDATE
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN) :: YPROGRAM  ! program calling surf. schemes
INTEGER,            INTENT(IN) :: KI
REAL,DIMENSION(KI), INTENT(IN) :: PTS_IN
REAL,DIMENSION(KI), INTENT(IN) :: PITM
CHARACTER(LEN=2),   INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
REAL(KIND=JPRB)             :: ZHOOK_HANDLE
CHARACTER(LEN=10)           :: YVAR    ! Name of the prognostic variable (in LFI file)
CHARACTER(LEN=100)          :: YPREFIX ! Prefix of the prognostic variable  (in LFI file)
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

#ifdef LFI
  CFILEIN_LFI = CPGDFILE        ! input PGD file (orography)
  CFILE_LFI=CFILEIN_LFI
#endif
  CALL INIT_IO_SURF_n(YPROGRAM,'SEA   ','SURF  ','READ ')
  !
  !  Read orography
  !
  CALL READ_SURF(YPROGRAM,'ZS',        ZALT,  IRESP)
  CALL END_IO_SURF_n(YPROGRAM)
  CALL IO_BUFF_CLEAN_n

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
CALL INIT_IO_SURF_n(YPROGRAM,'SEA   ','SURF  ','READ ')
!
!  Read prognostic variables
!
CALL READ_SURF(YPROGRAM,'TG1',       PTS,   IRESP)
CALL READ_SURF(YPROGRAM,'SST',       PSST,  IRESP)
CALL END_IO_SURF_n(YPROGRAM)
CALL IO_BUFF_CLEAN_n

! Read SST from file or set it to input SST
IF ( .NOT. LAESST ) THEN

  ! Set SST to input
  ZSST(:) = PTS_IN(:)
ELSE

  ! SST analysed in CANARI 
  ZSST(:)    = XUNDEF
  WHERE (PITM(:)< 0.5 .AND. PTS(:)==XUNDEF )
     ZSST(:) = PTS_IN(:)   ! set SST analysis from CANARI
  END WHERE
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
  IF (ZSST(I)/=XUNDEF) THEN
    PSST(I) = ZSST(I)
  ELSE
    IF ( LEXTRAP_SEA ) THEN
      OINTERP_SST(I) = .TRUE.
    PSST(I) = XUNDEF
  ENDIF
  ENDIF
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
  IF (LPRINT) THEN
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

! Write updated SST field
WRITE(*,*) 'WRITING UPDATED SST'
!
#ifdef LFI
CFILEOUT_LFI=CPREPFILE
#endif
CALL FLAG_UPDATE(.FALSE.,.TRUE.,.FALSE.,.FALSE.)
CALL INIT_IO_SURF_n(YPROGRAM,'SEA   ','SURF  ','WRITE')

YVAR='SST'
YPREFIX='X_Y_SST (K)                                       '
CALL WRITE_SURF(YPROGRAM,YVAR,PSST,IRESP,HCOMMENT=YPREFIX)

CALL END_IO_SURF_n(YPROGRAM)
CALL IO_BUFF_CLEAN_n

IF (LHOOK) CALL DR_HOOK('ASSIM_SEA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ASSIM_SEA_n
