!     ###############################################################################
SUBROUTINE ASSIM_ISBA_n(HPROGRAM,KI,                                   &
                        PCON_RAIN, PSTRAT_RAIN, PCON_SNOW, PSTRAT_SNOW,&
                        PCLOUDS,   PLSM,        PEVAPTR,   PEVAP,      &
                        PSWEC,     PTSC,                               &
                        PTS,       PT2M,        PHU2M,     PSWE,       &
                        HTEST )

!     ###############################################################################
!
!!****  *ASSIM_ISBA_n * - Chooses the surface assimilation schemes for ISBA
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
USE MODD_ASSIM,          ONLY : CASSIM_ISBA,LAESNM,LEXTRAP_NATURE,LPRINT
USE MODD_SURF_ATM_n,     ONLY : NR_NATURE
USE MODD_SURF_ATM_GRID_n,ONLY : XLAT, XLON
USE MODN_IO_OFFLINE,     ONLY : CPGDFILE,CPREPFILE
!
#ifdef LFI
USE MODD_IO_SURF_LFI,    ONLY : CFILEIN_LFI,CFILE_LFI,CFILEOUT_LFI
#endif
!
USE YOMHOOK,             ONLY : LHOOK,   DR_HOOK
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
USE MODI_ASSIM_ISBA_UPDATE_SNOW
USE MODI_ASSIM_NATURE_ISBA_EKF
USE MODI_ASSIM_NATURE_ISBA_OI
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),    INTENT(IN) :: HPROGRAM  ! program calling surf. schemes
INTEGER,             INTENT(IN) :: KI
REAL, DIMENSION(KI), INTENT(IN) :: PCON_RAIN
REAL, DIMENSION(KI), INTENT(IN) :: PSTRAT_RAIN
REAL, DIMENSION(KI), INTENT(IN) :: PCON_SNOW
REAL, DIMENSION(KI), INTENT(IN) :: PSTRAT_SNOW
REAL, DIMENSION(KI), INTENT(IN) :: PCLOUDS
REAL, DIMENSION(KI), INTENT(IN) :: PLSM
REAL, DIMENSION(KI), INTENT(IN) :: PEVAPTR
REAL, DIMENSION(KI), INTENT(IN) :: PEVAP
REAL, DIMENSION(KI), INTENT(IN) :: PSWEC
REAL, DIMENSION(KI), INTENT(IN) :: PTSC
REAL, DIMENSION(KI), INTENT(IN) :: PTS
REAL, DIMENSION(KI), INTENT(IN) :: PT2M
REAL, DIMENSION(KI), INTENT(IN) :: PHU2M
REAL, DIMENSION(KI), INTENT(IN) :: PSWE
CHARACTER(LEN=2),    INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
REAL(KIND=JPRB)                    :: ZHOOK_HANDLE
LOGICAL, DIMENSION(:), ALLOCATABLE :: OINTERP_NATURE
REAL,    DIMENSION(:), ALLOCATABLE :: PLON
REAL,    DIMENSION(:), ALLOCATABLE :: PLAT
REAL,    DIMENSION(:), ALLOCATABLE :: PTS_EP,ZTS_EP
REAL,    DIMENSION(:), ALLOCATABLE :: PTP_EP,ZTP_EP
REAL,    DIMENSION(:), ALLOCATABLE :: PWS_EP,ZWS_EP
REAL,    DIMENSION(:), ALLOCATABLE :: PWP_EP,ZWP_EP
REAL,    DIMENSION(:), ALLOCATABLE :: PTL_EP,ZTL_EP
REAL,    DIMENSION(:), ALLOCATABLE :: PSNS_EP,ZSNS_EP
REAL,    DIMENSION(:), ALLOCATABLE :: ZALT
INTEGER                            :: I,IRESP
CHARACTER(LEN=10)                  :: YVAR    ! Name of the prognostic variable (in LFI file)
CHARACTER(LEN=100)                 :: YPREFIX ! Prefix of the prognostic variable  (in LFI file)

IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_N',0,ZHOOK_HANDLE)

IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_ISBA_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF

IF ( LEXTRAP_NATURE ) THEN
  ALLOCATE(OINTERP_NATURE(KI))
  ALLOCATE(PLON(KI))
  ALLOCATE(PLAT(KI))
  ALLOCATE(PTS_EP(KI),ZTS_EP(KI))
  ALLOCATE(PTP_EP(KI),ZTP_EP(KI))
  ALLOCATE(PWS_EP(KI),ZWS_EP(KI))
  ALLOCATE(PWP_EP(KI),ZWP_EP(KI))
  ALLOCATE(PTL_EP(KI),ZTL_EP(KI))
  ALLOCATE(PSNS_EP(KI),ZSNS_EP(KI))
  ALLOCATE(ZALT(KI))

  ! Set longitudes/latitudes for water point
  DO I=1,KI
    PLON(I)=XLON(NR_NATURE(I))
    PLAT(I)=XLAT(NR_NATURE(I))
  ENDDO
  OINTERP_NATURE = .FALSE.

  DO I=1,KI
    IF ( PLSM(I) < 0.5 ) THEN
      OINTERP_NATURE(I) = .TRUE.
    ENDIF
  ENDDO

  !   Read orography
#ifdef LFI
  CFILEIN_LFI = CPGDFILE
  CFILE_LFI=CFILEIN_LFI
#endif
  CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','SURF  ','READ ')
  CALL READ_SURF(HPROGRAM,'ZS',        ZALT,     IRESP)
  !
  CALL END_IO_SURF_n(HPROGRAM)
  CALL IO_BUFF_CLEAN_n

  !   Read field to extrapolate from
#ifdef LFI
  CFILEIN_LFI = CPREPFILE
  CFILE_LFI=CFILEIN_LFI
#endif
  CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','SURF  ','READ ')

  CALL READ_SURF(HPROGRAM,'WG1',       PWS_EP,   IRESP)
  CALL READ_SURF(HPROGRAM,'WG2',       PWP_EP,   IRESP)
  CALL READ_SURF(HPROGRAM,'TG1',       PTS_EP,   IRESP)
  CALL READ_SURF(HPROGRAM,'TG2',       PTP_EP,   IRESP)
  CALL READ_SURF(HPROGRAM,'WGI2',      PTL_EP,   IRESP)
  CALL READ_SURF(HPROGRAM,'WSNOW_VEG1',PSNS_EP,  IRESP)
  !
  CALL END_IO_SURF_n(HPROGRAM)
  CALL IO_BUFF_CLEAN_n

  ! Search for the nearest grid point values for land surface fields
  ! at locations where the CANARI land fraction is less than 50%
  ! and therefore useless values MIGTH be given
  !
  ZTS_EP(:) = PTS_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,PLAT,PLON,ZTS_EP,PLAT,PLON,PTS_EP,OINTERP_NATURE,ZALT)
  ZTP_EP(:) = PTP_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,PLAT,PLON,ZTP_EP,PLAT,PLON,PTP_EP,OINTERP_NATURE,ZALT)
  ZWS_EP(:) = PWS_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,PLAT,PLON,ZWS_EP,PLAT,PLON,PWS_EP,OINTERP_NATURE)
  ZWP_EP(:) = PWP_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,PLAT,PLON,ZWP_EP,PLAT,PLON,PWP_EP,OINTERP_NATURE)
  ZTL_EP(:) = PTL_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,PLAT,PLON,ZTL_EP,PLAT,PLON,PTL_EP,OINTERP_NATURE)
  ZSNS_EP(:) = PSNS_EP(:)
  CALL OI_HOR_EXTRAPOL_SURF(KI,PLAT,PLON,ZSNS_EP,PLAT,PLON,PSNS_EP,OINTERP_NATURE)

  !
  ! PRINT values produced by OI_HO_EXTRAPOL_SURF for TS
  !
  IF (LPRINT) THEN
    DO I=1,KI
     IF (OINTERP_NATURE(I)) THEN
       PRINT *,'Surface temperature set to ',PTS_EP(I),'from nearest neighbour at I=',NR_NATURE(I)
     ENDIF
    ENDDO
  ENDIF

  !
  !   Write extrapolated fields to file
  !
#ifdef LFI
  CFILEOUT_LFI=CPREPFILE
#endif
  CALL FLAG_UPDATE(.FALSE.,.TRUE.,.FALSE.,.FALSE.)
  CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','SURF  ','WRITE')
  !
  YVAR='WG1'
  YPREFIX='X_Y_WG1 (m3/m3)                                   '
  CALL WRITE_SURF(HPROGRAM,YVAR,PWS_EP,IRESP,HCOMMENT=YPREFIX)
  YVAR='WG2'
  YPREFIX='X_Y_WG2 (m3/m3)                                   '
  CALL WRITE_SURF(HPROGRAM,YVAR,PWP_EP,IRESP,HCOMMENT=YPREFIX)
  YVAR='WGI2'
  YPREFIX='X_Y_WGI2 (m3/m3)                                  '
  CALL WRITE_SURF(HPROGRAM,YVAR,PTL_EP,IRESP,HCOMMENT=YPREFIX)
  YVAR='TG1'
  YPREFIX='X_Y_TG1 (K)                                       '
  CALL WRITE_SURF(HPROGRAM,YVAR,PTS_EP,IRESP,HCOMMENT=YPREFIX)
  YVAR='TG2'
  YPREFIX='X_Y_TG2 (K)                                       '
  CALL WRITE_SURF(HPROGRAM,YVAR,PTP_EP,IRESP,HCOMMENT=YPREFIX)
  YVAR='WSNOW_VEG1'
  YPREFIX='X_Y_WSNOW_VEG1 (kg/m2)                            '
  CALL WRITE_SURF(HPROGRAM,YVAR,PSNS_EP,IRESP,HCOMMENT=YPREFIX)
  !
  CALL END_IO_SURF_n(HPROGRAM)
  CALL IO_BUFF_CLEAN_n

  DEALLOCATE(OINTERP_NATURE)
  DEALLOCATE(PLON)
  DEALLOCATE(PLAT)
  DEALLOCATE(PTS_EP,ZTS_EP)
  DEALLOCATE(PTP_EP,ZTP_EP)
  DEALLOCATE(PWS_EP,ZWS_EP)
  DEALLOCATE(PWP_EP,ZWP_EP)
  DEALLOCATE(PTL_EP,ZTL_EP)
  DEALLOCATE(PSNS_EP,ZSNS_EP)
  DEALLOCATE(ZALT)
ENDIF

! Snow analysis/update
IF (LAESNM) THEN
  WRITE(*,*) 'UPDATE SNOW FROM ANALYSED VALUES'
  CALL ASSIM_ISBA_UPDATE_SNOW(HPROGRAM,KI,PSWE,HTEST)
ELSE
  WRITE(*,*) 'SNOW IS NOT UPDATED FROM ANALYSED VALUES'
ENDIF

! Soil assimilation
IF ( CASSIM_ISBA == 'EKF  ' ) THEN

  CALL ASSIM_NATURE_ISBA_EKF(HPROGRAM,KI,   &
                             PT2M,    PHU2M,&
                             HTEST )

ELSEIF ( CASSIM_ISBA == 'OI   ' ) THEN

  CALL ASSIM_NATURE_ISBA_OI(HPROGRAM, KI,                                  &
                            PCON_RAIN, PSTRAT_RAIN, PCON_SNOW, PSTRAT_SNOW,&
                            PCLOUDS,   PLSM,        PEVAPTR,   PEVAP,      &
                            PSWEC,     PTSC,                               &
                            PTS,       PT2M,        PHU2M,                 &
                            HTEST )
ELSE
  CALL ABOR1_SFX(CASSIM_ISBA//' is not a defined scheme for ASSIM_ISBA_N')
ENDIF

IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ASSIM_ISBA_n
