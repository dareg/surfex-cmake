SUBROUTINE ASSIM_NATURE_ISBA_OI(YPROGRAM, KI,                     &
                                PRRCL,    PRRSL,  PRRCN,   PRRSN, &
                                PATMNEB,  PITM,   PEVAPTR, PEVAP, &
                                PSNC,     PTSC,                   &
                                PTS_O,    PT2M_O, PHU2M_O, PSWE,  &
                                HTEST )

! ------------------------------------------------------------------------------------------
!  *****************************************************************************************
!
!  Routine to perform OI within SURFEX 
!  a soil analysis for water content and temperature 
!  using the Meteo-France optimum interpolation technique of Giard and Bazile (2000)
!
!  Derived from CANARI subroutines externalized by Lora Taseva (Dec. 2007)
!
!  Author : Jean-Francois Mahfouf (01/2008)
!
!  Modifications : 
!   (05/2008)  : The I/O of this version follow the newly available LFI format in SURFEX  
!   (01/2009)  : Read directly atmospheric FA files using XRD library instead of using "edf"
!   (06/2009)  : Modifications to allow the assimilation of ASCAT superficial soil moisture
!   (09/2010)  : More parameters to goto_surfex
!   (03/2011)  : Initialization of ZEVAPTR (F.Bouyssel)
!   (07/2011)  : Read pgd+prep (B. Decharme)
!   (04/2012)  : Made as a subroutine (T. Aspelien)
! ******************************************************************************************
! ------------------------------------------------------------------------------------------
 USE MODD_TYPE_DATE_SURF,  ONLY : DATE_TIME
 USE MODD_CSTS,            ONLY : XDAY, XPI, XG, XRHOLW, XLVTT, NDAYSEC
 USE MODD_SURF_PAR,        ONLY : XUNDEF 
!
 USE MODD_ASSIM,           ONLY : LOBSWG,ITRAD,LPRINT,NECHGU,RD1,RSCALDW,&
                                  RTHR_QC,SIGWGB,SIGWGO,SIGWGO_MAX
 USE MODN_IO_OFFLINE,     ONLY  : CPGDFILE,CPREPFILE
!
#ifdef LFI
 USE MODD_IO_SURF_LFI,     ONLY : CFILEIN_LFI, CFILEOUT_LFI
#endif
!
 USE MODD_SURF_ATM_n,      ONLY : NDIM_NATURE,NR_NATURE
 USE MODD_SURF_ATM_GRID_n, ONLY : XLAT, XLON
!
 USE YOMHOOK,              ONLY : LHOOK,   DR_HOOK
 USE PARKIND1,             ONLY : JPRB
!
 USE MODI_ABOR1_SFX
 USE MODI_INIT_IO_SURF_n
 USE MODI_READ_SURF
 USE MODI_END_IO_SURF_n
 USE MODI_TRANS_CHAINE
 USE MODI_IO_BUFF_CLEAN_n
 USE MODI_OI_BC_SOIL_MOISTURE
 USE MODI_OI_CACSTS
 USE MODI_FLAG_UPDATE
 USE MODI_WRITE_SURF
!
 IMPLICIT NONE
 CHARACTER(LEN=6),    INTENT(IN) :: YPROGRAM  ! program calling surf. schemes
 INTEGER,             INTENT(IN) :: KI
 REAL, DIMENSION(KI), INTENT(IN) :: PRRCL
 REAL, DIMENSION(KI), INTENT(IN) :: PRRSL
 REAL, DIMENSION(KI), INTENT(IN) :: PRRCN
 REAL, DIMENSION(KI), INTENT(IN) :: PRRSN
 REAL, DIMENSION(KI), INTENT(IN) :: PATMNEB
 REAL, DIMENSION(KI), INTENT(IN) :: PITM
 REAL, DIMENSION(KI), INTENT(IN) :: PEVAPTR
 REAL, DIMENSION(KI), INTENT(IN) :: PEVAP
 REAL, DIMENSION(KI), INTENT(IN) :: PSNC
 REAL, DIMENSION(KI), INTENT(IN) :: PTSC
 REAL, DIMENSION(KI), INTENT(IN) :: PTS_O
 REAL, DIMENSION(KI), INTENT(IN) :: PT2M_O
 REAL, DIMENSION(KI), INTENT(IN) :: PHU2M_O
 REAL, DIMENSION(KI), INTENT(OUT):: PSWE
 CHARACTER(LEN=2),    INTENT(IN) :: HTEST ! must be equal to 'OK'
!
 INTEGER                                  :: IDAT
!
!    Declarations of local variables
!
 CHARACTER(LEN=3)                         :: YREAD
 CHARACTER(LEN=2)                         :: CMONTH
 INTEGER                                  :: IYEAR                      ! current year (UTC)
 INTEGER                                  :: IMONTH                     ! current month (UTC)
 INTEGER                                  :: IDAY                       ! current day (UTC)
 INTEGER                                  :: NSSSSS                     ! current time since start of the run (s)
 INTEGER                                  :: IRESP                      ! return code
 TYPE (DATE_TIME)                         :: TTIME                      ! Current date and time  

!
! Arrays for soil OI analysis
!
 
 REAL, DIMENSION (KI) :: PWS, PWP, PTS, PTP, PTL, PLAI, PVEG, PRSMIN, PD2, PSAB, PARG,        &
                         PTSUN, PZENITH, PAZIMSOL, PTCLS, PHCLS, PRAIN, PSNOW,                &
                         PWIND, PSWD, PSWS, PUCLS, PVCLS, PSNS,                               &
                         ZT2INC, ZH2INC, ZWS, ZWP, ZTL, ZTS, ZTP, ZTCLS, ZHCLS, ZUCLS,        &
                         ZVCLS, PSSTC, PWPINC1, PWPINC2, PWPINC3, PT2MBIAS, PH2MBIAS,         &
                         PALBF, PEMISF, PZ0F, PIVEG, PZ0H,                                    &
                         PTPC, PWSC, PWPC, ZEVAP, ZEVAPTR, PGELAT, PGELAM, PGEMU,             &
                         ZWSINC, ZWPINC, ZSNS, ZTSINC, ZTPINC, ZTLINC,                        &
                         PSM_O,  PSIG_SMO, PLSM_O, PWS_O,  ZWGINC ,                           &
                         ZALT  
 
 REAL,DIMENSION(KI)   :: PLON,  PLAT
 INTEGER              :: I,J
 CHARACTER(LEN=10)    :: YVAR    ! Name of the prognostic variable (in LFI file)
 CHARACTER(LEN=100)   :: YPREFIX ! Prefix of the prognostic variable  (in LFI file)
 INTEGER              :: INOBS   ! number of observations
 REAL                 :: ZTHRES
 REAL(KIND=JPRB)      :: ZHOOK_HANDLE
!
! ----------------------------------------------------------------------------------
!
 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_OI',0,ZHOOK_HANDLE)

 IF (HTEST/='OK') THEN
   CALL ABOR1_SFX('ASSIM_ISBA_n: FATAL ERROR DURING ARGUMENT TRANSFER')
 END IF
 
 PRINT *,'--------------------------------------------------------------------------'
 PRINT *,'|                                                                        |'
 PRINT *,'|                             ENTER OI_ASSIM                             |'
 PRINT *,'|                                                                        |'
 PRINT *,'--------------------------------------------------------------------------'

!
!   Update some constants dependant from NACVEG
!
!  scaling of soil moisture increments when assimilation window is different
!  from 6 hours
 RSCALDW = REAL(NECHGU)/6.0
!  half assimilation window in sec
  ITRAD   = NECHGU*1800

!
!
 ! READ PGD FILE
!------------------------------------------------------------
!
!   File handling definition
!
#ifdef LFI
 CFILEIN_LFI = CPGDFILE        ! input PGD file (surface fields)
#endif
!
!   Read grid dimension for allocation
!
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','SURF  ','READ ')

 CALL READ_SURF(YPROGRAM,'SAND',      PSAB,  IRESP)
 CALL READ_SURF(YPROGRAM,'CLAY',      PARG,  IRESP)
 CALL READ_SURF(YPROGRAM,'ZS',        ZALT,  IRESP)
!
 CALL END_IO_SURF_n(YPROGRAM)
!
!------------------------------------------------------------
! READ PREP FILE
!------------------------------------------------------------
!
!   File handling definition
!
#ifdef LFI
 CFILEIN_LFI = CPREPFILE        ! input PREP file (surface fields)
#endif
!
!   Read grid dimension for allocation
!
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','SURF  ','READ ')
!
!  Read prognostic variables
!
 CALL READ_SURF(YPROGRAM,'WG1',       PWS,   IRESP)
 CALL READ_SURF(YPROGRAM,'WG2',       PWP,   IRESP)
 CALL READ_SURF(YPROGRAM,'TG1',       PTS,   IRESP)
 CALL READ_SURF(YPROGRAM,'TG2',       PTP,   IRESP)
 CALL READ_SURF(YPROGRAM,'WGI2',      PTL,   IRESP)
 CALL READ_SURF(YPROGRAM,'WSNOW_VEG1',PSNS,  IRESP)
!
 CALL READ_SURF(YPROGRAM,'STORAGETYPE', YREAD, IRESP)
 IF (YREAD=='ALL') THEN
   CALL READ_SURF(YPROGRAM,'T2M',       PTCLS, IRESP)
   CALL READ_SURF(YPROGRAM,'HU2M',      PHCLS, IRESP)
   CALL READ_SURF(YPROGRAM,'ZON10M',    PUCLS, IRESP)
   CALL READ_SURF(YPROGRAM,'MER10M',    PVCLS, IRESP)
 ENDIF
!
! Read diag of constant surface fields (present in PREP file)
!
 CALL READ_SURF(YPROGRAM,'RSMIN',     PRSMIN,IRESP)
 CALL READ_SURF(YPROGRAM,'DG2',       PD2,   IRESP)
 CALL READ_SURF(YPROGRAM,'LAI',       PLAI,  IRESP)
 CALL READ_SURF(YPROGRAM,'VEG',       PVEG,  IRESP)

 ! Set PIVEG (SURFIND.VEG.DOMI) since it is not available
 PIVEG = 0.0

 !   Find current time
 !
 CALL READ_SURF(YPROGRAM,'DTCUR',TTIME,IRESP)
 !
 !   Time initializations 
 !
 IYEAR  = TTIME%TDATE%YEAR
 IMONTH = TTIME%TDATE%MONTH
 IDAY   = TTIME%TDATE%DAY
 NSSSSS = TTIME%TIME
 IF (NSSSSS > NDAYSEC) NSSSSS = NSSSSS - NDAYSEC
 CALL TRANS_CHAINE(CMONTH,IMONTH,2)
 IDAT = IYEAR*10000. + IMONTH*100. + IDAY

!
! PRINT 
!
 IF (LPRINT) THEN
  PRINT *,'value in PREP file => WG1       ',SUM(PWS)/KI
  PRINT *,'value in PREP file => WG2       ',SUM(PWP)/KI
  PRINT *,'value in PREP file => TG1       ',SUM(PTS)/KI
  PRINT *,'value in PREP file => TG2       ',SUM(PTP)/KI
  PRINT *,'value in PREP file => WGI2      ',SUM(PTL)/KI
  PRINT *,'value in PREP file => WSNOW_VEG1',SUM(PSNS)/KI
  PRINT *,'value in PREP file => LAI       ',SUM(PLAI)/KI
  PRINT *,'value in PREP file => VEG       ',SUM(PVEG)/KI
  PRINT *,'value in PREP file => RSMIN     ',SUM(PRSMIN)/KI
  PRINT *,'value in PREP file => DATA_DG2  ',SUM(PD2)/KI
  PRINT *,'value in PREP file => SAND      ',SUM(PSAB)/KI
  PRINT *,'value in PREP file => CLAY      ',SUM(PARG)/KI
  PRINT *,'value in PREP file => ZS        ',SUM(ZALT)/KI
 ENDIF
!
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n
!
!
!  Read ASCAT SM observations (in percent)
!


 INOBS = 0
 IF (LOBSWG) THEN
   OPEN(UNIT=111,FILE='ASCAT_SM.DAT')
   DO I=1,KI
     READ(111,*) PSM_O(I),PSIG_SMO(I),PLSM_O(I)
     IF (PLSM_O(I) < 1.0)          PSM_O(I) = 999.0 ! data rejection if not on land
     IF (PSIG_SMO(I) > SIGWGO_MAX) PSM_O(I) = 999.0 ! data rejection of error too large
     IF (PSM_O(I) /= 999.0) INOBS = INOBS + 1
   ENDDO
   CLOSE(UNIT=111)
   PRINT *,'READ ASCAT SM OK'
 ELSE
   PSM_O(:)    = 999.0
   PSIG_SMO(:) = 999.0
   PLSM_O(:)   = 0.0
 ENDIF
 PRINT *,' NUMBER OF ASCAT OBSERVATIONS AFTER INITIAL CHECKS  :: ',INOBS
 INOBS = 0
!
! Perform bias correction of SM observations
!
 CALL OI_BC_SOIL_MOISTURE(KI,PSM_O,PSAB,PWS_O)

 !Set longitudes/latitudes
 DO I=1,NDIM_NATURE
   PLON(I)=XLON(NR_NATURE(I))
   PLAT(I)=XLAT(NR_NATURE(I))
 ENDDO

!
! Screen-level innovations
!
 ZT2INC(:) = PT2M_O(:) - PTCLS(:)
 ZH2INC(:) = PHU2M_O(:) - PHCLS(:)
!
! Threshold for background check
!
 ZTHRES=RTHR_QC*SQRT(SIGWGO**2 + SIGWGB**2)
!
! Superficial soil moisture innovations in (m3/m3)
!
 DO I=1,KI
   IF (PWS_O(I) /= 999.0) THEN
     ZWGINC(I) = PWS_O(I) - PWS(I)
     IF (ABS(ZWGINC(I)) > ZTHRES) THEN 
       ZWGINC(I) = 0.0 ! background check
     ELSE
       INOBS = INOBS + 1
     ENDIF
   ELSE
     ZWGINC(I) = 0.0
   ENDIF
 ENDDO
 PRINT *,' NUMBER OF ASCAT OBSERVATIONS AFTER BACKGROUND CHECK  :: ',INOBS
!
 PRINT *,'           '
 PRINT *,'Mean T2m increments  ',SUM(ZT2INC)/KI
 PRINT *,'Mean HU2m increments ',SUM(ZH2INC)/KI
 PRINT *,'           '
!
! Interface (define arrays and perform unit conversions)
!
 PARG(:)     = PARG(:)*100.0
 PSAB(:)     = PSAB(:)*100.0
!
 ZWS(:)      = PWS(:)*RD1*XRHOLW     ! conversion of m3/m3 -> mm
 ZWP(:)      = PWP(:)*PD2(:)*XRHOLW  ! conversion of m3/m3 -> mm
 ZTL(:)      = PTL(:)*PD2(:)*XRHOLW  ! conversion of m3/m3 -> mm
 ZTCLS(:)    = PTCLS(:)
 ZHCLS(:)    = PHCLS(:)
 ZUCLS(:)    = PUCLS(:)
 ZVCLS(:)    = PVCLS(:)

 ! SST not used in cacsts
 PSSTC(:)    = 0.

 PWPINC1(:)  = XUNDEF
 PWPINC2(:)  = XUNDEF
 PWPINC3(:)  = XUNDEF
 PT2MBIAS(:) = XUNDEF
 PH2MBIAS(:) = XUNDEF

!
! Sea-ice surface properties
!
 PALBF(:)    = XUNDEF
 PEMISF(:)   = XUNDEF
 PZ0F(:)     = XUNDEF
 PZ0H(:)     = XUNDEF
!
! Climatological arrays set to missing values
!
 PWSC(:)     =  XUNDEF
 PWPC(:)     =  XUNDEF
 PTPC(:)     =  XUNDEF
! 
 DO I=1,KI
   PGELAT(I)   = PLAT(I) 
   PGELAM(I)   = PLON(I) 
   PGEMU(I)    = SIN(PLAT(I)*XPI/180.)
 ENDDO
!
 ZEVAP(:)   =  (PEVAP(:)/XLVTT*XDAY)/(NECHGU*3600.) ! conversion W/m2 -> mm/day
 ZEVAPTR(:) =  PEVAPTR(:)*XDAY 
 ZSNS(:)    =  PSNS(:)
!
 DO I=1,KI
   ZTS(I) = PTS(I)
   ZTP(I) = PTP(I)
 ENDDO
!
!
!  Soil analysis based on optimal interpolation
!
 write(*,*) 'PERFORMING OI SOIL ANALYSIS'
 CALL OI_CACSTS(KI,ZT2INC,ZH2INC,ZWGINC,PWS_O,                         &
                IDAT,NSSSSS,                                           &
                ZTP,ZWP,ZTL,ZSNS,ZTS,ZWS,                              &
                ZTCLS,ZHCLS,ZUCLS,ZVCLS,PSSTC,PWPINC1,PWPINC2,PWPINC3, &
                PT2MBIAS,PH2MBIAS,                                     &
                PRRCL,PRRSL,PRRCN,PRRSN,PATMNEB,ZEVAP,ZEVAPTR,         &
                PITM,PVEG,PALBF,PEMISF,PZ0F,                           &
                PIVEG,PARG,PD2,PSAB,PLAI,PRSMIN,PZ0H,                  &
                PTSC,PTPC,PWSC,PWPC,PSNC,                              &
                PGELAT,PGELAM,PGEMU)  

! Update snow
 PSWE=ZSNS

!
!  Perform soil moiture analyses
!
 ZWSINC(:) = 0.0
 ZWPINC(:) = 0.0
 ZTLINC(:) = 0.0
!
 ZWSINC(:) = ZWS(:) - PWS(:)*(RD1*XRHOLW)    
 ZWPINC(:) = ZWP(:) - PWP(:)*(PD2(:)*XRHOLW) 
 ZTLINC(:) = ZTL(:) - PTL(:)*(PD2(:)*XRHOLW) 
!
 PWS(:)  = ZWS(:)/(RD1*XRHOLW)
 PWP(:)  = ZWP(:)/(PD2(:)*XRHOLW)
 PTL(:)  = ZTL(:)/(PD2(:)*XRHOLW)
!
!  Perform temperature analyses
!
!
 ZTSINC(:) = 0.0
 ZTPINC(:) = 0.0
!
 ZTSINC(:) = ZTS(:) - PTS(:)
 ZTPINC(:) = ZTP(:) - PTP(:)
!
 PTS(:)    = ZTS(:)
 PTP(:)    = ZTP(:)

!
!
! PRINT statistics of the soil analysis
!
 PRINT *,'---------------------------------------------------------------'
 PRINT *,'Mean WS increments over NATURE ',SUM(ZWSINC)/KI
 PRINT *,'Mean WP increments over NATURE ',SUM(ZWPINC)/KI
 PRINT *,'Mean TS increments over NATURE ',SUM(ZTSINC)/KI
 PRINT *,'Mean TP increments over NATURE ',SUM(ZTPINC)/KI
 PRINT *,'Mean TL increments over NATURE ',SUM(ZTLINC)/KI
 PRINT *,'---------------------------------------------------------------'
!
!   Write analysis in LFI file PREP
!
#ifdef LFI
 CFILEOUT_LFI=CPREPFILE
#endif
 CALL FLAG_UPDATE(.FALSE.,.TRUE.,.FALSE.,.FALSE.)
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','SURF  ','WRITE')
!
 YVAR='WG1'
 YPREFIX='X_Y_WG1 (m3/m3)                                   '
 CALL WRITE_SURF(YPROGRAM,YVAR,PWS,IRESP,HCOMMENT=YPREFIX)
 YVAR='WG2'
 YPREFIX='X_Y_WG2 (m3/m3)                                   '
 CALL WRITE_SURF(YPROGRAM,YVAR,PWP,IRESP,HCOMMENT=YPREFIX)
 YVAR='WGI2'
 YPREFIX='X_Y_WGI2 (m3/m3)                                  '
 CALL WRITE_SURF(YPROGRAM,YVAR,PTL,IRESP,HCOMMENT=YPREFIX)
 YVAR='TG1'
 YPREFIX='X_Y_TG1 (K)                                       '
 CALL WRITE_SURF(YPROGRAM,YVAR,PTS,IRESP,HCOMMENT=YPREFIX)
 YVAR='TG2'
 YPREFIX='X_Y_TG2 (K)                                       '
 CALL WRITE_SURF(YPROGRAM,YVAR,PTP,IRESP,HCOMMENT=YPREFIX)
!
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n

!
! -------------------------------------------------------------------------------------
 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_OI',1,ZHOOK_HANDLE)
 END SUBROUTINE ASSIM_NATURE_ISBA_OI
