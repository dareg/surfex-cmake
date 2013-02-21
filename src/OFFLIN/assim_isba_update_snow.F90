SUBROUTINE ASSIM_ISBA_UPDATE_SNOW(YPROGRAM, KI, PSWE, PSWE_ORIG, LINITSNOW, LINC, HTEST )

! ------------------------------------------------------------------------------------------
!  *****************************************************************************************
!
!  Routine to update snow field for ISBA
!
!
! ******************************************************************************************
! ------------------------------------------------------------------------------------------
#ifdef LFI
 USE MODD_IO_SURF_LFI, ONLY : CFILEIN_LFI,CFILEOUT_LFI
#endif
 USE MODN_IO_OFFLINE,  ONLY : CPREPFILE
 USE YOMHOOK,          ONLY : LHOOK,DR_HOOK
 USE PARKIND1,         ONLY : JPRB
 USE MODD_CSTS,        ONLY : XTT
 USE MODD_SURF_PAR,    ONLY : XUNDEF
 USE MODD_SNOW_PAR,    ONLY : XANSMIN, XANSMAX, XRHOSMIN, XRHOSMAX
!
 USE MODI_ABOR1_SFX
 USE MODI_INIT_IO_SURF_n 
 USE MODI_READ_SURF 
 USE MODI_END_IO_SURF_n
 USE MODI_IO_BUFF_CLEAN_n
 USE MODI_FLAG_UPDATE
 USE MODI_WRITE_SURF

!
 IMPLICIT NONE
!
 CHARACTER(LEN=6),    INTENT(IN)    :: YPROGRAM  ! program calling surf. schemes
 INTEGER,             INTENT(IN)    :: KI
 REAL, DIMENSION(KI), INTENT(IN)    :: PSWE
 REAL, DIMENSION(KI), INTENT(INOUT) :: PSWE_ORIG
 LOGICAL,             INTENT(IN)    :: LINITSNOW
 LOGICAL,             INTENT(IN)    :: LINC
 CHARACTER(LEN=2),    INTENT(IN)    :: HTEST     ! must be equal to 'OK'
!
!    Declarations of local variables
!
 REAL(KIND=JPRB)                   :: ZHOOK_HANDLE
 CHARACTER(LEN=10)                 :: YVAR    ! Name of the prognostic variable (in LFI file)
 CHARACTER(LEN=100)                :: YPREFIX ! Prefix of the prognostic variable  (in LFI file)
 INTEGER                           :: IRESP
 INTEGER                           :: IVERSION, IBUGFIX
 REAL, DIMENSION(KI)               :: ZSWE     ! Snow before update
 REAL, DIMENSION(KI)               :: ZSWEINC
 REAL, DIMENSION(KI)               :: PTS
!    Addtional snow fields with D95 snow scheme 
 REAL, DIMENSION(KI)               :: ZSNR     ! Snow density 
 REAL, DIMENSION(KI)               :: ZSNA     ! Snow albedo 
!
!
! ----------------------------------------------------------------------------------
!
 IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_UPDATE_SNOW',0,ZHOOK_HANDLE)
 
 IF (HTEST/='OK') THEN
   CALL ABOR1_SFX('ASSIM_ISBA_n: FATAL ERROR DURING ARGUMENT TRANSFER')
 END IF

! If we only do second step, we must set working SWE as input SWE
 ZSWE=PSWE

 IF ( LINITSNOW ) THEN
   !   read initial snow before update
#ifdef LFI
   CFILEIN_LFI = CPREPFILE        ! input PREP file (surface fields)
#endif
   print*,CFILEIN_LFI
   CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','SURF  ','READ ')
   CALL READ_SURF(YPROGRAM,'VERSION',IVERSION,IRESP)
   CALL READ_SURF(YPROGRAM,'BUG',IBUGFIX,IRESP)
   IF (IVERSION<7 .OR. IVERSION==7 .AND. IBUGFIX<3) THEN
     CALL READ_SURF(YPROGRAM,'WSNOW_VEG1',ZSWE,  IRESP)
   ELSE
     CALL READ_SURF(YPROGRAM,'WSN_VEG1',ZSWE,  IRESP)
   ENDIF
   CALL READ_SURF(YPROGRAM,'TG1',       PTS,   IRESP)
   CALL END_IO_SURF_n(YPROGRAM)
   CALL IO_BUFF_CLEAN_n

   ! Set snow=0 where 1. guess = 0 and Ts>0, to avoid that the snow analysis introduce snow where it is no snow.
   WHERE ( PSWE(:)/=XUNDEF .AND. PSWE(:)<1.0E-10 .AND. PTS(:)>XTT )
      ZSWE(:)   = 0.0
   END WHERE

   ZSWE=PSWE
   ! Set snow=0 where 1. guess = 0 and Ts>0, to avoid that the snow analysis introduce snow where it is no snow.
   WHERE ( PSWE(:)/=XUNDEF .AND. PSWE(:)<1.0E-10 .AND. PTS(:)>XTT )
      ZSWE(:)   = 0.0
   END WHERE
 ENDIF

 ! Update snow
 IF ( LINC ) THEN
 
   ! Calculate increments
   ZSWEINC(:) = ZSWE(:) - PSWE_ORIG(:)
   WRITE(*,'("  SURFRESERV.NEIGE - min, mean, max: ",3E13.4)') MINVAL(ZSWE),MAXVAL(ZSWE),SUM(ZSWE)/KI
   WRITE(*,*) 'Mean SN increments over NATURE ',SUM(ZSWEINC)/KI
 
   !   read additional snow fields 
#ifdef LFI 
   CFILEIN_LFI = CPREPFILE        ! input PREP file (surface fields) 
#endif 
print*,CFILEIN_LFI
   CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','SURF  ','READ ') 
   CALL READ_SURF(YPROGRAM,'VERSION',IVERSION,IRESP)
   CALL READ_SURF(YPROGRAM,'BUG',IBUGFIX,IRESP)
   IF (IVERSION<7 .OR. IVERSION==7 .AND. IBUGFIX<3) THEN
     CALL READ_SURF(YPROGRAM,'RSNOW_VEG1',ZSNR,  IRESP)
     CALL READ_SURF(YPROGRAM,'ASNOW_VEG',ZSNA,  IRESP)
   ELSE
     CALL READ_SURF(YPROGRAM,'RSN_VEG1',ZSNR,  IRESP)
     CALL READ_SURF(YPROGRAM,'ASN_VEG',ZSNA,  IRESP)
   ENDIF   
   CALL END_IO_SURF_n(YPROGRAM) 
   CALL IO_BUFF_CLEAN_n 
 
   ! Snow albedo and density are given initial values in points  
   ! which get initial snow in the snow analysis 
   WHERE ( PSWE_ORIG(:) < 1.0E-10 .AND. ZSWE(:)>= 1.0E-10 ) 
   ZSNA(:)    = 0.5 * ( XANSMIN + XANSMAX ) 
   ZSNR(:)    = 0.5 * ( XRHOSMIN + XRHOSMAX ) 
   END WHERE 
   
 ENDIF

 ! Write updated snow
 !
#ifdef LFI
 CFILEOUT_LFI=CPREPFILE
#endif
 CALL FLAG_UPDATE(.FALSE.,.TRUE.,.FALSE.,.FALSE.)
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','SURF  ','WRITE')

 YVAR='WSN_VEG1'
 YPREFIX='X_Y_WSN_VEG1 (kg/m2)'
 CALL WRITE_SURF(YPROGRAM,YVAR,ZSWE,IRESP,HCOMMENT=YPREFIX)
 
 IF ( LINC ) THEN          
   YVAR='RSN_VEG1' 
   YPREFIX='X_Y_RSNOW_VEG1 (kg/m3)                            ' 
   CALL WRITE_SURF(YPROGRAM,YVAR,ZSNR,IRESP,HCOMMENT=YPREFIX) 
   YVAR='ASN_VEG' 
   YPREFIX='X_Y_ASNOW_VEG1 (%)                            ' 
   CALL WRITE_SURF(YPROGRAM,YVAR,ZSNA,IRESP,HCOMMENT=YPREFIX) 
 ENDIF 
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n

!
! -------------------------------------------------------------------------------------
 IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_UPDATE_SNOW',1,ZHOOK_HANDLE)
 END SUBROUTINE ASSIM_ISBA_UPDATE_SNOW
