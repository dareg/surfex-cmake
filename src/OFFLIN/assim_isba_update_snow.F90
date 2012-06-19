SUBROUTINE ASSIM_ISBA_UPDATE_SNOW(HPROGRAM, KI, PSWE, HTEST )

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
 CHARACTER(LEN=6),    INTENT(IN) :: HPROGRAM  ! program calling surf. schemes
 INTEGER,             INTENT(IN) :: KI
 REAL, DIMENSION(KI), INTENT(IN) :: PSWE
 CHARACTER(LEN=2),    INTENT(IN) :: HTEST     ! must be equal to 'OK'
!
!    Declarations of local variables
!
 CHARACTER(LEN=6)                :: YPROGRAM  = 'LFI   '
 REAL(KIND=JPRB)                 :: ZHOOK_HANDLE
 CHARACTER(LEN=10)               :: YVAR    ! Name of the prognostic variable (in LFI file)
 CHARACTER(LEN=100)              :: YPREFIX ! Prefix of the prognostic variable  (in LFI file)
 INTEGER                         :: IRESP
 REAL, DIMENSION(KI)             :: ZSWE     ! Snow before update
 REAL, DIMENSION(KI)             :: ZSWEINC  ! Snow increment
!
!
! ----------------------------------------------------------------------------------
!
 IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_UPDATE_SNOW',0,ZHOOK_HANDLE)
 
 IF (HTEST/='OK') THEN
   CALL ABOR1_SFX('ASSIM_ISBA_n: FATAL ERROR DURING ARGUMENT TRANSFER')
 END IF

 !   read initial snow before update
#ifdef LFI
 CFILEIN_LFI = CPREPFILE        ! input PREP file (surface fields)
#endif
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','SURF  ','READ ')
 CALL READ_SURF(YPROGRAM,'WSNOW_VEG1',ZSWE,  IRESP)
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n

 ! Update snow
 ZSWEINC=0.
 ZSWEINC(:) = PSWE(:) - ZSWE(:)

 WRITE(*,'("  SURFRESERV.NEIGE - min, mean, max: ",3E13.4)') MINVAL(PSWE),MAXVAL(PSWE),SUM(PSWE)/KI

 WRITE(*,*) 'Mean SN increments over NATURE ',SUM(ZSWEINC)/KI

 ! Write updated snow
 !
#ifdef LFI
 CFILEOUT_LFI=CPREPFILE
#endif
 CALL FLAG_UPDATE(.FALSE.,.TRUE.,.FALSE.,.FALSE.)
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','SURF  ','WRITE')

 YVAR='WSNOW_VEG1'
 YPREFIX='X_Y_WSNOW_VEG1 (kg/m2)                            '
 CALL WRITE_SURF(YPROGRAM,YVAR,PSWE,IRESP,HCOMMENT=YPREFIX)

 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n

!
! -------------------------------------------------------------------------------------
 IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_UPDATE_SNOW',1,ZHOOK_HANDLE)
 END SUBROUTINE ASSIM_ISBA_UPDATE_SNOW
