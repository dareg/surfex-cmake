!     #########
      SUBROUTINE INIT_IO_SURF_LFI_n(HMASK,HACTION)
!     ######################
!
!!****  *INIT_IO_SURF_LFI* Keep in memory the output files
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      P. Le Moigne 04/2004: distinguish in and out file name
!!      P. Le Moigne 04/2006: special HACTION='GTMSK' to initialize
!!                            a mask different of 'FULL ' in order 
!!                            to read dimensions only.
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE MODD_SURF_PAR,   ONLY: NUNDEF
USE MODD_IO_SURF_LFI,ONLY: CFILE_LFI, CFILEIN_LFI,CFILEOUT_LFI,   &
                           NMASK,CLUOUT_LFI,NFULL,CMASK, NLUOUT,  &
                           NFULL_SURF,                            &
                           NIB, NIE, NJB, NJE, NIU, NJU,          &
                           NIB_SURF, NIE_SURF, NJB_SURF, NJE_SURF,&
                           NIU_SURF, NJU_SURF  
!
USE MODI_GET_LUOUT
USE MODI_READ_SURF
USE MODI_INIT_IO_SURF_MASK_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_DIM_FULL_n
!
USE MODI_GET_TYPE_DIM_n
!
IMPLICIT NONE
!
CHARACTER(LEN=6),  INTENT(IN)  :: HMASK    
CHARACTER(LEN=5),  INTENT(IN)  :: HACTION    
!
INTEGER                        :: ILU,IRET, IL
!
INTEGER                :: INB ! number of articles in the file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_LFI_N',0,ZHOOK_HANDLE)
CALL GET_LUOUT('LFI   ',NLUOUT)
!
IF (HACTION=='GTMSK') THEN
   CALL FMOPEN(CFILEIN_LFI,'OLD',CLUOUT_LFI,0,1,1,INB,IRET)
   CFILE_LFI = CFILEIN_LFI
   CMASK = HMASK
   IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_LFI_N',1,ZHOOK_HANDLE)
   RETURN
ENDIF
!
IF (HACTION == 'READ ') THEN
   CALL FMOPEN(CFILEIN_LFI,'OLD',CLUOUT_LFI,0,1,1,INB,IRET)
   CFILE_LFI = CFILEIN_LFI
   IF (HMASK == 'FULL  ') THEN
      CMASK = HMASK
      CALL READ_SURF('LFI   ','DIM_FULL',ILU,IRET)
      NFULL_SURF = ILU
      NIB_SURF = NIB
      NIE_SURF = NIE
      NJB_SURF = NJB
      NJE_SURF = NJE
      NIU_SURF = NIU
      NJU_SURF = NJU
   ENDIF
   NFULL = NFULL_SURF
ENDIF
!
!
IF (HACTION=='WRITE') THEN
   CALL FMOPEN(CFILEOUT_LFI,'UNKNOWN',CLUOUT_LFI,0,1,1,INB,IRET)
   CFILE_LFI = CFILEOUT_LFI
   CMASK = HMASK
   CALL GET_DIM_FULL_n(NFULL)
ENDIF
!
!*       initialisation of 2D arrays
! 
IF (NIB_SURF/=NUNDEF) THEN
NIB = NIB_SURF
NIE = NIE_SURF
NJB = NJB_SURF
NJE = NJE_SURF
NIU = NIU_SURF
NJU = NJU_SURF
END IF
!------------------------------------------------------------------------------
!
IL = NFULL
CALL GET_TYPE_DIM_n(HMASK,IL)
!
CALL INIT_IO_SURF_MASK_n(HMASK, IL, NLUOUT, NFULL, NMASK)
!
!------------------------------------------------------------------------------
CMASK = HMASK
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_LFI_N',1,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
!
END SUBROUTINE INIT_IO_SURF_LFI_n
