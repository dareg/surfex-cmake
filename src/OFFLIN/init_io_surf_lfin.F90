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
!!      V. Masson   *Meteo France*
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
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NINDEX, NPIO, NSIZE
!
USE MODD_SURF_PAR,   ONLY: NUNDEF
!
USE MODD_IO_SURF_LFI,ONLY: CFILE_LFI, CFILEIN_LFI,CFILEOUT_LFI,   &
                           NMASK,CLUOUT_LFI,NFULL,CMASK, NLUOUT,  &
                           NFULL_SURF,                            &
                           NIB, NIE, NJB, NJE, NIU, NJU,          &
                           NIB_SURF, NIE_SURF, NJB_SURF, NJE_SURF,&
                           NIU_SURF, NJU_SURF  
!
USE MODI_GET_LUOUT
USE MODI_READ_SURF
USE MODI_GET_DIM_FULL_n
USE MODI_GET_SIZE_FULL_n
USE MODI_GET_TYPE_DIM_n
USE MODI_INIT_IO_SURF_MASK_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HMASK    
 CHARACTER(LEN=5),  INTENT(IN)  :: HACTION    
!
INTEGER                        :: ILU,IRET, IL
INTEGER                :: INB ! number of articles in the file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_LFI_N',0,ZHOOK_HANDLE)
!
 CALL GET_LUOUT('LFI   ',NLUOUT)
!
!$OMP BARRIER
!
IF (HACTION=='GTMSK') THEN
  IF (NRANK==NPIO) THEN 
!$OMP SINGLE          
    CALL FMOPEN(CFILEIN_LFI,'OLD',CLUOUT_LFI,0,1,1,INB,IRET)
!$OMP END SINGLE    
    CFILE_LFI = CFILEIN_LFI
  ENDIF
  CMASK = HMASK
  IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_LFI_N',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
IF (HACTION == 'READ ') THEN
  IF (NRANK==NPIO) THEN
!$OMP SINGLE            
    CALL FMOPEN(CFILEIN_LFI,'OLD',CLUOUT_LFI,0,1,1,INB,IRET)
!$OMP END SINGLE    
    CFILE_LFI = CFILEIN_LFI
  ENDIF
  CALL READ_SURF(IOB, &
                 'LFI   ','DIM_FULL',NFULL,IRET,HDIR='A')
  IF (HMASK == 'FULL  ') THEN
    NFULL_SURF = NFULL
    NIB_SURF = NIB
    NIE_SURF = NIE
    NJB_SURF = NJB
    NJE_SURF = NJE
    NIU_SURF = NIU
    NJU_SURF = NJU
   ENDIF
ELSE
  CALL GET_DIM_FULL_n(U, &
                      NFULL)
ENDIF
!
!
IF (HACTION=='WRITE' .AND. NRANK==NPIO) THEN
!$OMP SINGLE        
  CALL FMOPEN(CFILEOUT_LFI,'UNKNOWN',CLUOUT_LFI,0,1,1,INB,IRET)
!$OMP END SINGLE   
  CFILE_LFI = CFILEOUT_LFI
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
!
! nindex is needed for call to get_size_full_n. In init_index_mpi, 
! it's not initialized for first readings.
IF (.NOT.ALLOCATED(NINDEX)) THEN
  ALLOCATE(NINDEX(NFULL))
  NINDEX(:) = 0
ENDIF  
!
!------------------------------------------------------------------------------
!
! MASK is sized according to the mpi task running
 CALL GET_SIZE_FULL_n(U, &
                      'LFI   ',NFULL,ILU)
IF (ILU>NSIZE) NSIZE = ILU
!
IL = ILU
 CALL GET_TYPE_DIM_n(DTCO, U, &
                     HMASK,IL)
 CALL INIT_IO_SURF_MASK_n(DTCO, U, &
                          HMASK, IL, NLUOUT, ILU, NMASK)
!
!------------------------------------------------------------------------------
CMASK = HMASK
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_LFI_N',1,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
!
END SUBROUTINE INIT_IO_SURF_LFI_n
