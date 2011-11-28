!     #########
      SUBROUTINE INIT_IO_SURF_ASC_n(HMASK,HACTION)
!     ######################
!
!!****  *INIT_IO_SURF_ASC* Keep in memory the output files
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
USE MODD_IO_SURF_ASC,ONLY:NUNIT,CFILEIN,CFILEOUT,NMASK,NLUOUT,NFULL,CMASK, &
                          LOPEN_READ
USE MODI_GET_LUOUT
USE MODI_READ_SURF
USE MODI_INIT_IO_SURF_MASK_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_TYPE_DIM_n
!
IMPLICIT NONE
!
CHARACTER(LEN=6),  INTENT(IN)  :: HMASK    
CHARACTER(LEN=5),  INTENT(IN)  :: HACTION    
!
INTEGER                        :: ILU,IRET, IL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_ASC_N',0,ZHOOK_HANDLE)
!
CALL GET_LUOUT('ASCII ',NLUOUT)
!
LOPEN_READ=.FALSE.
!
NUNIT=20
!
IF (HACTION=='GTMSK') THEN
   OPEN(UNIT=NUNIT,FILE=CFILEIN,FORM='FORMATTED')
   CMASK = HMASK
   IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_ASC_N',1,ZHOOK_HANDLE)
   RETURN
ENDIF
!
IF (HACTION == 'READ ') THEN
   OPEN(UNIT=NUNIT,FILE=CFILEIN,FORM='FORMATTED')
   LOPEN_READ=.TRUE. 
   IF (HMASK == 'FULL  ') THEN
      CMASK = HMASK
      CALL READ_SURF('ASCII ','DIM_FULL',ILU,IRET)
      NFULL = ILU
   ENDIF
ELSE
   OPEN(UNIT=NUNIT,FILE=CFILEOUT,FORM='FORMATTED')
ENDIF
!
!------------------------------------------------------------------------------
!
IL = NFULL
CALL GET_TYPE_DIM_n(HMASK,IL)
!
CALL INIT_IO_SURF_MASK_n(HMASK, IL, NLUOUT, NFULL, NMASK)
!
!------------------------------------------------------------------------------
CMASK = HMASK
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_ASC_N',1,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
END SUBROUTINE INIT_IO_SURF_ASC_n
