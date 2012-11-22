!     ###############################################################################
SUBROUTINE ASSIM_TEB_n(YPROGRAM,KI,PT2M_O,HTEST)

!     ###############################################################################
!
!!****  *ASSIM_TOWN_n * - Chooses the surface schemes for TOWN parts  
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
USE MODD_CSTS,          ONLY : XPI
USE MODD_SURF_PAR,      ONLY : XUNDEF
USE MODD_SURF_ATM_n,    ONLY : CTOWN
USE MODN_IO_OFFLINE,    ONLY : CPREPFILE
!
#ifdef LFI
USE MODD_IO_SURF_LFI,   ONLY : CFILEIN_LFI, CFILE_LFI,CFILEOUT_LFI
#endif
!
USE YOMHOOK,            ONLY : LHOOK,   DR_HOOK
USE PARKIND1,           ONLY : JPRB
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
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN) :: YPROGRAM  ! program calling surf. schemes
INTEGER,            INTENT(IN) :: KI
REAL,DIMENSION(:), INTENT(IN) :: PT2M_O
CHARACTER(LEN=2),   INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
REAL, DIMENSION (SIZE(PT2M_O)) :: ZTRD3
REAL, DIMENSION (SIZE(PT2M_O)) :: ZT2INC
REAL, DIMENSION (SIZE(PT2M_O)) :: ZTCLS
CHARACTER(LEN=10)    :: YVAR    ! Name of the prognostic variable (in LFI file)
CHARACTER(LEN=100)   :: YPREFIX ! Prefix of the prognostic variable  (in LFI file)
CHARACTER(LEN=3)     :: YREAD
INTEGER              :: IRESP
REAL(KIND=JPRB)      :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('ASSIM_TEB_N',0,ZHOOK_HANDLE)

IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_TEB_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF

WRITE(*,*) 'UPDATING TOWN FOR SCHEME: ',TRIM(CTOWN)

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
CALL INIT_IO_SURF_n(YPROGRAM,'TOWN  ','SURF  ','READ ')
!
!  Read prognostic variables
!
CALL READ_SURF(YPROGRAM,'T_ROAD3',   ZTRD3,  IRESP)
CALL READ_SURF(YPROGRAM,'STORAGETYPE',YREAD,  IRESP)
IF (YREAD=='ALL')THEN
CALL READ_SURF(YPROGRAM,'T2M',       ZTCLS, IRESP)
ENDIF
!
CALL END_IO_SURF_n(YPROGRAM)
CALL IO_BUFF_CLEAN_n

!
! Screen-level innovations
!
ZT2INC(:) = PT2M_O(:) - ZTCLS(:)

PRINT *,'Mean T2m increments  ',SUM(ZT2INC)/KI

!
! c) Temperature analysis of TOWN points
!
WHERE (ZTRD3(:)/=XUNDEF)
  ZTRD3(:) = ZTRD3(:) + ZT2INC(:)/(2.0*XPI)
END WHERE
!

WRITE(*,*) 'Mean T_ROAD3 increments over TOWN ',SUM(ZT2INC)/KI

CFILEOUT_LFI=CPREPFILE
CALL FLAG_UPDATE(.FALSE.,.TRUE.,.FALSE.,.FALSE.)
CALL INIT_IO_SURF_n(YPROGRAM,'TOWN  ','SURF  ','WRITE')

YVAR='T_ROAD3'
YPREFIX='X_Y_T_ROAD3 (K)                                   '
CALL WRITE_SURF(YPROGRAM,YVAR,ZTRD3,IRESP,HCOMMENT=YPREFIX)

CALL END_IO_SURF_n(YPROGRAM)
CALL IO_BUFF_CLEAN_n

IF (LHOOK) CALL DR_HOOK('ASSIM_TEB_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ASSIM_TEB_n
