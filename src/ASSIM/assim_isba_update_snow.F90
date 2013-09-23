SUBROUTINE ASSIM_ISBA_UPDATE_SNOW(YPROGRAM, KI, PSWE, PSWE_ORIG, LINITSNOW, LINC, HTEST )

! ------------------------------------------------------------------------------------------
!  *****************************************************************************************
!
!  Routine to update snow field for ISBA
!  Trygve Aspelien, Separating IO  06/2013
!
!
! ******************************************************************************************
! ------------------------------------------------------------------------------------------
USE MODN_IO_OFFLINE,  ONLY : CPREPFILE
USE YOMHOOK,          ONLY : LHOOK,DR_HOOK
USE PARKIND1,         ONLY : JPRB
USE MODD_CSTS,        ONLY : XTT
USE MODD_ISBA_n,      ONLY : XTG,TSNOW
USE MODD_SURF_PAR,    ONLY : XUNDEF
USE MODD_SNOW_PAR,    ONLY : XANSMIN, XANSMAX, XRHOSMIN, XRHOSMAX

USE MODI_ABOR1_SFX


IMPLICIT NONE
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
INTEGER                           :: NLAYER,NPATCH
REAL, DIMENSION(KI)               :: ZSWE     ! Snow before update
REAL, DIMENSION(KI)               :: ZSWEINC
REAL, DIMENSION(KI)               :: ZTS
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

IF ( TSNOW%SCHEME=='D95' ) THEN
  NLAYER=1
  NPATCH=1
ELSE
  CALL ABOR1_SFX("Update of snow is only implemented for D95")
ENDIF

IF ( LINITSNOW ) THEN
  ZSWE(:)=TSNOW%WSNOW(:,NLAYER,NPATCH)
  ZTS(:)=XTG(:,1,NPATCH)
  PSWE_ORIG(:)=ZSWE(:)

  ZSWE=PSWE
  ! Set snow=0 where 1. guess = 0 and Ts>0, to avoid that the snow analysis introduce snow where it is no snow.
  WHERE ( PSWE(:)/=XUNDEF .AND. PSWE(:)<1.0E-10 .AND. ZTS(:)>XTT )
     ZSWE(:)   = 0.0
  END WHERE
  TSNOW%WSNOW(:,NLAYER,NPATCH)=ZSWE
ENDIF


! Update snow
IF ( LINC ) THEN

  ZSWE(:)=TSNOW%WSNOW(:,NLAYER,NPATCH)  
  ZSNA(:)=TSNOW%ALB(:,NPATCH)
  ZSNR(:)=TSNOW%RHO(:,NLAYER,NPATCH)

  ! If we only do second step, we must set working SWE as input SWE
  IF ( .NOT. LINITSNOW ) THEN
    ZSWE(:)=PSWE(:)
  ENDIF
 
  ! Calculate increments
  ZSWEINC(:) = ZSWE(:) - PSWE_ORIG(:)
  WRITE(*,'("  SURFRESERV.NEIGE - min, mean, max: ",3E13.4)') MINVAL(ZSWE),MAXVAL(ZSWE),SUM(ZSWE)/KI
  WRITE(*,*) 'Mean SN increments over NATURE ',SUM(ZSWEINC)/KI

  ! Snow albedo and density are given initial values in points  
  ! which get initial snow in the snow analysis 
  WHERE ( PSWE_ORIG(:) < 1.0E-10 .AND. ZSWE(:)>= 1.0E-10 ) 
    ZSNA(:)    = 0.5 * ( XANSMIN + XANSMAX ) 
    ZSNR(:)    = 0.5 * ( XRHOSMIN + XRHOSMAX ) 
  END WHERE 
  TSNOW%ALB(:,NPATCH)=ZSNA(:)
  TSNOW%RHO(:,NLAYER,NPATCH)=ZSNR(:)
  TSNOW%WSNOW(:,NLAYER,NPATCH)=ZSWE(:)
ENDIF

!
! -------------------------------------------------------------------------------------
 IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_UPDATE_SNOW',1,ZHOOK_HANDLE)
 END SUBROUTINE ASSIM_ISBA_UPDATE_SNOW
