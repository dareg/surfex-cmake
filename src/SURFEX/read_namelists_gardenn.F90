!     #########
SUBROUTINE READ_NAMELISTS_GARDEN_n(HPROGRAM, HINIT)
!     #######################################################
!
!------------------------------------
!
USE MODN_TEB_GARDEN_n   
!
USE MODI_DEFAULT_ISBA
USE MODI_DEFAULT_CH_DEP
USE MODI_DEFAULT_CH_BIO_FLUX
!
USE MODI_READ_DEFAULT_TEB_GARDEN_n
USE MODI_READ_TEB_GARDEN_CONF_n
!
USE MODI_READ_NAM_PREP_GARDEN_n
!
!------------------------------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=3),   INTENT(IN)  :: HINIT     ! choice of fields to initialize
!
CHARACTER(LEN=3)                  :: HRAIN 
LOGICAL                           :: GCANOPY_DRAG
LOGICAL                           :: GGLACIER
LOGICAL                           :: GTRIP
LOGICAL                           :: GFLOOD
LOGICAL                           :: GVEGUPD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                      
IF (LHOOK) CALL DR_HOOK('READ_NAMELISTS_GARDEN_N',0,ZHOOK_HANDLE)
!
CALL DEFAULT_ISBA(XTSTEP, XOUT_TSTEP,                            &
                     CROUGH,CRUNOFF,CALBEDO,CSCOND,              &
                     CC1DRY, CSOILFRZ, CDIFSFCOND, CSNOWRES,     &
                     CCPSURF, XCGMAX, XCDRAG, CKSAT, CSOM,       &
                     CTOPREG, HRAIN, CHORT, GFLOOD, GTRIP,       &
                     GGLACIER, GCANOPY_DRAG, GVEGUPD             )      
!
CALL DEFAULT_CH_DEP(CCH_DRY_DEP)
CALL DEFAULT_CH_BIO_FLUX(LCH_BIO_FLUX)
!
CALL READ_DEFAULT_TEB_GARDEN_n(HPROGRAM)
!
CALL READ_TEB_GARDEN_CONF_n(HPROGRAM)
!
IF (HINIT=='PRE') CALL READ_NAM_PREP_GARDEN_n(HPROGRAM)
!
IF (LHOOK) CALL DR_HOOK('READ_NAMELISTS_GARDEN_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_NAMELISTS_GARDEN_n
