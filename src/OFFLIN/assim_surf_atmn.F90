!     #################################################################################
SUBROUTINE ASSIM_SURF_ATM_n(HPROGRAM,KI,                                                &
                            PCON_RAIN, PSTRAT_RAIN, PCON_SNOW, PSTRAT_SNOW,             &
                            PCLOUDS,   PLSM,        PEVAPTR,   PEVAP,                   &
                            PSWEC,     PTSC,                                            &
                            PTS,       PT2M,        PHU2M,     PSWE,                    &
                            HTEST )
!     #################################################################################
!
!
!!****  *ASSIM_SURF_ATM_n * - Driver to call the schemes for the 
!!       four surface types (SEA, WATER, NATURE, TOWN)
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
!!-------------------------------------------------------------
!
USE MODD_SURF_CONF,      ONLY : CPROGNAME
USE MODD_SURF_ATM_n,     ONLY : NDIM_FULL,                             &
                                NSIZE_SEA, NSIZE_WATER, NSIZE_TOWN, NSIZE_NATURE, &
                                NR_SEA,    NR_WATER,    NR_TOWN,    NR_NATURE
!
USE MODD_ASSIM,          ONLY : LREAD_SST_FROM_FILE
!
USE YOMHOOK,             ONLY : LHOOK,   DR_HOOK
USE PARKIND1,            ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_ASSIM_READ_SST_FROM_FILE
USE MODI_ASSIM_SEA_n
USE MODI_ASSIM_INLAND_WATER_n
USE MODI_ASSIM_NATURE_n
USE MODI_ASSIM_TOWN_n
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=6),    INTENT(IN) :: HPROGRAM     ! program calling surf. schemes
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
REAL, DIMENSION(KI), INTENT(INOUT) :: PTS
REAL, DIMENSION(KI), INTENT(IN) :: PT2M
REAL, DIMENSION(KI), INTENT(IN) :: PHU2M
REAL, DIMENSION(KI), INTENT(IN) :: PSWE
 CHARACTER(LEN=2),   INTENT(IN)  :: HTEST        ! must be equal to 'OK'

!
!*      0.2    declarations of local variables
!
INTEGER                              :: JTILE                        ! loop on type of surface
LOGICAL                              :: GNATURE, GTOWN, GWATER, GSEA ! .T. if the corresponding surface is represented
REAL(KIND=JPRB)                      :: ZHOOK_HANDLE

!
!-------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ASSIM_SURF_ATM_N',0,ZHOOK_HANDLE)
CPROGNAME=HPROGRAM
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_SURF_ATMN: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
!
!-------------------------------------------------------------------------------------
! Preliminaries: Tile related operations
!-------------------------------------------------------------------------------------

! FLAGS for the various surfaces:
!
GSEA      = NSIZE_SEA    >0
GWATER    = NSIZE_WATER  >0
GTOWN     = NSIZE_TOWN   >0
GNATURE   = NSIZE_NATURE >0
!
! Tile counter:
!
JTILE     = 0 
!
!
!--------------------------------------------------------------------------------------
! Call interfaces for sea, water, nature and town here...
!--------------------------------------------------------------------------------------
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! SEA Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE = JTILE + 1
!
IF(GSEA)THEN
!
  IF ( LREAD_SST_FROM_FILE ) CALL ASSIM_READ_SST_FROM_FILE(HPROGRAM,NDIM_FULL,PTS,PLSM,HTEST)
  CALL ASSIM_TREAT_SURF(JTILE,NSIZE_SEA,NR_SEA)
!
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! INLAND WATER Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE = JTILE + 1
!
IF(GWATER)THEN
!
  CALL ASSIM_TREAT_SURF(JTILE,NSIZE_WATER,NR_WATER)
!
ENDIF 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! NATURAL SURFACE Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE = JTILE + 1
!
IF(GNATURE)THEN
!
  CALL ASSIM_TREAT_SURF(JTILE,NSIZE_NATURE,NR_NATURE)
!
ENDIF 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! URBAN Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE = JTILE + 1
!
IF(GTOWN)THEN
!
  CALL ASSIM_TREAT_SURF(JTILE,NSIZE_TOWN,NR_TOWN)
!
ENDIF

IF (LHOOK) CALL DR_HOOK('ASSIM_SURF_ATM_N',1,ZHOOK_HANDLE)
!
!=======================================================================================
CONTAINS
!=======================================================================================
SUBROUTINE ASSIM_TREAT_SURF(KTILE,KSIZE,KMASK)
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                   :: KTILE
INTEGER, INTENT(IN)                   :: KSIZE
INTEGER, INTENT(IN), DIMENSION(KSIZE) :: KMASK
INTEGER                               :: JJ,JI
REAL,DIMENSION(KSIZE)                 :: ZP_PCON_RAIN
REAL,DIMENSION(KSIZE)                 :: ZP_PSTRAT_RAIN
REAL,DIMENSION(KSIZE)                 :: ZP_PCON_SNOW
REAL,DIMENSION(KSIZE)                 :: ZP_PSTRAT_SNOW
REAL,DIMENSION(KSIZE)                 :: ZP_PCLOUDS
REAL,DIMENSION(KSIZE)                 :: ZP_PLSM
REAL,DIMENSION(KSIZE)                 :: ZP_PEVAPTR
REAL,DIMENSION(KSIZE)                 :: ZP_PEVAP
REAL,DIMENSION(KSIZE)                 :: ZP_PSWEC
REAL,DIMENSION(KSIZE)                 :: ZP_PTSC
REAL,DIMENSION(KSIZE)                 :: ZP_PTS
REAL,DIMENSION(KSIZE)                 :: ZP_PT2M
REAL,DIMENSION(KSIZE)                 :: ZP_PHU2M
REAL,DIMENSION(KSIZE)                 :: ZP_PSWE

DO JJ=1,KSIZE
  JI=KMASK(JJ)
  ZP_PCON_RAIN(JJ)   = PCON_RAIN(JI)
  ZP_PSTRAT_RAIN(JJ) = PSTRAT_RAIN(JI)
  ZP_PCON_SNOW(JJ)   = PCON_SNOW(JI)
  ZP_PSTRAT_SNOW(JJ) = PSTRAT_SNOW(JI)
  ZP_PCLOUDS(JJ)     = PCLOUDS(JI)
  ZP_PLSM(JJ)        = PLSM(JI)
  ZP_PEVAPTR(JJ)     = PEVAPTR(JI)
  ZP_PEVAP(JJ)       = PEVAP(JI)
  ZP_PSWEC(JJ)       = PSWEC(JI)
  ZP_PTSC(JJ)        = PTSC(JI)
  ZP_PTS(JJ)         = PTS(JI) 
  ZP_PT2M(JJ)        = PT2M(JI)
  ZP_PHU2M(JJ)       = PHU2M(JI)
  ZP_PSWE(JJ)        = PSWE(JI)
ENDDO

IF (KTILE==1) THEN
 
  WRITE(*,*) '*********************************************'
  WRITE(*,*) '*      ASSIMILATIONS FOR SEA POINTS         *'
  WRITE(*,*) '*********************************************'
 
  CALL ASSIM_SEA_n(HPROGRAM,KSIZE,ZP_PTS,ZP_PLSM,HTEST)

ELSEIF (KTILE==2) THEN
  
  WRITE(*,*) '*********************************************'  
  WRITE(*,*) '*      ASSIMILATIONS FOR WATER POINTS       *'
  WRITE(*,*) '*********************************************'
  CALL ASSIM_INLAND_WATER_n(HPROGRAM,KSIZE,ZP_PTS,ZP_PLSM,HTEST)

ELSEIF (KTILE==3) THEN
  
  WRITE(*,*) '*********************************************'  
  WRITE(*,*) '*      ASSIMILATIONS FOR NATURE POINTS      *'
  WRITE(*,*) '*********************************************'
  CALL ASSIM_NATURE_n(HPROGRAM,KSIZE,                                             &
                      ZP_PCON_RAIN, ZP_PSTRAT_RAIN, ZP_PCON_SNOW, ZP_PSTRAT_SNOW, &
                      ZP_PCLOUDS,   ZP_PLSM,        ZP_PEVAPTR,   ZP_PEVAP,       & 
                      ZP_PSWEC,     ZP_PTSC,                                      &
                      ZP_PTS,       ZP_PT2M,        ZP_PHU2M,     ZP_PSWE,        & 
                      HTEST )
  
ELSEIF (KTILE==4) THEN
  
  WRITE(*,*) '*********************************************'  
  WRITE(*,*) '*      ASSIMILATIONS FOR URBAN POINTS       *'
  WRITE(*,*) '*********************************************'
  CALL ASSIM_TOWN_n(HPROGRAM,KSIZE,ZP_PT2M,HTEST)
  
ENDIF

END SUBROUTINE ASSIM_TREAT_SURF
!=======================================================================================
END SUBROUTINE ASSIM_SURF_ATM_n
!=======================================================================================

