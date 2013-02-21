!     #########################
      SUBROUTINE ARRANGE_COVER
!     #########################
!
!!**** *ARRANGE_COVER*
!!
!!    PURPOSE
!!    -------
!!
!!    change water and intertidal (not lake) to nature and/or town to rock : arrange cover properly
!!
!!    METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    B. Decharme        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    03/2009
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODD_SURF_ATM_n,     ONLY : LWATER_TO_NATURE, LTOWN_TO_ROCK
!
USE MODD_DATA_COVER,     ONLY : XDATA_TOWN, XDATA_NATURE, XDATA_WATER, XDATA_SEA, &
                                XDATA_ROOT_DEPTH, XDATA_GROUND_DEPTH, XDATA_DICE, &
                                XDATA_VEGTYPE, XDATA_LAI, XDATA_LAI_ALL_YEARS,    &
                                XDATA_GARDEN, XDATA_ALB_VEG_NIR, XDATA_ALB_VEG_VIS, &
                                XDATA_ALB_SOIL_NIR, XDATA_ALB_SOIL_VIS   
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE, JPCOVER, NROCK, NWATER, NVT_ROCK
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL     :: ZWORK
!
INTEGER  :: JCOVER, JVEGTYPE, JL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
! Change water (not lake) to nature
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ARRANGE_COVER',0,ZHOOK_HANDLE)
IF(LWATER_TO_NATURE)THEN
  DO JCOVER=1,JPCOVER
     IF(JCOVER/=NWATER(1).AND.JCOVER/=NWATER(2).AND.JCOVER/=NWATER(3).AND.XDATA_WATER(JCOVER)>0.0)THEN
       XDATA_NATURE(JCOVER)=XDATA_NATURE(JCOVER)+XDATA_WATER(JCOVER)
       XDATA_WATER (JCOVER)=0.0
     ENDIF
     !Only cover 242
     IF(XDATA_SEA(JCOVER)>0.0.AND.XDATA_SEA(JCOVER)<1.0)THEN
           XDATA_NATURE(JCOVER)=XDATA_NATURE(JCOVER)+XDATA_SEA(JCOVER)
           XDATA_SEA(JCOVER)=0.0
     ENDIF
  ENDDO
ENDIF
!
!-------------------------------------------------------------------------------
! Change town to rock but keep other natural fraction
!-------------------------------------------------------------------------------
!
IF(LTOWN_TO_ROCK)THEN
!        
  DO JCOVER=1,JPCOVER
     IF(XDATA_TOWN(JCOVER)>0.0.OR.XDATA_GARDEN(JCOVER)>0.0)THEN
!
       XDATA_NATURE(JCOVER) = XDATA_NATURE(JCOVER) + XDATA_GARDEN(JCOVER) * XDATA_TOWN(JCOVER)
       XDATA_TOWN  (JCOVER) = XDATA_TOWN  (JCOVER) * ( 1. - XDATA_GARDEN(JCOVER))
       XDATA_GARDEN(JCOVER) = 0.0
!
       ZWORK=XDATA_NATURE(JCOVER)+XDATA_TOWN(JCOVER)
!       
       DO JVEGTYPE=1,NVEGTYPE
             XDATA_VEGTYPE(JCOVER,JVEGTYPE)=XDATA_VEGTYPE(JCOVER,JVEGTYPE)*XDATA_NATURE(JCOVER)/ZWORK
       ENDDO
!      
       XDATA_VEGTYPE(JCOVER,NVT_ROCK) = XDATA_VEGTYPE(JCOVER,NVT_ROCK)+XDATA_TOWN(JCOVER)/ZWORK
!
       XDATA_NATURE(JCOVER)=XDATA_NATURE(JCOVER)+XDATA_TOWN(JCOVER)
!       
       XDATA_TOWN  (JCOVER)=0.0
!
!      Initialise some variables
       XDATA_LAI          (JCOVER,:,NVT_ROCK) = 0.0
       XDATA_LAI_ALL_YEARS(JCOVER,:,NVT_ROCK) = 0.0
       XDATA_ROOT_DEPTH   (JCOVER,  NVT_ROCK) = 0.2
       XDATA_GROUND_DEPTH (JCOVER,  NVT_ROCK) = 0.2
       XDATA_DICE         (JCOVER,  NVT_ROCK) = 0.2
       XDATA_ALB_VEG_NIR  (JCOVER,:,NVT_ROCK) = 0.3
       XDATA_ALB_VEG_VIS  (JCOVER,:,NVT_ROCK) = 0.1
       XDATA_ALB_SOIL_NIR (JCOVER,:,NVT_ROCK) = 0.3
       XDATA_ALB_SOIL_VIS (JCOVER,:,NVT_ROCK) = 0.1
!
     ENDIF
  ENDDO
!
ENDIF        
IF (LHOOK) CALL DR_HOOK('ARRANGE_COVER',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ARRANGE_COVER
