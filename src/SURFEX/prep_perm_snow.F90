!     #########
SUBROUTINE PREP_PERM_SNOW (IO, IP, IR)
!          ################################################
!
!
!!****  *PREP_PERM_SNOW* - takes into account permanent snow into prognostic snow
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      B. Decharme 03/2009: Consistency with Arpege permanent
!!                                          snow/ice treatment
!!      B. Decharme 07/2012: 3-L or Crocus adjustments
!!      M. Lafaysse 09/2012: adaptation with new snow age in Crocus
!!------------------------------------------------------------------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM, NPROC
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE MODD_TYPE_SNOW
USE MODD_CSTS,           ONLY : XTT
USE MODD_DATA_COVER_PAR, ONLY : NVT_SNOW
USE MODD_SNOW_PAR,       ONLY : XRHOSMAX, XANSMAX, XANSMIN, &
                                XAGLAMAX, XAGLAMIN, XHGLA,  &
                                XRHOSMAX_ES
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODD_ISBA_PAR,       ONLY : XWGMIN
!
USE MODI_VEGTYPE_TO_PATCH
USE MODI_SNOW_HEAT_TO_T_WLIQ
USE MODI_SNOW_T_WLIQ_TO_HEAT
USE MODI_MKFLAG_SNOW
USE MODE_SURF_SNOW_FRAC
USE MODE_SNOW3L
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
!*      0.1    declarations of arguments
!
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
!*      0.2    declarations of local parameter
!
REAL, PARAMETER :: ZRHOL1 = 150.
!
!*      0.3    declarations of local variables
!
INTEGER                             :: JLAYER      ! loop counter on snow layers
REAL, DIMENSION(:),   ALLOCATABLE   :: ZWSNOW_PERM ! snow total reservoir due to perm. snow
REAL, DIMENSION(:),   ALLOCATABLE   :: ZWSNOW      ! initial snow total reservoir
REAL, DIMENSION(:),   ALLOCATABLE   :: ZD          ! new snow total depth
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZDEPTH      ! depth of each layer
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZT          ! new snow temperature profile
REAL, DIMENSION(:),   ALLOCATABLE   :: ZPSN        ! permanent snow fraction
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZWAT        ! 
!
LOGICAL, DIMENSION(:,:), ALLOCATABLE :: GWORK
INTEGER                              :: IWORK
!
#ifdef SFX_MPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
INTEGER :: INFOMPI
INTEGER :: ISNOW          ! patch number where permanent snow is
!
REAL, DIMENSION(0:NPROC-1) :: ZPSN0
REAL :: ZSUM_PSN
REAL              ::ZRHOSMAX
REAL              ::ZAGE_NOW
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!*       1.    Snow where permanent snow is
!              ----------------------------
!
!* snow fraction must be at least equal to permanent snow fraction
!  The snow fraction is computed as Wsnow/(Wsnow+XWCRN)
!
!
IF (LHOOK) CALL DR_HOOK('PREP_PERM_SNOW',0,ZHOOK_HANDLE)
!
ISNOW = VEGTYPE_TO_PATCH(NVT_SNOW,IO%NPATCH)
!
ZRHOSMAX=XRHOSMAX
IF(IR%TSNOW%SCHEME=='3-L'.OR.IR%TSNOW%SCHEME=='CRO')THEN
  ZRHOSMAX=XRHOSMAX_ES
ENDIF
!
ALLOCATE(ZPSN(SIZE(IR%XTG,1)))
ZPSN(:) = MIN ( IP%XVEGTYPE_PATCH(:,NVT_SNOW,ISNOW) , 0.9999 )
!
!* if no permanent snow present
!
ZSUM_PSN = SUM(ZPSN(:))
IF (NPROC>1) THEN
#ifdef SFX_MPI
  CALL MPI_ALLGATHER(ZSUM_PSN,KIND(ZSUM_PSN)/4,MPI_REAL,&
                     ZPSN0,KIND(ZPSN0)/4,MPI_REAL,NCOMM,INFOMPI)
#endif
ELSE
  ZPSN0(:) = ZSUM_PSN
ENDIF

IF (ALL(ZPSN0(:)==0.)) THEN
  DEALLOCATE(ZPSN) 
  IF (LHOOK) CALL DR_HOOK('PREP_PERM_SNOW',1,ZHOOK_HANDLE)
  RETURN
END IF
!
!* total snow amount due to permanent snow
!
ALLOCATE(ZWSNOW_PERM(SIZE(IR%XTG,1)))
ZWSNOW_PERM(:) = WSNOW_FROM_SNOW_FRAC_GROUND(ZPSN)
!
!* limitation of maximum snow amount
!
IF(IO%LGLACIER)THEN
!  limited to 33.3 meters of aged snow
   ZWSNOW_PERM(:) = MIN(ZWSNOW_PERM(:),XHGLA * ZRHOSMAX )
ELSE
!  limited to 2. meters of aged snow
   ZWSNOW_PERM(:) = MIN(ZWSNOW_PERM(:),2.0 * ZRHOSMAX )
ENDIF
!
!* permanent snow can be added only if deep soil temperature is below 5 C
!  (glaciers on subgrid mountains tops that are contained in the grid mesh are neglected)
!
IF(.NOT.IO%LGLACIER)THEN
  WHERE(IR%XTG(:,SIZE(IR%XTG,2),ISNOW)>XTT+5.) ZWSNOW_PERM(:) = 0.
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*       2.    Other parameters of new snow, except temperature
!              ------------------------------------------------
!
!* rho must be defined for snow 3-L before temperature and heat computations
!
ALLOCATE(GWORK(SIZE(IR%XTG,1),IR%TSNOW%NLAYER))
!
DO JLAYER=1,IR%TSNOW%NLAYER
!
  GWORK(:,JLAYER)=.FALSE.
!
  IF(IO%LGLACIER)THEN
      WHERE(ZWSNOW_PERM(:)>0.)GWORK(:,JLAYER)=.TRUE.
  ELSE
      WHERE(ZWSNOW_PERM(:)>0..AND.IR%TSNOW%WSNOW(:,JLAYER,ISNOW)==0.)GWORK(:,JLAYER)=.TRUE.
  ENDIF
!
!* rho
!
  WHERE(GWORK(:,JLAYER))
    IR%TSNOW%RHO(:,JLAYER,ISNOW) = ZRHOSMAX
  END WHERE
!
!* albedo
!
  IF(IO%LGLACIER)THEN
    WHERE(GWORK(:,JLAYER))
         IR%TSNOW%ALB(:,ISNOW) = (XAGLAMAX+XAGLAMIN)/2.0
    END WHERE
  ELSE
    WHERE(GWORK(:,JLAYER))
         IR%TSNOW%ALB(:,ISNOW) = (XANSMAX+XANSMIN)/2.0
    END WHERE
  ENDIF
!
END DO
!
IF (IR%TSNOW%SCHEME=='3-L'.OR.IR%TSNOW%SCHEME=='CRO') THEN
!
! * optimized rho perm snow profile
!
  IF(IO%LGLACIER.AND.IR%TSNOW%NLAYER>=6)THEN
    WHERE(GWORK(:,1))
      IR%TSNOW%RHO(:,1,ISNOW) = ZRHOL1
    END WHERE 
    IF(IR%TSNOW%NLAYER>=6.AND.IR%TSNOW%NLAYER<12)THEN
      WHERE(GWORK(:,2))
       IR%TSNOW%RHO(:,2,ISNOW) = ZRHOL1 + 100.
      END WHERE 
      WHERE(GWORK(:,3))
       IR%TSNOW%RHO(:,3,ISNOW) = ZRHOL1 + 250.
      END WHERE 
    ELSE
      DO JLAYER=2,IR%TSNOW%NLAYER
         WHERE(GWORK(:,JLAYER))
              IR%TSNOW%RHO(:,JLAYER,ISNOW) = MIN(ZRHOSMAX,IR%TSNOW%RHO(:,JLAYER-1,ISNOW)+100.)
         END WHERE     
      ENDDO
    ENDIF
  ENDIF
!
! * Snow age profile
!
  DO JLAYER=1,IR%TSNOW%NLAYER/4
    WHERE(GWORK(:,JLAYER))
      IR%TSNOW%AGE(:,JLAYER,ISNOW) = 365.0*FLOAT(JLAYER-1)/ FLOAT(IR%TSNOW%NLAYER)
    END WHERE
  END DO
  DO JLAYER=1+IR%TSNOW%NLAYER/4,IR%TSNOW%NLAYER
    WHERE(GWORK(:,JLAYER))
      IR%TSNOW%AGE(:,JLAYER,ISNOW) = 3650.*FLOAT(JLAYER-1)/ FLOAT(IR%TSNOW%NLAYER) 
    END WHERE
  END DO
!
  IF(IO%LGLACIER)THEN
    WHERE(GWORK(:,:))IR%TSNOW%AGE(:,:,ISNOW) = 0.0
  ENDIF
!
END IF
!
IF (IR%TSNOW%SCHEME=='CRO') THEN
DO JLAYER=1,IR%TSNOW%NLAYER/4
  WHERE(GWORK(:,JLAYER))
    IR%TSNOW%GRAN1(:,JLAYER,ISNOW) = MIN(-1.,-99.* (1.-4*FLOAT(JLAYER)/FLOAT(IR%TSNOW%NLAYER))) 
    IR%TSNOW%GRAN2(:,JLAYER,ISNOW) = 50. 
    IR%TSNOW%HIST(:,JLAYER,ISNOW) = 0 
  END WHERE
END DO
DO JLAYER=1+IR%TSNOW%NLAYER/4,IR%TSNOW%NLAYER
  WHERE(GWORK(:,JLAYER))
    IR%TSNOW%GRAN1(:,JLAYER,ISNOW) = 99. 
    IR%TSNOW%GRAN2(:,JLAYER,ISNOW) = 0.0003 
    IR%TSNOW%HIST(:,JLAYER,ISNOW) = 0 
  END WHERE
END DO
END IF
!
!-------------------------------------------------------------------------------------
!
!*       3.    Modification of snow reservoir profile
!              --------------------------------------
!
!* initial snow content
!
ALLOCATE(ZWSNOW(SIZE(IR%XTG,1)))
ZWSNOW(:) = 0.
DO JLAYER=1,IR%TSNOW%NLAYER
  ZWSNOW(:) = ZWSNOW(:) + IR%TSNOW%WSNOW(:,JLAYER,ISNOW) 
END DO
!
!* new total snow content
!
ZWSNOW_PERM(:) = MAX(ZWSNOW_PERM(:),ZWSNOW(:))
!
!* new total snow depth
!
ALLOCATE(ZD(SIZE(IR%XTG,1)))
ZD(:) = 0.
DO JLAYER=1,IR%TSNOW%NLAYER
  ZD(:) = ZD(:) + IR%TSNOW%WSNOW(:,JLAYER,ISNOW)/IR%TSNOW%RHO(:,JLAYER,ISNOW)
END DO
ZD(:) = ZD(:) + (ZWSNOW_PERM(:)-ZWSNOW(:))/ZRHOSMAX
!
!* modified snow content profile
!
SELECT CASE(IR%TSNOW%SCHEME)
  CASE('D95','1-L','EBA')
    GWORK(:,1)=.FALSE.
    IF(IO%LGLACIER)THEN
       WHERE(ZWSNOW(:)>=0..AND.IR%TSNOW%WSNOW(:,1,ISNOW)/=XUNDEF)GWORK(:,1)=.TRUE.
    ELSE
       WHERE(ZWSNOW(:)==0..AND.IR%TSNOW%WSNOW(:,1,ISNOW)/=XUNDEF)GWORK(:,1)=.TRUE.
    ENDIF
    WHERE(GWORK(:,1))
      IR%TSNOW%WSNOW(:,1,ISNOW) = ZWSNOW_PERM(:)
    END WHERE
  CASE('3-L','CRO')
    !* grid
    ALLOCATE(ZDEPTH(SIZE(IR%XTG,1),IR%TSNOW%NLAYER))
    CALL SNOW3LGRID(ZDEPTH,ZD)
    DO JLAYER=1,IR%TSNOW%NLAYER
      WHERE(ZWSNOW(:)>= 0. .AND. IR%TSNOW%WSNOW(:,JLAYER,ISNOW)/=XUNDEF)
        IR%TSNOW%WSNOW(:,JLAYER,ISNOW) = ZDEPTH(:,JLAYER) * IR%TSNOW%RHO(:,JLAYER,ISNOW)
      END WHERE
   END DO
   DEALLOCATE(ZDEPTH)

END SELECT
!
DEALLOCATE(ZD)
!-------------------------------------------------------------------------------------
!
!*       4.    Temperature of new snow
!              -----------------------
!
ALLOCATE(ZT(SIZE(IR%TSNOW%WSNOW,1),SIZE(IR%TSNOW%WSNOW,2),SIZE(IR%TSNOW%WSNOW,3)))
!       
SELECT CASE(IR%TSNOW%SCHEME)
  CASE('1-L')
    ZT(:,:,:) = IR%TSNOW%T (:,:,:)
  CASE('3-L','CRO')
    CALL SNOW_HEAT_TO_T_WLIQ(IR%TSNOW%HEAT,IR%TSNOW%RHO,ZT)
END SELECT
!
!* new snow is set to deep ground temperature
!
DO JLAYER=1,IR%TSNOW%NLAYER
!
  GWORK(:,JLAYER)=.FALSE.
!
  IF(IO%LGLACIER)THEN
      WHERE(ZWSNOW_PERM(:)>0.)GWORK(:,JLAYER)=.TRUE.
  ELSE
      WHERE(ZWSNOW_PERM(:)>0. .AND. ZWSNOW(:)==0)GWORK(:,JLAYER)=.TRUE.
  ENDIF
!  
  WHERE(GWORK(:,JLAYER))
      ZT(:,JLAYER,ISNOW) = MIN(IR%XTG(:,SIZE(IR%XTG,2),ISNOW),XTT)
  END WHERE
!
END DO
!
!
SELECT CASE(IR%TSNOW%SCHEME)
  CASE('1-L')
    IR%TSNOW%T (:,:,:) = ZT(:,:,:)
  CASE('3-L','CRO')
    CALL SNOW_T_WLIQ_TO_HEAT(IR%TSNOW%HEAT,IR%TSNOW%RHO,ZT)
END SELECT
!
DEALLOCATE(ZT   )
DEALLOCATE(GWORK)
!
!
!-------------------------------------------------------------------------------------
!
!*       5.    Soil ice initialization for LGLACIER
!              -----------------------
!
ALLOCATE(ZWAT(SIZE(IR%XTG,1),SIZE(IR%XTG,2)))
!
IF(IO%LGLACIER)THEN
!
  IF (IO%CISBA == 'DIF') THEN
      IWORK=IO%NGROUND_LAYER
      ZWAT(:,:)=IP%XWFC(:,:)
  ELSE
      IWORK=2
      ZWAT(:,:)=IP%XWSAT(:,:)
  ENDIF
!
  DO JLAYER=1,IWORK
     WHERE(IP%XVEGTYPE_PATCH(:,NVT_SNOW,ISNOW)>0.0)
           IR%XWGI(:,JLAYER,ISNOW) = MAX(IR%XWGI(:,JLAYER,ISNOW),ZWAT(:,JLAYER)*ZPSN(:))
           IR%XWG (:,JLAYER,ISNOW) = MIN(IR%XWG (:,JLAYER,ISNOW), &
                        MAX(IP%XWSAT(:,JLAYER)-IR%XWGI(:,JLAYER,ISNOW),XWGMIN))
     END WHERE
     WHERE(IR%XWG(:,JLAYER,ISNOW) /= XUNDEF .AND. (IR%XWG(:,JLAYER,ISNOW) + &
                        IR%XWGI(:,JLAYER,ISNOW)) > IP%XWSAT(:,JLAYER) )
           IR%XWGI(:,JLAYER,ISNOW) = IP%XWSAT(:,JLAYER)-IR%XWG (:,JLAYER,ISNOW) !WGT<=WSAT
     END WHERE
  ENDDO
!
ENDIF
!
DEALLOCATE(ZWAT)
DEALLOCATE(ZPSN)
!
!-------------------------------------------------------------------------------------
!
!*       6.    Masking where there is no snow
!              ------------------------------
!
 CALL MKFLAG_SNOW(IR%TSNOW)
IF (LHOOK) CALL DR_HOOK('PREP_PERM_SNOW',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_PERM_SNOW
