!###################
MODULE MODE_COUPLING_VAR_SFX_TRIP
!###################
!
!!****  *MODE_COUPLING_VAR_SFX_TRIP*
!!
!!    PURPOSE
!!    -------
!    
!      The purpose of this routine is to store here all routines 
!      used to get or to put each variable used in the coupling of SFX - TRIP.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!       NONE          
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	B. Decharme       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/05/08
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  INTERFACE GET_COUPLING_VAR_TRIP_n
      MODULE PROCEDURE GET_COUPLING_VAR_TRIP_n
  END INTERFACE

  INTERFACE GET_COUPLING_VAR_SFX_n
      MODULE PROCEDURE GET_COUPLING_VAR_SFX_n
  END INTERFACE
  
  INTERFACE PUT_COUPLING_VAR_SFX_n
      MODULE PROCEDURE PUT_COUPLING_VAR_SFX_n
  END INTERFACE
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
!     ####################################################
      SUBROUTINE GET_COUPLING_VAR_TRIP_n(PFFLOOD,PPIFLOOD)
!     ####################################################
!
!!    PURPOSE
!!    -------
!    
!     Get TRIP - ISBA coupling variables
!
USE MODD_TRIP_GRID_n, ONLY : GMASK
USE MODD_TRIP_n,      ONLY : XFFLOOD,XPIFLOOD
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      declarations of arguments
!
REAL, DIMENSION(:), INTENT(OUT)  :: PFFLOOD
REAL, DIMENSION(:), INTENT(OUT)  :: PPIFLOOD
!
!*      declarations of local variables
!
INTEGER I,J,ICOUNT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      procedure
!
IF (LHOOK) CALL DR_HOOK('MODE_COUPLING_VAR_SFX_TRIP:GET_COUPLING_VAR_TRIP_n',0,ZHOOK_HANDLE)
ICOUNT=0
!
DO J=1,SIZE(GMASK,2)
   DO I=1,SIZE(GMASK,1)
      ICOUNT=ICOUNT+1
      IF(.NOT.GMASK(I,J))CYCLE
      PFFLOOD (ICOUNT)=XFFLOOD (I,J)
      PPIFLOOD(ICOUNT)=XPIFLOOD(I,J)
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_COUPLING_VAR_SFX_TRIP:GET_COUPLING_VAR_TRIP_n',1,ZHOOK_HANDLE)
END SUBROUTINE GET_COUPLING_VAR_TRIP_n
!
!-------------------------------------------------------------------------------
!
!     ###################################################################
      SUBROUTINE GET_COUPLING_VAR_SFX_n(OFLOOD,PDRAIN,PRUNOFF,PSRC_FLOOD)
!     ###################################################################
!
!!    PURPOSE
!!    -------
!    
!     Get TRIP - ISBA coupling variables
!     Put ISBA variables in TRIP dimension (kg/m² --> kg)
!
USE MODD_SURF_PAR,  ONLY : XUNDEF
!
USE MODD_ISBA_n, ONLY : XCPL_DRAIN,XCPL_RUNOFF,XCPL_EFLOOD,  &
                        XCPL_PFLOOD,XCPL_IFLOOD,XCPL_ICEFLUX,&
                        LGLACIER, TSNOW
!
USE MODD_SURF_ATM_n,      ONLY : NR_NATURE, XNATURE
USE MODD_SURF_ATM_GRID_n, ONLY : XMESH_SIZE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      declarations of arguments
!
LOGICAL,            INTENT(IN )           :: OFLOOD
!
REAL, DIMENSION(:), INTENT(OUT)           :: PDRAIN
REAL, DIMENSION(:), INTENT(OUT)           :: PRUNOFF
REAL, DIMENSION(:), INTENT(OUT)           :: PSRC_FLOOD
!
!*      declarations of local variables
!
INTEGER :: II, JI, INATURE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! Initialize
!
IF (LHOOK) CALL DR_HOOK('MODE_COUPLING_VAR_SFX_TRIP:GET_COUPLING_VAR_SFX_n',0,ZHOOK_HANDLE)
INATURE=SIZE(XCPL_DRAIN)
!
PDRAIN    (:) = 0.0
PRUNOFF   (:) = 0.0
PSRC_FLOOD(:) = 0.0
!
! Some calculation
!
IF(LGLACIER.AND.TSNOW%SCHEME=='D95')THEN
   XCPL_RUNOFF(:) = XCPL_RUNOFF(:) + XCPL_ICEFLUX(:)     
ENDIF
!
! Get variable over nature tile to global field
!
DO JI=1,INATURE
   II = NR_NATURE(JI)
   PDRAIN    (II) = XCPL_DRAIN (JI)
   PRUNOFF   (II) = XCPL_RUNOFF(JI)
ENDDO
!
! kg/m2 -> kg
!
PDRAIN    (:) = PDRAIN    (:) * XMESH_SIZE(:) * XNATURE(:)
PRUNOFF   (:) = PRUNOFF   (:) * XMESH_SIZE(:) * XNATURE(:)
!
! re-initialize cumulative field
!
XCPL_DRAIN  (:) = 0.0
XCPL_RUNOFF (:) = 0.0
XCPL_ICEFLUX(:) = 0.0
!
! Floodplains
!
IF(OFLOOD)THEN
!
  DO JI=1,INATURE
     II = NR_NATURE(JI)
     PSRC_FLOOD(II) = XCPL_PFLOOD(JI)-XCPL_IFLOOD(JI)-XCPL_EFLOOD(JI)
  ENDDO 
!        
  PSRC_FLOOD (:) = PSRC_FLOOD(:) * XMESH_SIZE(:) * XNATURE(:)
!        
  XCPL_PFLOOD(:) = 0.0
  XCPL_IFLOOD(:) = 0.0
  XCPL_EFLOOD(:) = 0.0
!  
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_COUPLING_VAR_SFX_TRIP:GET_COUPLING_VAR_SFX_n',1,ZHOOK_HANDLE)
END SUBROUTINE GET_COUPLING_VAR_SFX_n
!
!-------------------------------------------------------------------------------
!
!     ##########################################################
      SUBROUTINE PUT_COUPLING_VAR_SFX_n(PFFLOOD,PPIFLOOD,PTSTEP)
!     ##########################################################
!
!!    PURPOSE
!!    -------
!    
!     Get TRIP - ISBA coupling variables
!
USE MODD_SURF_ATM_GRID_n, ONLY : XMESH_SIZE
!
USE MODD_SURF_ATM_n,  ONLY : NR_NATURE, XNATURE
!
USE MODD_ISBA_n,      ONLY : XFFLOOD,XPIFLOOD,XTSTEP_COUPLING
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN) :: PFFLOOD
REAL, DIMENSION(:), INTENT(IN) :: PPIFLOOD
!
REAL, INTENT(IN), OPTIONAL     :: PTSTEP
!
INTEGER :: II, JI, INATURE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      declarations of local variables
!
!*      procedure
!
! Initialize
!
IF (LHOOK) CALL DR_HOOK('MODE_COUPLING_VAR_SFX_TRIP:PUT_COUPLING_VAR_SFX_n',0,ZHOOK_HANDLE)
IF(PRESENT(PTSTEP)) XTSTEP_COUPLING = PTSTEP
!
INATURE=SIZE(XFFLOOD)
!
! Get variable over global field to nature tile kg -> kg/m²/s
!
DO JI=1,INATURE
  II = NR_NATURE(JI)  
  XFFLOOD (JI) = PFFLOOD (II) * XNATURE(II)
  XPIFLOOD(JI) = PPIFLOOD(II) /(XTSTEP_COUPLING*XMESH_SIZE(II)*XNATURE(II))
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_COUPLING_VAR_SFX_TRIP:PUT_COUPLING_VAR_SFX_n',1,ZHOOK_HANDLE)
END SUBROUTINE PUT_COUPLING_VAR_SFX_n
!
!-------------------------------------------------------------------------------
!
END MODULE MODE_COUPLING_VAR_SFX_TRIP      
