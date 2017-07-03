!     #########
      SUBROUTINE TRIP_DIAG_RUN (TPDG, TPG, &
                                KLISTING,KLON,KLAT,PRUNTIME)  
!     #####################################################
!
!!****  *TRIP_DIAG_RUN*  
!!
!!    PURPOSE
!!    -------
!
!     TRIP river routing run mean outputs.
!     
!!      
!!    AUTHOR
!!    ------
!!      B. Decharme
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/05/05 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_TRIP_DIAG, ONLY : TRIP_DIAG_t
USE MODD_TRIP_GRID, ONLY : TRIP_GRID_t
!
USE MODN_TRIP,       ONLY : CGROUNDW, CVIT, LFLOOD
USE MODN_TRIP_RUN,   ONLY : LDIAG_MISC
USE MODD_TRIP_OASIS, ONLY : LCPL_LAND, LCPL_FLOOD
!                           
!
USE MODE_RW_TRIP
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(TRIP_DIAG_t), INTENT(INOUT) :: TPDG
TYPE(TRIP_GRID_t), INTENT(INOUT) :: TPG
!
INTEGER, INTENT(IN)             :: KLISTING
INTEGER, INTENT(IN)             :: KLON
INTEGER, INTENT(IN)             :: KLAT
!
REAL, INTENT(IN)                :: PRUNTIME 
!
!*      0.2    declarations of local variables
!
 CHARACTER(LEN=18), PARAMETER         :: YRUN         ='TRIP_DIAG_RUN.nc'
 CHARACTER(LEN=50)                    :: YFILE
 CHARACTER(LEN=10)                    :: YVNAME
!
REAL,   DIMENSION(KLON,KLAT)         :: ZWRITE
LOGICAL,DIMENSION(KLON,KLAT)         :: LMASK
LOGICAL,DIMENSION(KLON,KLAT)         :: LMASK_GW
LOGICAL,DIMENSION(KLON,KLAT)         :: LMASK_FLD
!
INTEGER :: ITNUM, ITVAL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG_RUN',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
! * Trip mask
!-------------------------------------------------------------------------------
!
LMASK(:,:) = TPG%GMASK(:,:)
!
! * Groundwater specific mask
!
IF(CGROUNDW/='DEF')THEN
  LMASK_GW(:,:) = TPG%GMASK_GW(:,:)        
ENDIF
!
IF(LFLOOD)THEN
  LMASK_FLD(:,:) = TPG%GMASK_FLD(:,:)
ENDIF
!
!-------------------------------------------------------------------------------
! * outputs attributes
!-------------------------------------------------------------------------------
!
ITNUM = 1
ITVAL = 0
YFILE =YRUN
!
!-------------------------------------------------------------------------------
! * River mass, fluxes, and velocity
!-------------------------------------------------------------------------------
!
YVNAME = 'rivw'
ZWRITE = TPDG%TDIAG_RUN%XSURF_STO / PRUNTIME
CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL)        
!        
YVNAME = 'rivdis'
ZWRITE = TPDG%TDIAG_RUN%XQDIS / PRUNTIME 
CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL)        
!
YVNAME = 'rivin'
ZWRITE = TPDG%TDIAG_RUN%XQIN / PRUNTIME
CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL) 
!
IF(CVIT=='VAR')THEN
!
  YVNAME = 'rivh'
  ZWRITE = TPDG%TDIAG_RUN%XVEL / PRUNTIME
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL)        
!
  YVNAME = 'rivv'
  ZWRITE = TPDG%TDIAG_RUN%XHS / PRUNTIME
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL)        
!
ENDIF
!
!-------------------------------------------------------------------------------
! * Groundwater mass, fluxes, and depth
!-------------------------------------------------------------------------------
!
IF(CGROUNDW/='DEF')THEN
!
  YVNAME = 'gw'
  ZWRITE = TPDG%TDIAG_RUN%XGROUND_STO / PRUNTIME
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL)        
!
  YVNAME = 'gwToriv'
  ZWRITE = TPDG%TDIAG_RUN%XQGF / PRUNTIME 
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL)
!
ENDIF
!
IF(CGROUNDW=='DIF')THEN
!
  YVNAME = 'gwh'
  ZWRITE = TPDG%TDIAG_RUN%XHGROUND / PRUNTIME
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL)
!
  YVNAME = 'wtd'
  ZWRITE = TPDG%TDIAG_RUN%XWTD / PRUNTIME
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL)
!
ENDIF
!
!-------------------------------------------------------------------------------
! * Floodplains
!-------------------------------------------------------------------------------
!
IF(LFLOOD)THEN
!
  YVNAME = 'fldw'
  ZWRITE = TPDG%TDIAG_RUN%XFLOOD_STO / PRUNTIME
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL)
!
  YVNAME = 'fldf'
  ZWRITE = TPDG%TDIAG_RUN%XFF / PRUNTIME
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL)
!
  YVNAME = 'fldh'
  ZWRITE = TPDG%TDIAG_RUN%XHF / PRUNTIME
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL)
!
ENDIF
!
IF(LCPL_LAND)THEN                
  YVNAME = 'RUNOFF'
  ZWRITE = TPDG%TDIAG_RUN%XRUNOFF !mm of water
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL)
  YVNAME = 'DRAIN'
  ZWRITE = TPDG%TDIAG_RUN%XDRAIN !mm of water
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL)
ENDIF
!
IF(LCPL_FLOOD)THEN   
  YVNAME = 'FSOURCE'
  ZWRITE = TPDG%TDIAG_RUN%XSOURCE !mm of water
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL)
ENDIF
!
IF(LDIAG_MISC)THEN
!
  IF(CGROUNDW=='DIF')THEN                
!
    YVNAME = 'gwdf'
    ZWRITE = TPDG%TDIAG_RUN%XFWTD / PRUNTIME
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL)
!
    YVNAME = 'gwfTocell'
    ZWRITE = TPDG%TDIAG_RUN%XQGCELL / PRUNTIME
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL)
!
    YVNAME = 'gwh_less_rivh'
    ZWRITE = TPDG%TDIAG_RUN%XHGHS / PRUNTIME
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL)
!
  ENDIF
!
  IF(LFLOOD)THEN
!
    YVNAME = 'fldh_less_rivh'
    ZWRITE = TPDG%TDIAG_RUN%XHSF / PRUNTIME
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_FLD,ZWRITE,ITNUM,ITVAL)
!
    YVNAME = 'fld_width'
    ZWRITE = TPDG%TDIAG_RUN%XWF / PRUNTIME
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_FLD,ZWRITE,ITNUM,ITVAL)
!
    YVNAME = 'fld_length'
    ZWRITE = TPDG%TDIAG_RUN%XLF / PRUNTIME
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_FLD,ZWRITE,ITNUM,ITVAL)       
!
  ENDIF

ENDIF
!
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG_RUN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE TRIP_DIAG_RUN
