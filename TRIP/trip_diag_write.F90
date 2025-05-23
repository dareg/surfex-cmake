!     #########
      SUBROUTINE TRIP_DIAG_WRITE (TPDG, TPG, &
                                  KLISTING,KLON,KLAT,KDIAG,PTSTEP_DIAG,OXIOS)  
!     ############################################################
!
!!****  *TRIP_DIAG_WRITE*  
!!
!!    PURPOSE
!!    -------
!
!     TRIP river routing outputs.
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
USE MODD_TRIP_OASIS, ONLY : LCPL_LAND
!
USE MODD_TRIP_PAR,   ONLY : XTIME_DIAG
!
#ifdef WXIOS
use XIOS
#endif
!
USE MODE_RW_TRIP
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
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
INTEGER, INTENT(IN)             :: KDIAG
REAL,    INTENT(IN)             :: PTSTEP_DIAG
LOGICAL, INTENT(IN)             :: OXIOS
!
!*      0.2    declarations of local variables
!
 CHARACTER(LEN=19), PARAMETER         :: YDIAG  = 'TRIP_DIAG.nc'
 CHARACTER(LEN=50)                    :: YFILE
 CHARACTER(LEN=10)                    :: YVNAME
!
REAL,   DIMENSION(KLON,KLAT)         :: ZWRITE
LOGICAL,DIMENSION(KLON,KLAT)         :: LMASK
LOGICAL,DIMENSION(KLON,KLAT)         :: LMASK_GW
!
INTEGER :: ITNUM, ITVAL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG_WRITE',0,ZHOOK_HANDLE)
!
! * Trip mask
!
LMASK(:,:) = TPG%GMASK(:,:)
!
! * Groundwater specific mask
!
IF(CGROUNDW/='DEF')THEN
  LMASK_GW(:,:) = TPG%GMASK_GW(:,:)        
ENDIF
!
!-------------------------------------------------------------------------------
!outputs
!-------------------------------------------------------------------------------
!
! * Recup diag file
!
YFILE =YDIAG
!
! * Time attribut
!
ITNUM = KDIAG
!
IF (OXIOS) THEN
#ifdef WXIOS
  CALL XIOS_UPDATE_CALENDAR(KDIAG)
#endif
  ITVAL = 0
ELSE
  ITVAL = (KDIAG-1)*INT(PTSTEP_DIAG/XTIME_DIAG)
ENDIF
!
! * Store output in diag file
!
YVNAME = 'SURF_STO'
ZWRITE = TPDG%TDIAG%XSURF_STO / PTSTEP_DIAG
 CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)   
 CALL CUM_DIAG(TPDG%TDIAG%XSURF_STO,TPDG%TDIAG_RUN%XSURF_STO)
!                
YVNAME = 'QDIS'
ZWRITE = TPDG%TDIAG%XQDIS / PTSTEP_DIAG
 CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS) 
 CALL CUM_DIAG(TPDG%TDIAG%XQDIS,TPDG%TDIAG_RUN%XQDIS)
!
IF(LDIAG_MISC)THEN                
  YVNAME = 'QSIN'
  ZWRITE = TPDG%TDIAG%XQIN / PTSTEP_DIAG
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)
  CALL CUM_DIAG(TPDG%TDIAG%XQIN,TPDG%TDIAG_RUN%XQIN)
ENDIF
!
IF(LCPL_LAND.AND.LDIAG_MISC)THEN                
! kg/m2 (can be used to force the model offline) 
  YVNAME = 'RUNOFF'
  ZWRITE = TPDG%TDIAG%XRUNOFF
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)
  CALL CUM_DIAG(TPDG%TDIAG%XRUNOFF,TPDG%TDIAG_RUN%XRUNOFF)
! kg/m2 (can be used to force the model offline)
  YVNAME = 'DRAIN'
  ZWRITE = TPDG%TDIAG%XDRAIN
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)
  CALL CUM_DIAG(TPDG%TDIAG%XDRAIN,TPDG%TDIAG_RUN%XDRAIN)
!          
ENDIF
!
IF(CGROUNDW/='DEF')THEN
!
  YVNAME = 'QGF'
  ZWRITE = TPDG%TDIAG%XQGF / PTSTEP_DIAG 
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)
  CALL CUM_DIAG(TPDG%TDIAG%XQGF,TPDG%TDIAG_RUN%XQGF)
!          
  YVNAME = 'GROUND_STO'
  ZWRITE = TPDG%TDIAG%XGROUND_STO / PTSTEP_DIAG
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS) 
  CALL CUM_DIAG(TPDG%TDIAG%XGROUND_STO,TPDG%TDIAG_RUN%XGROUND_STO)
!  
ENDIF
!
IF(CGROUNDW=='DIF')THEN
!
  YVNAME = 'HGROUND'
  ZWRITE = TPDG%TDIAG%XHGROUND / PTSTEP_DIAG
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)
  CALL CUM_DIAG(TPDG%TDIAG%XHGROUND,TPDG%TDIAG_RUN%XHGROUND)
!       
  YVNAME = 'FWTD'
  ZWRITE = TPDG%TDIAG%XFWTD / PTSTEP_DIAG
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)  
  CALL CUM_DIAG(TPDG%TDIAG%XFWTD,TPDG%TDIAG_RUN%XFWTD)
!       
  YVNAME = 'WTD'
  ZWRITE = TPDG%TDIAG%XWTD / PTSTEP_DIAG 
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)
  CALL CUM_DIAG(TPDG%TDIAG%XWTD,TPDG%TDIAG_RUN%XWTD)
!   
  IF(LDIAG_MISC)THEN      
!      
      YVNAME = 'QGCELL'
      ZWRITE = TPDG%TDIAG%XQGCELL / PTSTEP_DIAG 
      CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)
      CALL CUM_DIAG(TPDG%TDIAG%XQGCELL,TPDG%TDIAG_RUN%XQGCELL)
!       
      YVNAME = 'HGHRIV'
      ZWRITE = TPDG%TDIAG%XHGHS / PTSTEP_DIAG
      CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK_GW,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)
      CALL CUM_DIAG(TPDG%TDIAG%XHGHS,TPDG%TDIAG_RUN%XHGHS)
!
  ENDIF
! 
ENDIF
!
IF(CVIT=='VAR')THEN
!
  YVNAME = 'VEL'
  ZWRITE = TPDG%TDIAG%XVEL / PTSTEP_DIAG
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)
  CALL CUM_DIAG(TPDG%TDIAG%XVEL,TPDG%TDIAG_RUN%XVEL)
!
  YVNAME = 'HSTREAM'
  ZWRITE = TPDG%TDIAG%XHS / PTSTEP_DIAG
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS) 
  CALL CUM_DIAG(TPDG%TDIAG%XHS,TPDG%TDIAG_RUN%XHS)
!
ENDIF
!
IF(LFLOOD)THEN
!
  YVNAME = 'FFLOOD'
  ZWRITE = TPDG%TDIAG%XFF / PTSTEP_DIAG
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)   
  CALL CUM_DIAG(TPDG%TDIAG%XFF,TPDG%TDIAG_RUN%XFF)
!
  YVNAME = 'FLOOD_STO'
  ZWRITE = TPDG%TDIAG%XFLOOD_STO / PTSTEP_DIAG
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)   
  CALL CUM_DIAG(TPDG%TDIAG%XFLOOD_STO,TPDG%TDIAG_RUN%XFLOOD_STO)
!
  YVNAME = 'HFLOOD'
  ZWRITE = TPDG%TDIAG%XHF / PTSTEP_DIAG
  CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)  
  CALL CUM_DIAG(TPDG%TDIAG%XHF,TPDG%TDIAG_RUN%XHF)
!
  IF(LDIAG_MISC)THEN
!
! kg/m2 (can be used to force the model offline)
    YVNAME = 'FSOURCE'
    ZWRITE = TPDG%TDIAG%XSOURCE
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS) 
    CALL CUM_DIAG(TPDG%TDIAG%XSOURCE,TPDG%TDIAG_RUN%XSOURCE)
!
! kg/m2/s
    YVNAME = 'QFR'
    ZWRITE = TPDG%TDIAG%XQFR / PTSTEP_DIAG
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)  
    CALL CUM_DIAG(TPDG%TDIAG%XQFR,TPDG%TDIAG_RUN%XQFR)
!                
    YVNAME = 'QRF'
    ZWRITE = TPDG%TDIAG%XQRF / PTSTEP_DIAG
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS) 
    CALL CUM_DIAG(TPDG%TDIAG%XQRF,TPDG%TDIAG_RUN%XQRF)
!
    YVNAME = 'VFIN'
    ZWRITE = TPDG%TDIAG%XVFIN / PTSTEP_DIAG
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS) 
    CALL CUM_DIAG(TPDG%TDIAG%XVFIN,TPDG%TDIAG_RUN%XVFIN) 
!
    YVNAME = 'VFOUT'
    ZWRITE = TPDG%TDIAG%XVFOUT / PTSTEP_DIAG
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS) 
    CALL CUM_DIAG(TPDG%TDIAG%XVFOUT,TPDG%TDIAG_RUN%XVFOUT)
!
    YVNAME = 'HSF'
    ZWRITE = TPDG%TDIAG%XHSF / PTSTEP_DIAG
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)
    CALL CUM_DIAG(TPDG%TDIAG%XHSF,TPDG%TDIAG_RUN%XHSF)
!
    YVNAME = 'WF'
    ZWRITE = TPDG%TDIAG%XWF / PTSTEP_DIAG
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)  
    CALL CUM_DIAG(TPDG%TDIAG%XWF,TPDG%TDIAG_RUN%XWF)
!
    YVNAME = 'LF'
    ZWRITE = TPDG%TDIAG%XLF / PTSTEP_DIAG
    CALL WRITE_TRIP(KLISTING,YFILE,YVNAME,LMASK,ZWRITE,ITNUM,ITVAL,OXIOS=OXIOS)
    CALL CUM_DIAG(TPDG%TDIAG%XLF,TPDG%TDIAG_RUN%XLF)
!      
  ENDIF
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG_WRITE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
 CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE CUM_DIAG(PDIAG,PRUN)
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:),INTENT(INOUT) :: PDIAG
REAL, DIMENSION(:,:),INTENT(OUT  ) :: PRUN
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG_WRITE:CUM_DIAG',0,ZHOOK_HANDLE)
!
PRUN(:,:) = PRUN(:,:) + PDIAG(:,:)
!
PDIAG(:,:) = 0.0
!
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG_WRITE:CUM_DIAG',1,ZHOOK_HANDLE)
!
END SUBROUTINE CUM_DIAG
!
!-------------------------------------------------------------------------------
END SUBROUTINE TRIP_DIAG_WRITE
