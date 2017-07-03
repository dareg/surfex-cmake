!     #########
      SUBROUTINE INIT_TRIP_DIAG (TPDG, TPG, &
                                  KLISTING,HFILE,KLON,KLAT,HTITLE,HTIMEUNIT,OTIME)
!     #######################################################################
!
!!****  *INIT_TRIP_DIAG*  
!!
!!    PURPOSE
!!    -------
!
!     Define the name and unit of each trip output variable.
!     
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!      B. Decharme
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    27/05/08 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
!
USE MODD_TRIP_DIAG, ONLY : TRIP_DIAG_t
USE MODD_TRIP_GRID, ONLY : TRIP_GRID_t
!
USE MODE_TRIP_NETCDF
!
USE MODN_TRIP_RUN,   ONLY : LDIAG_MISC
USE MODN_TRIP,       ONLY : CGROUNDW, CVIT, LFLOOD
!
USE MODD_TRIP_OASIS, ONLY : LCPL_LAND, LCPL_FLOOD
!
USE MODD_TRIP_PAR, ONLY : XUNDEF, LNCPRINT
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE NETCDF
!
USE MODI_GET_LONLAT_TRIP
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
!
!
TYPE(TRIP_DIAG_t), INTENT(INOUT) :: TPDG
TYPE(TRIP_GRID_t), INTENT(INOUT) :: TPG
!
 CHARACTER(LEN=*), INTENT(IN) :: HFILE, HTITLE, HTIMEUNIT
!
INTEGER, INTENT(IN)          :: KLISTING, KLON, KLAT
!
LOGICAL, INTENT(IN)          :: OTIME
!
!*      0.2    declarations of output variables
!
INTEGER, PARAMETER :: INDIAG = 100
!
 CHARACTER(LEN=NF90_MAX_NAME), DIMENSION(INDIAG) :: YVNAME  !Name of each output variable
 CHARACTER(LEN=NF90_MAX_NAME), DIMENSION(INDIAG) :: YVLNAME !Long name of each output variables
 CHARACTER(LEN=NF90_MAX_NAME), DIMENSION(INDIAG) :: YUNIT   !Unit of each output variable
!
 CHARACTER(LEN=NF90_MAX_NAME) :: YFILE,YTITLE,YTIMEUNIT
!
REAL, DIMENSION(KLON) ::  ZLON
REAL, DIMENSION(KLAT) ::  ZLAT
!
INTEGER :: INCID, INUM
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
! * Number of output variable
!
IF (LHOOK) CALL DR_HOOK('INIT_TRIP_DIAG',0,ZHOOK_HANDLE)
INUM   = 0
!
!-------------------------------------------------------------------------------
! * River mass, fluxes, and velocity
!-------------------------------------------------------------------------------
!
INUM = INUM + 1
YVNAME (INUM) = 'rivw                      '
YVLNAME(INUM) = 'River storage             '
YUNIT  (INUM) = 'kg m-2                    '
!
INUM = INUM + 1
YVNAME (INUM) = 'rivdis                    '
YVLNAME(INUM) = 'River discharge           '
YUNIT  (INUM) = 'm3 s-1                    '
!
INUM = INUM + 1
YVNAME (INUM) = 'rivin                     '
YVLNAME(INUM) = 'River Inflow              '
YUNIT  (INUM) = 'm3 s-1                    '
!
IF(CVIT=='VAR')THEN
!
  INUM = INUM + 1
  YVNAME (INUM) = 'rivh                      '
  YVLNAME(INUM) = 'River height              '
  YUNIT  (INUM) = 'm                         '
!
  INUM = INUM + 1
  YVNAME (INUM) = 'rivv                      '
  YVLNAME(INUM) = 'River velocity            '
  YUNIT  (INUM) = 'm s-1                     '
!
ENDIF
!
!-------------------------------------------------------------------------------
! * Groundwater mass, fluxes, and depth
!-------------------------------------------------------------------------------
!
IF(CGROUNDW/='DEF')THEN
!
  INUM = INUM + 1
  YVNAME (INUM) = 'gw                '
  YVLNAME(INUM) = 'Groundwater storage       '
  YUNIT  (INUM) = 'kg m-2                    '
!
  INUM = INUM + 1
  YVNAME (INUM) = 'gwToriv                   '
  YVLNAME(INUM) = 'Groundwater-river exchange'
  YUNIT  (INUM) = 'kg m-2 s-1                '
!
ENDIF
!
IF(CGROUNDW=='DIF')THEN
!
  INUM = INUM + 1
  YVNAME (INUM) = 'gwh                               '
  YVLNAME(INUM) = 'Groundwater height from 0 altitude'
  YUNIT  (INUM) = 'm                                 '
!
  INUM = INUM + 1
  YVNAME (INUM) = 'wtd                                  '
  YVLNAME(INUM) = 'Water Table Depth (positive downward)'
  YUNIT  (INUM) = 'm                                    '
!
ENDIF
!
!-------------------------------------------------------------------------------
! * Floodplains
!-------------------------------------------------------------------------------
!
IF(LFLOOD)THEN
!        
  INUM = INUM + 1
  YVNAME (INUM) = 'fldw                      '
  YVLNAME(INUM) = 'Floodplain storage        '
  YUNIT  (INUM) = 'kg m-2                    '
!
  INUM = INUM + 1
  YVNAME (INUM) = 'fldf                      '
  YVLNAME(INUM) = 'Floodplain fraction       '
  YUNIT  (INUM) = '-                         '
!
  INUM = INUM + 1
  YVNAME (INUM) = 'fldh                      '
  YVLNAME(INUM) = 'Floodplain depth          '
  YUNIT  (INUM) = 'm                         '
!
ENDIF
!
!-------------------------------------------------------------------------------
! * Forcing variables can be used to force TRIP offline
!-------------------------------------------------------------------------------
!
IF(LCPL_LAND)THEN
! 
  INUM = INUM + 1
  YVNAME (INUM) = 'RUNOFF                    '
  YVLNAME(INUM) = 'Input surface runoff (can be used to force TRIP offline)'
  YUNIT  (INUM) = 'mm of water               '
! 
  INUM = INUM + 1
  YVNAME (INUM) = 'DRAIN                     '
  YVLNAME(INUM) = 'Input drainage or recharge (can be used to force TRIP offline)'
  YUNIT  (INUM) = 'mm of water               '  
! 
ENDIF
!
IF(LCPL_FLOOD)THEN
  INUM = INUM + 1
  YVNAME (INUM)= 'FSOURCE                   '
  YVLNAME(INUM)= 'Floodplains source (Pf-Ef-If) (can be used to force TRIP offline)'
  YUNIT  (INUM)= 'mm of water               '
ENDIF
!
!-------------------------------------------------------------------------------
! * MISC fields
!-------------------------------------------------------------------------------
!
IF(LDIAG_MISC)THEN
!
  IF(CGROUNDW=='DIF')THEN
!
    INUM = INUM + 1
    YVNAME (INUM) = 'gwdf                                   '
    YVLNAME(INUM) = 'grid-cell fraction of Groundwater depth'
    YUNIT  (INUM) = '-                                      '
!
    INUM = INUM + 1
    YVNAME (INUM) = 'gwfTocell                               ' 
    YVLNAME(INUM) = 'Groundwater fluxes between adjacent cell'
    YUNIT  (INUM) = 'kg m-2 s-1                              '
!
    INUM = INUM + 1
    YVNAME (INUM)= 'gwh_less_rivh              '
    YVLNAME(INUM)= 'Hground - Hriver           '
    YUNIT  (INUM)= 'm                          '
!
  ENDIF
!
  IF(LFLOOD)THEN
!
    INUM = INUM + 1
    YVNAME (INUM)= 'fldh_less_rivh              '
    YVLNAME(INUM)= 'Hflood - Hriver             '
    YUNIT  (INUM)= 'm                           '
!
    INUM = INUM + 1
    YVNAME (INUM)= 'fld_width                   '
    YVLNAME(INUM)= 'Flood width                 '
    YUNIT  (INUM)= 'm                           '
!
    INUM = INUM + 1
    YVNAME (INUM)= 'fld_length                  '
    YVLNAME(INUM)= 'Flood length                '
    YUNIT  (INUM)= 'm                           '
!
  ENDIF
!
ENDIF     
!
!-------------------------------------------------------------------------------
! * Create netcdf file
!-------------------------------------------------------------------------------
!
! * Create netcdf file
!
YFILE     = HFILE(1:LEN_TRIM(HFILE))
YTITLE    = HTITLE(1:LEN_TRIM(HTITLE))
YTIMEUNIT = HTIMEUNIT(1:LEN_TRIM(HTIMEUNIT))
!
 CALL GET_LONLAT_TRIP(TPG,KLON,KLAT,ZLON,ZLAT)
!
 CALL NCCREATE(KLISTING,YFILE,YTITLE,YTIMEUNIT,YVNAME,YVLNAME,YUNIT,ZLON,ZLAT,XUNDEF,LNCPRINT,INCID,OTIME)
!
 CALL NCCLOSE(KLISTING,LNCPRINT,YFILE,INCID)
!
IF (LHOOK) CALL DR_HOOK('INIT_TRIP_DIAG',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
END SUBROUTINE INIT_TRIP_DIAG
