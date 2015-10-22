!     #########
SUBROUTINE OL_READ_ATM_CONF (YSC, &
                             HSURF_FILETYPE, HFORCING_FILETYPE,  &
                              ODELAYEDSTART_NC, KDATESTOP,      &
                              PDURATION,        &
                              PTSTEP_FORC, KNI, KYEAR,KMONTH,     &
                              KDAY, PTIME, PLAT, PLON,            &
                              PZS, PZREF, PUREF,KTIMESTARTINDEX)  
!
!==================================================================
!!****  *OL_READ_ATM_CONF* - Initialization routine
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      F. Habets   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified by P. Le Moigne (04/2005): cleaning and checking
!!      Modified by P. Le Moigne (04/2006): init_io_surf for nature
!!                  with GTMSK to read dimensions.
!!      Modified by M. Lafaysse 04/2015 : option ODELAYEDSTART_NC
!==================================================================
!
!
USE MODD_SURFEX_n, ONLY : SURFEX_t
!
USE MODI_OL_READ_ATM_CONF_NETCDF
USE MODI_OL_READ_ATM_CONF_ASCII
USE MODD_SURF_CONF,      ONLY : CPROGNAME
!==================================================================
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE(SURFEX_t), INTENT(INOUT) :: YSC
!
 CHARACTER(LEN=6), INTENT(IN)  :: HSURF_FILETYPE
 CHARACTER(LEN=6), INTENT(IN)  :: HFORCING_FILETYPE
 LOGICAL, INTENT(IN)           :: ODELAYEDSTART_NC !Allow the simulation to start from a different time step than the first record of a netcdf file
INTEGER,DIMENSION(4),INTENT(IN) :: KDATESTOP !Allow the simulation to end at a different time step than the last record of a netcdf file
INTEGER,          INTENT(OUT) :: KNI
INTEGER,          INTENT(OUT) :: KYEAR, KMONTH, KDAY
REAL,             INTENT(OUT) :: PDURATION,PTSTEP_FORC
REAL,             INTENT(OUT) :: PTIME
INTEGER,          INTENT(OUT) :: KTIMESTARTINDEX ! index from which we start reading FORCING.nc
REAL, DIMENSION(:),  POINTER  :: PLAT, PLON
REAL, DIMENSION(:),  POINTER  :: PZS
REAL, DIMENSION(:),  POINTER  :: PZREF, PUREF
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!==================================================================
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_CONF',0,ZHOOK_HANDLE)
CPROGNAME = HSURF_FILETYPE
!
IF (HFORCING_FILETYPE == 'NETCDF') THEN
!
 CALL OL_READ_ATM_CONF_NETCDF(YSC, &
                              HSURF_FILETYPE,                     &
                                ODELAYEDSTART_NC,KDATESTOP,       &
                                PDURATION,         &
                                PTSTEP_FORC, KNI, KYEAR,KMONTH,     &
                                KDAY, PTIME, PLAT, PLON,            &
                                PZS, PZREF, PUREF ,KTIMESTARTINDEX)
!
ELSE IF (HFORCING_FILETYPE == 'ASCII ' .OR. HFORCING_FILETYPE == 'BINARY') THEN
!
 CALL OL_READ_ATM_CONF_ASCII(YSC, &
                              HSURF_FILETYPE,HFORCING_FILETYPE,   &
                                PDURATION,                          &
                                PTSTEP_FORC, KNI, KYEAR,KMONTH,     &
                                KDAY, PTIME, PLAT, PLON,            &
                                PZS, PZREF, PUREF                   )
 KTIMESTARTINDEX=1
!                              
ENDIF
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_CONF',1,ZHOOK_HANDLE)
!
!==================================================================
!
END SUBROUTINE OL_READ_ATM_CONF
