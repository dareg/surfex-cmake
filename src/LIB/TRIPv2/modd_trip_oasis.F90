!###############
MODULE MODD_TRIP_OASIS
!###############
!
!!****  *MODD_TRIP_OASIS - declaration of variable for TRIP-OASIS coupling
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       10/13
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
! * Trip - Oasis coupling general key :
!
!-------------------------------------------------------------------------------
!
LOGICAL             :: LCPL_LAND    = .FALSE. ! Fields from trip to land surface model (e.g. SURFEX, ...)
LOGICAL             :: LCPL_CALVING = .FALSE. ! Calving flux from land surface glacier
LOGICAL             :: LCPL_GW      = .FALSE. ! Fields for/from trip groundwater scheme
LOGICAL             :: LCPL_FLOOD   = .FALSE. ! Fields for/from trip floodplains scheme
LOGICAL             :: LCPL_SEA     = .FALSE. ! Fields from trip to ocean model (e.g. NEMO, ...)
LOGICAL             :: LCPL_CALVSEA = .FALSE. ! Calving fluxes from trip to ocean model (e.g. NEMO, ...)
!
!-------------------------------------------------------------------------------
!
! *  Land surface variables for Trip - Oasis coupling
!
!-------------------------------------------------------------------------------
!
! Input variables
!
INTEGER             :: NRUNOFF_ID    ! Surface runoff id
INTEGER             :: NDRAIN_ID     ! Drainage id
INTEGER             :: NRECHARGE_ID  ! Groundwater recharge id
INTEGER             :: NCALVING_ID   ! Calving flux id
INTEGER             :: NPFLOOD_ID    ! Floodplains precipitation interception id
INTEGER             :: NEFLOOD_ID    ! Floodplains evaporation id
INTEGER             :: NIFLOOD_ID    ! Floodplains infiltration id
!
! Output variables
!
INTEGER             :: NFFLOOD_ID    ! Floodplains fraction id
INTEGER             :: NPIFLOOD_ID   ! Potential flood infiltration id
INTEGER             :: NWTD_ID       ! Water table depth id
INTEGER             :: NFWTD_ID      ! grid-cell fraction of WTD to rise id
!
!-------------------------------------------------------------------------------
!
! * Sea variables for Trip - Oasis coupling
!
!-------------------------------------------------------------------------------
!
! Output variables
!
INTEGER             :: NRIVDIS_ID    ! River discharge id
INTEGER             :: NCALVGRE_ID   ! Calving flux over greenland id
INTEGER             :: NCALVANT_ID   ! Calving flux over antarctica id
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_TRIP_OASIS
