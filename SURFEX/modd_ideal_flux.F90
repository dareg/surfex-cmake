!     ######################
      MODULE MODD_IDEAL_FLUX
!     ######################
!
!!****  *MODD_IDEAL_FLUX * - Defines the quantities for ideal surface fluxes.
!!
!!    PURPOSE
!!    -------
!
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
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!	    V. Masson   * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2003 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
REAL                              :: XFREQ      ! frequency of surface forcing
INTEGER, PARAMETER                :: NFORC =48  ! number of surface forcing instants
!
REAL, DIMENSION(0:NFORC)          :: XSFTH      ! hourly data of heat surface flux        (W/m2)
REAL, DIMENSION(0:NFORC)          :: XSFTQ      ! hourly data of water vapor surface flux (kg/m2/s) or (W/m2)
REAL, DIMENSION(:,:), ALLOCATABLE :: XSFTS      ! hourly data of scalar surface flux      (kg/m2/s)
REAL, DIMENSION(0:NFORC)          :: XSFCO2     ! hourly data of CO2 surface flux         (kg/m2/s)
CHARACTER(LEN=5)                  :: CUSTARTYPE ! type of computation for friction
                                                ! 'USTAR'
                                                ! 'Z0   '
REAL, DIMENSION(0:NFORC)          :: XUSTAR     ! hourly data of friction                 (m2/s2)
REAL                              :: XZ0        ! roughness length (m)
REAL                              :: XALB       ! albedo (-)
REAL                              :: XEMIS      ! emissivity (-)
REAL, DIMENSION(0:NFORC)          :: XTSRAD     ! radiative temperature (K)
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_IDEAL_FLUX
