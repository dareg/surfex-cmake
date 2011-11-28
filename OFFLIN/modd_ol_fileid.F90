
!     ######################
      MODULE MODD_OL_FILEID
!     ######################
!
!!****  *MODD_OL_FILEID* Keep in memory the netcdf ID of the output files
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
!!	F. Habets   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!------------------------------------------------------------------------------
!
!* variables for each patch
!
CHARACTER(LEN=200), DIMENSION(7)  :: XNETCDF_FILENAME_IN= &
                                       (/'PARAMS.nc                  ',&
                                         'FORCING.nc                 ',&
                                         'SXPOST.nc                  ', &
                                         'ISBA_VEG_EVOLUTION_P.OUT.nc', &
                                         'ISBA_VEG_EVOLUTION_A.OUT.nc',&
                                         'ISBA_PROGNOSTIC.OUT.nc     ',&
                                         'ISBA_DIAGNOSTICS.OUT.nc    '/)
CHARACTER(LEN=200), DIMENSION(21) :: XNETCDF_FILENAME_OUT= &
                                       (/'ISBA_VEG_EVOLUTION.OUT.nc  ',&
                                         'ISBA_VEG_EVOLUTION_P.OUT.nc', &
                                         'ISBA_VEG_EVOLUTION_A.OUT.nc', &
                                         'ISBA_PROGNOSTIC.OUT.nc     ',&
                                         'ISBA_DIAGNOSTICS.OUT.nc    ',&
                                         'ISBA_DIAG_CUMUL.OUT.nc     ',&
                                         'SEAFLUX_PROGNOSTIC.OUT.nc  ',&
                                         'SEAFLUX_DIAGNOSTICS.OUT.nc ',&
                                         'SEAFLUX_DIAG_CUMUL.OUT.nc  ',&
                                         'WATFLUX_PROGNOSTIC.OUT.nc  ',&
                                         'WATFLUX_DIAGNOSTICS.OUT.nc ',&
                                         'WATFLUX_DIAG_CUMUL.OUT.nc  ',&
                                         'FLAKE_PROGNOSTIC.OUT.nc    ',&
                                         'FLAKE_DIAGNOSTICS.OUT.nc   ',&
                                         'FLAKE_DIAG_CUMUL.OUT.nc    ',&
                                         'TEB_PROGNOSTIC.OUT.nc      ',&
                                         'TEB_DIAGNOSTICS.OUT.nc     ',&
                                         'TEB_CANOPY.OUT.nc          ',&
                                         'TEB_DIAG_CUMUL.OUT.nc      ',&
                                         'SURF_ATM.OUT.nc            ',&
                                         'SURF_ATM_DIAGNOSTICS.OUT.nc'/)  
CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE :: XVAR_TO_FILEOUT, XVAR_TO_FILEIN !contains names
                                                                                !of variables to write
INTEGER*4, DIMENSION(:), ALLOCATABLE :: XID, XID_IN, XID_VARIN  !contains ids of
                                                                !opened files for each 
                                                                !variable to write

CHARACTER(LEN=20), DIMENSION(:), POINTER :: XVAR_SURF, XVAR_NATURE, &
                                             XVAR_SEA, XVAR_WATER, XVAR_TOWN  
INTEGER*4, DIMENSION(:), POINTER :: XID_SURF, XID_NATURE, XID_SEA,  &
                                     XID_WATER, XID_TOWN  
INTEGER         :: XOUT, XIN

INTEGER :: XCOUNT
!------------------------------------------------------------------------------
!
END MODULE MODD_OL_FILEID

