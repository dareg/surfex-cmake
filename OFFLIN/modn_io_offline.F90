!     #######################
      MODULE  MODN_IO_OFFLINE
!     #######################
!
!!****  *MODN_IO_OFFLINE* define the variables and namelist for SURFEX
!                         offline programs (pgd, prep, offline)
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2008
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*    Types of files
!     --------------
!
CHARACTER(LEN=6) :: CSURF_FILETYPE       = 'ASCII ' ! type of SURFEX surface files
!                                                   ! 'NETDF '
!                                                   ! 'FA    '
!                                                   ! 'ASCII '
!                                                   ! 'LFI   '
CHARACTER(LEN=6) :: CTIMESERIES_FILETYPE = 'NONE  ' ! type of the files contining the
!                                                   ! output diagnostic time series
!                                                   ! 'NETCDF ', 'TEXTE '
CHARACTER(LEN=6) :: CFORCING_FILETYPE    = 'NETCDF' ! type of atmospheric FORCING files
!                                                   ! 'NETDF', 'BINARY', or 'ASCII '
!
!
!*    Names of files
!     --------------
!
CHARACTER(LEN=28):: CPGDFILE  ='PGD'          ! name of the PGD file
CHARACTER(LEN=28):: CPREPFILE ='PREP'         ! name of the INITIAL file
CHARACTER(LEN=28):: CSURFFILE ='SURFOUT'      ! name of the final output CSURFEX file
CHARACTER(LEN=28):: CNAMELIST ='OPTIONS.nam'  ! name of namelist file
!
!
!*    General flags defining output options
!     -------------------------------------
!
LOGICAL          :: LPRINT   = .FALSE.  ! write some information on screen 
LOGICAL          :: LRESTART = .FALSE.  ! write restart file
LOGICAL          :: LINQUIRE = .FALSE.  ! inquiry mode
!      
LOGICAL          :: LWRITE_COORD = .FALSE. ! write lat/lon of the target grid
!
LOGICAL          :: LOUT_TIMENAME = .FALSE.! change the name of output file at the end of a day
                                           ! (ex: 19860502_00h00 -> 19860501_24h00)
!
!*    Time steps
!     ----------
!
REAL             :: XTSTEP_SURF   = 300.   ! time step of the surface 
REAL             :: XTSTEP_OUTPUT = 1800.  ! time step of the output time-series
INTEGER          :: NB_READ_FORC  = 0      ! subdivisions of the reading of forcings
!
LOGICAL          :: LLAND_USE = .FALSE.
!
!*    General flag for coherence between forcing file orography and surface file orography
!     ----------
!
LOGICAL          :: LSET_FORC_ZS =.FALSE.  ! .T. : the orography of the
!                                          !  forcing file is
!                                          !  automatically set to the same
!                                          !  value as in the surface file
!                                          ! .F. : the orography of the
!                                          !  forcing file is kept as it is
!
!*    General flag for coherence between forcing Qair and calculated Qsat(Tair)
!     ----------
!
LOGICAL          :: LLIMIT_QAIR = .FALSE. ! .T. : Qair always <= Qsat(Tair)
                                          ! .F. : No limitation
!
!-------------------------------------------------------------------------------
!
!*       1.    NAMELISTS
!              ---------
!
NAMELIST/NAM_IO_OFFLINE/CSURF_FILETYPE, CTIMESERIES_FILETYPE, CFORCING_FILETYPE, &
                          CPGDFILE, CPREPFILE, CSURFFILE,                          &
                          LPRINT, LRESTART, LINQUIRE,                              &
                          XTSTEP_SURF, XTSTEP_OUTPUT,                              &
                          LSET_FORC_ZS, LWRITE_COORD, LOUT_TIMENAME, LLIMIT_QAIR,  &
                          NB_READ_FORC, LLAND_USE  
!
!-------------------------------------------------------------------------------
END MODULE MODN_IO_OFFLINE
