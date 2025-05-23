!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###################
      MODULE MODD_PGDWORK
!     ###################
!
!!****  *MODD_PGDWORK* - declaration of work arrays and variables
!!
!!    PURPOSE
!!    -------  
!!
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/09/95   
!!                  15/03/96 orographic filter parameter                   
!!                  25/07/97 directional z0 computations
!!                  15/03/99 add XSUMCOVER
!!                  03/2004  externalization
!!                  01/2012  add aggregation with the MAJORITY rule
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
! to get same results in parallel mode (if results differ, change XPREC to a
! lower value)
REAL, PARAMETER :: XPREC = 1.0E+8
!
!*        0.1    summation variables
!                -------------------
!
REAL, DIMENSION(:,:), ALLOCATABLE   :: XSUMVAL  
                            ! Sum of data in each mesh
!
REAL, DIMENSION(:,:), ALLOCATABLE   :: XEXT_ALL
                            ! Sum of square data in each mesh 
!        
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XALL
                            ! Sum of each cover type data in each mesh                           
!  
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSSO_ALL
                            ! max of orography in a SSO subgrid square
!
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: NSSO_ALL
                            ! presence of data in a SSO subgrid square
!  
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XFSSO_ALL
                            ! mean of orography in a fractional SSO subgrid square
!
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: NFSSO_ALL
                            ! number of data in a fractional SSO subgrid square      
!
INTEGER, DIMENSION(:,:), ALLOCATABLE:: NSIZE
!                          ! Number of points inside each mesh of the domain
INTEGER, DIMENSION(:,:), ALLOCATABLE:: NSIZE_ALL
!                          ! Number of points inside each mesh of the domain
!
 CHARACTER(LEN=3)                  :: CATYPE = 'ARI'
!                          ! Type of averaging:
!                          ! 'ARI' : arithmetic
!                          ! 'INV' : inverse
!                          ! 'CDN' : neutral CD
!
!*        0.2    variables for SSO computations
!                ------------------------------
!
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: XSSQO ! max of orography in a
!                                               ! SSO subgrid square from
!                                               ! ZMAXSSQ averaged values
!                                               ! 1st dim: NSSO (x direc.)
!                                               ! 2st dim: NSSO (y direc.)
!                                               ! 3nd dim: number of grid meshes
LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: LSSQO ! presence of data in a SSO
!                                               ! subgrid square
INTEGER :: NSSO                                 ! number of SSO subgrid squares
!                                               ! in each direction in grid mesh
!
! Fractional SSO
!
LOGICAL                                :: LORORAD     ! activate orographic radiation parameters
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: XFSSQO      ! mean of orography in a
!                                                     ! SSO subgrid square from
!                                                     ! NFSSQO averaged values
!                                                     ! 1st dim: NFSSO (x direc.)
!                                                     ! 2st dim: NFSSO (y direc.)
!                                                     ! 3nd dim: number of grid meshes
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: NFSSQO      ! number of contributors to a
!                                                     ! NFSSO subgrid square
INTEGER                                :: NFSSO       ! number of FSSO subgrid squares
!                                                     ! in each direction in grid mesh
REAL                                   :: XRFSSO      ! reduction factor for computing NFSSO
INTEGER                                :: NFSSOMAX    ! max for NFSSO (limit for memory reasons)
REAL                                   :: XHALORADIUS ! radius of the halo in which the horizon is computed (m)
REAL, PARAMETER                        :: XFLATRAD=0.0! slope angle (rad) under which mesh (or submesh) is considered flat
!
!
!*        0.3    variables for topographic index statistics computations
!                -------------------------------------------------------
!
REAL, DIMENSION(:), ALLOCATABLE   :: XMIN_WORK
REAL, DIMENSION(:), ALLOCATABLE   :: XMAX_WORK
REAL, DIMENSION(:), ALLOCATABLE   :: XMEAN_WORK
REAL, DIMENSION(:), ALLOCATABLE   :: XSTD_WORK
REAL, DIMENSION(:), ALLOCATABLE   :: XSKEW_WORK
!
!
!*        0.4    Variables for the Majority aggregation rule
!                -------------------------------------------
!
INTEGER, PARAMETER :: JPVALMAX=20  ! Maximum number of different values 
!                                  ! in each grid mesh
INTEGER, DIMENSION(:,:),   ALLOCATABLE :: NVALNBR
!                                  ! number of different values 
!                                  ! in each grid mesh
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: NVALCOUNT
!                                  ! Number of times each value has been 
!                                  ! counted in each grid mesh
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: XVALLIST
!                                  ! List of Values encountered in each grid mesh
!-------------------------------------------------------------------------------
!
!ek_beg
!
!*        0.5    Variables for the lake database
!                ------------------------------------------- 
!
REAL, DIMENSION(:), ALLOCATABLE :: XFRAC_LDB   ! Fraction of lakes calculated from the lake database
!-------------------------------------------------------------------------------
!
END MODULE MODD_PGDWORK
