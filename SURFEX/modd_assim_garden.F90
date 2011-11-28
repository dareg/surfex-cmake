!     ##################
      MODULE MODD_ASSIM_GARDEN
!     ##################
!
!!****  *MODD_ASSIM_GARDEN - declaration of keys for ISBA assimilation scheme (2DVAR, Bouyssel et al.)
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
!!	L. Jarlan   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       23/02/05
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE

!-------------------------------------------------------------------------------
!
! ISBA assimilation Scheme Options:
!
  LOGICAL                        :: LASSIM   ! Assimilation or not
                                             ! '.TRUE.'
                                             ! '.FALSE.'
  CHARACTER(LEN=5)               :: CASSIM   ! type of correction
!                                            ! 'PLUS ' (default)
!                                            ! 'AVERA'            
!                                            ! '2DVAR'
!
! Constants and options of the soil OI analysis
!
 LOGICAL ::  LHUMID,  LIMVEG, LISSEW,  L_SM_WP, LFGEL, LCLIM, LLDHMT,      &
              LOBSWG,  LOBS2M, LPRINT,  LAROME  
 INTEGER ::  MINDJ,   NNEBUL, NNEIGT,  NNEIGW,  NR_SM_WP, NECHGU, NTVGLA,  &
              NSEAICE, NLISSEW,         IDJ,     ITRAD  
 REAL    ::  ANEBUL, RCLIMN, RCLIMTP,  RCLIMTS, RCLIMV,  RCLIMWP, RCLIMWS, &
              SCOEFH, SCOEFT, SEVAP,    SIGH2MO, SIGT2MO, SNEIGT,  SNEIGW,  &
              SPRECIP, SWFC,  V10MX,    RD1,     RTINER,  WCRIN,   WPMX,    &
              WSMX,   TMERGL, RZHZ0G,   RCLIMCA, RCLISST, RWPIA,   RWPIB,   &
              RSNSA,  RSNSB,  SALBM,    SALBB,   SEMIB,   SZZ0B,   SMU0,    &
              SICE,   SEMIM,  RA_SM_WP, RSCALDW, SPRECIP2,                  &
              REPSM,  RCDTR,  SIGHP1,   SIGT2MR, SIGH2MR, RSABR,            &
              RARGR,  GWFC,   EWFC,     GWWILT,  EWWILT,  G1WSAT,  G2WSAT,  &
              REPS1,  REPS2,  REPS3,    ADWR,    SODELX(0:9),               &
              SIGWGO, SIGWGB, SIGW2B,   RTHR_QC, SIGWGO_MAX, RSCAL_JAC  
!
END MODULE MODD_ASSIM_GARDEN
