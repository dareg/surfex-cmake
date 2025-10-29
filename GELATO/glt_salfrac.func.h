FUNCTION glt_salfrac( pv )
! 
! If S is the ocean salinity, as sea ice forms at a rate V (in m.s-1),
! the salinity of the new formed sea ice is : Si = keff.Sw
! Here we follow Cox and Weeks (1988) to compute keff as a function of V.
!
  IMPLICIT NONE
  REAL ::  &
    & glt_salfrac 
  REAL, INTENT(in) ::  &
    & pv 
  REAL ::  &
    & zalpha 

!!  IF ( pv > 3.8e-7 ) THEN
!!      glt_salfrac = 0.26 / ( 0.26+0.74*exp( -724300.*pv ) )
!!    ELSE IF ( pv > 3.4e-7 ) THEN
!!      zalpha = ( pv-3.4e-7 ) / 0.4e-7 
!!      glt_salfrac = zalpha * 0.26 / ( 0.26+0.74*exp( -724300.*pv ) ) +  &
!!        ( 1.-zalpha ) * 0.8925 + 0.0568*log( 100.*pv )
!!    ELSE IF ( pv > 2.2e-8 ) THEN
!!      glt_salfrac = 0.8925 + 0.0568*log( 100.*pv )
!!    ELSE IF ( pv > 1.8e-8 ) THEN
!!      zalpha = ( pv-1.8e-8 ) / 0.4e-8 
!!      glt_salfrac = zalpha * ( 0.8925 + 0.0568*log( 100.*pv ) ) +  &
!!        ( 1.-zalpha ) * 0.12
!!    ELSE 
!!      glt_salfrac = 0.12
!!  ENDIF
  IF ( pv > 3.6e-7 ) THEN
      glt_salfrac = 0.26 / ( 0.26+0.74*EXP( -724300.*pv ) )
    ELSE IF ( pv > 1.24e-8 ) THEN
      glt_salfrac = 0.8925 + 0.0568*LOG( 100.*pv )
    ELSE 
      glt_salfrac = 0.12
  ENDIF
!        
END FUNCTION glt_salfrac
