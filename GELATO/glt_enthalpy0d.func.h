FUNCTION glt_enthalpy0d(pt,ps)
!
!   The input arguments are temperature profile, in Celsius
! and salinity (g.kg-1). Note that both arguments are compulsory.
!
  USE modd_glt_const_thm
!
  IMPLICIT NONE
!
  REAL, INTENT(in) ::  &
    & pt 
  REAL, INTENT(in) ::  &
    & ps 
  REAL ::  &
    & glt_enthalpy0d 
!
  REAL ::  &
    & ztice_m 
!
!
! 1. Initializations
! ===================
!
! ..  Compute sea ice melting point as a function of salinity
!
  ztice_m = -mu * ps
!
! 
! 2. If the slab is salty ice
! ============================
!
!* Compute the amount of energy needed to raise sea ice temperature to
! melting point and melt sea ice completely
!
  IF ( ps>0. ) THEN
! If temperature is lower than melting point
    IF ( pt<ztice_m ) THEN
        glt_enthalpy0d = cpice0*( pt-ztice_m ) - xmhofusn0*( 1.-ztice_m/pt )
      ELSE
! If temperature is melting point
        glt_enthalpy0d = 0.
    ENDIF
!
!* Add a term for the energy needed to increase the meltwater temperature to 0
! Celsius
!
    glt_enthalpy0d = glt_enthalpy0d + cpsw*ztice_m
!
  ELSE
!
! 
! 3. If the slab is pure ice
! ===========================
!
    IF ( pt<0. ) THEN
      glt_enthalpy0d = cpice0*pt - xmhofusn0
    ELSE
      glt_enthalpy0d = 0.
    ENDIF
!
  ENDIF   
!
END FUNCTION glt_enthalpy0d
