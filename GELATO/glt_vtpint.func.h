! * If pvtpo (old vertical tracer profile) has dimension (n), 
! plevo (old, normalized levels: bottom=0, top=1) should have 
! dimension (n+1)
! * Output: integral between z=0 and z=pz of vertical tracer.
!
FUNCTION glt_vtpint(pz,pvtpo,plevo)
!
  IMPLICIT NONE
!
  REAL, INTENT(in) ::  &
    & pz 
  REAL, DIMENSION(:), INTENT(in) ::  &
    & pvtpo 
  REAL, DIMENSION(:), INTENT(in) ::  &
    & plevo 
  REAL ::  &
    & glt_vtpint 
!
  INTEGER :: jl
!
  glt_vtpint = 0.
!
  DO jl=1,SIZE(pvtpo)
    glt_vtpint = glt_vtpint +  &
      & pvtpo(jl)*  &
      & ( ( AMIN1(pz,plevo(jl+1))-plevo(jl)) *  &
      & ( 1.+SIGN(1.,pz-plevo(jl)) )/2. ) 
  END DO
!
END FUNCTION glt_vtpint
