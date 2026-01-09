FUNCTION glt_interpz(plevn,pvtpo,plevo,nilay) RESULT(tab_interp)
!
! Goal: Interpolate a vertical tracer profile, pvtpo (dimension n),
! defined on a vertical, normalized grid :
!   [ plevo(1)=1, plevo(2), ..., plevo(n), plevo(n+1)=1 ],
! where plevo(jl),plevo(jl+1) define the height (from ice/water bottom
! interface) of respectively the lower and upper boundaries of layer jl.
! Note that n can be any number.
! 
! The glt_output is delivered on the model's standard vertical levels.
!
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: nilay
  REAL, DIMENSION(nilay+1), INTENT(in) ::  &
    & plevn 
  REAL, DIMENSION(:), INTENT(in) ::  &
    & pvtpo 
  REAL, DIMENSION(:), INTENT(in) ::  &
    & plevo 
  REAL, DIMENSION(nilay) ::  &
    & tab_interp 
!
  INTEGER :: jl
!
  DO jl=1,nilay
    tab_interp(jl) =  &
      & ( glt_vtpint(plevn(jl+1),pvtpo,plevo)-  &
        & glt_vtpint(plevn(jl),pvtpo,plevo) ) /  &
      & ( plevn(jl+1)-plevn(jl) ) 
  END DO
!
END FUNCTION glt_interpz
