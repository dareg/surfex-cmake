!     ######################################################################
      SUBROUTINE SUBSCALE_Z0EFF_1D(PAOSIP,PAOSIM,PAOSJP,PAOSJM,            &
                                PHO2IP,PHO2IM,PHO2JP,PHO2JM,PZ0VEG,        &
                                PZ0EFFIP,PZ0EFFIM,PZ0EFFJP,PZ0EFFJM,       &
                                OMASK                                      )
!     ######################################################################
!
!!*SUBSCALE_Z0EFF  computes an effective roughness lenght deduced
!!                 from the subgrid-scale orography.
!!
!!
!!    METHOD
!!    ------
!!    See M.Georgelin and al. July 1994, Monthly Weather Review.
!!   
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    M. Georgelin      Laboratoire d'Aerologie
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    18/12/95
!!                22/12/97 (V Masson) call with dummy arguments
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_CSTS,     ONLY : XKARMAN
USE MODD_ISBA_PAR, ONLY : XCDZ0EFF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
REAL, DIMENSION(:), INTENT(IN)  :: PAOSIP  ! A/S for increasing x
REAL, DIMENSION(:), INTENT(IN)  :: PAOSIM  ! A/S for decreasing x
REAL, DIMENSION(:), INTENT(IN)  :: PAOSJP  ! A/S for increasing y
REAL, DIMENSION(:), INTENT(IN)  :: PAOSJM  ! A/S for decreasing y
REAL, DIMENSION(:), INTENT(IN)  :: PHO2IP  ! h/2 for increasing x
REAL, DIMENSION(:), INTENT(IN)  :: PHO2IM  ! h/2 for decreasing x
REAL, DIMENSION(:), INTENT(IN)  :: PHO2JP  ! h/2 for increasing y
REAL, DIMENSION(:), INTENT(IN)  :: PHO2JM  ! h/2 for decreasing y
REAL, DIMENSION(:), INTENT(IN)  :: PZ0VEG  ! vegetation roughness length
!
REAL, DIMENSION(:), INTENT(INOUT) :: PZ0EFFIP! roughness length for increasing x
REAL, DIMENSION(:), INTENT(INOUT) :: PZ0EFFIM! roughness length for decreasing x
REAL, DIMENSION(:), INTENT(INOUT) :: PZ0EFFJP! roughness length for increasing y
REAL, DIMENSION(:), INTENT(INOUT) :: PZ0EFFJM! roughness length for decreasing y
!
LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: OMASK ! mask where computations
                                                       ! are done
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
LOGICAL, DIMENSION(SIZE(PZ0EFFIM)) :: GMASK
!
INTEGER         :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SUBSCALE_Z0EFF_1D',0,ZHOOK_HANDLE)
IF (PRESENT(OMASK)) THEN
  GMASK=OMASK
ELSE
  GMASK=(PAOSIP/=XUNDEF)    ! computations always performed where SSO data exist
  PZ0EFFIP = XUNDEF
  PZ0EFFIM = XUNDEF
  PZ0EFFJP = XUNDEF
  PZ0EFFJM = XUNDEF
END IF
!
!*    1.     Computations from A/S and h/2
!            -----------------------------
!
DO JJ=1,SIZE(PHO2JP) 
  IF (GMASK(JJ)) THEN        
    CALL GET_Z0EFF(PZ0VEG(JJ),PHO2JP(JJ),PAOSJP(JJ),PZ0EFFJP(JJ))
    CALL GET_Z0EFF(PZ0VEG(JJ),PHO2JM(JJ),PAOSJM(JJ),PZ0EFFJM(JJ))
    CALL GET_Z0EFF(PZ0VEG(JJ),PHO2IM(JJ),PAOSIM(JJ),PZ0EFFIM(JJ))
    CALL GET_Z0EFF(PZ0VEG(JJ),PHO2IP(JJ),PAOSIP(JJ),PZ0EFFIP(JJ))
  ENDIF
!
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SUBSCALE_Z0EFF_1D',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
CONTAINS
!
SUBROUTINE GET_Z0EFF(PZ0,PHO,PAO,PZ0EFF)
!
IMPLICIT NONE
!
REAL, INTENT(IN) :: PZ0
REAL, INTENT(IN) :: PHO
REAL, INTENT(IN) :: PAO
REAL, INTENT(OUT):: PZ0EFF
!
REAL :: ZLOC1,ZLOC2,ZLOC3
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SUBSCALE_Z0EFF_1D:GET_ZOEFF',0,ZHOOK_HANDLE)
!
IF ( PHO > PZ0 .AND. (PZ0.NE.0..OR.PAO.NE.0.)) THEN 
  ZLOC1  = XCDZ0EFF/(2.*XKARMAN**2)*PAO
  IF ( PZ0 > 0. ) THEN
    ZLOC2 = 1./(ALOG(PHO/PZ0))**2
  ELSE
    ZLOC2 = 0.
  ENDIF 
  ZLOC3  = SQRT(1./(ZLOC1+ZLOC2))
  PZ0EFF = PHO * EXP(-ZLOC3)
ELSE
  PZ0EFF = PZ0 
ENDIF
!
IF (LHOOK) CALL DR_HOOK('SUBSCALE_Z0EFF_1D:GET_ZOEFF',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_Z0EFF
! 
END SUBROUTINE SUBSCALE_Z0EFF_1D
