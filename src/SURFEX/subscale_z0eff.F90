!     ######spl
      SUBROUTINE SUBSCALE_Z0EFF(MSS,PZ0VEG,OZ0REL,OMASK  )
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
USE MODD_CSTS,       ONLY : XKARMAN
USE MODD_ISBA_PAR,   ONLY : XCDZ0EFF
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODD_SSO_n, ONLY : SSO_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
TYPE(SSO_t), INTENT(INOUT) :: MSS
REAL, DIMENSION(:,:), INTENT(IN)  :: PZ0VEG  ! vegetation roughness length
!
LOGICAL, INTENT(IN) :: OZ0REL
LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: OMASK ! mask where computations
                                                       ! are done
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
REAL,    DIMENSION(SIZE(MSS%XAOSIP)) :: ZLOC
LOGICAL, DIMENSION(SIZE(MSS%XZ0EFFIM,1)) :: GMASK
!
INTEGER :: IPATCH  ! number of patches
INTEGER :: JPATCH  ! loop counter on number of patches
INTEGER :: JJ      ! loop counter on points
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_SUBSCALE_Z0EFF:SUBSCALE_Z0EFF_1D_PATCH',0,ZHOOK_HANDLE)
!
IF (.NOT.PRESENT(OMASK)) THEN
  MSS%XZ0EFFIP = XUNDEF
  MSS%XZ0EFFIM = XUNDEF
  MSS%XZ0EFFJP = XUNDEF
  MSS%XZ0EFFJM = XUNDEF
ENDIF
!
IPATCH = SIZE(PZ0VEG,2)
!----------------------------------------------------------------------------
DO JPATCH=1,IPATCH
!----------------------------------------------------------------------------
!
  IF (PRESENT(OMASK)) THEN
    GMASK=OMASK
  ELSEIF (ALL(PZ0VEG(:,:)==0.)) THEN
    GMASK = (MSS%XAOSIP/=XUNDEF)    ! computations always performed where SSO data exist
  ELSE
    GMASK=PZ0VEG(:,JPATCH) /= XUNDEF    ! computations always performed where defined
  END IF
!
!*    1.     Computations from A/S and h/2
!            -----------------------------
!   
!print*,'MASK ',GMASK(1)
!print*,'Z0VEG ',PZ0VEG(1,JPATCH)
!print*,'HO2JP ',MSS%XHO2JP(1)
!print*,'AOSJP ',MSS%XAOSJP(1)
! print*,'Z0EFFJP ',MSS%XZ0EFFJP(1,JPATCH)
 CALL GET_Z0EFF(GMASK(:),PZ0VEG(:,JPATCH),MSS%XHO2JP(:),MSS%XAOSJP(:),MSS%XZ0EFFJP(:,JPATCH))
! print*,'Z0EFFJP ',MSS%XZ0EFFJP(1,JPATCH)
 CALL GET_Z0EFF(GMASK(:),PZ0VEG(:,JPATCH),MSS%XHO2JM(:),MSS%XAOSJM(:),MSS%XZ0EFFJM(:,JPATCH))
 CALL GET_Z0EFF(GMASK(:),PZ0VEG(:,JPATCH),MSS%XHO2IM(:),MSS%XAOSIM(:),MSS%XZ0EFFIM(:,JPATCH))
 CALL GET_Z0EFF(GMASK(:),PZ0VEG(:,JPATCH),MSS%XHO2IP(:),MSS%XAOSIP(:),MSS%XZ0EFFIP(:,JPATCH))
!
END DO
!
IF (PRESENT(OMASK)) THEN
  GMASK=OMASK
ELSE
  GMASK=(MSS%XAOSIP/=XUNDEF)
END IF
!
IF (OZ0REL) THEN
        
  MSS%XZ0REL=XUNDEF
  !
  ZLOC(:) = 0.
  !
  WHERE (GMASK(:))
    ZLOC  (:) = 0.25 * XCDZ0EFF/(2.*XKARMAN**2)                  &
                     * (MSS%XAOSIP(:) + MSS%XAOSIM(:) + MSS%XAOSJP(:) + MSS%XAOSJM(:))        
    WHERE ( ZLOC(:) > 0. )
      MSS%XZ0REL(:) = 0.25 * (MSS%XHO2IP(:) + MSS%XHO2IM(:) + MSS%XHO2JP(:) + MSS%XHO2JM(:)) &
                       * EXP(-SQRT(1./ZLOC(:)))
      MSS%XZ0REL(:) = MAX(MSS%XZ0REL(:),1E-10)
    ELSEWHERE
      MSS%XZ0REL(:) = 0.
    END WHERE
  END WHERE      
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('SUBSCALE_Z0EFF',1,ZHOOK_HANDLE)
!
CONTAINS
!
SUBROUTINE GET_Z0EFF(OCOMPUT,PZ0,PHO,PAO,PZ0EFF)
!
IMPLICIT NONE
!
LOGICAL, DIMENSION(:), INTENT(IN) :: OCOMPUT
REAL,    DIMENSION(:), INTENT(IN) :: PZ0
REAL,    DIMENSION(:), INTENT(IN) :: PHO
REAL,    DIMENSION(:), INTENT(IN) :: PAO
REAL,    DIMENSION(:), INTENT(INOUT):: PZ0EFF
!
LOGICAL, DIMENSION(SIZE(PZ0)) :: LWORK1
!
REAL    :: ZLOC1,ZLOC2,ZLOC3
INTEGER :: JJ, INI
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SUBSCALE_Z0EFF:GET_ZOEFF',0,ZHOOK_HANDLE)
!
INI=SIZE(PZ0)
!
LWORK1(:)=(PHO(:)>PZ0(:).AND.(PZ0(:)/=0.0.OR.PAO(:)/=0.0))
!
DO JJ=1,INI
  IF (OCOMPUT(JJ)) THEN
    IF (LWORK1(JJ)) THEN 
      ZLOC1  = (XCDZ0EFF/(2.*XKARMAN**2))*PAO(JJ)
      IF ( PZ0(JJ) > 0. ) THEN
        ZLOC2 = 1./(ALOG(PHO(JJ)/PZ0(JJ)))**2
      ELSE
        ZLOC2 = 0.
      ENDIF 
      ZLOC3  = SQRT(1./(ZLOC1+ZLOC2))
      PZ0EFF(JJ) = PHO(JJ) * EXP(-ZLOC3)
    ELSE
      PZ0EFF(JJ) = PZ0(JJ) 
    ENDIF
  ENDIF
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SUBSCALE_Z0EFF:GET_ZOEFF',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_Z0EFF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SUBSCALE_Z0EFF
