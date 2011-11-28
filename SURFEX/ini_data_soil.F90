!     #########
      SUBROUTINE INI_DATA_SOIL(HISBA,PSURF,PSURF2,&
        PROOT_DEPTH,PGROUND_DEPTH,PROOT_EXTINCTION,PROOT_LIN,&
        PDG_IN,PDG_OUT,PCUM_ROOT_FRAC)
!     #########################
!
!!**** *INI_DATA_SOIL* initializes soil depth and root fraction for a given
!!                     number of soil layers
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    01/04/2003
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODI_SOILGRID
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=*), INTENT(IN) :: HISBA   ! type of soil (Force-Restore OR Diffusion)
REAL, DIMENSION(:), OPTIONAL, INTENT(IN) :: PSURF
REAL, DIMENSION(:), OPTIONAL, INTENT(IN) :: PSURF2
!
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PROOT_DEPTH
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PGROUND_DEPTH
REAL, DIMENSION(:,:), INTENT(IN) :: PROOT_EXTINCTION
REAL, DIMENSION(:,:), INTENT(IN) :: PROOT_LIN
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(IN) :: PDG_IN
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT) :: PDG_OUT
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PCUM_ROOT_FRAC


!
!*    0.2    Declaration of local variables
!      ------------------------------
!
!
LOGICAL            :: LSURF
INTEGER            :: JLOOP    ! class loop counter
INTEGER            :: JLAYER   ! soil layer loop counter
INTEGER            :: JVEG     ! vegetation types loop counter
!
REAL,DIMENSION(SIZE(PCUM_ROOT_FRAC,1),SIZE(PCUM_ROOT_FRAC,2),SIZE(PCUM_ROOT_FRAC,3)) :: ZDG
REAL               :: ZJACKSON ! Jackson (1996) formulation for cumulative root fraction
REAL               :: ZUNIF    ! linear formulation for cumulative root fraction
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                              ! (i.e. uniform root fraction)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*    1.     Allocations
!            -----------
!
IF (LHOOK) CALL DR_HOOK('INI_DATA_SOIL',0,ZHOOK_HANDLE)
!
ZDG(:,:,:) = XUNDEF
!
PCUM_ROOT_FRAC = XUNDEF
IF (PRESENT(PDG_OUT)) PDG_OUT = XUNDEF
!
!-------------------------------------------------------------------------------
!
!*    2.     loop on cover types
!            -------------------
!
LSURF=.FALSE. 
!
DO JLOOP=1,SIZE(PCUM_ROOT_FRAC,1)
! 
  IF (PRESENT(PSURF2) .AND. PRESENT(PSURF)) THEN
    LSURF=(PSURF(JLOOP)==0. .AND. PSURF2(JLOOP)==0.)
  ELSEIF (PRESENT(PSURF)) THEN
    LSURF=(PSURF(JLOOP)==0.)
  ENDIF
!
!-------------------------------------------------------------------------------
!
!* nature exists for this cover type
!
  IF (LSURF) CYCLE
!
!-------------------------------------------------------------------------------
!
!*    3.     soil depth
!            ----------
!
!*    3.1    force-restore case (2 layers)
!            ------------------
  IF (PRESENT(PDG_IN)) THEN

    ZDG(JLOOP,:,:)=PDG_IN(JLOOP,:,:)

  ELSE
!
    IF (HISBA=='2-L') THEN
      WHERE(PROOT_DEPTH(JLOOP,:) /= XUNDEF)
        ZDG(JLOOP,1,:) = 0.01
        ZDG(JLOOP,2,:) = PROOT_DEPTH(JLOOP,:)
      END WHERE
!
!
!*    3.2    force-restore case (3 layers)
!            ------------------
!
    ELSEIF (HISBA=='3-L') THEN
      WHERE(PGROUND_DEPTH(JLOOP,:) /= XUNDEF)
        ZDG(JLOOP,1,:) = 0.01
        ZDG(JLOOP,2,:) = PROOT_DEPTH  (JLOOP,:)
        ZDG(JLOOP,3,:) = PGROUND_DEPTH(JLOOP,:)
      END WHERE
!
!
!*    3.3    Diffusion case (at least 4 soil layers)
!            --------------
!
    ELSE
      DO JVEG=1,SIZE(PGROUND_DEPTH,2)
        IF (PGROUND_DEPTH(JLOOP,JVEG) ==XUNDEF) CYCLE
        CALL SOILGRID(PGROUND_DEPTH(JLOOP:JLOOP,JVEG),ZDG(JLOOP:JLOOP,:,JVEG))
      END DO
    END IF

    IF (PRESENT(PDG_OUT)) PDG_OUT(JLOOP,:,:)=ZDG(JLOOP,:,:)

  ENDIF
!
!-------------------------------------------------------------------------------
!
!*    4.     Cumulative root fraction (diffusion case only)
!            ------------------------
!
    IF (HISBA=='DIF') THEN

      DO JVEG=1,SIZE(PGROUND_DEPTH,2)
        IF (PGROUND_DEPTH(JLOOP,JVEG) ==XUNDEF) CYCLE
        !
        DO JLAYER=1,SIZE(PCUM_ROOT_FRAC,2)
          !
          ZJACKSON = (1.-PROOT_EXTINCTION(JLOOP,JVEG)**(ZDG (JLOOP,JLAYER,JVEG)*100.)) &
                     / (1.-PROOT_EXTINCTION(JLOOP,JVEG)**(PROOT_DEPTH(JLOOP,JVEG)*100.))  
          ZJACKSON = MIN(ZJACKSON, 1.)
          !
          ZUNIF    = ZDG(JLOOP,JLAYER,JVEG)/PROOT_DEPTH(JLOOP,JVEG)
          ZUNIF    = MIN(ZUNIF, 1.)
          !
          PCUM_ROOT_FRAC(JLOOP,JLAYER,JVEG) = &
                          PROOT_LIN(JLOOP,JVEG)  * ZUNIF    &
                          +(1.-PROOT_LIN(JLOOP,JVEG)) * ZJACKSON
                
      END DO
    END DO
  END IF
!
!-------------------------------------------------------------------------------
!
END DO
IF (LHOOK) CALL DR_HOOK('INI_DATA_SOIL',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_DATA_SOIL
