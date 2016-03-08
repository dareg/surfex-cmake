!     ##########
      SUBROUTINE EXP_DECAY_SOIL_FR (HISBA, PF, IP, MX    )  
!     ##########################################################################
!
!!****  *EXP_DECAY_SOIL_FR*  
!!
!!    PURPOSE
!!    -------
!
!     We caculate the hydraulic coductivity decay factor for each FR-coefficients.
!     Also, we redefine the surface hydraulic coductivity at saturation for
!     convective precipitation parametrisation.
!     
!!**  METHOD
!!    ------
!
!     Direct calculation
!
!!    EXTERNAL
!!    --------
!
!     None
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
!!      B. Decharme     
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/11/03 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t
!
USE MODD_SURF_PAR,ONLY : XUNDEF
USE MODD_SGH_PAR, ONLY : X2                                
USE MODD_CSTS,    ONLY : XDAY
#ifdef TOPD
USE MODD_DUMMY_EXP_PROFILE,ONLY : XC_DEPTH_RATIO
#endif
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=*)                  :: HISBA   ! hydrology/soil:
!                                            ! '2-L'  = single column
!                                            ! '3-L'  = root zone/baseflow layer
!                                            ! 'DIF'  = N-layer diffusion: Richard's Eq.
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PF
!                                    PF = exponential decay factor (1/m)
!
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: MX
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PF,1))         :: ZD_G_TOT, ZC_DEPTH, ZKSAT_NOEXP, ZC_DEPTH_RATIO
!                                    ZD_G_TOT = depth of the soil column (m)
!                                    ZC_DEPTH = assumed as the depth where the vertical 
!                                               satured hydraulic conductivities reach
!                                               the compacted value given in Clapp and
!                                               Hornberger. (m)
!                                               For ISBA-FR, we take the root depth.
!
INTEGER :: JP
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('EXP_DECAY_SOIL_FR',0,ZHOOK_HANDLE)
!
DO JP=1,SIZE(MX%XDG,3)
  !
  IF (IP%NSIZE_NATURE_P(JP) == 0 ) CYCLE
  !
  ZD_G_TOT(:) = MX%XDG(:,2,JP)
  IF(HISBA=='3-L')ZD_G_TOT(:) = MX%XDG(:,3,JP)
  !
  ZKSAT_NOEXP(:) = IP%XCONDSAT(:,2,JP)
  !
  ZC_DEPTH_RATIO(:) = 1.
  !
#ifdef TOPD
  IF (ALLOCATED(XC_DEPTH_RATIO)) ZC_DEPTH_RATIO(:) = XC_DEPTH_RATIO(:)
#endif
  !
  WHERE(ZD_G_TOT(:)/=XUNDEF)
    !
    !compacted depth
    !
    ZC_DEPTH(:) = MX%XDG(:,2,JP)*ZC_DEPTH_RATIO(:)
    !ZC_DEPTH(:) = MX%XDG(:,2)
    !
    !surface hydraulic conductivity at saturation
    !
    IP%XCONDSAT(:,1,JP) = IP%XCONDSAT(:,1,JP)*EXP(PF(:,JP)*ZC_DEPTH(:))
    !
    !mean hydraulic conductivity at saturation over the root zone
    !   
    IP%XCONDSAT(:,2,JP) = ZKSAT_NOEXP(:)*( EXP(PF(:,JP)*ZC_DEPTH)-EXP(PF(:,JP)*(ZC_DEPTH(:)-MX%XDG(:,2,JP))) )   &
                      /(PF(:,JP)*MX%XDG(:,2,JP))
    !   
    !mean hydraulic conductivity at saturation over the first soil centimeters
    !   
    IP%XKSAT_ICE(:,JP) = ZKSAT_NOEXP(:)*( EXP(PF(:,JP)*ZC_DEPTH)-EXP(PF(:,JP)*(ZC_DEPTH(:)-MX%XD_ICE(:,JP))) )   &
                     /(PF(:,JP)*MX%XD_ICE(:,JP))  
    !
    !decay factor for C1 coef
    !   
    IP%XC1SAT(:,JP) = IP%XC1SAT(:,JP)*SQRT( EXP(-PF(:,JP)*ZC_DEPTH(:)) )
    !
    !decay factor for C2 coef 
    !
    IP%XC2REF(:,JP)=IP%XC2REF(:,JP)+( IP%XCONDSAT(:,2,JP)-ZKSAT_NOEXP(:) ) * XDAY/MX%XDG(:,2,JP) 
    !
    !C3 coef with exponential decay in root soil layer 
    !
    IP%XC3(:,1,JP)=IP%XC3(:,1,JP)*( EXP(PF(:,JP)*ZC_DEPTH(:))-EXP(PF(:,JP)*(ZC_DEPTH(:)-MX%XDG(:,2,JP))) ) / &
            (PF(:,JP)*MX%XDG(:,2,JP))
    !
  ENDWHERE
  !
  IF(HISBA=='3-L')THEN
    ! 
    WHERE(MX%XDG(:,2,JP)< ZD_G_TOT(:).AND.MX%XDG(:,2,JP)/=XUNDEF)
      !           
      !  C3 coef with exponential decay in deep soil layer 
      !
      IP%XC3(:,2,JP)=IP%XC3(:,2,JP)*( EXP(PF(:,JP)*(ZC_DEPTH(:)-MX%XDG(:,2,JP)))-EXP(PF(:,JP)*(ZC_DEPTH(:)-ZD_G_TOT(:))) )      &
                       / (PF(:,JP)*(ZD_G_TOT(:)-MX%XDG(:,2,JP)))  
      ! 
      !  decay factor for C4 coef
      !      
      IP%XC4REF(:,JP)=IP%XC4REF(:,JP)*( EXP(PF(:,JP)*(ZC_DEPTH(:)-MX%XDG(:,2,JP)/X2))-EXP(PF(:,JP)*(ZC_DEPTH(:)&
                           -((MX%XDG(:,2,JP)+ZD_G_TOT(:))/2.))) ) * X2/(PF(:,JP)*ZD_G_TOT(:))        
      !
    ENDWHERE
    !
  ENDIF
  !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('EXP_DECAY_SOIL_FR',1,ZHOOK_HANDLE)
!
END SUBROUTINE EXP_DECAY_SOIL_FR
