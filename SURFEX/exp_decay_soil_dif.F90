!     #########################
SUBROUTINE EXP_DECAY_SOIL_DIF (PD_G,PROOTFRAC,PF,PCONDSAT_EXP,PEXP_DIF)
!     ##########################################################################
!
!!****  *EXP_DECAY_SOIL_DIF*  
!!
!!    PURPOSE
!!    -------
!
!     We caculate the hydraulic coductivity decay factor for each interfacial conductivity (for diffusion option).
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
!!	B. Decharme     
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/11/03 
!-------------------------------------------------------------------------------
!
USE MODD_SGH_PAR, ONLY : X2,XF_DECAY
USE MODD_SURF_PAR,ONLY : XUNDEF
!
!*      0.1    declarations of arguments
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:),INTENT(IN   ) :: PD_G          !layer depth
REAL, DIMENSION(:,:),INTENT(IN   ) :: PROOTFRAC     !root fraction
REAL, DIMENSION(:),  INTENT(INOUT) :: PF            !Exponential decay factor
REAL, DIMENSION(:,:),INTENT(INOUT) :: PCONDSAT_EXP  !exponential hydraulic conductivity at saturation (m s-1)
REAL, DIMENSION(:,:),INTENT(INOUT) :: PEXP_DIF      !interfacial exponential profile coeff
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PF))                   :: ZD_INTER, ZC_DEPTH
REAL, DIMENSION(SIZE(PF),SIZE(PD_G(:,:),2)) :: ZD_MID
!
INTEGER                       :: I, INLVLD, IDEPTH, INI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('EXP_DECAY_SOIL_DIF',0,ZHOOK_HANDLE)
INLVLD=SIZE(PD_G(:,:),2)
INI=SIZE(PF)
!
!-------------------------------------------------------------------------------
!
!Mid point depth
!
ZD_MID(:,:)=XUNDEF
!
WHERE(PD_G(:,1)/=XUNDEF)ZD_MID(:,1)=PD_G(:,1)/X2
!
DO IDEPTH=2,INLVLD
  WHERE(PD_G(:,IDEPTH)/=XUNDEF)   
    ZD_MID(:,IDEPTH)=(PD_G(:,IDEPTH-1)+PD_G(:,IDEPTH))/X2
  ENDWHERE
ENDDO
!
!-------------------------------------------------------------------------------
!
!depth where the vertical satured hydraulic conductivities reach
!the compacted value given in Clapp and Hornberger (root depth)
!
ZC_DEPTH=0.0
!
!For instance, We use PROOTFRAC but after the real root depth "d2" will have to be specified
DO I=1,INI
   IF(PD_G(I,1)/=XUNDEF)THEN
     DO IDEPTH=1,INLVLD
        IF(PROOTFRAC(I,IDEPTH)==1.0)THEN
          ZC_DEPTH(I)= PD_G(I,IDEPTH)
          EXIT
        ENDIF
     ENDDO
   ENDIF
ENDDO
!
!No vegetation specific case
WHERE(ZC_DEPTH(:)==0.0.AND.PD_G(:,1)/=XUNDEF) ZC_DEPTH(:) = PD_G(:,INLVLD)   
!   
!-------------------------------------------------------------------------------
!
WHERE(ZC_DEPTH(:)/=0.0)PF(:)=4.0/ZC_DEPTH(:)
PF(:)=MIN(PF(:),XF_DECAY)
!
!-------------------------------------------------------------------------------
! Exponential conductivity of heach mid point layer 
!-------------------------------------------------------------------------------
!
DO IDEPTH=1,INLVLD
   WHERE(ZD_MID(:,IDEPTH)/=XUNDEF)             
         PCONDSAT_EXP (:,IDEPTH) = PCONDSAT_EXP (:,IDEPTH) * EXP(-(ZD_MID(:,IDEPTH)-ZC_DEPTH(:))*PF(:))   
   ENDWHERE
ENDDO
!
!-------------------------------------------------------------------------------
! Interfacial conductivity coef 
! 1      -> First interblock
! 2      -> Second interblock 
! ...
! INLVLD -> Ending mind point layer (if no water table)
!-------------------------------------------------------------------------------
!   
!Interfacial layers
!
DO IDEPTH=1,INLVLD-1
!
  WHERE(ZD_MID(:,IDEPTH)/=XUNDEF)   
!           
!  here, ZD_INTER is the depth to layer interface (interblock depth)
!
    ZD_INTER (:) = (ZD_MID(:,IDEPTH)+ZD_MID(:,IDEPTH+1))/X2
!
    PEXP_DIF (:,IDEPTH) = EXP(-(ZD_INTER(:)-ZC_DEPTH(:))*PF(:))
!   
  ENDWHERE
!
ENDDO
!
!Last layer (fluxes are calculated at the Mid point depth)
!
PEXP_DIF (:,INLVLD) = EXP(-(ZD_MID(:,INLVLD)-ZC_DEPTH(:))*PF(:))
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('EXP_DECAY_SOIL_DIF',1,ZHOOK_HANDLE)
!
END SUBROUTINE EXP_DECAY_SOIL_DIF





