!     #########
      SUBROUTINE SOIL_HEATDIF(PTSTEP,PD_G,PSOILCONDZ,PSOILHCAPZ,PCT,             &
                               PTERM1,PTERM2,PTDEEP,PTG                            )  
!     ############################################################################
!
!!****  *SOIL_HEATDIF*  
!
!!    PURPOSE
!!    -------
!     This subroutine solves the 1-d diffusion of 'PTG' using a
!     layer averaged set of equations which are time differenced
!     using the backward-difference scheme (implicit). 
!     The eqs are solved rapidly by taking advantage of the
!     fact that the matrix is tridiagonal. This is a very
!     general routine and can be used for the 1-d diffusion of any 
!     quantity as long as the diffusity is not a function of the
!     quantity being diffused. Aaron Boone 8-98. Soln to the eq:
!
!                    dQ    d    dQ       
!                 c  -- =  -- k -- 
!                    dt    dx   dx    
!
!     where k = k(x) (thermal conductivity), c = c(x) (heat capacity)
!     Diffusivity is k/c
!
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
!!    USE MODD_PARAMETERS
!!    USE MODI_TRIDIAG_GROUND
!!      
!!    REFERENCE
!!    ---------
!!
!!    Boone (2000)
!!    Boone et al. (2000)
!!      
!!    AUTHOR
!!    ------
!!	A. Boone          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    16/02/00   Boone
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODI_TRIDIAG_GROUND
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
REAL, INTENT(IN)                    :: PTSTEP ! Model time step (s)
!
REAL, DIMENSION(:), INTENT(IN)      :: PCT, PTERM1, PTERM2, PTDEEP
!                                      PCT    = thermal inertia [(m2 K)/J]
!                                      PTERM1 = coefficient of linearization
!                                               of surface energy budget 
!                                      PTERM2 = coefficient of linearization
!                                               of surface energy budget 
!                                      PTDEEP = Deep soil temperature (prescribed)
!                                               which models heating/cooling from
!                                               below the diurnal wave penetration
!                                               (surface temperature) depth. If it
!                                               is FLAGGED as undefined, then the zero
!                                               flux lower BC is applied.
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSOILCONDZ, PSOILHCAPZ, PD_G
!                                      PSOILCONDZ = soil thermal conductivity [W/(K m)]
!                                      PSOILHCAPZ = soil heat capacity [J/(m3 K)]
!                                      PD_G       = depth of bottom of soil layer (m)
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTG
!                                      PTG    = soil temperature (K)
!
!
!*      0.2    declarations of local variables
!
INTEGER JJ
!
INTEGER                                        :: INLVLD ! Number of grid layers
!
REAL, DIMENSION(SIZE(PTG,1),SIZE(PTG,2)) :: ZTGM, ZDTERM, ZCTERM,   &
                                                    ZFRCV, ZAMTRX, ZBMTRX,     &
                                                    ZCMTRX  
!
REAL, DIMENSION(SIZE(PTG,1),SIZE(PTG,2)) :: ZWORK1, ZWORK2, ZDZDIF, ZDZG
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!-------------------------------------------------------------------------------
!
! Initialize local variables:
!
IF (LHOOK) CALL DR_HOOK('SOIL_HEATDIF',0,ZHOOK_HANDLE)
ZDTERM(:,:)    = 0.0
ZCTERM(:,:)    = 0.0
ZFRCV(:,:)     = 0.0
ZAMTRX(:,:)    = 0.0
ZBMTRX(:,:)    = 0.0
ZCMTRX(:,:)    = 0.0
ZTGM(:,:)   = PTG(:,:)
!
INLVLD         = SIZE(PTG(:,:),2)
!
!-------------------------------------------------------------------------------
!
!* soil layers thicknesses
!
ZDZG(:,1) = PD_G(:,1)
DO JJ=2,INLVLD
  ZDZG(:,JJ) = PD_G(:,JJ) - PD_G(:,JJ-1)
END DO
!
!-------------------------------------------------------------------------------
!
! Calculate tri-diagonal matrix coefficients:
!
DO JJ=1,INLVLD-1
   ZDZDIF(:,JJ)  = ZDZG(:,JJ) + ZDZG(:,JJ+1)
ENDDO
ZDZDIF(:,INLVLD) = ZDZG(:,INLVLD) 
!
ZWORK1(:,:)      =  ZDZG(:,:)  *PSOILCONDZ(:,:)
DO JJ=1,INLVLD-1
   ZWORK2(:,JJ)  =  ZDZG(:,JJ+1)*PSOILCONDZ(:,JJ+1)
ENDDO
ZWORK2(:,INLVLD) = 0.0
!
ZDTERM(:,:)      = 2.0*(ZWORK1(:,:)+ZWORK2(:,:))/(ZDZDIF(:,:)*ZDZDIF(:,:))
!
ZCTERM(:,:)      = PSOILHCAPZ(:,:)*ZDZG(:,:)/PTSTEP
!
! - - -------------------------------------------------
!
! Upper BC
!
ZAMTRX(:,1) =  0.0
ZBMTRX(:,1) =  1./(PCT(:)*PTSTEP)
ZCMTRX(:,1) = -PTERM2(:)*ZBMTRX(:,1)
ZFRCV(:,1)  =  PTERM1(:)*ZBMTRX(:,1)
!
!
! Interior Grid
!
DO JJ=2,INLVLD-1
   ZAMTRX(:,JJ) = -ZDTERM(:,JJ-1) 
   ZBMTRX(:,JJ) =  ZCTERM(:,JJ) + ZDTERM(:,JJ-1) + ZDTERM(:,JJ)
   ZCMTRX(:,JJ) = -ZDTERM(:,JJ)
   ZFRCV(:,JJ)  =  ZCTERM(:,JJ)*ZTGM(:,JJ) 
ENDDO
!
! Lower BC: 2 currently accounted for, Either zero flux
! or a fixed temperature 'TDEEP' 
!
ZAMTRX(:,INLVLD) = -ZDTERM(:,INLVLD-1) 
ZCMTRX(:,INLVLD) =  0.0                
!
WHERE(PTDEEP(:) /= XUNDEF)
   ZBMTRX(:,INLVLD) =  ZCTERM(:,INLVLD) + ZDTERM(:,INLVLD-1) + ZDTERM(:,INLVLD)
   ZFRCV(:,INLVLD)  =  ZCTERM(:,INLVLD)*ZTGM(:,INLVLD) + ZDTERM(:,INLVLD)*PTDEEP(:)
ELSEWHERE
   ZBMTRX(:,INLVLD) =  ZCTERM(:,INLVLD) + ZDTERM(:,INLVLD-1) 
   ZFRCV(:,INLVLD)  =  ZCTERM(:,INLVLD)*ZTGM(:,INLVLD) 
END WHERE
!
! - - -------------------------------------------------
!
! Compute ZTGM (solution vector) 
! used for systems of equations involving tridiagonal 
! matricies.
!
CALL TRIDIAG_GROUND(ZAMTRX,ZBMTRX,ZCMTRX,ZFRCV,ZTGM)
!
!
! Update values in time:
!
PTG(:,:) = ZTGM(:,:)
IF (LHOOK) CALL DR_HOOK('SOIL_HEATDIF',1,ZHOOK_HANDLE)
!
!
!
END SUBROUTINE SOIL_HEATDIF
