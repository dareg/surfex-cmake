!     #########
      SUBROUTINE ICE_SOILDIF(HDIF, PTSTEP, PTAUICE, PKSFC_IVEG,                  &
                            PSOILHCAPZ, PWSATZ, PMPOTSATZ, PBCOEFZ,              &
                            PTG, PWGI, PWG, PDELTAT, PALPHA, PN, PM,             &
                            PWRES                                                )  
!     ##########################################################################
!
!!****  *ICE_SOILDIF*  
!
!!    PURPOSE
!!    -------
!     This subroutine calculates soil water phase changes using the
!     available/excess energy approach. Soil temperature and volumetric
!     ice content are adjusted due to phase changes. See the references
!     Boone et al., 39, JAM, 2000 and Giard and Bazile, 128, MWR, 2000.
!     NOTE that more recently a modification was made: freeze/thaw follows
!     a relationship between liquid water and temperature derriving from
!     the Clausius Clapeyron Eq. This results in little to no freezing for
!     sufficiently dry but cold (below freezing) soils. Scatter about this
!     curve results due to 'phase change efficiencies' and the surface insolation
!     coefficient.
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
!!    Boone (2000)
!!    Boone et al., (2000)
!!      
!!    AUTHOR
!!    ------
!!	A. Boone          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/02/00   Boone
!!      Modified    24/11/09   Boone
!!                             Limit energy available for phase change by
!                              local amount owing to diffusion. Has almost
!                              no impact except under rare circumstances
!                              (avoids rare but possible oscillatory behavior)
!                              Also, add minimum (numerical) melt/freeze efficieny to prevent
!                              prolonged periods of small ice amounts approaching zero.
!      Modified    05/2010     Decharme
!                              Possibility to use Brook and Corey or Van Genuchten
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,     ONLY : XLMTT, XTT, XG, XCI, XRHOLI, XRHOLW
USE MODD_ISBA_PAR, ONLY : XWGMIN
!
USE MODE_HYDRO_DIF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP  ! Model time step (s)
!
CHARACTER(LEN=*),  INTENT(IN)      :: HDIF  ! NOTE: Only used when HISBA = DIF
!                                           ! 'BC' = Brook and Corey
!                                           ! 'VG' = Van Genuchten
!
REAL, DIMENSION(:), INTENT(IN)      :: PTAUICE, PKSFC_IVEG
!                                      PKSFC_IVEG = effect of surface layer insolation on phase changes
!                                                    Giard and Bazile (2000): non-dimensional
!                                      PTAUICE    = soil phase change characteristic time scale (s)
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSOILHCAPZ, PWSATZ, PDELTAT
!                                      PSOILHCAPZ = soil heat capacity [J/(m3 K)]
!                                      PWSATZ     = soil porosity (m3/m3)
!                                      PDELTAT    = change in temperature over the time
!                                                   step before adjustment owing to phase 
!                                                   changes (K)
REAL, DIMENSION(:,:), INTENT(IN)    :: PMPOTSATZ, PBCOEFZ
!                                      PMPOTSATZ  = matric potential at saturation (m)
!                                      PBCOEFZ    = slope of the water retention curve (-)
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTG, PWGI, PWG
!                                      PTG    = soil temperature (K)
!                                      PWGI   = soil volumetric ice content (m3/m3)
!                                      PWGI   = soil volumetric liquid water content (m3/m3)
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PALPHA      ! Van Genuchten parameter
REAL, DIMENSION(:,:), INTENT(IN)    :: PN          ! Van Genuchten parameter
REAL, DIMENSION(:,:), INTENT(IN)    :: PM          ! Van Genuchten parameter
REAL, DIMENSION(:,:), INTENT(IN)    :: PWRES       ! Van Genuchten parameter
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ   ! loop control
!
INTEGER                             :: INLVLD  ! Number of explicit soil layers
!
REAL, DIMENSION(SIZE(PTG,1),SIZE(PTG,2)) :: ZDELTAT, ZPHASE, ZTGM,            &
                                              ZWGM, ZWGIM, ZEFFIC, ZK, ZTAUICE, &
                                              ZDELTATI  
!
REAL, DIMENSION(SIZE(PTG,1),SIZE(PTG,2)) :: ZWGMAX, ZPSIMAX
!
REAL, DIMENSION(SIZE(PTG,1),SIZE(PTG,2)) :: ZPHASEM, ZPHASEF, ZTGMAX, ZPSI
!
!*      0.3    declarations of local parameters
!
REAL, PARAMETER                          :: ZTWGHT     = 0.50 ! (-)   (0 < ZTWGHT <= 1/2)
!                                                             ! Weight for averaging the actual and flux corrected
!                                                             ! temperature depressions. Default is 1/2
!                                                             ! If ZTWGHT=0, then flux correction is OFF.
REAL, PARAMETER                          :: ZEFFIC_MIN = 0.01 ! (-)   (0 <= ZEFFIC_MIN << 1)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                                             ! This parameter ensures
!                                                             ! a small minimum melt or freeze efficiency...
!                                                             ! It is numerical. When it is small, it has
!                                                             ! a only small impact on results, except
!                                                             ! that it keeps very small values of ice from persisting
!                                                             ! over long periods of time as they approach zero.
!                                                             ! If it is zero, then this effect off.
!-------------------------------------------------------------------------------
!
! Initialization:
! ---------------
!
IF (LHOOK) CALL DR_HOOK('ICE_SOILDIF',0,ZHOOK_HANDLE)
ZWGIM(:,:) = PWGI(:,:)
ZWGM(:,:)  = PWG(:,:)
ZTGM(:,:)  = PTG(:,:)
!
INLVLD     = SIZE(PTG(:,:),2)
!
!-------------------------------------------------------------------------------
!
! Eventually, could make PTAUICE a function of
! texture (therefore possibly soil depth):
!
DO JJ=1,INLVLD
   ZTAUICE(:,JJ) = PTAUICE(:)
ENDDO
!
! Surface layer vegetation insulation coefficient (-)
! and effective heat capacity [J/(m3 K)]:
!
ZK(:,1)        = PKSFC_IVEG(:)
ZK(:,2:INLVLD) = 1.0
!
!
! Begin soil ice evolution computations:
! --------------------------------------
!
! The maximum liquid water content as
! as function of temperature (sub-freezing)
! based on Gibbs free energy (Fuchs et al., 1978):
!
IF(HDIF=='BC')THEN
  ZPSIMAX(:,:)=PMPOTSATZ(:,:)
ELSE
  ZPSIMAX(:,:)=0.0
ENDIF
!
ZPSIMAX(:,:) = MIN(ZPSIMAX(:,:), XLMTT*(PTG(:,:)-XTT)/(XG*PTG(:,:)) )
ZWGMAX(:,:)  = W_FUNC(HDIF,ZPSIMAX,PWSATZ,PMPOTSATZ,PBCOEFZ,PALPHA,PN,PM,PWRES)
!
! Calculate maximum temperature for ice based on Gibbs free energy: first
! compute soil water potential using Clapp and Hornberger (1978) model:
!
ZPSI(:,:)     = PSI_FUNC(HDIF,PWG,PWSATZ,PMPOTSATZ,PBCOEFZ,PALPHA,PN,PM,PWRES)
ZTGMAX(:,:)   = XLMTT*XTT/(XLMTT - XG*ZPSI(:,:))
ZDELTATI(:,:) = PTG(:,:) - ZTGMAX(:,:)           ! initial temperature depression 
!
!
! Limit phase change energy by energy actually available for phase changes:
! (this has little effect generally but was found to prevent oscillations
! between soil T and water/ice in rare circumstances, notably for large time steps).
! We take an average of this flux-limited value and actual
! (one could eventually put in a time dependence on the weights used in averaging,
! perhaps giving larger weight to limited flux for larger time steps...for now
! just a constant factor assumed)
!
ZDELTAT(:,:)  = ZDELTATI(:,:)
ZDELTAT(:,:)  = MIN(ZDELTAT(:,:), MAX(0.0,PDELTAT(:,:)))
ZDELTAT(:,:)  = MAX(ZDELTAT(:,:), MIN(0.0,PDELTAT(:,:)))
ZDELTAT(:,:)  = ZTWGHT*ZDELTAT(:,:) + (1.0-ZTWGHT)*ZDELTATI(:,:)
WHERE(PDELTAT(:,:)*ZDELTAT(:,:) < 0.0)ZDELTAT(:,:) = 0.0
!
! NOTE that making the terms ZEFFIC(:,:) ==> 1 (and ZK also) (below) has the effect of
! reducing scatter about Fuch's T-depression curve. 
!
! *Melt* ice if energy and ice available:
!
ZEFFIC(:,:)   = MAX(ZEFFIC_MIN, MIN(1.0, ZWGIM(:,:)/(PWSATZ(:,:)-XWGMIN)))
ZPHASEM(:,:)  = (PTSTEP/ZTAUICE(:,:))*MIN(ZK(:,:)*ZEFFIC(:,:)*XCI*XRHOLI*MAX(0.,ZDELTAT(:,:)),   &
                    ZWGIM(:,:)*XLMTT*XRHOLW)  
!
! *Freeze* liquid water if energy and water available:
!
ZEFFIC(:,:)   = MAX(ZEFFIC_MIN, ZWGM(:,:)-XWGMIN)/PWSATZ(:,:)
ZPHASEF(:,:)  = (PTSTEP/ZTAUICE(:,:))*MIN(ZK(:,:)*ZEFFIC(:,:)*XCI*XRHOLI*MAX(0.0,-ZDELTAT(:,:)), &
                   MAX(0.0,ZWGM(:,:)-ZWGMAX(:,:))*XLMTT*XRHOLW)  
!
!
! Update ice, liquid and heat content if melting or freezing
!
ZPHASE(:,:) = ZPHASEF(:,:)  - ZPHASEM(:,:)                   ! (J m-3)
PWGI(:,:)   = ZWGIM(:,:)    + ZPHASE(:,:)/(XLMTT*XRHOLW)     ! (m3/m3)
PWG(:,:)    = ZWGM(:,:)     - ZPHASE(:,:)/(XLMTT*XRHOLW)     ! (m3/m3)
PTG(:,:)    = ZTGM(:,:)     + ZPHASE(:,:)/PSOILHCAPZ(:,:)    ! (K)
!
!    
! Prevent keeping track of very small numbers for ice content: (melt it)
! and conserve energy:
!
WHERE(PWGI < 1.0E-10)
   PWG(:,:)  = PWG(:,:) + PWGI(:,:)
   PTG(:,:)  = PTG(:,:) - PWGI(:,:)*XLMTT*XRHOLW/PSOILHCAPZ(:,:)
   PWGI(:,:) = 0.0
END WHERE
IF (LHOOK) CALL DR_HOOK('ICE_SOILDIF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ICE_SOILDIF
