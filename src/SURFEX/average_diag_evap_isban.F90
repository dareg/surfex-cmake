!#############################
SUBROUTINE AVERAGE_DIAG_EVAP_ISBA_n (OSURF_BUDGETC, DE, DEC, DEP, DEPC, NP, KNPATCH, &
                                     OGLACIER, OMEB_PATCH, PTSTEP, PRAIN, PSNOW)
!#############################
!
!
!!****  *AVERAGE_DIAG_EVAP_ISBA_n*  
!!
!!    PURPOSE
!!    -------
!      Average the cumulated diagnostics from all ISBA tiles
!     
!!**  METHOD
!!    ------
!
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
!!      P. Le Moigne           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/03
!!      B. Decharme 2008     New diag for the water budget
!!      B. Decharme 2012     New diag for snow 
!!                                        carbon
!!                                        isab water budget
!!                  2013                  Sublimation
!!                                        Subsurface runoff if SGH (DIF option only)
!!      P. Samuelsson 10/2014: MEB
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODE_DIAG
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t, DIAG_EVAP_ISBA_PATCH_t
USE MODD_ISBA_n, ONLY : ISBA_NP_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
LOGICAL, INTENT(IN) :: OSURF_BUDGETC
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DE
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEC 
TYPE(DIAG_EVAP_ISBA_PATCH_t), INTENT(INOUT) :: DEP
TYPE(DIAG_EVAP_ISBA_PATCH_t), INTENT(INOUT) :: DEPC
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
INTEGER, INTENT(IN) :: KNPATCH
!
LOGICAL, INTENT(IN) :: OGLACIER
LOGICAL, DIMENSION(:), INTENT(IN) :: OMEB_PATCH
!
REAL,                  INTENT(IN) :: PTSTEP        ! time step (s)
REAL,    DIMENSION(:), INTENT(IN) :: PRAIN         ! rainfall rate
REAL,    DIMENSION(:), INTENT(IN) :: PSNOW         ! snowfall rate
!
!
!*      0.2    declarations of local variables
!
INTEGER :: JP ! tile loop counter
INTEGER :: JI, IMASK
INTEGER           :: ISIZE_LMEB_PATCH   ! Number of patches where multi-energy balance should be applied
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!       0.     Initialization
!              --------------
!
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_EVAP_ISBA_N',0,ZHOOK_HANDLE)
!
ISIZE_LMEB_PATCH=COUNT(OMEB_PATCH(:))
!
!       1.     Surface Energy fluxes
!              -----------------------
!
IF (DE%LSURF_EVAP_BUDGET) THEN
  !
  CALL INIT_EVAP_BUD(DE)
  IF (ISIZE_LMEB_PATCH>0) CALL INIT_MEB_BUD(DE)
  !
  IF(DE%LWATER_BUDGET)THEN
    !  
    CALL INIT_WATER_BUD(DE)
    !
    DE%XRAINFALL  (:) = PRAIN(:) * PTSTEP
    DE%XSNOWFALL  (:) = PSNOW(:) * PTSTEP
    !
  ENDIF
  !
  CALL MAKE_AVERAGE_EVAP(DE,DEP)  
  ! 
  ! Ice calving flux
  !  
  IF(OGLACIER)THEN 

    DE%XICEFLUX(:)= 0.
    DO JP=1,KNPATCH
      DO JI=1,NP%AL(JP)%NSIZE_P
        IMASK = NP%AL(JP)%NR_P(JI)
        DE%XICEFLUX(IMASK) = DE%XICEFLUX(IMASK) + NP%AL(JP)%XPATCH(JI) * DEP%AL(JP)%XICEFLUX(JI)      
      END DO
    END DO

  END IF
  !  
ENDIF
!
!       2.     Surface Cumulated Energy fluxes
!              -------------------------------
!
IF (OSURF_BUDGETC) THEN
  !
  !
  CALL INIT_EVAP_BUD(DEC)
  !
  IF (ISIZE_LMEB_PATCH>0) THEN
    CALL INIT_MEB_BUD(DEC)
  ENDIF
  !
  ! Isba water budget and reservoir time tendencies
  !
  IF(DE%LWATER_BUDGET)THEN
    !  
    CALL INIT_WATER_BUD(DEC)
    !
    DEC%XRAINFALL  (:) = DEC%XRAINFALL (:) + PRAIN(:) * PTSTEP
    DEC%XSNOWFALL  (:) = DEC%XSNOWFALL (:) + PSNOW(:) * PTSTEP
    !
  ENDIF
  !
  CALL MAKE_AVERAGE_EVAP(DEC,DEPC)
  !
  ! Ice calving flux
  !  
  IF(OGLACIER)THEN 

    DEC%XICEFLUX(:)= 0.
    DO JP=1,KNPATCH
      DO JI=1,NP%AL(JP)%NSIZE_P
        IMASK = NP%AL(JP)%NR_P(JI)
        DEC%XICEFLUX(IMASK) = DEC%XICEFLUX(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPC%AL(JP)%XICEFLUX(JI)      
      END DO
    END DO

  END IF
!  
END IF
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_EVAP_ISBA_N',1,ZHOOK_HANDLE)
!
CONTAINS
!
SUBROUTINE MAKE_AVERAGE_EVAP(DEA,DEPA)
!
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEA
TYPE(DIAG_EVAP_ISBA_PATCH_t), INTENT(INOUT) :: DEPA
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_EVAP_ISBA_N:MAKE_AVERAGE_EVAP',0,ZHOOK_HANDLE)
!
DO JP=1,KNPATCH
  DO JI=1,NP%AL(JP)%NSIZE_P
    IMASK = NP%AL(JP)%NR_P(JI)
    !
    ! Latent heat of evaporation over the ground
    DEA%XLEG (IMASK) = DEA%XLEG (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLEG(JI)
    !
    ! Surface soil ice sublimation
    DEA%XLEGI(IMASK) = DEA%XLEGI(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLEGI(JI)
    !
    ! Latent heat of evaporation over vegetation
    DEA%XLEV (IMASK) = DEA%XLEV (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLEV(JI)
    !
    ! Latent heat of sublimation over snow
    DEA%XLES (IMASK) = DEA%XLES (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLES(JI)
    !
    ! Latent heat of evaporation of liquid water over snow
    DEA%XLESL(IMASK) = DEA%XLESL(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLESL(JI)
    !
    ! Evaporation from canopy water interception
    DEA%XLER (IMASK) = DEA%XLER (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLER(JI)
    !
    ! Evapotranspiration of the vegetation
    DEA%XLETR(IMASK) = DEA%XLETR(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLETR(JI)
    !
    ! Blowing snow sublimation (ES or Crocus)
    DEA%XSNDRIFT(IMASK) = DEA%XSNDRIFT(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XSNDRIFT(JI)
    !
    ! Soil drainage flux
    DEA%XDRAIN (IMASK) = DEA%XDRAIN (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XDRAIN(JI)
    !
    ! Soil lateral subsurface flux
    DEA%XQSB   (IMASK) = DEA%XQSB   (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XQSB(JI)        
    !
    ! Supersaturation runoff
    DEA%XRUNOFF(IMASK) = DEA%XRUNOFF(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XRUNOFF(JI)
    !
    ! Horton runoff
    DEA%XHORT  (IMASK) = DEA%XHORT  (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XHORT(JI)
    !
    ! Vegetation dripping
    DEA%XDRIP  (IMASK) = DEA%XDRIP  (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XDRIP(JI)
    !
    ! Precipitation intercepted by the vegetation
    DEA%XRRVEG (IMASK) = DEA%XRRVEG (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XRRVEG(JI)
    !      
    ! Snow melt
    DEA%XMELT  (IMASK) = DEA%XMELT  (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XMELT(JI)
    !      
    ! Flood infiltartion
    DEA%XIFLOOD(IMASK) = DEA%XIFLOOD(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XIFLOOD(JI)
    !      
    ! Precipitation intercepted by the floodplains   
    DEA%XPFLOOD(IMASK) = DEA%XPFLOOD(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XPFLOOD(JI)
    !      
    ! Floodplains evaporation  
    DEA%XLE_FLOOD (IMASK) = DEA%XLE_FLOOD (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLE_FLOOD (JI)
    DEA%XLEI_FLOOD(IMASK) = DEA%XLEI_FLOOD(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLEI_FLOOD(JI)
    !      
    ! irrigation rate (as soil input)
    DEA%XIRRIG_FLUX(IMASK) = DEA%XIRRIG_FLUX(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XIRRIG_FLUX(JI)
    !
    ! Gross primary production
    DEA%XGPP      (IMASK) = DEA%XGPP      (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XGPP(JI)
    !
    ! Autotrophic respiration 
    DEA%XRESP_AUTO(IMASK) = DEA%XRESP_AUTO(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XRESP_AUTO(JI)
    !
    ! Ecosystem respiration
    DEA%XRESP_ECO (IMASK) = DEA%XRESP_ECO (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XRESP_ECO(JI)  
    !        
    IF (ISIZE_LMEB_PATCH>0) THEN
      DEA%XLEVCV    (IMASK) = DEA%XLEVCV    (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLEVCV(JI)
      DEA%XLESC     (IMASK) = DEA%XLESC     (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLESC(JI)
      DEA%XLETRCV   (IMASK) = DEA%XLETRCV   (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLETRCV(JI)
      DEA%XLELITTER (IMASK) = DEA%XLELITTER (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLELITTER(JI)
      DEA%XLELITTERI(IMASK) = DEA%XLELITTERI(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLELITTERI(JI)
      DEA%XDRIPLIT  (IMASK) = DEA%XDRIPLIT  (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XDRIPLIT(JI)
      DEA%XRRLIT    (IMASK) = DEA%XRRLIT    (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XRRLIT(JI)
      DEA%XLERCV    (IMASK) = DEA%XLERCV    (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLERCV(JI)
      DEA%XLE_C_A   (IMASK) = DEA%XLE_C_A   (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLE_C_A(JI)
      DEA%XLE_V_C   (IMASK) = DEA%XLE_V_C   (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLE_V_C(JI)
      DEA%XLE_G_C   (IMASK) = DEA%XLE_G_C   (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLE_G_C(JI)
      DEA%XLE_N_C   (IMASK) = DEA%XLE_N_C   (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLE_N_C(JI)
      DEA%XSWNET_V  (IMASK) = DEA%XSWNET_V  (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XSWNET_V(JI)
      DEA%XSWNET_G  (IMASK) = DEA%XSWNET_G  (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XSWNET_G(JI)
      DEA%XSWNET_N  (IMASK) = DEA%XSWNET_N  (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XSWNET_N(JI)
      DEA%XSWNET_NS (IMASK) = DEA%XSWNET_NS (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XSWNET_NS(JI)
      DEA%XLWNET_V  (IMASK) = DEA%XLWNET_V  (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLWNET_V(JI)
      DEA%XLWNET_G  (IMASK) = DEA%XLWNET_G  (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLWNET_G(JI)
      DEA%XLWNET_N  (IMASK) = DEA%XLWNET_N  (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLWNET_N(JI)
      DEA%XSWDOWN_GN(IMASK) = DEA%XSWDOWN_GN(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XSWDOWN_GN(JI)
      DEA%XLWDOWN_GN(IMASK) = DEA%XLWDOWN_GN(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XLWDOWN_GN(JI)
      DEA%XH_V_C    (IMASK) = DEA%XH_V_C    (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XH_V_C(JI)
      DEA%XH_G_C    (IMASK) = DEA%XH_G_C    (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XH_G_C(JI)
      DEA%XH_C_A    (IMASK) = DEA%XH_C_A    (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XH_C_A(JI)
      DEA%XH_N_C    (IMASK) = DEA%XH_N_C    (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XH_N_C(JI)
      DEA%XSR_GN    (IMASK) = DEA%XSR_GN    (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XSR_GN(JI)
      DEA%XMELTCV   (IMASK) = DEA%XMELTCV   (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XMELTCV(JI)
      DEA%XFRZCV    (IMASK) = DEA%XFRZCV    (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XFRZCV(JI)
    ENDIF
  END DO
ENDDO
!
! Isba water budget and reservoir time tendencies
!
IF(DE%LWATER_BUDGET)THEN
  !  
  DO JP=1,KNPATCH
  DO JI=1,NP%AL(JP)%NSIZE_P
    IMASK = NP%AL(JP)%NR_P(JI)
      !
      DEA%XDWG   (IMASK) = DEA%XDWG   (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XDWG   (JI)
      DEA%XDWGI  (IMASK) = DEA%XDWGI  (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XDWGI  (JI)
      DEA%XDWR   (IMASK) = DEA%XDWR   (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XDWR   (JI)
      DEA%XDSWE  (IMASK) = DEA%XDSWE  (IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XDSWE  (JI)
      DEA%XWATBUD(IMASK) = DEA%XWATBUD(IMASK) + NP%AL(JP)%XPATCH(JI) * DEPA%AL(JP)%XWATBUD(JI)

    ENDDO
  ENDDO
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_EVAP_ISBA_N:MAKE_AVERAGE_EVAP',1,ZHOOK_HANDLE)
!
END SUBROUTINE MAKE_AVERAGE_EVAP
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE_DIAG_EVAP_ISBA_n
