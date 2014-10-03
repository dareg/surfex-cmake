!     ########################################################################## 
      SUBROUTINE ISBA_SWNET_MEB(PLAI,PF_TALLVEG,                                             &
           PALBVIS_SNOW,PALBNIR_SNOW,                                                        &
           PALBNIR_VEG, PALBVIS_VEG,PALBNIR_SOIL, PALBVIS_SOIL, PFALB,                       &
           PVEG,PFF,PPSN,PPSNA,PTAU_N,                                                       &
           PZENITH,PSCA_SW,PSW_RAD,PSWUP,PSWNET_N,PSWNET_V,PSWNET_G,PSWNET_NS,PALBEDO,PALBG, &
           PSWDOWN_GN                                                                        )

!
!!****  *ISBA_SWNET_MEB*
!!
!!    PURPOSE
!!    -------
!
!     Calculates the net shortwave radiation budget terms for fully
!     coupled snow, soil-understory vegetation and canopy vegetation.
!         
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    Noilhan and Planton (1989)
!!    Belair (1995)
!!    * to be done * (2011)
!!      
!!    AUTHOR
!!    ------
!!
!!	A. Boone           * Meteo-France *
!!      P. Samuelsson      * SMHI *
!!      S. Gollvik         * SMHI * 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    22/01/11 
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
USE MODD_MEB_PAR, ONLY : XTAU_V_CF
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

!
!*      0.1    Declaration of Arguments
!
REAL, DIMENSION(:),   INTENT(IN)  :: PLAI, PF_TALLVEG
!
REAL, DIMENSION(:),   INTENT(IN)  :: PZENITH, PSW_RAD, PSCA_SW,                             &
                                     PPSN, PPSNA, PTAU_N
!
REAL, DIMENSION(:),   INTENT(IN)  :: PVEG, PFF
!
REAL, DIMENSION(:),   INTENT(IN)  :: PFALB
!
REAL, DIMENSION(:),   INTENT(IN)  :: PALBNIR_SNOW, PALBVIS_SNOW
!
REAL, DIMENSION(:),   INTENT(IN)  :: PALBNIR_VEG, PALBVIS_VEG,PALBNIR_SOIL, PALBVIS_SOIL
!
REAL, DIMENSION(:),   INTENT(OUT) :: PSWUP, PSWNET_N, PSWNET_V, PSWNET_G, PSWNET_NS,        &
                                     PALBEDO, PALBG, PSWDOWN_GN
!                                    PALBEDO    = effective surface all-wavelength albedo (-)
!                                    PSWDOWN_GN = all wavelength SWdown transmitted through
!                                                 canopy to ground-snow surface (W m-2) 
!
!
!*      0.2    Declaration of local Arguments
!
INTEGER, PARAMETER             :: NSPB = 2   ! spectral bands considered...NOTE! coded herein for 
!                                            ! this number (and for arguments!)
!
INTEGER                        :: JZENITHA, JBAND, JJ ! loop control
!
INTEGER                        :: INI                 ! Number of grid points
!
REAL, DIMENSION(SIZE(PZENITH)) :: ZWORK
!
REAL, DIMENSION(SIZE(PZENITH)) :: ZSWUP, ZSWNET_NS,           &
                                  ZSWNET_N,                   &
                                  ZSWNET_G, ZSWNET_V

REAL, DIMENSION(SIZE(PZENITH)) :: ZSWUP_SUM, ZSWNET_N_SUM,    &
                                  ZSWNET_G_SUM, ZSWNET_V_SUM

REAL, DIMENSION(SIZE(PZENITH)) :: ZCHI_DIF, ZCHI_DIR, ZF_DIF, &
                                  ZSIGMA_F, ZSIGMA_FN,        &
                                  ZSIGMA_F_DIF, ZSIGMA_FN_DIF,&
                                  ZPSNN

REAL, DIMENSION(SIZE(PZENITH)) :: ZALPHA_VS_AVG, ZDELTA,      &
                                  ZALPHA_V_O_S, ZALPHA_VS,    &
                                  ZALPHA_VS_NIR, ZSWDOWN,     &
                                  ZALPHA_V, ZALPHA_G,         &
                                  ZALPHA_N, ZALPHA_S

REAL, DIMENSION(SIZE(PZENITH)) :: ZSWDOWN_A,  ZSWDOWN_B,      &
                                  ZSWDOWN_C,  ZSWDOWN_J,      &
                                  ZSWDOWN_E,  ZSWDOWN_H,      &
                                  ZSWDOWN_K

REAL, DIMENSION(SIZE(PZENITH)) :: ZSWDOWNN_A, ZSWDOWNN_B,     &
                                  ZSWDOWNN_C, ZSWDOWNN_J,     &
                                  ZSWDOWNN_E, ZSWDOWNN_H,     &
                                  ZSWDOWNN_K

REAL, DIMENSION(SIZE(PZENITH)) :: ZLAIN, ZSWNET,              &
                                  ZINV_COSZENITH, ZCOSZENITH
!
REAL, DIMENSION(NSPB)               :: ZF_BND
REAL, DIMENSION(SIZE(PZENITH),NSPB) :: ZALPHA_VS_BND, ZALPHA_VT_BND, ZALPHA_G_BND, ZALPHA_N_BND
!
!*      0.3    Declaration of local Parameters
!
INTEGER,                   PARAMETER :: IZENITHA     = 3        
REAL, DIMENSION(IZENITHA), PARAMETER :: ZZENITHAWGHT = (/ 0.308, 0.514, 0.178/) 
REAL, DIMENSION(IZENITHA), PARAMETER :: ZCOSZENITHA  = (/0.96592581, 0.70710677, 0.25881907 /)
!                                                      ! NOTE for zenith = 15, 45, 75 (deg)
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       0.     Initialization:
!               ---------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_SWNET_MEB',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
INI = SIZE(PZENITH)
!
! Assume that input vegetation albedo corresponds to tall vegetation and
! set short-vegetation albedo values equal to tall vegetation albedo 
!
ZALPHA_VT_BND(:,1) = PALBVIS_VEG(:)
ZALPHA_VT_BND(:,2) = PALBNIR_VEG(:)
ZALPHA_VS_BND(:,:) = ZALPHA_VT_BND(:,:)
ZALPHA_G_BND(:,1)  = PALBVIS_SOIL(:)
ZALPHA_G_BND(:,2)  = PALBNIR_SOIL(:)
ZALPHA_N_BND(:,1)  = PALBVIS_SNOW(:)
ZALPHA_N_BND(:,2)  = PALBNIR_SNOW(:)
!
ZF_BND(1)          = 0.48              ! VIS
ZF_BND(2)          = 0.52              ! NIR
!
!*       1.a    View factors: transmission: non-snow buried canopy
!               --------------------------------------------------
!
ZCHI_DIF(:)  = 0.0
DO JZENITHA=1,IZENITHA
   DO JJ=1,INI
      ZCHI_DIF(JJ) = ZCHI_DIF(JJ) +                                                        &
                     ZZENITHAWGHT(JZENITHA)*EXP(-XTAU_V_CF*PLAI(JJ)/ZCOSZENITHA(JZENITHA)) 
   ENDDO
ENDDO

ZCOSZENITH(:)     = MAX(0.0, COS(PZENITH(:)) )

ZINV_COSZENITH(:) = 1.0/MAX(0.01, ZCOSZENITH(:)) ! inverse of cos(zenith)

ZCHI_DIR(:)   = EXP(-XTAU_V_CF*PLAI(:)*ZINV_COSZENITH(:))

! - compute diffuse fraction

ZF_DIF(:)     = PSCA_SW(:)/MAX(0.01, PSW_RAD(:))

! - view factors/transmission factors:

ZSIGMA_F(:)    = 1.0 - ( (1.0 - ZF_DIF(:))*ZCHI_DIR(:) + ZF_DIF(:)*ZCHI_DIF(:) )
!
ZSIGMA_F_DIF(:)= 1.0 - ZCHI_DIF(:) ! save diffuse only value
!
!
!*       1.b    View factors: transmission: partially snow-buried canopy
!               --------------------------------------------------------
!
! First compute effective snow-free LAI:
!
ZLAIN(:) = PLAI(:)*(1.0 - PPSNA(:))
!
ZCHI_DIF(:)  = 0.0
DO JZENITHA=1,IZENITHA
   DO JJ=1,INI
      ZCHI_DIF(JJ)  = ZCHI_DIF(JJ) +                                                         &
                      ZZENITHAWGHT(JZENITHA)*EXP(-XTAU_V_CF*ZLAIN(JJ)/ZCOSZENITHA(JZENITHA)) 
   ENDDO
ENDDO

ZCHI_DIR(:)  = EXP(-XTAU_V_CF*ZLAIN(:)*ZINV_COSZENITH(:) )

! - view factors/transmission factors:

ZSIGMA_FN(:)      = 1.0 - ( (1.0 - ZF_DIF(:))*ZCHI_DIR(:) + ZF_DIF(:)*ZCHI_DIF(:) )
!
ZSIGMA_FN_DIF(:)  = 1.0 - ZCHI_DIF(:) ! save diffuse only value
!
!*       2.     Albedo
!               ------
! 
! Prelim albedo: Compute short veg NIR:

! - all wavelength short veg:

ZALPHA_VS_AVG(:)  = SUM(ZALPHA_VS_BND,2)/SIZE(ZALPHA_VS_BND,2)

! - clear sky short veg albedo

ZDELTA(:)         = 1.0
WHERE(ZCOSZENITH(:)>0.5)ZDELTA(:) = 0.0
ZALPHA_V_O_S(:)   = ZALPHA_VS_AVG(:)*( (ZDELTA(:)/(0.5+ZCOSZENITH(:))) + &
                      (1.-ZDELTA(:))*(0.5 + (1/(1.+2*ZCOSZENITH(:)))) )

!  - assume cloudy = all wavelength:

ZALPHA_VS(:)      = (1.-ZF_DIF(:))*ZALPHA_V_O_S(:) + ZF_DIF(:)*ZALPHA_VS_AVG(:)

ZALPHA_VS_NIR(:)  = (ZALPHA_VS(:)*PSW_RAD(:) -                                 &
                       ZALPHA_VS_AVG(:)*PSW_RAD(:)*ZF_BND(1) )/                &
                      (PSW_RAD(:)*ZF_BND(2) + 1.E-4)

!
!*       3.     Spectral Band Integration
!               -------------------------
!
! Initialize sums:
!
ZSWUP_SUM(:)     = 0.0
ZSWNET_G_SUM(:)  = 0.0
ZSWNET_V_SUM(:)  = 0.0
ZSWNET_N_SUM(:)  = 0.0
!
PALBG(:)         = 0.0
PSWDOWN_GN(:)    = 0.0
!
!
! Effective snow fraction for reflection from snow-free ground back down
! in the presence of snow-buried vegetation:
!
ZPSNN(:) = PPSN(:)*(1. - PPSNA(:))
!
!
SPECTRAL_BAND_INT: DO JBAND=1,SIZE(ZALPHA_VS_BND,2)
         
   DO JJ=1,INI

! Get shortwave radiation for this part of spectrum:

      ZSWDOWN(JJ)        = PSW_RAD(JJ)*ZF_BND(JBAND)

! Set up albedos for each band

! - veg is a mix of tall and short note, short veg for NIR get's special treatment

      IF(JBAND==2)THEN ! NIR
         ZALPHA_VS(JJ)   = ZALPHA_VS_NIR(JJ)
      ELSE
         ZALPHA_VS(JJ)   = ZALPHA_VS_BND(JJ,JBAND)
      ENDIF

      ZALPHA_V(JJ)       =     PF_TALLVEG(JJ)* ZALPHA_VT_BND(JJ,JBAND) +                   &
                           (1.-PF_TALLVEG(JJ))*ZALPHA_VS(JJ)

! - soilJJ

      ZALPHA_G(JJ)       = ZALPHA_G_BND(JJ,JBAND)

! - snow

      ZALPHA_N(JJ)       = ZALPHA_N_BND(JJ,JBAND)

! - net (below canopy) surface
!   incorporate any understory vegetation...we assume it has the same albedo as the overlying vegetation.
!   Also, consider flooding....for now, just use constant water albedo

      ZALPHA_S(JJ)     = (1.-PFF(JJ))*((1.-PVEG(JJ))*ZALPHA_G(JJ) + PVEG(JJ)*ZALPHA_V(JJ)) + PFF(JJ)*PFALB(JJ)
!
!*       3.a    Radiation passing through part of canopy with no snow below
!               ----------------------------------------------------------- 

      ZSWDOWN_A(JJ)  = ZSWDOWN(JJ)  *(1.-PPSN(JJ))
      ZSWDOWN_B(JJ)  = ZSWDOWN_A(JJ)*(1.-ZSIGMA_F(JJ))
      ZSWDOWN_C(JJ)  = ZSWDOWN_B(JJ)*                ZALPHA_S(JJ)*(1.-PPSN(JJ))
      ZWORK(JJ)      = ZSWDOWN_C(JJ)* ZSIGMA_F(JJ) * ZALPHA_V(JJ)
      ZSWDOWN_E(JJ)  = ZWORK(JJ)    * (1.-ZPSNN(JJ))
      ZSWDOWN_H(JJ)  = ZWORK(JJ)    *     ZPSNN(JJ)
      ZSWDOWN_J(JJ)  = ZSWDOWN_C(JJ)*(1.-ZSIGMA_F_DIF(JJ))
      ZSWDOWN_K(JJ)  = ZSWDOWN_A(JJ)*    ZSIGMA_F(JJ) *ZALPHA_V(JJ)
!
!
!*       3.b    Radiation passing through part of canopy with snow below
!               --------------------------------------------------------
!
      ZSWDOWNN_A(JJ) = ZSWDOWN(JJ)   * PPSN(JJ)
      ZSWDOWNN_B(JJ) = ZSWDOWNN_A(JJ)*(1.-ZSIGMA_FN(JJ))
      ZSWDOWNN_C(JJ) = ZSWDOWNN_B(JJ)*                   ZALPHA_N(JJ) * PPSN(JJ)
      ZWORK(JJ)      = ZSWDOWNN_C(JJ)* ZSIGMA_FN(JJ)   * ZALPHA_V(JJ)
      ZSWDOWNN_E(JJ) = ZWORK(JJ)     * (1.-PPSN(JJ))
      ZSWDOWNN_H(JJ) = ZWORK(JJ)     *     PPSN(JJ)
      ZSWDOWNN_J(JJ) = ZSWDOWNN_C(JJ)*(1.-ZSIGMA_FN_DIF(JJ))
      ZSWDOWNN_K(JJ) = ZSWDOWNN_A(JJ)*    ZSIGMA_FN(JJ) *ZALPHA_V(JJ)
!
!*       3.c    Diagnostics
!               -----------
!
! Downwelling SW radiation arriving at ground/snow surfaceJJ
!
      PSWDOWN_GN(JJ) = PSWDOWN_GN(JJ) + ZSWDOWNN_B(JJ) + ZSWDOWN_B(JJ)  

! Compute all-wavelength effective ground albedoJJ

      PALBG(JJ)      = PALBG(JJ) + ZALPHA_S(JJ)*ZF_BND(JBAND)

   ENDDO

!
!*       3.d    Net budgets
!               -----------

   ZSWUP(:)      =  ZSWDOWN_J(:)  + ZSWDOWN_K(:)  +   &
                    ZSWDOWNN_J(:) + ZSWDOWNN_K(:) 

   ZSWNET_N(:)   =  ZSWDOWNN_B(:) - ZSWDOWNN_C(:) +   & 
                    ZSWDOWNN_H(:) + ZSWDOWN_H(:) 
      
   ZSWNET_G(:)   =  ZSWDOWN_B(:)  - ZSWDOWN_C(:)  +   &
                    ZSWDOWNN_E(:) + ZSWDOWN_E(:) 

! Here, use more simple computation for net vegetation (to get same result) and to reduce
! any small numerical errors:
!
   ZSWNET(:)     = ZSWDOWN(:) - ZSWUP(:)
   ZSWNET_V(:)   = ZSWNET(:)  - ZSWNET_G(:) - ZSWNET_N(:)
!
!
!*       3.e    Spectral integration (sums)
!               ---------------------------
!
   ZSWUP_SUM(:)     = ZSWUP_SUM(:)     + ZSWUP(:)
   ZSWNET_N_SUM(:)  = ZSWNET_N_SUM(:)  + ZSWNET_N(:)
   ZSWNET_G_SUM(:)  = ZSWNET_G_SUM(:)  + ZSWNET_G(:) 
   ZSWNET_V_SUM(:)  = ZSWNET_V_SUM(:)  + ZSWNET_V(:) 

ENDDO SPECTRAL_BAND_INT

! Save fluxes integrated over spectral bands:

PSWUP(:)     = ZSWUP_SUM(:)
PSWNET_G(:)  = ZSWNET_G_SUM(:)
PSWNET_V(:)  = ZSWNET_V_SUM(:)
PSWNET_N(:)  = ZSWNET_N_SUM(:)

! Quantity of net shortwave radiation *absorbed* in surface snow layer
! (rest goes into snowpack sub-sfc layers = ZSWNET_NS*PTAU_N in snow scheme)
! This quantity is used in snow surface energy budget computation.

PSWNET_NS(:)  =  PSWNET_N(:)*(1.0 - PTAU_N(:))
   
! All wavelength effective albedo (-):

PALBEDO(:)   = PSWUP(:)/MAX(1.E-5, PSW_RAD(:))

IF (LHOOK) CALL DR_HOOK('ISBA_SWNET_MEB',1,ZHOOK_HANDLE)

end SUBROUTINE ISBA_SWNET_MEB
