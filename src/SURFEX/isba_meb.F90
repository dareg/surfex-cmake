!     #########
      SUBROUTINE ISBA_MEB(OMEB, OFORC_MEASURE, TPTIME,                         &  
        PPS, PZENITH, PSCA_SW, PSW_RAD,                                        &
        PALBNIR_TVEG, PALBVIS_TVEG,PALBNIR_TSOIL, PALBVIS_TSOIL, PFALB,        &
        PSNOWALB, PSNOWALBVIS, PSNOWALBNIR, PSNOWALBFIR,                       &
        PVEG, PFF, PPSN, PPALPHAN, PZF_TALLVEG, PLAIV,                         &
        PSNOWRHO, PSNOWSWE, PSNOWHEAT, PSNOWTEMP, PSNOWDZ, PEMISNOW,           &
        PSWUP, PSWNET_N, PSWNET_V, PSWNET_G, PSWNET_NS, PALBT, PSWDOWN_GN      )
!     ##########################################################################
!
!!****  *ISBA_MEB*  
!!
!!    PURPOSE
!!    -------
!       Monitor for the calculation of the surface fluxes and of the
!     prognostic variables of the surface over natural areas 
!     with an explicit vegetation layer
!
!     NOTE...currently MEB can be coupled with 
!     HISBA='DIF' or '3-L' soil options
!     HSNOW='3-L' snow scheme
!     Soon, HSNOW=CRO and HPHOTO/=NON (i.e. Ags will be added)
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
!!    Noilhan and Planton (1989)
!!      
!!    AUTHOR
!!    ------
!!	A. Boone           * Meteo-France *
!!      P. Samuelsson      * SMHI *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2014
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_CSTS,       ONLY : XLVTT, XLSTT, XCPD, XRHOLW 
!
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE MODE_THERMOS
USE MODE_MEB,    ONLY     : SNOW_INTERCEPT_EFF
!
USE MODI_SOILSTRESS
USE MODI_WET_LEAVES_FRAC
USE MODI_VEG
USE MODI_SNOW_LEAVES_FRAC_MEB
USE MODI_PREPS_FOR_MEB_EBUD_RAD
USE MODI_ISBA_SWNET_MEB
USE MODI_ISBA_LWNET_MEB
USE MODI_DRAG_MEB
USE MODI_E_BUDGET_MEB
USE MODI_ISBA_FLUXES_MEB
USE MODI_SNOW_LOAD_MEB
USE MODI_HYDRO_VEG
USE MODI_SNOW3L_ISBA

!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!              -------------------------
!
!
!* general variables
!  -----------------
!
TYPE(DATE_TIME), INTENT(IN)         :: TPTIME        ! current date and time
!
LOGICAL, INTENT(IN)                 :: OMEB          ! True = patch with multi-energy balance 
!                                                    ! False = patch with classical ISBA 
LOGICAL, INTENT(IN)                 :: OFORC_MEASURE ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PPS           ! Pressure
REAL, DIMENSION(:),   INTENT(IN)    :: PZENITH       ! solar zenith angle
REAL, DIMENSION(:),   INTENT(IN)    :: PSW_RAD       ! solar   incoming radiation
REAL, DIMENSION(:),   INTENT(IN)    :: PSCA_SW       ! solar diffuse incoming radiation
!
REAL, DIMENSION(:),   INTENT(IN)    :: PZF_TALLVEG   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLAIV         ! explicit canopy overstory LAI..."PLAI" is the 
!                                                    ! composite surface  LAI 
REAL, DIMENSION(:),   INTENT(IN)    :: PVEG          ! fraction of vegetation of the
!                                                    ! mesh covered by natural or
!                                                    ! agricultural areas
!                                                    ! 1-PVEG --> bare soil
REAL, DIMENSION(:),   INTENT(IN)    :: PFF           ! Floodplain fraction at the surface
REAL, DIMENSION(:),   INTENT(IN)    :: PPSN          ! fraction of the grid covered
!                                                    ! by snow
REAL, DIMENSION(:),   INTENT(IN)    :: PPALPHAN      ! snow/canopy transition coefficient
REAL, DIMENSION(:),   INTENT(IN)    :: PFALB         ! Floodplain albedo
REAL, DIMENSION(:),   INTENT(IN)    :: PALBNIR_TVEG  ! albedo of vegetation in NIR 
!                                                    ! (needed for LM_TR or MEB)
REAL, DIMENSION(:),   INTENT(IN)    :: PALBVIS_TVEG  ! albedo of vegetation in VIS 
!                                                    ! (needed for LM_TR or MEB)
REAL, DIMENSION(:),   INTENT(IN)    :: PALBNIR_TSOIL ! albedo of bare soil in NIR 
!                                                    ! (needed for LM_TR or MEB)
REAL, DIMENSION(:),   INTENT(IN)    :: PALBVIS_TSOIL ! albedo of bare soil in VIS 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALB      ! Snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALBVIS   ! Snow VIS albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALBNIR   ! Snow NIR albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALBFIR   ! Snow FIR albedo
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWSWE      ! Snow model layer liquid water equivalent or 
!                                                    ! SWE (kg m-2)  
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWHEAT     ! Snow layer heat content (J/m3) 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO      ! Snow layer average density (kg/m3)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWTEMP     ! Snow layer average temperature (K)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ       ! Snow layer thickness (m)
!
REAL, DIMENSION(:),   INTENT(OUT)   :: PEMISNOW      ! Snow surface emissivity (-)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWNET_N      ! net snow shortwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWNET_NS     ! net snow shortwave radiation for 
!                                                    ! the *surface* snow layer 
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWNET_V      ! net vegetation canopy shortwave radiation 
!                                                    ! [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWNET_G      ! net surface (ground) shortwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWUP         ! net upwelling shortwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PALBT         ! total surface albedo
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWDOWN_GN    ! total shortwave radiation transmitted through 
                                                     ! the vegetation canopy
!
!
!*      0.2    declarations of local variables
!
!
REAL, PARAMETER                                    :: ZTSTEP_EB     = 300. ! s Minimum time tstep required 
!                                                                          !   to time-split MEB energy budget
INTEGER                                            :: KTSPLIT_EB
!
REAL, DIMENSION(SIZE(PSNOWSWE,1),SIZE(PSNOWSWE,2)) :: ZSNOWCOND            ! snow thermal conductivity  [W/(m K)] 
REAL, DIMENSION(SIZE(PSNOWSWE,1),SIZE(PSNOWSWE,2)) :: ZSNOWHCAP            ! snow heat capacity [J/(m3 K)]
REAL, DIMENSION(SIZE(PSNOWSWE,1))                  :: ZSNOWALBVIS          ! Snow VIS albedo
REAL, DIMENSION(SIZE(PSNOWSWE,1))                  :: ZSNOWALBNIR          ! Snow NIR albedo
REAL, DIMENSION(SIZE(PPS))                         :: ZCHIP                ! 
REAL, DIMENSION(SIZE(PPS))                         :: ZTAU_N               ! 
REAL, DIMENSION(SIZE(PPS))                         :: ZSIGMA_F             ! 
REAL, DIMENSION(SIZE(PPS))                         :: ZALBG                ! effective ground albedo
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      1.0    Preliminaries
!              -------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB',0,ZHOOK_HANDLE)
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      1.0    Preliminaries for energy and radiation budget
!              ---------------------------------------------
!
CALL PREPS_FOR_MEB_EBUD_RAD(PPS,                                     &
        PLAIV,PSNOWRHO,PSNOWSWE,PSNOWHEAT,                           &
        PSNOWTEMP,PSNOWDZ,ZSNOWCOND,ZSNOWHCAP,PEMISNOW,ZTAU_N,       &
        ZSIGMA_F,ZCHIP)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      2.0    Shortwave radiative transfer
!              ----------------------------  
!
!
! NOTE, currently MEB only uses 2 of 3 potential snow albedo spectral bands
!
CALL ISBA_SWNET_MEB(PLAIV,PZF_TALLVEG,                                       &
        PSNOWALBVIS,PSNOWALBNIR,                                             &
        PALBNIR_TVEG, PALBVIS_TVEG,PALBNIR_TSOIL, PALBVIS_TSOIL, PFALB,      &
        PVEG,PFF,PPSN,PPALPHAN,ZTAU_N,PZENITH,PSCA_SW,PSW_RAD,               &
        PSWUP,PSWNET_N,PSWNET_V,PSWNET_G,PSWNET_NS,                          &
        PALBT,ZALBG,PSWDOWN_GN                                               )

!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ISBA_MEB
