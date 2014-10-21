!     #########
      SUBROUTINE ISBA_MEB(OMEB, OFORC_MEASURE, TPTIME,               &                
        PPS,                                                         &
        PLAIV,PSNOWRHO,PSNOWSWE,PSNOWHEAT,                           &
        PSNOWTEMP,PSNOWDZ,PEMISNOW                                   )
!     ##########################################################################
!
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
USE MODE_MEB,  ONLY       : SNOW_INTERCEPT_EFF

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
TYPE(DATE_TIME), INTENT(IN)         :: TPTIME     ! current date and time
!
LOGICAL, INTENT(IN)                 :: OMEB       ! True = patch with multi-energy balance 
!                                                 ! False = patch with classical ISBA 
LOGICAL, INTENT(IN)                 :: OFORC_MEASURE
!
REAL, DIMENSION(:), INTENT(IN)      :: PPS        ! Pressure
REAL, DIMENSION(:), INTENT(IN)      :: PLAIV      ! explicit canopy overstory LAI..."PLAI" is the 
!                                                 ! composite surface  LAI 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWSWE   ! Snow model layer liquid water equivalent or SWE (kg m-2)  
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWHEAT  ! Snow layer heat content (J/m3) 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO   ! Snow layer average density (kg/m3)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWTEMP  ! Snow layer average temperature (K)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ    ! Snow layer thickness (m)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PEMISNOW   ! Snow surface emissivity (-)
!
!*      0.2    declarations of local variables
!
!
REAL, PARAMETER                                    :: ZTSTEP_EB     = 300. ! s Minimum time tstep required 
!                                                                          !   to time-split MEB energy budget
INTEGER                                            :: KTSPLIT_EB
!
REAL, DIMENSION(SIZE(PSNOWSWE,1),SIZE(PSNOWSWE,2)) :: ZSNOWCOND  ! snow thermal conductivity  [W/(m K)] 
REAL, DIMENSION(SIZE(PSNOWSWE,1),SIZE(PSNOWSWE,2)) :: ZSNOWHCAP  ! snow heat capacity [J/(m3 K)]
REAL, DIMENSION(SIZE(PPS))                         :: ZCHIP                ! 
REAL, DIMENSION(SIZE(PPS))                         :: ZTAU_N               ! 
REAL, DIMENSION(SIZE(PPS))                         :: ZSIGMA_F             ! 
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
IF (LHOOK) CALL DR_HOOK('ISBA_MEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ISBA_MEB
