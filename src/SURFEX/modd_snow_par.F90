!     ######################
      MODULE MODD_SNOW_PAR
!     ######################
!
!!****  *MODD_SNOW_PAR* - declaration of parameters related
!!                          to the snow parameterization
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     parameters related to the surface parameterization of snow.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004                    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!--------------------------------------------------------------------------------
! Snow on the ground: Given in ini_surf_csts and/or in NAM_SURF_CSTS
!--------------------------------------------------------------------------------
!
! Snow emissivity:
!
REAL, SAVE       :: XEMISSN
!
! Minimum and maximum values of the albedo of snow:
!
REAL, SAVE       :: XANSMIN
REAL, SAVE       :: XANSMAX 
!
! Minimum and maximum values of the albedo of permanet snow/ice:
!
REAL, SAVE       :: XAGLAMIN
REAL, SAVE       :: XAGLAMAX
! 
! Height (m) of aged snow in glacier case (allows Pn=1)
!
REAL, SAVE       :: XHGLA
! 
! Coefficient for calculation of snow fraction over vegetation
!
REAL, SAVE       :: XWSNV
!
! Roughness length of pure snow surface (m)
!
REAL, SAVE       :: XZ0SN  
!
! Roughness length for heat of pure snow surface (m)
!
REAL, SAVE       :: XZ0HSN
!
!--------------------------------------------------------------------------------
! Snow on the ground: PARAMETER
!--------------------------------------------------------------------------------
!
! Critical value of the equivalent water content
! of the snow reservoir for snow fractional coverage and albedo computations
!
REAL, PARAMETER       :: XWCRN = 10.0   ! (kg m-2) Veg (default value)
REAL, PARAMETER       :: XWCRN_EXPL =  1.0   ! (kg m-2) Veg explicit
REAL, PARAMETER       :: XWCRN_ROOF =  1.0   ! (kg m-2)  Roofs 
REAL, PARAMETER       :: XWCRN_ROAD =  1.0   ! (kg m-2)  Roads
REAL, PARAMETER       :: XWCRN_VEG  =  1.0   ! (kg m-2)  Urban veg
!
!
! Critical value of snow emissivity
!
REAL, PARAMETER       :: XEMCRIN = 0.98
!
! Minimum and maximum values of the albedo of snow:
!
REAL, PARAMETER       :: XANSMIN_ROOF = 0.30 ! (-)   Roofs
REAL, PARAMETER       :: XANSMIN_ROAD = 0.15 ! (-)   Roads
!
REAL, PARAMETER       :: XANSMAX_ROOF = 0.85 ! (-)   Roofs
REAL, PARAMETER       :: XANSMAX_ROAD = 0.85 ! (-)   Roads
!
! Snow aging coefficients (albedo and Force-Restore density):
!
REAL, PARAMETER       :: XANS_TODRY    = 0.008     ! (-) Veg (default value)
REAL, PARAMETER       :: XANS_TODRY_ROOF = 0.008   ! (-)  Roofs
REAL, PARAMETER       :: XANS_TODRY_ROAD = 0.008   ! (-)  Roads
!
REAL, PARAMETER       :: XANS_T        = 0.240     ! (-) Veg (default value)
REAL, PARAMETER       :: XANS_T_ROOF     = 0.174   ! (-)  Roofs
REAL, PARAMETER       :: XANS_T_ROAD     = 0.174   ! (-)  Roads (alley simul)
!
! Minimum and maximum values of the density of snow 
! for Force-Restore snow option
!
REAL, PARAMETER       :: XRHOSMIN = 100.       ! (kg m-3)   Veg (Default value)
REAL, PARAMETER       :: XRHOSMIN_ROOF = 100.  ! (kg m-3)   Roofs
REAL, PARAMETER       :: XRHOSMIN_ROAD = 100.  ! (kg m-3)   Roads
!
REAL, PARAMETER       :: XRHOSMAX = 300.       ! (kg m-3)   Veg (Default value)
REAL, PARAMETER       :: XRHOSMAX_ROOF = 300.  ! (kg m-3)   Roofs
REAL, PARAMETER       :: XRHOSMAX_ROAD = 350.  ! (kg m-3)   Roads
!
! Minimum and maximum values of the density of snow 
! for ISBA-ES snow option
!
REAL, PARAMETER       :: XRHOSMIN_ES =  50.  ! (kg m-3)
REAL, PARAMETER       :: XRHOSMAX_ES = 750.  ! (kg m-3)
!
! ISBA-ES Critical snow depth at which snow grid thicknesses constant
!
REAL, PARAMETER                      :: XSNOWCRITD = 0.03  ! (m)
!                                       
! ISBA-ES Minimum total snow depth for model 
!
REAL, PARAMETER                      :: XSNOWDMIN = 0.000001  ! (m)
!                                       
! Maximum Richardson number limit for very stable conditions using the ISBA-ES 'RIL' option
!
REAL, PARAMETER                      :: X_RI_MAX = 0.20
!                                       
! ISBA-ES Maximum snow liquid water holding capacity (fraction by mass) parameters:
!
REAL, PARAMETER                      :: XWSNOWHOLDMAX2   = 0.10  ! (-) 
REAL, PARAMETER                      :: XWSNOWHOLDMAX1   = 0.03  ! (-)
REAL, PARAMETER                      :: XSNOWRHOHOLD     = 200.0 ! (kg/m3)
!
END MODULE MODD_SNOW_PAR












