!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ####################
      MODULE MODD_SURF_ATM
!     ####################
!
!!****  *MODD_SURF_ATM - declaration of surface ATM
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      P. Le Moigne *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       10/2007
!
!*       0.   DECLARATIONS
!             ------------
!
!!      B.Decharme     04/2009 Add flag used to Read/Write precipitation forcing from/into the restart file for ARPEGE/ALADIN run
!!      B.Decharme     08/2009 Add flag used to know if you use SURFEX in the Earth System Model
!!      B.Decharme     09/2012 New wind implicitation key option
!!      B.Decharme     04/2013 Flag used to Read/Write some field from/into the restart file for coupling with ARPEGE/ALADIN
!!                             Delete LRW_PRECIP, LSAVE_PRECIP
!!                             Vertical shift for LW and Precip
!!      R. Séférian    03/2014 Adding key for decouple CO2 for photosynthesis (XCO2UNCPL) 
!!      Y. Seity       08/2023 Add XACMAX
!!      J. Masek       08/2023 New variables XZ0_OFFSET, XWNEW, XWCRN, XFACZ0,
!!                             XTAU_ICE, LZ0_AVG_EXACT, LZ0_EFF; removed LALDZ0H
!
IMPLICIT NONE
!
!-----------------------------------------------------------------------------------------------------
REAL       :: XZ0_OFFSET ! flag for approximate/exact roughness averaging (0/1)
REAL       :: XCISMIN  ! minimum wind shear
REAL       :: XVMODMIN ! minimum wind speed
REAL       :: XWNEW    ! factor in snow albedo calculation for melting/no-melting case
REAL       :: XWCRN    ! critical value of snow reservoir in calculation of snow fraction
REAL       :: XFACZ0   ! multiplicative factor for orographic roughness length
REAL       :: XTAU_ICE ! characteristic time for ice in force-restore (s)
LOGICAL    :: LALDTHRES! activate aladin threshold for wind

REAL       :: XVMODFAC ! factor for threshold of wind speed in drag calculations
!
LOGICAL    :: LDRAG_COEF_ARP ! activate aladin formulation for Cd and Ch
LOGICAL    :: LZ0_AVG_EXACT  ! activate exact formula for z0 averaging
                             ! (via unapproximated neutral drag coefficient)
LOGICAL    :: LZ0_EFF        ! include orog. component in mechanical roughness
!
LOGICAL    :: LNOSOF   ! No parameterization of Subgrid Orography effects on atmospheric Forcing  
LOGICAL    :: LSLOPE   ! Parameterization of homogeneous slope effect on atmospheric Forcing
LOGICAL    :: LVERTSHIFT ! vertical shift from atmospheric orography to surface orography
LOGICAL    :: LVSHIFT_LW ! vertical shift for LW
LOGICAL    :: LVSHIFT_PRCP ! vertical shift for Precip
!
LOGICAL    :: LVZIUSTAR0_ARP   ! activate aladin formulation for zoh over sea
LOGICAL    :: LRRGUST_ARP      ! activate aladin formulation for CD CH, CDN, correction due to moist gustiness
LOGICAL    :: LCPL_ARP         ! activate aladin formulation for Cp and L
LOGICAL    :: LQVNPLUS         ! An option for the resolution of the surface temperature equation
LOGICAL    :: LZ0SNOWH_ARP     ! Snow effect on roughness length via snow height
!
LOGICAL    :: LCPL_GCM         ! Flag used to Read/Write some field from/into the restart file for coupling with ARPEGE/ALADIN
!
REAL       :: XEDB
REAL       :: XEDC
REAL       :: XEDD
REAL       :: XEDK
REAL       :: XUSURIC
REAL       :: XUSURID
REAL       :: XUSURICL
REAL       :: XVCHRNK
REAL       :: XVZ0CM
REAL       :: XRIMAX
REAL       :: XRISHIFT   ! shift of Ri numbers in drag  calculations
REAL       :: XCD_COEFF1
REAL       :: XCD_COEFF2
REAL       :: XCH_COEFF1
REAL       :: XACMAX
REAL       :: XDELTA_MAX ! Maximum fraction of the foliage covered by intercepted water for high vegetation
!
REAL       :: XWINDMIN ! minimum wind speed (canopy)
!
REAL       :: XRZHZ0M
REAL       :: XVZIUSTAR0
REAL       :: XRRSCALE
REAL       :: XRRGAMMA
REAL       :: XUTILGUST
REAL       :: XRZ0_TO_HEIGHT  ! roughness to obstacle height ratio
!
REAL       :: XCO2UNCPL ! uncoupled CO2 values (ppmv)
!
!-----------------------------------------------------------------------------------------------------
!
END MODULE MODD_SURF_ATM
