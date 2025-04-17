!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SFX_LIC for details. version 1.
!!
!!    #####################
      MODULE MODN_SURF_ATM
!!    #####################
!!
!!*** *MODN_DUST*
!!
!!    PURPOSE
!!    -------
!       Namelist for wind threshold
!!
!!**  AUTHOR
!!    ------
!!    P. Le Moigne      *CNRM*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 10/2007
!!    J. Masek 08/2023 new variables LZ0_AVG_EXACT, LZ0_EFF, XWNEW, XWCRN,
!!                     XFACZ0, XTAU_ICE; removed LALDZ0H, LARP_PN
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_SURF_ATM, ONLY :   XZ0_OFFSET,                   &
                            XCISMIN, XVMODMIN, XWNEW,     &
                            XWCRN, XFACZ0, XTAU_ICE,      &
                            LALDTHRES, LDRAG_COEF_ARP,    &
                            LZ0_AVG_EXACT, LZ0_EFF,       &
                            LNOSOF, LSLOPE, LCPL_GCM,     &
                            XEDB, XEDC, XEDD, XEDK,       &
                            XUSURIC, XUSURID, XUSURICL,   &
                            XVCHRNK, XVZ0CM, XDELTA_MAX,  &
                            XRIMAX, XRISHIFT, XACMAX,     &
                            LVERTSHIFT,                   &
                            LVZIUSTAR0_ARP, LRRGUST_ARP,  &
                            XVZIUSTAR0,XRZHZ0M,           &
                            XRRSCALE, XRRGAMMA,           &
                            XUTILGUST, LCPL_ARP, LQVNPLUS,&
                            LVSHIFT_LW, LVSHIFT_PRCP,     &
                            XCO2UNCPL, XWINDMIN,          &
                            LZ0SNOWH_ARP, XRZ0_TO_HEIGHT, &
                            XCD_COEFF1, XCD_COEFF2,       &
                            XCH_COEFF1, XVMODFAC
!
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
SAVE
NAMELIST /NAM_SURF_ATM/  LZ0_AVG_EXACT, LZ0_EFF,        &
                         XWNEW, XWCRN, XFACZ0, XTAU_ICE,&
                         XCISMIN, XVMODMIN, LALDTHRES,  &
                         LDRAG_COEF_ARP,                &
                         LNOSOF, LSLOPE, LCPL_GCM,      &
                         XEDB, XEDC, XEDD, XEDK,        &
                         XUSURIC, XUSURID, XUSURICL,    &
                         XVCHRNK, XVZ0CM, XDELTA_MAX,   &
                         XRIMAX, XRISHIFT, XACMAX,      &
                         LVERTSHIFT,                    &
                         LVZIUSTAR0_ARP, LRRGUST_ARP,   &
                         XVZIUSTAR0,XRZHZ0M,            &
                         XRRSCALE, XRRGAMMA,            &
                         XUTILGUST, LCPL_ARP, LQVNPLUS, &
                         LVSHIFT_LW, LVSHIFT_PRCP,      &
                         XCO2UNCPL, XWINDMIN,           &
                         LZ0SNOWH_ARP, XRZ0_TO_HEIGHT,  &
                         XCD_COEFF1, XCD_COEFF2,        &
                         XCH_COEFF1, XVMODFAC
!
END MODULE MODN_SURF_ATM
