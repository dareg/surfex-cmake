!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!!
!!    #####################
      MODULE MODD_TREEDRAG
!!    #####################
!!
!!*** *MODD_TREEDRAG*
!!
!!    PURPOSE
!!    -------
!       Declaration to take into account tree drag in the atmospheric model
!              instead of SURFEX. The Z0 forest is therefore reduced to
!              the Z0 grass
!!
!!**  AUTHOR
!!    ------
!!    C.Lac                   *CNRM*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 30/06/11
!!    P. Samuelsson, SMHI  02/2020: Added XALLEN_TERM and XGRASS_H_DNM
!!    S. Viana, AEMET      06/2020: Added LFAKETREE,XHFAKETREE,XFFAKETREE
!!    P. Samuelsson, SMHI  05/2022: Extended LFAKETREE to a vector
!!
!!    08/2023 J. Masek added array XSCALE_H_TREE_ECOFG
!-----------------------------------------------------------------------------


USE MODD_DATA_COVER_PAR, ONLY: NVEGTYPE_MAX

!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
LOGICAL    ::     LTREEDRAG=.FALSE.    ! flag used to  take into account tree drag in 
!                                      ! the atmospheric model instead of SURFEX.
REAL       ::     XSCALE_H_TREE_ECOSG=1.0 ! Scale the tree height with this factor
REAL       ::     XALLEN_TERM=3.5      ! The term is in the expression for height of crops
REAL       ::     XGRASS_H_DNM=6.0     ! The denominater value is in the expression for height of grass
REAL       ::     XHFAKETREE=10.       ! Height of fake trees 
REAL       ::     XFFAKETREE=0.1       ! Fraction of fake trees 
REAL       ::     XZ0_MIN_LIMIT=0.     ! Limit Z0 (m) in Z0V_FROM_LAI to this minimum value.
REAL       ::     XZ0_MAX_LIMIT=999.   ! Limit Z0 (m) in Z0V_FROM_LAI to this maximum value.
REAL       ::     XFORLAT1=1.0         ! 
REAL       ::     XFORLAT2=0.0         ! 
REAL       ::     XFORFRAC1=1.0        ! 
REAL       ::     XFORFRAC2=0.0        ! 
REAL       :: XSCALE_H_TREE_ECOFG(1:NVEGTYPE_MAX)=1.  ! scaling factor for the tree height
!
! Flag to activate fake trees option for low vegetation ECOSG Vegtypes.
! The logical vector has 7 positions where the positions represent:
! 1 NVT_BOGR  boreal grass
! 2 NVT_GRAS  grassland
! 3 NVT_TROG  tropical grassland
! 4 NVT_C3W   C3W cultures types
! 5 NVT_C3S   C3S cultures types
! 6 NVT_C4    C4 cultures types
! 7 NVT_FLGR  flooded grassland
LOGICAL, DIMENSION(7) :: LFAKETREE=.FALSE.    ! Flag to activate fake trees option in low vegetation areas.
!
END MODULE MODD_TREEDRAG
