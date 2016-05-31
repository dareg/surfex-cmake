!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ############################
      SUBROUTINE INIT_TYPES_PARAM
!     ############################
!
!!**** *INIT_TYPES_PARAM* initializes cover-field correspondance arrays
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    S.Faroux        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    23/03/11
!!
!!    R. Alkama    05/2012 : read 19 vegtypes rather than 12
!     10/2014 : add status='old' for ecoclimap.bin files E. Martin
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------

USE MODD_TYPE_DATE_SURF
!
USE MODD_DATA_COVER,     ONLY : XDATA_TOWN, XDATA_NATURE, XDATA_SEA, XDATA_WATER,   &
                                XDATA_VEGTYPE, XDATA_GARDEN, XDATA_Z0_TOWN, &
                                XDATA_BLD, XDATA_BLD_HEIGHT, XDATA_WALL_O_HOR
!
USE MODD_DATA_COVER_PAR, ONLY : NTYPE, NUT_DENS, NUT_SUB1, NUT_SUB2,   &
                                NUT_ZIC, NUT_ROAD, NUT_PORT, NUT_AIR, NUT_MINE, &
                                NUT_PARK, NUT_SPOR, NVT_NO, NVT_TEBD, NVT_PARK, &
                                NVT_BOBD, NVT_GRAS, IDX_TWN_ECOSG
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER               :: JTYPE
!
!*    0.3    Declaration of namelists
!            ------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_TYPES_PARAM',0,ZHOOK_HANDLE)
!
! NATURE
!
DO JTYPE = 1,NTYPE(1)
  XDATA_NATURE(JTYPE) = 1.
  XDATA_VEGTYPE(JTYPE,JTYPE) = 1.
ENDDO
!
! SEA
!
XDATA_SEA(NTYPE(1)+1:NTYPE(1)+NTYPE(2)) = 1.
!
! WATER
!
XDATA_WATER(SUM(NTYPE(1:2))+1:SUM(NTYPE(1:3))) = 1.
!
! TOWN
!
XDATA_TOWN(SUM(NTYPE(1:3))+1:SUM(NTYPE(1:4))) = 1.
!
! index of town in ECOSG
IDX_TWN_ECOSG = SUM(NTYPE(1:3)) + 2
!
XDATA_GARDEN    (NUT_DENS) = 0.10            
XDATA_Z0_TOWN   (NUT_DENS) = 3.0000
XDATA_BLD       (NUT_DENS) = 0.5  
XDATA_BLD_HEIGHT(NUT_DENS) = 30.
XDATA_WALL_O_HOR(NUT_DENS) = 1.
XDATA_VEGTYPE   (NUT_DENS,NVT_TEBD) = 0.5
XDATA_VEGTYPE   (NUT_DENS,NVT_PARK) = 0.5 
!
DO JTYPE = NUT_SUB1,NUT_SUB2
  XDATA_GARDEN    (JTYPE) = 0.10          
  XDATA_Z0_TOWN   (JTYPE) = 3.0000
  XDATA_BLD       (JTYPE) = 0.5  
  XDATA_BLD_HEIGHT(JTYPE) = 30.
  XDATA_WALL_O_HOR(JTYPE) = 1.
ENDDO
!
XDATA_VEGTYPE   (NUT_SUB1,NVT_TEBD) = 0.5
XDATA_VEGTYPE   (NUT_SUB1,NVT_PARK) = 0.5 
!
XDATA_VEGTYPE   (NUT_SUB2,NVT_BOBD) = 0.5 
XDATA_VEGTYPE   (NUT_SUB2,NVT_PARK) = 0.5
!
XDATA_GARDEN    (NUT_ZIC) = 0.10       
XDATA_Z0_TOWN   (NUT_ZIC) = 2.0000
XDATA_BLD       (NUT_ZIC) = 0.5  
XDATA_BLD_HEIGHT(NUT_ZIC) = 5.
XDATA_WALL_O_HOR(NUT_ZIC) = 0.5
XDATA_VEGTYPE   (NUT_ZIC,NVT_NO) = 1.0
!
XDATA_GARDEN    (NUT_ROAD) = 0.10       
XDATA_Z0_TOWN   (NUT_ROAD) = 0.5
XDATA_BLD       (NUT_ROAD) = 0.1  
XDATA_BLD_HEIGHT(NUT_ROAD) = 5.
XDATA_WALL_O_HOR(NUT_ROAD) = 0.5
XDATA_VEGTYPE   (NUT_ROAD,NVT_GRAS) = 1.0
!
XDATA_GARDEN    (NUT_PORT) = 0.10       
XDATA_Z0_TOWN   (NUT_PORT) = 2.0
XDATA_BLD       (NUT_PORT) = 0.5  
XDATA_BLD_HEIGHT(NUT_PORT) = 20.
XDATA_WALL_O_HOR(NUT_PORT) = 1.
XDATA_VEGTYPE   (NUT_PORT,NVT_NO) = 1.0
!
XDATA_GARDEN    (NUT_AIR) = 0.70       
XDATA_Z0_TOWN   (NUT_AIR) = 0.01
XDATA_BLD       (NUT_AIR) = 0.1 
XDATA_BLD_HEIGHT(NUT_AIR) = 10.
XDATA_WALL_O_HOR(NUT_AIR) = 0.5
XDATA_VEGTYPE   (NUT_AIR,NVT_GRAS) = 1.0
!
XDATA_GARDEN    (NUT_MINE) = 0.90      
XDATA_Z0_TOWN   (NUT_MINE) = 0.1
XDATA_BLD       (NUT_MINE) = 0.1  
XDATA_BLD_HEIGHT(NUT_MINE) = 5.
XDATA_WALL_O_HOR(NUT_MINE) = 0.5
XDATA_VEGTYPE   (NUT_MINE,NVT_NO) = 1.0
!
XDATA_GARDEN    (NUT_PARK) = 0.90     
XDATA_Z0_TOWN   (NUT_PARK) = 0.5
XDATA_BLD       (NUT_PARK) = 0.1 
XDATA_BLD_HEIGHT(NUT_PARK) = 5.
XDATA_WALL_O_HOR(NUT_PARK) = 0.5
XDATA_VEGTYPE   (NUT_PARK,NVT_TEBD) = 0.5
XDATA_VEGTYPE   (NUT_PARK,NVT_PARK) = 0.5
!
XDATA_GARDEN    (NUT_SPOR) = 0.80  
XDATA_Z0_TOWN   (NUT_SPOR) = 1.0
XDATA_BLD       (NUT_SPOR) = 0.5  
XDATA_BLD_HEIGHT(NUT_SPOR) = 10.
XDATA_WALL_O_HOR(NUT_SPOR) = 1.
XDATA_VEGTYPE   (NUT_SPOR,NVT_TEBD) = 0.2
XDATA_VEGTYPE   (NUT_SPOR,NVT_PARK) = 0.8
!
IF (LHOOK) CALL DR_HOOK('INIT_TYPES_PARAM',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_TYPES_PARAM
