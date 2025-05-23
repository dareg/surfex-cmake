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

USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE MODD_DATA_COVER,     ONLY : XDATA_TOWN, XDATA_NATURE, XDATA_SEA, XDATA_WATER,   &
                                XDATA_VEGTYPE, XDATA_GARDEN, XDATA_Z0_TOWN, &
                                XDATA_BLD, XDATA_BLD_HEIGHT, XDATA_WALL_O_HOR, &
                                XDATA_H_TRAFFIC, XDATA_H_INDUSTRY
!
USE MODD_DATA_COVER_PAR, ONLY : NTYPE, NUT_CPHR, NUT_CPMR, NUT_CPLR, NUT_OPHR, &
                                NUT_OPMR, NUT_OPLR, NUT_LWLR, NUT_LALR, NUT_SPAR, &
                                NUT_INDU, NVT_NO, NVT_TEBD, NVT_BOBD, NVT_GRAS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_TYPES_PARAM',0,ZHOOK_HANDLE)
!
! SEA
!
XDATA_SEA(1:NTYPE(1)) = 1.
!
! WATER
!
XDATA_WATER(NTYPE(1)+1:SUM(NTYPE(1:2))) = 1.
!
! NATURE
!
DO JTYPE = 1,NTYPE(3)
  XDATA_NATURE (SUM(NTYPE(1:2))+JTYPE) = 1.
  XDATA_VEGTYPE(SUM(NTYPE(1:2))+JTYPE,JTYPE) = 1.
ENDDO
!
! TOWN
!
XDATA_TOWN(SUM(NTYPE(1:3))+1:SUM(NTYPE(1:4))) = 1.
!
XDATA_GARDEN(NUT_CPHR) = 0.10     
XDATA_GARDEN(NUT_CPMR) = 0.10
XDATA_GARDEN(NUT_CPLR) = 0.10
XDATA_GARDEN(NUT_OPHR) = 0.30     
XDATA_GARDEN(NUT_OPMR) = 0.30
XDATA_GARDEN(NUT_OPLR) = 0.40
XDATA_GARDEN(NUT_LWLR) = 0.10
XDATA_GARDEN(NUT_LALR) = 0.10
XDATA_GARDEN(NUT_SPAR) = 0.70
XDATA_GARDEN(NUT_INDU) = 0.40
!
XDATA_Z0_TOWN(NUT_CPHR) = 3.
XDATA_Z0_TOWN(NUT_CPMR) = 1.
XDATA_Z0_TOWN(NUT_CPLR) = 0.5
XDATA_Z0_TOWN(NUT_OPHR) = 2.0
XDATA_Z0_TOWN(NUT_OPMR) = 0.5
XDATA_Z0_TOWN(NUT_OPLR) = 0.5
XDATA_Z0_TOWN(NUT_LWLR) = 0.25
XDATA_Z0_TOWN(NUT_LALR) = 0.25
XDATA_Z0_TOWN(NUT_SPAR) = 0.5
XDATA_Z0_TOWN(NUT_INDU) = 0.5
!
XDATA_BLD(NUT_CPHR) = 0.5
XDATA_BLD(NUT_CPMR) = 0.55
XDATA_BLD(NUT_CPLR) = 0.55
XDATA_BLD(NUT_OPHR) = 0.30
XDATA_BLD(NUT_OPMR) = 0.30
XDATA_BLD(NUT_OPLR) = 0.30
XDATA_BLD(NUT_LWLR) = 0.75
XDATA_BLD(NUT_LALR) = 0.40
XDATA_BLD(NUT_SPAR) = 0.10
XDATA_BLD(NUT_INDU) = 0.25
!
XDATA_BLD_HEIGHT(NUT_CPHR) = 75
XDATA_BLD_HEIGHT(NUT_CPMR) = 20
XDATA_BLD_HEIGHT(NUT_CPLR) = 5
XDATA_BLD_HEIGHT(NUT_OPHR) = 75  
XDATA_BLD_HEIGHT(NUT_OPMR) = 20
XDATA_BLD_HEIGHT(NUT_OPLR) = 5
XDATA_BLD_HEIGHT(NUT_LWLR) = 3
XDATA_BLD_HEIGHT(NUT_LALR) = 5
XDATA_BLD_HEIGHT(NUT_SPAR) = 5
XDATA_BLD_HEIGHT(NUT_INDU) = 10
!
XDATA_WALL_O_HOR(NUT_CPHR) = 4.0
XDATA_WALL_O_HOR(NUT_CPMR) = 1.3
XDATA_WALL_O_HOR(NUT_CPLR) = 0.9
XDATA_WALL_O_HOR(NUT_OPHR) = 1.4  
XDATA_WALL_O_HOR(NUT_OPMR) = 0.7
XDATA_WALL_O_HOR(NUT_OPLR) = 0.7
XDATA_WALL_O_HOR(NUT_LWLR) = 0.75
XDATA_WALL_O_HOR(NUT_LALR) = 0.24
XDATA_WALL_O_HOR(NUT_SPAR) = 0.36
XDATA_WALL_O_HOR(NUT_INDU) = 0.45
!
XDATA_H_TRAFFIC(NUT_CPHR) = 20
XDATA_H_TRAFFIC(NUT_CPMR) = 10
XDATA_H_TRAFFIC(NUT_CPLR) = 5
XDATA_H_TRAFFIC(NUT_OPHR) =  20 
XDATA_H_TRAFFIC(NUT_OPMR) = 10
XDATA_H_TRAFFIC(NUT_OPLR) = 5
XDATA_H_TRAFFIC(NUT_LWLR) = 5
XDATA_H_TRAFFIC(NUT_LALR) = 5
XDATA_H_TRAFFIC(NUT_SPAR) = 5
XDATA_H_TRAFFIC(NUT_INDU) = 5
!
XDATA_H_INDUSTRY(NUT_CPHR) = 0.
XDATA_H_INDUSTRY(NUT_CPMR) = 0.
XDATA_H_INDUSTRY(NUT_CPLR) = 0.
XDATA_H_INDUSTRY(NUT_OPHR) = 0.
XDATA_H_INDUSTRY(NUT_OPMR) = 0.
XDATA_H_INDUSTRY(NUT_OPLR) = 0.
XDATA_H_INDUSTRY(NUT_LWLR) = 0.
XDATA_H_INDUSTRY(NUT_LALR) = 0.
XDATA_H_INDUSTRY(NUT_SPAR) = 0.
XDATA_H_INDUSTRY(NUT_INDU) = 0.
!
XDATA_VEGTYPE(NUT_CPHR,NVT_GRAS) = 0.0
XDATA_VEGTYPE(NUT_CPMR,NVT_GRAS) = 0.0
XDATA_VEGTYPE(NUT_CPLR,NVT_GRAS) = 0.0 
XDATA_VEGTYPE(NUT_OPHR,NVT_GRAS) = 0.4  
XDATA_VEGTYPE(NUT_OPMR,NVT_GRAS) = 0.4
XDATA_VEGTYPE(NUT_OPLR,NVT_GRAS) = 0.4
XDATA_VEGTYPE(NUT_LWLR,NVT_GRAS) = 0.5
XDATA_VEGTYPE(NUT_LALR,NVT_GRAS) = 0.5
XDATA_VEGTYPE(NUT_SPAR,NVT_GRAS) = 0.5
XDATA_VEGTYPE(NUT_INDU,NVT_GRAS) = 0.4

!
XDATA_VEGTYPE(NUT_CPHR,NVT_NO) = 0.0
XDATA_VEGTYPE(NUT_CPMR,NVT_NO) = 0.0
XDATA_VEGTYPE(NUT_CPLR,NVT_NO) = 0.0 
XDATA_VEGTYPE(NUT_OPHR,NVT_NO) = 0.2  
XDATA_VEGTYPE(NUT_OPMR,NVT_NO) = 0.2
XDATA_VEGTYPE(NUT_OPLR,NVT_NO) = 0.2
XDATA_VEGTYPE(NUT_LWLR,NVT_NO) = 0.5
XDATA_VEGTYPE(NUT_LALR,NVT_NO) = 0.5
XDATA_VEGTYPE(NUT_SPAR,NVT_NO) = 0.2
XDATA_VEGTYPE(NUT_INDU,NVT_NO) = 0.6
!
XDATA_VEGTYPE(NUT_CPHR,NVT_TEBD) = 1.0 
XDATA_VEGTYPE(NUT_CPMR,NVT_TEBD) = 1.0
XDATA_VEGTYPE(NUT_CPLR,NVT_TEBD) = 1.0
XDATA_VEGTYPE(NUT_OPHR,NVT_TEBD) = 0.4  
XDATA_VEGTYPE(NUT_OPMR,NVT_TEBD) = 0.4
XDATA_VEGTYPE(NUT_OPLR,NVT_TEBD) = 0.4
XDATA_VEGTYPE(NUT_LWLR,NVT_TEBD) = 0.0
XDATA_VEGTYPE(NUT_LALR,NVT_TEBD) = 0.0
XDATA_VEGTYPE(NUT_SPAR,NVT_TEBD) = 0.3
XDATA_VEGTYPE(NUT_INDU,NVT_TEBD) = 0.0


!
IF (LHOOK) CALL DR_HOOK('INIT_TYPES_PARAM',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_TYPES_PARAM
