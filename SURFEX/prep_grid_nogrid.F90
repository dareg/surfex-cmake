!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE PREP_GRID_NOGRID(HFILETYPE,HINTERP_TYPE,KNI)
!     ##########################################################################
!
!**** *PREP_GRID_NOGRID*  - Support for a 1D geometry

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014

USE MODI_READ_SURF

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE
CHARACTER(LEN=6),  INTENT(IN)    :: HFILETYPE    ! file type
CHARACTER(LEN=6),  INTENT(OUT)   :: HINTERP_TYPE ! Grid type
INTEGER,           INTENT(OUT)   :: KNI          ! number of points

CHARACTER(LEN=12) :: YRECFM    
INTEGER           :: IRESP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('PREP_GRID_NOGRID',0,ZHOOK_HANDLE)

YRECFM = 'NOGRID_SIZE'
CALL READ_SURF(HFILETYPE,YRECFM,KNI,IRESP)

HINTERP_TYPE='BUFFER'

IF (LHOOK) CALL DR_HOOK('PREP_GRID_NOGRID',1,ZHOOK_HANDLE)

END SUBROUTINE PREP_GRID_NOGRID
