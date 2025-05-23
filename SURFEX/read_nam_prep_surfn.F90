!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE READ_NAM_PREP_SURF_n(HPROGRAM)
!     #######################################################
!
!---------------------------------------
!
USE MODD_SURF_PAR, ONLY : XUNDEF, NUNDEF
USE MODN_PREP_SURF_ATM
!
USE MODI_TEST_NAM_VAR_SURF
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODE_POS_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
!
INTEGER :: ILUNAM         ! logical unit of namelist file
INTEGER :: ILUOUT
LOGICAL :: GFOUND         ! Return code when searching namelist
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!---------------------------------------
!
IF (LHOOK) CALL DR_HOOK('READ_NAM_PREP_SURF_N',0,ZHOOK_HANDLE)
!
NHALO_PREP = 0
!
NYEAR=NUNDEF
NMONTH=NUNDEF
NDAY=NUNDEF
XTIME=XUNDEF
CFILE     = '                         '
CFILETYPE = '      '
CFILEPGD     = '                         '
CFILEPGDTYPE = '      '
LWRITE_EXTERN = .FALSE.
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!
 CALL OPEN_NAMELIST(HPROGRAM,ILUNAM)
 CALL POSNAM(ILUNAM,'NAM_PREP_SURF_ATM',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_PREP_SURF_ATM)
 CALL CLOSE_NAMELIST(HPROGRAM,ILUNAM)
!
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CFILETYPE',   CFILETYPE,   '      ','GRIB  ','MESONH','ASCII ','LFI   ',&
         'FA    ','NC    ')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CFILEPGDTYPE',   CFILEPGDTYPE,   '      ','GRIB  ','MESONH','ASCII ',&
 'LFI   ','FA    ','NC    ')
IF (LHOOK) CALL DR_HOOK('READ_NAM_PREP_SURF_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_NAM_PREP_SURF_n
