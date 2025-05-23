!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE CREATE_FILE(HFILE,KDIMS,HNAM_DIM,KFILE_ID,KDIM_ID)
!
USE MODI_HANDLE_ERR
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
USE NETCDF
!
IMPLICIT NONE
!
!
! dummy arguments
!
 CHARACTER(LEN=*),            INTENT(IN)          :: HFILE
INTEGER, DIMENSION(:),        INTENT(IN)          :: KDIMS
 CHARACTER(LEN=*),DIMENSION(:),INTENT(IN)         :: HNAM_DIM
INTEGER,                      INTENT(INOUT)       :: KFILE_ID
INTEGER, DIMENSION(:),        INTENT(INOUT)       :: KDIM_ID
!
!* local variables
!
INTEGER :: IRET,ILEN,JNBDIM,INBDIM
 CHARACTER(LEN=50) :: YFILE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!----------------------------------------------------
!
! 1.0  create files
!-------------------
IF (LHOOK) CALL DR_HOOK('CREATE_FILE',0,ZHOOK_HANDLE)
ILEN  = LEN_TRIM(HFILE)
YFILE = HFILE(:LEN_TRIM(HFILE))
IRET  = NF90_CREATE(YFILE, NF90_64BIT_OFFSET, KFILE_ID)

IF (IRET.NE.NF90_NOERR) CALL HANDLE_ERR(IRET,'CREATE_FILE')

! 2.0  define dimensions
!-----------------------
INBDIM=SIZE(KDIMS)
DO JNBDIM=1,INBDIM
   IRET = NF90_DEF_DIM(KFILE_ID,HNAM_DIM(JNBDIM),KDIMS(JNBDIM),KDIM_ID(JNBDIM))
ENDDO
IF (LHOOK) CALL DR_HOOK('CREATE_FILE',1,ZHOOK_HANDLE)

END SUBROUTINE CREATE_FILE
