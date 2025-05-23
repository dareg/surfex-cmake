!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE OPEN_NAMELIST_OL(HPROGRAM,KLUNAM,HFILE)
!     #######################################################
!
!!****  *OPEN_NAMELIST_OL* - opens namelists files for surface (OFFLINE universe)
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004 
!       10/2014 : add status='old'  E. Martin
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER,           INTENT(OUT) :: KLUNAM   ! logical unit of namelist
 CHARACTER(LEN=28), INTENT(IN)  :: HFILE ! ASCII file to open
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
 CHARACTER(LEN=28) :: YNAM
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
INTEGER            :: IERR
!
!-------------------------------------------------------------------------------
!
!* reading of namelist
!  -------------------
!
IF (LHOOK) CALL DR_HOOK('OPEN_NAMELIST_OL',0,ZHOOK_HANDLE)
IF (LEN_TRIM(HFILE)>0) THEN
  YNAM = HFILE
ELSE
  YNAM='OPTIONS.nam'
END IF
!
KLUNAM=11
OPEN(KLUNAM,FILE=YNAM,ACTION='READ',FORM="FORMATTED",POSITION="REWIND", &
        STATUS='OLD',IOSTAT=IERR)
 IF (IERR /= 0 ) THEN
    CALL ABOR1_SFX ('ERROR WHILE OPENING '//YNAM//' THIS FILE IS MISSING'// &
                  ' IN THE RUN DIRECTORY')
  ENDIF


IF (LHOOK) CALL DR_HOOK('OPEN_NAMELIST_OL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE OPEN_NAMELIST_OL
