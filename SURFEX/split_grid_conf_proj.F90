!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###########################################################
      SUBROUTINE SPLIT_GRID_CONF_PROJ(HPROGRAM,KDIM_FULL,KSIZE_FULL,KGRID_PAR,PGRID_PAR,KHALO)
!     ###########################################################
!!
!!    PURPOSE
!!    -------
!!   This program splits a PGD grid on several processors (according to host program)
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson                   Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     08/11
!!      M.Moge     02/15  using PGRID_PAR(11) instead of KDIM_FULL
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODE_GRIDTYPE_CONF_PROJ
USE MODE_SPLIT_GRID_PARAMETER
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
 CHARACTER(LEN=6),   INTENT(IN)    :: HPROGRAM  ! host program 
INTEGER,            INTENT(IN)    :: KDIM_FULL ! total number of points
INTEGER,            INTENT(OUT)   :: KSIZE_FULL! number of points on this processor
INTEGER,            INTENT(INOUT) :: KGRID_PAR ! size of PGRID_PAR pointer
REAL, DIMENSION(:), POINTER       :: PGRID_PAR ! parameters defining this grid
INTEGER, OPTIONAL,  INTENT(IN)    :: KHALO ! size of the Halo
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!* original grid
REAL                            :: ZLAT0, ZLON0, ZRPK, ZBETA, ZLATOR, ZLONOR
INTEGER                         :: IIMAX, IJMAX
REAL, DIMENSION(INT(PGRID_PAR(11)))   :: ZX, ZY, ZDX, ZDY
!
!* splitted grid on processor
INTEGER                         :: IIMAX_SPLIT, IJMAX_SPLIT
INTEGER                         :: ILONE, ILATE
INTEGER :: IWIDTH_I_X ! width of I zone
INTEGER :: IWIDTH_I_Y ! width of I zone
REAL    :: ZTRUNC     ! spectral truncation factor
REAL, DIMENSION(:), ALLOCATABLE :: ZX_SPLIT, ZY_SPLIT, ZDX_SPLIT, ZDY_SPLIT
!
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPLIT_GRID_CONF_PROJ',0,ZHOOK_HANDLE)
!
!*    1.      Gets Parameters of the Grid
!
 CALL GET_GRIDTYPE_CONF_PROJ(PGRID_PAR,ZLAT0,ZLON0,ZRPK,ZBETA,&
                            ZLATOR,ZLONOR,IIMAX,IJMAX,       &
                            ZX,ZY,ZDX,ZDY,KLATE=ILATE,KLONE=ILONE,&
                            KWIDTH_I_X=IWIDTH_I_X,KWIDTH_I_Y=IWIDTH_I_Y, &
                            PTRUNC=ZTRUNC)
!
!
!*    2.      Splits the (pertinent) parameters of the grid
!
IF (PRESENT(KHALO)) THEN
  CALL SPLIT_GRID_PARAMETERN0(HPROGRAM,'CONF PROJ ','IMAX  ',IIMAX,IIMAX_SPLIT,KHALO)
  CALL SPLIT_GRID_PARAMETERN0(HPROGRAM,'CONF PROJ ','JMAX  ',IJMAX,IJMAX_SPLIT,KHALO)
ELSE
  CALL SPLIT_GRID_PARAMETERN0(HPROGRAM,'CONF PROJ ','IMAX  ',IIMAX,IIMAX_SPLIT)
  CALL SPLIT_GRID_PARAMETERN0(HPROGRAM,'CONF PROJ ','JMAX  ',IJMAX,IJMAX_SPLIT)
ENDIF
!
KSIZE_FULL = IIMAX_SPLIT * IJMAX_SPLIT
!
ALLOCATE(ZX_SPLIT (KSIZE_FULL))
ALLOCATE(ZY_SPLIT (KSIZE_FULL))
ALLOCATE(ZDX_SPLIT(KSIZE_FULL))
ALLOCATE(ZDY_SPLIT(KSIZE_FULL))
!
IF (PRESENT(KHALO)) THEN
  CALL SPLIT_GRID_PARAMETERX1(HPROGRAM,'CONF PROJ ','XX    ',SIZE(ZX),KSIZE_FULL,ZX,ZX_SPLIT,IIMAX,IJMAX,KHALO)
  CALL SPLIT_GRID_PARAMETERX1(HPROGRAM,'CONF PROJ ','YY    ',SIZE(ZY),KSIZE_FULL,ZY,ZY_SPLIT,IIMAX,IJMAX,KHALO)
  CALL SPLIT_GRID_PARAMETERX1(HPROGRAM,'CONF PROJ ','DX    ',SIZE(ZDX),KSIZE_FULL,ZDX,ZDX_SPLIT,IIMAX,IJMAX,KHALO)
  CALL SPLIT_GRID_PARAMETERX1(HPROGRAM,'CONF PROJ ','DY    ',SIZE(ZDY),KSIZE_FULL,ZDY,ZDY_SPLIT,IIMAX,IJMAX,KHALO)
ELSE
  CALL SPLIT_GRID_PARAMETERX1(HPROGRAM,'CONF PROJ ','XX    ',KDIM_FULL,KSIZE_FULL,ZX,ZX_SPLIT)
  CALL SPLIT_GRID_PARAMETERX1(HPROGRAM,'CONF PROJ ','YY    ',KDIM_FULL,KSIZE_FULL,ZY,ZY_SPLIT)
  CALL SPLIT_GRID_PARAMETERX1(HPROGRAM,'CONF PROJ ','DX    ',KDIM_FULL,KSIZE_FULL,ZDX,ZDX_SPLIT)
  CALL SPLIT_GRID_PARAMETERX1(HPROGRAM,'CONF PROJ ','DY    ',KDIM_FULL,KSIZE_FULL,ZDY,ZDY_SPLIT)
ENDIF
!
!
!*    3.      Stores Parameters of the Grid in grid pointer
!
NULLIFY(PGRID_PAR)
 CALL PUT_GRIDTYPE_CONF_PROJ(PGRID_PAR,ZLAT0,ZLON0,ZRPK,ZBETA,       &
                            ZLATOR,ZLONOR,IIMAX_SPLIT,IJMAX_SPLIT,  &
                            ZX_SPLIT,ZY_SPLIT,ZDX_SPLIT,ZDY_SPLIT, &
                            ILATE,ILONE,IWIDTH_I_X,IWIDTH_I_Y,ZTRUNC)
                            !
!
KGRID_PAR = SIZE(PGRID_PAR)
!
DEALLOCATE(ZX_SPLIT )
DEALLOCATE(ZY_SPLIT )
DEALLOCATE(ZDX_SPLIT)
DEALLOCATE(ZDY_SPLIT)
!
IF (LHOOK) CALL DR_HOOK('SPLIT_GRID_CONF_PROJ',1,ZHOOK_HANDLE)
!_______________________________________________________________________________
!
END SUBROUTINE SPLIT_GRID_CONF_PROJ

