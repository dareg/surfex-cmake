!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_SURF_ATM_CONF(HPROGRAM)
!     #######################################################
!
!!****  *READ_SURF_ATM_CONF* - reads the general configuration for surface
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
!!      Original    01/2003 
!!      B. Decharme 06/2013 CIMPLICIT_WIND is now in NAM_SURF_REPROD_OPER
!!      J. Masek    08/2023 Setting of XZ0_OFFSET; writing final values of
!!                          important NAM_SURF_ATM variables
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODE_POS_SURF
!
USE MODN_CHS_ORILAM
USE MODN_SURF_ATM
USE MODN_WRITE_SURF_ATM
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling GROUND
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
LOGICAL           :: GFOUND         ! Return code when searching namelist
INTEGER           :: ILUOUT         ! logical unit of output file
INTEGER           :: INAM           ! logical unit of namelist file
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!* get output listing file logical unit
!
IF (LHOOK) CALL DR_HOOK('READ_SURF_ATM_CONF',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!* open namelist file
!
 CALL OPEN_NAMELIST(HPROGRAM,INAM)
!
!* reading of namelist
!  -------------------
!
 CALL POSNAM(INAM,'NAM_CHS_ORILAM',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=INAM,NML=NAM_CHS_ORILAM)
!
 CALL POSNAM(INAM,'NAM_SURF_ATM',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=INAM,NML=NAM_SURF_ATM)
!
IF (LSLOPE.AND.(.NOT.LNOSOF)) THEN
  CALL ABOR1_SFX(' if LSLOPE=T, LNOSOF must be TRUE ')
ENDIF
!
 CALL POSNAM(INAM,'NAM_WRITE_SURF_ATM',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=INAM,NML=NAM_WRITE_SURF_ATM)
!
!
!* close namelist file
!
 CALL CLOSE_NAMELIST(HPROGRAM,INAM)

IF (LZ0_AVG_EXACT) THEN
  XZ0_OFFSET=1.
ELSE
  XZ0_OFFSET=0.
ENDIF

WRITE(ILUOUT,'(A)')              'IMPORTANT SURFEX NAM_SURF_ATM SETTINGS:'
WRITE(ILUOUT,'(A,L)')            '  LDRAG_COEF_ARP = ',LDRAG_COEF_ARP
WRITE(ILUOUT,'(A,L)')            '  LALDTHRES      = ',LALDTHRES
WRITE(ILUOUT,'(A,L)')            '  LNOSOF         = ',LNOSOF
WRITE(ILUOUT,'(A,L)')            '  LSLOPE         = ',LSLOPE
WRITE(ILUOUT,'(A,L)')            '  LVERTSHIFT     = ',LVERTSHIFT
WRITE(ILUOUT,'(A,L)')            '  LVSHIFT_LW     = ',LVSHIFT_LW
WRITE(ILUOUT,'(A,L)')            '  LVSHIFT_PRCP   = ',LVSHIFT_PRCP
WRITE(ILUOUT,'(A,L)')            '  LVZIUSTAR0_ARP = ',LVZIUSTAR0_ARP
WRITE(ILUOUT,'(A,L)')            '  LRRGUST_ARP    = ',LRRGUST_ARP
WRITE(ILUOUT,'(A,L)')            '  LCPL_ARP       = ',LCPL_ARP
WRITE(ILUOUT,'(A,L)')            '  LQVNPLUS       = ',LQVNPLUS
WRITE(ILUOUT,'(A,L)')            '  LCPL_GCM       = ',LCPL_GCM
WRITE(ILUOUT,'(A,L)')            '  LZ0_AVG_EXACT  = ',LZ0_AVG_EXACT
WRITE(ILUOUT,'(A,L)')            '  LZ0_EFF        = ',LZ0_EFF
WRITE(ILUOUT,'(A,ES11.4)')       '  XZ0_OFFSET     = ',XZ0_OFFSET
WRITE(ILUOUT,'(A,L)')            '  LZ0SNOWH_ARP   = ',LZ0SNOWH_ARP
WRITE(ILUOUT,'(A,ES11.4)')       '  XRZ0_TO_HEIGHT = ',XRZ0_TO_HEIGHT
WRITE(ILUOUT,'(A,ES11.4)')       '  XFACZ0         = ',XFACZ0
WRITE(ILUOUT,'(A,ES11.4)')       '  XTAU_ICE       = ',XTAU_ICE
WRITE(ILUOUT,'(A,ES11.4)')       '  XCISMIN        = ',XCISMIN
WRITE(ILUOUT,'(A,ES11.4)')       '  XVMODMIN       = ',XVMODMIN
WRITE(ILUOUT,'(A,ES11.4)')       '  XWNEW          = ',XWNEW
WRITE(ILUOUT,'(A,ES11.4)')       '  XWCRN          = ',XWCRN
WRITE(ILUOUT,'(A,ES11.4)')       '  XEDB           = ',XEDB
WRITE(ILUOUT,'(A,ES11.4)')       '  XEDC           = ',XEDC
WRITE(ILUOUT,'(A,ES11.4)')       '  XEDD           = ',XEDD
WRITE(ILUOUT,'(A,ES11.4)')       '  XEDK           = ',XEDK
WRITE(ILUOUT,'(A,ES11.4)')       '  XUSURIC        = ',XUSURIC
WRITE(ILUOUT,'(A,ES11.4)')       '  XUSURID        = ',XUSURID
WRITE(ILUOUT,'(A,ES11.4)')       '  XUSURICL       = ',XUSURICL
WRITE(ILUOUT,'(A,ES11.4)')       '  XVCHRNK        = ',XVCHRNK
WRITE(ILUOUT,'(A,ES11.4)')       '  XVZ0CM         = ',XVZ0CM
WRITE(ILUOUT,'(A,ES11.4)')       '  XRIMAX         = ',XRIMAX
WRITE(ILUOUT,'(A,ES11.4)')       '  XDELTA_MAX     = ',XDELTA_MAX
WRITE(ILUOUT,'(A,ES11.4)')       '  XWINDMIN       = ',XWINDMIN
WRITE(ILUOUT,'(A,ES11.4)')       '  XRZHZ0M        = ',XRZHZ0M
WRITE(ILUOUT,'(A,ES11.4)')       '  XVZIUSTAR0     = ',XVZIUSTAR0
WRITE(ILUOUT,'(A,ES11.4)')       '  XRRSCALE       = ',XRRSCALE
WRITE(ILUOUT,'(A,ES11.4)')       '  XRRGAMMA       = ',XRRGAMMA
WRITE(ILUOUT,'(A,ES11.4)')       '  XUTILGUST      = ',XUTILGUST
WRITE(ILUOUT,'(A,ES11.4)')       '  XCO2UNCPL      = ',XCO2UNCPL

IF (LHOOK) CALL DR_HOOK('READ_SURF_ATM_CONF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_SURF_ATM_CONF
