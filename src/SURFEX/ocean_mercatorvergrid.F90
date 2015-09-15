!     #########
    SUBROUTINE OCEAN_MERCATORVERGRID
!   ######################################################################
!
!!****  *OCEAN_MERCATORVERGRID*  
!!
!!    PURPOSE
!!    -------
!
!     Define the vertical ocean grid
!         
!     
!!**  METHOD
!!    ------
!
!
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_OCEAN_CST
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      C. Lebeaupin Brossier  * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2008
!       D.BARBARY 11/2014 : HPROGRAM,OUNIF in Calling OCEAN_MERCATORVERGRID
!                           Reading oceanic level and depth
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_OCEAN_GRID
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.2    declarations of local variables
!
!
REAL                   :: ZUP,ZDOWN
INTEGER            :: JLOOP
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!       1.     Allocations
!              -----------
IF (LHOOK) CALL DR_HOOK('OCEAN_MERCATORVERGRID',0,ZHOOK_HANDLE)

!       1.     Principal grid
!              --------------
!
!-------------------------------------------------------------------------------
!       2.     Ocean Secondary Grids 
!              ---------------------
!
ALLOCATE(XK1        (NOCKMIN:NOCKMAX))
ALLOCATE(XK2        (NOCKMIN:NOCKMAX))
ALLOCATE(XK3        (NOCKMIN:NOCKMAX))
ALLOCATE(XK4        (NOCKMIN:NOCKMAX))
!
ALLOCATE(XZ2        (NOCKMIN:NOCKMAX))
ALLOCATE(XDZ1       (NOCKMIN:NOCKMAX))
ALLOCATE(XDZ2       (NOCKMIN:NOCKMAX))
!
ALLOCATE(XRAY       (NOCKMIN:NOCKMAX))
!
!ADD STEPH TO GET SAME RESULTS THAN BEFORE
!XZ2(NOCKMIN) = 0.
DO JLOOP = NOCKMIN,NOCKMAX-1
!DO JLOOP = NOCKMIN+1,NOCKMAX-1
!END ADD STEPH
  XZ2(JLOOP) = ((XZHOC(JLOOP+1)-XZHOC(JLOOP))/2.) + XZHOC(JLOOP)
ENDDO
!CLB
XZ2(NOCKMAX) = XZHOC(NOCKMAX) + ((XZHOC(NOCKMAX)-XZHOC(NOCKMAX-1))/2.) !dernier niveau+delta/2
!
!
XDZ1(0) = XUNDEF
XDZ2(0) = XUNDEF
DO JLOOP = NOCKMIN+1,NOCKMAX-1
  XDZ1(JLOOP) = XZHOC(JLOOP) - XZHOC(JLOOP+1)
  XDZ2(JLOOP) = XZ2(JLOOP-1) - XZ2(JLOOP)
ENDDO
!!
XDZ1(NOCKMAX) = XDZ1(NOCKMAX-1)
XDZ2(NOCKMAX) = XZ2(NOCKMAX-1) -XZ2(NOCKMAX) 
!
!!       3.     Grid Parameters
!!              ---------------
!
ZUP=1.
!ADD STEPH TO GET SAME RESULTS THAN BEFORE
ZDOWN = RAYO(XZ2(0))
XRAY(0) = ZUP - ZDOWN
ZUP = ZDOWN
!ADD STEPH END
!
XK1(0) = XUNDEF
XK4(0) = XUNDEF
DO JLOOP = NOCKMIN+1,NOCKMAX
  XK1(JLOOP) = -1. / ( XDZ2(JLOOP) * XDZ1(JLOOP-1) )
  XK4(JLOOP) =  1. / ( XDZ1(JLOOP) * XDZ1(JLOOP) )
  ZDOWN = RAYO(XZ2(JLOOP))
  XRAY(JLOOP) = ZUP - ZDOWN
  ZUP = ZDOWN
ENDDO
!
XK2(0) = XUNDEF
XK3(0) = XUNDEF
DO JLOOP = NOCKMIN+1,NOCKMAX-1
  XK2(JLOOP) = -1. / (XDZ2(JLOOP)*XDZ1(JLOOP))
  XK3(JLOOP) = -1. / (XDZ1(JLOOP)*XDZ2(JLOOP+1))
ENDDO
!
XK2(NOCKMAX) = XK2(NOCKMAX-1)
XK3(NOCKMAX) = 0.
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('OCEAN_MERCATORVERGRID',1,ZHOOK_HANDLE)
CONTAINS
!rayo
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!!              #########################################
             FUNCTION RAYO(Z) RESULT(RR)
!              #########################################
!
!
!!****  *RAYOFCTX* 
!!
!!    PURPOSE
!!    -------
!compute solar penetration coefficient
!    
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!     
!!    REFERENCE
!!    ---------
!!    Paulson and Simpson 1977
!!     
!!    AUTHOR
!!    ------
!!     C. Lebeaupin  *Meteo-France* (adapted from S. Belamari's code)
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     27/02/2006
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!
USE MODD_OCEAN_CSTS,ONLY : XR,XD1,XD2
!
!*      0.1    declarations of arguments
!
REAL :: RR,Z      
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE RAYO
!               ------------
!
IF (LHOOK) CALL DR_HOOK('RAYO',0,ZHOOK_HANDLE)
RR = XR*EXP(Z/XD1) + (1-XR)*EXP(Z/XD2)
IF (LHOOK) CALL DR_HOOK('RAYO',1,ZHOOK_HANDLE)
!
END FUNCTION RAYO
!
!-------------------------------------------------------------------------------
END SUBROUTINE OCEAN_MERCATORVERGRID
