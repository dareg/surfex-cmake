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
!!	C. Lebeaupin Brossier  * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2008
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XPI
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_OCEAN_GRID_n
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!*      0.2    declarations of local variables
!
!
REAL                   :: RUP,RDOWN
INTEGER            :: IK1,IK2,IKK1,IMIN,IMAX
INTEGER            :: JLOOP
REAL :: DELTAZ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!       1.     Allocations
!              -----------
IF (LHOOK) CALL DR_HOOK('OCEAN_MERCATORVERGRID',0,ZHOOK_HANDLE)
ALLOCATE(XK1        (NOCKMIN:NOCKMAX))
ALLOCATE(XK2        (NOCKMIN:NOCKMAX))
ALLOCATE(XK3        (NOCKMIN:NOCKMAX))
ALLOCATE(XK4        (NOCKMIN:NOCKMAX))
!
ALLOCATE(XZHOC      (NOCKMIN:NOCKMAX))
ALLOCATE(XZ2        (NOCKMIN:NOCKMAX))
ALLOCATE(XDZ1       (NOCKMIN:NOCKMAX))
ALLOCATE(XDZ2       (NOCKMIN:NOCKMAX))
!
ALLOCATE(XRAY       (NOCKMIN:NOCKMAX))
!-------------------------------------------------------------------------------
!       2.     Ocean Grid Levels
!              -----------------
!
DELTAZ=5.
IK1=NOCKMIN+1
IK2=NOCKMIN+2
IKK1=NOCKMAX-1
!NOCKMAX=40
!print*,SIZE(Z1)
!
XZHOC(0)=0.
XZ2(0)=0.
!print*,SIZE(ZDZ1)
XDZ1(0)=XUNDEF
XDZ2(0)=0.
!
XZHOC(IK1)    = -1.
XZHOC(2)      = -5.
XZHOC(3)      = -10.
XZHOC(4)      = -15.
XZHOC(5)      = -20.
XZHOC(6)      = -25.
XZHOC(7)      = -30.
XZHOC(8)      = -40.
XZHOC(9)      = -50.
XZHOC(10)     = -60.
XZHOC(11)     = -75.
XZHOC(12)     = -100.
XZHOC(13)     = -125.
XZHOC(14)     = -150.
XZHOC(15)     = -175.
XZHOC(16)     = -200.
XZHOC(17)     = -225.
XZHOC(18)     = -250.
XZHOC(19)     = -300.
XZHOC(20)     = -400.
XZHOC(21)     = -500.
XZHOC(22)     = -600.
XZHOC(23)     = -700.
XZHOC(24)     = -800.
XZHOC(25)     = -900.
XZHOC(26)     = -1000.
XZHOC(27)     = -1100.
XZHOC(28)     = -1200.
XZHOC(29)     = -1300.
XZHOC(30)     = -1400.
XZHOC(31)     = -1500.
XZHOC(32)     = -1750.
XZHOC(33)     = -2000.
XZHOC(34)     = -2250.
XZHOC(35)     = -2500.
XZHOC(36)     = -2750.
XZHOC(37)     = -3000.
XZHOC(38)     = -3250.
XZHOC(IKK1)   = -3500.
XZHOC(NOCKMAX)  = -4000.

DO JLOOP=IK1,IKK1
  XZ2(JLOOP)   = ((XZHOC(JLOOP+1)-XZHOC(JLOOP))/2.)+XZHOC(JLOOP)
ENDDO
XZ2(NOCKMAX)  = -4250.

DO JLOOP=IK1,IKK1
  XDZ1(JLOOP)  = -XZHOC(JLOOP+1)+XZHOC(JLOOP)
  XDZ2(JLOOP)  = -XZ2(JLOOP)+XZ2(JLOOP-1)
ENDDO

XDZ2(NOCKMAX) = -XZ2(NOCKMAX)+XZ2(IKK1)
XDZ1(NOCKMAX) = XDZ1(IKK1)
!
!!
!!
!!       3.     Grid Parameters
!!              ---------------
RUP=1.
DO JLOOP=IK1,NOCKMAX
  XK1(JLOOP)=-1. / (XDZ2(JLOOP)*XDZ1(JLOOP-1))
  XK4(JLOOP)= 1. / (XDZ1(JLOOP)*XDZ1(JLOOP))
  RDOWN= RAYO(XZ2(JLOOP))
  XRAY(JLOOP)= RUP-RDOWN
  RUP=RDOWN
ENDDO
DO JLOOP=IK1,IKK1
  XK2(JLOOP)=-1. / (XDZ2(JLOOP)*XDZ1(JLOOP))
  XK3(JLOOP)=-1. / (XDZ1(JLOOP)*XDZ2(JLOOP+1))
ENDDO
XK2(NOCKMAX)=XK2(IKK1)
XK3(NOCKMAX)=0.
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
