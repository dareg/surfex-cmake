!     #########
SUBROUTINE SLT_VELGRAV1D(PSIG, PRG, PTA, PRHODREF, PRHOP, PMU, PVGK,PDPK, PVGG, PDPG)
!!   #######################################
!!
!!   PURPOSE
!!   -------
!!
!!   REFERENCE
!!   ---------
!!   none
!!
!!   AUTHOR
!!    ------
!!   P. Tulet (meteo france)
!!
!!   MODIFICATIONS
!!    -------------
!!
! Entry variables:
!
! PM(IN)       -Array of moments
!
!*************************************************************
! Exit variables:
!
! PFSED(IN)  -Array of moment variation due to dry deposition
!
!*************************************************************
! Variables used during the deposition velocity calculation
! 
! PDPK       -Polydisperse diffusivity (m2/s)
! PVGK       -Polydisperse settling velocity of the kth moment (m/s)
!************************************************************
!!
!!   IMPLICIT ARGUMENTS
USE MODD_SLT_SURF
!!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

!! Declarations d'arguments
!
REAL,                     INTENT(IN) :: PRHOP
REAL, DIMENSION(:,:), INTENT(IN) :: PSIG, PRG
REAL, DIMENSION(:),   INTENT(IN) :: PTA, PRHODREF
REAL, DIMENSION(:,:), INTENT(OUT) :: PVGK,PDPK
REAL, DIMENSION(:),   INTENT(OUT) :: PMU
REAL, DIMENSION(:,:),   INTENT(OUT) :: PVGG, PDPG
!
!!!! Declarations de variables internes
!
REAL, DIMENSION(size(PSIG,1)) :: ZLAMBDA

REAL, DIMENSION(size(PSIG,1)) :: ZVG,ZRG,ZLN2S

REAL, DIMENSION(size(PSIG,1)) :: ZKNG



REAL, PARAMETER :: gasmw=28.9644d0
REAL, PARAMETER :: rgas=8314.
REAL :: ZK, ZRD, ZCPD, ZP00, ZAVOGADRO, ZBOLTZ, ZMD, ZPI, ZG

INTEGER :: IJ, II
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SLT_VELGRAV1D',0,ZHOOK_HANDLE)
ZPI = 2.*ASIN(1.)
ZP00 = 1.E5 !REFERENCE PRESSURE
ZAVOGADRO =  6.0221367E+23
ZBOLTZ = 1.380658E-23 
ZMD   = 28.9644E-3
ZRD  = ZAVOGADRO * ZBOLTZ / ZMD
ZCPD = 7.* ZRD /2.
ZG  = 9.80665


! Sutherland's equation for viscosity
PMU(:)=1.8325d-5*416.16/(PTA(:)+120)*(PTA(:)/296.16)*SQRT(PTA(:)/296.16)

! Mean free path (Seinfeld and Pandis p455)
ZLAMBDA(:)=PMU(:)/PRHODREF(:)*sqrt(1.89d-4*gasmw/PTA(:))*1.e6


DO II=1,JPMODE_SLT
  ZRG(:)=PRG(:,II) * 1E-6 
  ZLN2S(:)=LOG(PSIG(:,II))**2 
  !
  ZKNG(:)=ZLAMBDA(:) / PRG(:,II) 
  !
  PVGG(:,II)= 2.*ZG*PRHOP*ZRG(:)**2 /(9.*PMU(:))
  PDPG(:,II)=ZBOLTZ*PTA(:)/ (6.*ZPI* ZRG(:)*PMU(:))


 
  do IJ=0,2
!
    ZK=real(3*IJ)
   

    PDPK(:,3*II+IJ-2)=PDPG(:,II)*(exp((-2.*ZK+1.)/2.*ZLN2S(:))+1.246*ZKNG(:)*&
                exp((-4.*ZK+4)/2.*ZLN2S(:)))  

    PVGK(:,3*II+IJ-2)=PVGG(:,II)*&
      (exp((4.*ZK+4.)/2.*ZLN2S(:)) + 1.246*ZKNG(:)* exp((2.*ZK+1.)/2.*ZLN2S(:)))  
  enddo

ENDDO
IF (LHOOK) CALL DR_HOOK('SLT_VELGRAV1D',1,ZHOOK_HANDLE)
 

END SUBROUTINE SLT_VELGRAV1D
