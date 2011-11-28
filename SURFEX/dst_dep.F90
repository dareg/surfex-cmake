!     #########
     SUBROUTINE DST_DEP (PSVT, PFSVT,  PUSTAR, &
                      PRESA, PTA, PRHODREF)  
!###########################################################
  !
  !!                   
  !!                       
  !!
  !!    PURPOSE
  !!    -------
  !!      
  !!    Compute dry deposition velocity for dust species 
  !!
  !!    AUTHOR
  !!    ------
  !!      P.Tulet      * CNRM *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original      20/02/05 
  !!
  !-------------------------------------------------------------------------------
  !
  !*       0.    DECLARATIONS
  !              ------------
  !
  USE MODE_DST_SURF
  USE MODI_DST_VELGRAV1D
  !
!
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
  USE PARKIND1  ,ONLY : JPRB
!
  IMPLICIT NONE
  !
  !*       0.1   Declarations of dummy arguments :
  !
       REAL, DIMENSION(:,:),   INTENT(IN)    :: PSVT       ! friction velocity
       REAL, DIMENSION(:,:),   INTENT(INOUT) :: PFSVT      ! flux
       REAL, DIMENSION(:),     INTENT(IN)    :: PUSTAR       ! friction velocity
       REAL, DIMENSION(:),     INTENT(IN)    :: PRESA        ! aerodynamical resistance
       REAL, DIMENSION(:),     INTENT(IN)    :: PTA          ! ait temperature
       REAL, DIMENSION(:),     INTENT(IN)    :: PRHODREF     ! air density
  !
  !
  !*       0.2   Declarations of local variables :
  !
  REAL , DIMENSION(SIZE(PSVT,1), JPMODE_DST*3) :: ZRD ! surface  resistance
  REAL , DIMENSION(SIZE(PSVT,1), JPMODE_DST*3) :: ZVD ! [m/s] dry deposition velocity 
  REAL , DIMENSION(SIZE(PSVT,1), JPMODE_DST*3) :: Stn ! Stockes number
  REAL , DIMENSION(SIZE(PSVT,1), JPMODE_DST*3) :: Sc  ! Schmidt number
  REAL , DIMENSION(SIZE(PSVT,1))      :: ZUSTAR, ZRESA
  REAL , DIMENSION(SIZE(PSVT,1), JPMODE_DST*3):: ZWORK
  REAL, DIMENSION(SIZE(PSVT,1)) :: ZNU
  REAL, DIMENSION(SIZE(PSVT,1),JPMODE_DST*3)     :: Dg,zvs,zvsg, zdsg
  REAL, DIMENSION(SIZE(PSVT,1))        :: ZMU
  REAL, DIMENSION(SIZE(PSVT,1),JPMODE_DST*3)   :: ZVGK, ZDPK
  REAL, DIMENSION(SIZE(PSVT,1),JPMODE_DST) :: ZSIG, ZRG
  REAL, DIMENSION(SIZE(PSVT,1),JPMODE_DST) :: ZVG, ZDG
  REAL, DIMENSION(SIZE(PSVT,1), SIZE(PSVT,2)) :: ZSVT
  INTEGER :: JJ, JSV, JN
  REAL :: ZDEN2MOL, ZG, ZMD, ZAVOGADRO,ZMI, ZFAC, ZPI, ZRHOP
  INTEGER, DIMENSION(JPMODE_DST*3) :: NM0, NM3, NM6
  REAL(KIND=JPRB) :: ZHOOK_HANDLE


  !
  !============================================================================
  !
  !            Primilary
  !            ---------
  !Default values
  !--------------
! Cf Ackermann (all to black carbon except water)
IF (LHOOK) CALL DR_HOOK('DST_DEP',0,ZHOOK_HANDLE)
ZRHOP = XDENSITY_DUST
ZPI   = 2.*ASIN(1.)
ZFAC  = (4./3.)*ZPI*ZRHOP*1.e-9
ZMI   = XMOLARWEIGHT_DUST*1E3 ! molecular mass in g/mol
ZAVOGADRO = 6.0221367E+23
ZMD = 28.9644E-3
ZDEN2MOL = 1E-6 * XAVOGADRO / XMD
ZG       = 9.80665
ZMU(:) = 0.
ZVGK(:,:) = 0.
ZVG (:,:) = 0.
ZDPK(:,:) = 0.
ZUSTAR(:) = MAX(PUSTAR(:), 1.E-20)
ZRESA(:) = MIN(MAX(PRESA(:),   1.E-20), 9999.)

!Set counter for number, M3 and M6
IF (LVARSIG) THEN
DO JN=1,JPMODE_DST
    NM0(JN) = 1+(JN-1)*3
    NM3(JN) = 2+(JN-1)*3
    NM6(JN) = 3+(JN-1)*3
END DO
ELSE IF (LRGFIX_DST) THEN
DO JN=1,JPMODE_DST
    NM3(JN) = JN
END DO
ELSE
DO JN=1,JPMODE_DST
    NM0(JN) = 1+(JN-1)*2
    NM3(JN) = 2+(JN-1)*2
END DO
END IF

! Save scalars in local array
ZSVT(:,:)=MAX(PSVT(:,:),1E-60)

CALL DUSTMOMENT2SIZE(ZSVT, PRHODREF, PSIG1D=ZSIG, PRG1D=ZRG ) !don't need ZN, PN1D=ZN)
!DO JN=1,JPMODE_DST
!print*,'DST_DEP ZRG =', MINVAL(ZRG(:,JN)), MAXVAL(ZRG(:,JN))
!print*,'DST_DEP ZSIGG =', MINVAL(ZSIG(:,JN)), MAXVAL(ZSIG(:,JN))
!ENDDO

CALL DST_VELGRAV1D(ZSIG, ZRG, PTA, PRHODREF, ZRHOP, ZMU, ZVGK,ZDPK, ZVG, ZDG)

Dg(:,:)  = MAX(ZDPK(:,:),1.E-40)
zvs(:,:) = MAX(ZVGK(:,:),1.E-20)
ZNU(:)   = ZMU(:)/PRHODREF(:)
DO JN=1,JPMODE_DST
DO JJ= 0,2
zvsg(:,3*JN+JJ-2) =  MAX(ZVG(:,JN),1.E-20)
zdsg(:,3*JN+JJ-2) =  MAX(ZDG(:,JN),1.E-40)
END DO
END DO

!     compute Schmidt number
!     ----------------------
  DO  JN=1,JPMODE_DST*3
     !Sc(:,JN)= ZNU(:)/Dg(:,JN)
      Sc(:,JN)= ZNU(:)/zdsg(:,JN)
  END DO

Stn(:,:) =0.
  DO  JN=1,JPMODE_DST*3
   ZVD(:,JN) = 0.
   ZWORK(:,JN) = 0.
   !Stoke's number, Seinfeld & Pandis, pp 965
   Stn(:,JN)= zvsg(:,JN)*ZUSTAR(:)**2/(ZG*ZNU(:))
   Stn(:,JN)= MAX(Stn(:,JN), 0.05)

   !Get nominator of equation 19.18 Seinfeld & Pandis==> 1/rd
   ZRD(:,JN) = ZUSTAR(:) * (Sc(:,JN)**(-2./3.)+ &
          10**(-3./Stn(:,JN)))   
   !Limit to reasonable values
   !ZRD(:,JN) = MAX(ZRD(:,JN),1.E-10)

   !Get rd
   ZRD(:,JN) = 1. / ZRD(:,JN) 

   !Get ra + rd + ra*rd*vg which is equal to 
   !getting nominator of equation 19.7 Seinfeld & Pandis 
   ZWORK(:,JN)= ZRESA(:) + ZRD(:,JN) +&
       ZRESA(:)*ZRD(:,JN)*zvs(:,JN)  

   !Limit to reasonable values
   ZWORK(:,JN)= MAX(ZWORK(:,JN), 1.E-10)

   !Get the total dry dep velocity (Seinfeld & Pandis, eqn 19.7)
   IF (CVERMOD=='CMDVER') THEN 
     ZWORK(:,JN)= 10.0 * zvs(:,JN)  &  !Gravitation term 
                  + 1./ZWORK(:,JN)     !turbulence and surface resistance term
   else    
!     ZWORK(:,JN)=zvs(:,JN)  &  !Gravitation term 
!         + 1./ZWORK(:,JN)     !turbulence and surface resistance term
      ! The gravitation term as been computed by MesoNH (see sedim_dust.f90)
     ZWORK(:,JN) =  1./ZWORK(:,JN) ! turbulence and surface resistance term
   END IF    

   !         deposition velocity for each cover type
   !         ----------------------------------------
   ZVD(:,JN)=ZWORK(:,JN)
 END DO

! Only M3 flux (mass) has been used ; flux for over moment has been made after
  DO JSV=1,JPMODE_DST
     PFSVT(:,NM3(JSV)) =  PFSVT(:,NM3(JSV)) - PSVT(:,NM3(JSV))  * ZVD(:,2+(JSV-1)*3)
  ENDDO
IF (LHOOK) CALL DR_HOOK('DST_DEP',1,ZHOOK_HANDLE)


!
!							       	       
!---------------------------------------------------------------------
!
END SUBROUTINE DST_DEP
