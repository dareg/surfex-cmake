!     #########
       SUBROUTINE CLS_TQ_DIAN( PZONA, PMERA, PTA, PQA,  &
                            PPA, PPS, PHT,              &
                            PCD, PCH, PRI,              &
                            PTS, PHU, PQS, PZ0H, PH,         &
                            PTNM, PQNM, PHUNM,K2M      )  
!     #####################################################################
!
!!****  *CLS_TQ_DIAN*  
!!
!!    PURPOSE
!!    -------
!     modify T2M/RH2M diagnostics by DIAN 2016 final version N2M=3 recommended         
!     modify T2m/RH2M diagnostics unstable case accordingly N2M=4 
!!**  METHOD
!!    ------
!!    If N2M switch in surfex production namelist is set to 3 or 4
!!    this routine is called instead cls_tq.F90 but only for 2m values
!!    For 10m wind the Geleyn scheme in cls_tq.F90 is further used and called
!!    N2M=3: Dian 2016 code final version (lace report on RC-LACE.eu) 
!!    This reduces negative T2m bias in stable nighttime conditions as described
!!    by M. Dian. The tuning constants ZA/ZB are set as Dian 2016 DF, but
!!    can be tuned freely without changing the theory behind (Monin-Obuchov theory)
!!    See also Dian 2016 experiments.
!!
!!    Attention:
!!    In Alpine valleys and high resolution there is already warm bias
!!    in stable conditions (unresolved katabatic effects?),
!!    which can become worse with N2M=3 than N2M=2
!!    N2M=4: The modified Psi function is also applied to unstable case to
!!    reduce bias. This sometimes leads to quite warm patches in T2m and needs
!!    careful testing. The diagnostics has no direct impact on prognostic
!!    variables, but can affect DA results (soil+atmosphere) for 2m observations (modified FG
!!    values).   
!!
!! 
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    USE MODD_CST
!!    USE MODD_GROUND_PAR
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!    F. Meier 11/2020 based on M. Dian LACE report (implementation in old ISBA
!!    ALARO)
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    11/2020 F. Meier based on Dian(2016) and CLS_TQ code
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,     ONLY : XG, XCPD, XKARMAN,XRD,XRV
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODE_THERMOS
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
REAL, DIMENSION(:), INTENT(IN)       :: PZONA  ! atmospheric wind zonal
REAL, DIMENSION(:), INTENT(IN)       :: PMERA  ! atmospheric wind meridional
REAL, DIMENSION(:), INTENT(IN)       :: PTA    ! atmospheric temperature
REAL, DIMENSION(:), INTENT(IN)       :: PQA    ! atmospheric humidity (kg/kg)
REAL, DIMENSION(:), INTENT(IN)       :: PPA    ! atmospheric level pressure
REAL, DIMENSION(:), INTENT(IN)       :: PPS    ! surface pressure
REAL, DIMENSION(:), INTENT(IN)       :: PHT    ! atmospheric level height (temp)
REAL, DIMENSION(:), INTENT(IN)       :: PCD    ! drag coefficient for momentum
REAL, DIMENSION(:), INTENT(IN)       :: PCH    ! drag coefficient for heat
REAL, DIMENSION(:), INTENT(IN)       :: PRI    ! Richardson number
REAL, DIMENSION(:), INTENT(IN)       :: PTS    ! surface temperature
REAL, DIMENSION(:), INTENT(IN)       :: PHU    ! near-surface humidity (%)
REAL, DIMENSION(:), INTENT(IN)       :: PQS    ! near-surface specific humidity (kg/kg)
REAL, DIMENSION(:), INTENT(IN)       :: PZ0H   ! roughness length for heat
REAL, DIMENSION(:), INTENT(IN)       :: PH     ! height of diagnostic
!
REAL, DIMENSION(:), INTENT(INOUT)      :: PTNM   ! temperature at n meters
REAL, DIMENSION(:), INTENT(INOUT)      :: PQNM   ! specific humidity at n meters
REAL, DIMENSION(:), INTENT(INOUT)      :: PHUNM  ! relative humidity at n meters
INTEGER, INTENT(IN) :: K2M ! version of 2m diagnostics
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTA)) :: ZBNH,ZBH,ZRS,ZBD,ZHELP
REAL, DIMENSION(SIZE(PTA)) :: ZLOGS,ZCORS,ZIV,ZFAC
REAL, DIMENSION(SIZE(PTA)) :: ZQSATA, ZHUA
REAL, DIMENSION(SIZE(PTA)) :: ZQSATNM, ZPNM, ZQS, ZQSATS
REAL, DIMENSION(SIZE(PTA)) :: ZSIGMA,ZL,ZBHELP
 CHARACTER(LEN=2)           :: YHUMIDITY
REAL(KIND=JPRB) :: ZEPS1,ZRDIFF,ZA,ZB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('CLS_TQ_DIAN',0,ZHOOK_HANDLE)
PTNM (:) = XUNDEF
PQNM(:)  = XUNDEF
PHUNM(:) = XUNDEF
!
ZBNH   (:) = 0.
ZBH    (:) = 0.
ZBD    (:) = 0.
ZRS    (:) = 0.
ZLOGS  (:) = 0.
ZCORS  (:) = 0.
ZIV    (:) = 0.
ZQSATA (:) = 0.
ZHUA   (:) = 0.
ZQSATS (:) = 0.
ZPNM   (:) = 0.
ZQSATNM(:) = 0.
ZQS    (:) = 0.
ZL     (:) = 0.
ZHELP  (:) = 0.
ZFAC   (:) = 0.
ZBHELP (:) = 0.
ZEPS1=1.0E-8
ZA=1.0
ZB=0.8
!
!*      1.     preparatory calculations
!              ------------------------
!

ZBNH(:)=LOG( PHT(:)/PZ0H(:))
!
ZBH(:)=XKARMAN*SQRT( PCD(:) )/PCH(:) 
ZBD(:)=XKARMAN/SQRT(PCD(:))
!
ZRS(:)=MIN(PH(:)/PHT(:),1.)
!
ZLOGS(:)=LOG(1.+ZRS(:)*(EXP(ZBNH(:)) -1.))
!
!*      2.     Stability effects
!              -----------------
!
! use Dian 2016 in stable conditions N2M=3 or 4
WHERE (PRI(:)>=0.) ! stable case
  ZL(:)=(PMERA(:)**2.+PZONA(:)**2.)*ZBH(:)*PTS(:)/ &
    & MAX(ZBD(:)**2.,ZEPS1)/XG/MAX(ZEPS1,(PTA(:)-PTS(:)))
!ZCORS(:)=ZRS(:)*(ZBNH(:)-ZBH(:)) !Geleyn
  ZCORS(:)=LOG(1.+PH(:)/(ZL(:)/MAX(ZA,ZEPS1)+PHT(:)/(MAX(ZEPS1,EXP(ZBNH(:))-1.))))/ &
    & (LOG(1.+PHT(:)/(ZL(:)/MAX(ZA,ZEPS1)+PHT(:)/(MAX(ZEPS1,EXP(ZBNH(:))-1.))))) &
    & *(ZBNH(:)-ZBH(:))
  WHERE(PH(:)<2.)
    ZCORS(:)=ZRS(:)*(ZBNH(:)-ZBH(:)) !Geleyn
  END WHERE
END WHERE

IF(K2M == 3)THEN ! Dian 2016 use Geleyn scheme in unstable case - recommended

  WHERE (PRI(:)< 0.) ! unstable case
    ZCORS(:)=LOG(1.+ZRS(:)*(EXP(MAX(0.,ZBNH(:)-ZBH(:)))-1.)) ! Geleyn
  END WHERE

ELSE ! use Dian coefficients also in unstable case to reduce negative bias

  WHERE (PRI(:)< 0.) ! unstable case
    ZL(:)=(PMERA(:)**2.+PZONA(:)**2.)*ZBH(:)*PTS(:)/ &
      & (ZBD(:)**2.)/XG/(MIN(-ZEPS1,PTA(:)-PTS(:)))

    ZHELP(:)=(EXP(MAX(0.,ZBNH(:)-ZBH(:)))-1.)
    WHERE(ZHELP(:)==(EXP(ZBNH(:))-1.))
      ZCORS(:)=LOG(1.+ZRS(:)*(EXP(MAX(0.,ZBNH(:)-ZBH(:)))-1.)) ! Geleyn
    ELSEWHERE
      ZBHELP(:)=MAX(-0.999,MIN(0.999,ZB*ZL(:)/PHT(:)*ZHELP(:)* &
        & (EXP(ZBNH(:))-1.)/(ZHELP(:)-(EXP(ZBNH(:))-1.))))    

      ZFAC(:)=ZL(:)/PHT(:)*(EXP(ZBNH(:))-1.)*ZHELP(:)/(ZHELP(:)-(EXP(ZBNH(:))-1.))
      WHERE(ZFAC(:)==ZHELP.or.abs(ZBNH(:))<ZEPS1) ! ZCORS would become something devided by 0.
        ZCORS(:)=LOG(1.+ZRS(:)*(EXP(MAX(0.,ZBNH(:)-ZBH(:)))-1.)) ! Geleyn
      ELSEWHERE
        ZFAC(:)=ZBHELP(:)*PHT(:)/ZL(:)/(EXP(ZBNH(:))-1.)/ZHELP(:)*(ZHELP(:)-(EXP(ZBNH(:))-1.))
        WHERE(ZFAC(:)==1.)
          ZFAC(:)=0.001
        !ZCORS(:)=LOG(1.+ZRS(:)*(EXP(MAX(0.,ZBNH(:)-ZBH(:)))-1.)) ! Geleyn
        ELSEWHERE
          ZFAC(:)=1./(ZFAC(:)-1.)
        END WHERE
        ZCORS(:)=-ZFAC(:)*LOG(1.+(ZBHELP(:)*PH/ZL(:)*(ZHELP(:)-(EXP(ZBNH(:))-1.))-ZRS(:)* & 
        & ZHELP(:)*(EXP(ZBNH(:))-1.))/(ZBHELP(:)*PHT(:)/ZL(:)*(ZHELP(:)/(EXP(ZBNH(:))-1.)- &
        & 1.)-(EXP(ZBNH(:))-1.)))
      END WHERE
    END WHERE
  END WHERE
ENDIF
!
!*      3.     Interpolation of thermodynamical variables
!              ------------------------------------------
!
ZIV=MAX(0.,MIN(1.,(ZLOGS(:)-ZCORS(:))/ZBH(:)))
PTNM(:)=PTS(:)+ZIV(:)*(PTA(:)-PTS(:))

!
!*      4.     Interpolation of relative humidity
!              ----------------------------------
!
!* choice of interpolated variable
!
YHUMIDITY='Q '
!
ZPNM(:) = PPS(:) + PH(:)/PHT(:) * (PPA(:)-PPS(:))
ZQSATNM(:) = QSAT(PTNM(:),ZPNM(:))
!
IF (YHUMIDITY=='Q ') THEN
!        
  PQNM(:)   = PQS(:)+ZIV(:)*(PQA(:)-PQS(:))
  PQNM(:)   = MIN (ZQSATNM(:),PQNM(:)) !must be below saturation
  PHUNM(:)  = PQNM(:) / ZQSATNM(:)
!
ELSE IF (YHUMIDITY=='HU') THEN
!        
  ZQSATA(:) = QSAT(PTA(:),PPA(:))
  ZHUA(:)   = PQA(:) / ZQSATA(:)
  PHUNM(:)  = PHU(:)+ZIV(:)*(ZHUA(:)-PHU(:))
  PQNM(:)   = PHUNM(:) * ZQSATNM(:)
!
END IF
IF (LHOOK) CALL DR_HOOK('CLS_TQ_DIAN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CLS_TQ_DIAN
