!     #########
    SUBROUTINE MIXTL_n(PFSOL,PFNSOL,PSFTEAU,PSFU,PSFV,PSEATEMP)
!     #######################################################################
!
!!****  *MIXTLN (1D MODEL)*  
!!
!!    PURPOSE
!!    -------
!     Oceanic 1D model in TKE Closure scheme
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
!     Gaspar et al, 1990 : A simple eddy kinetic energy model for simulations of 
!     the oceanic vertical mixing : Tests at station Papa and Long-term upper ocean 
!     study site, JGR, 95,C9, 16,179--16,193. 
!!      
!!    AUTHOR
!!    ------
!!     C. Lebeaupin  *Météo-France* 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     02/2008
!!     01/2012 : H. Giordani, P. Peyrille
!!                       Add FLAGS: 
!!                       LREL_TS: relaxation on T, S
!!                       LREL_CUR: damping on current
!!                       FOR LREL_CUR: implicit and explicit codinf is made
!!                       for conveniency, DTREL=Ucur term, DSREL=VCURterm
!!                Corrections:
!!                   coriolis terms in current equation, 
!!                   richardson nb in diapycnal mixing, 
!!                   remove threshold value for mixing tendency
!!    07/2012, P. Le Moigne : CMO1D phasing
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS
USE MODD_OCEAN_CSTS
USE MODD_SEAFLUX_n, ONLY : XSEABATHY
USE MODD_OCEAN_n
USE MODD_OCEAN_GRID_n
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
! Module containing relaxation fields
USE MODD_OCEAN_REL_n , ONLY : XSEAT_REL, XSEAS_REL, XTAU_REL, &
                            LREL_CUR, LREL_TS, &
                            XQCORR, LFLX_CORR, LDIAPYCNAL,&
                            XSEAU_REL, XSEAV_REL
!
USE MODD_SEAFLUX_GRID_n, ONLY :  XLAT
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:)  ,INTENT(IN)       :: PFSOL   ! solar flux (W/m2)
REAL, DIMENSION(:)  ,INTENT(IN)       :: PFNSOL  ! non solar flux (W/m2)
REAL, DIMENSION(:)  ,INTENT(IN)       :: PSFTEAU ! fresh water flux(kg/m2/s)
REAL, DIMENSION(:)  ,INTENT(IN)       :: PSFU    ! zonal stress (Pa)
REAL, DIMENSION(:)  ,INTENT(IN)       :: PSFV    ! meridian stress (Pa)
!
REAL, DIMENSION(:)  ,INTENT(OUT)    :: PSEATEMP! sea surface temperature (K)
!*      0.2    declarations of local variables
!
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ADVT,ADVS !advection horiz. temperature and salinity
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ADVU,ADVV !advection horiz. of current
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ADVE      !advection of turbulent kinetic energy
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZSEAT,ZSEAS,ZSEAE,ZSEAV,ZSEAU
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZSEAT_REL,ZSEAS_REL,ZSEAV_REL,ZSEAU_REL
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZTDTREL  ! Tendancy derived from  relxation (K/s)
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZSDTREL  ! ---- salinity   ---------- (%/s)
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZUDTREL  ! Tendancy term derived from relaxation on U current (m/s2)
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZVDTREL  !  ------------------------------------    V        (m/s2)
!
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZLE,ZKMEL,ZKMELM
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZKMES,ZKMED,ZKMEWM,ZKMEWS
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZA, ZB, ZC, ZA2, ZB2, ZC2, ZYT, ZYS, ZYE        !matrices pour resolution numérique
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZPTH,ZPDY
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZTENDE,ZDIFFV
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZOMT, ZOMS, ZOME !vector for matrix inversion
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZWT, ZWS, ZWE
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZU,ZV,ZT,ZS,ZE
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZZDRHO
!
COMPLEX, DIMENSION(NOCKMIN:NOCKMAX,2) :: ZUC           ! vecteur vent en ecriture complexe
COMPLEX, DIMENSION(NOCKMIN:NOCKMAX)   :: ZAI,ZBI,ZCI     !matrices pour resolution numérique
COMPLEX, DIMENSION(NOCKMIN:NOCKMAX)   :: ZAU,ZBU,ZCU,ZYU !matrices pour resolution numérique
COMPLEX, DIMENSION(NOCKMIN:NOCKMAX)   :: ZOMU,ZWU
!
REAL, DIMENSION(NOCKMIN:NOCKMAX) :: ZDTFSOL
REAL :: ZDTFNSOL
!
REAL :: ZSFU, ZSFV, ZFNSOL, ZFSOL, ZSFTEAU, ZLAT
REAL :: ZSEAHMO, ZSEATEMP
!
REAL :: ZSEUIL,ZEMIN,ZEMAX,ZTEST
!
REAL :: ZF, ZEWS
REAL :: ZALG
REAL :: ZEE,ZPOT,ZXLME,ZXLPE,ZXROD,ZAUX,ZXDL
REAL :: ZDU,ZDV,ZRICH,ZDRHODZ
REAL :: ZT1, ZT2, ZT3, ZS1, ZS2
!
INTEGER :: BINF,BSUP
INTEGER :: J,JJ,JPT,JIN,IKHML
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
!       1.     Initializations
!__________________________________________________________________________
!Take account of an horizontale advection

!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MIXTL_N',0,ZHOOK_HANDLE)
BINF=NOCKMIN+1
BSUP=NOCKMAX
!
ZSEUIL=2.23E-3
ZTEST=0.
IKHML=1
ZEMIN=1.E-3
ZEMAX=0.1
!
ZT1=13.5
ZT2=-0.19494
ZT3=-0.49038E-2
ZS1=32.6
ZS2=0.77475
!
ZALG = XG / XRHOSW
!
!-------------------------------------------------------------------------------
!
!iterations on grid points
DO JPT=1,SIZE(PFSOL)
!
  !simplified variables inside this loop
  ZLAT=XLAT(JPT)
  ZFSOL=PFSOL(JPT)
  ZFNSOL=PFNSOL(JPT)
  ZSFTEAU=PSFTEAU(JPT)
  ZSFU=PSFU(JPT)
  ZSFV=PSFV(JPT)
  ZEWS=SQRT(ZSFU**2+ZSFV**2)/XRHOSW
  ZF=4.*XPI*SIN(ZLAT*XPI/180.)/86400.

  ZSEAHMO=0.
  DO J=BINF-1,BSUP
    ZSEAT(J)=XSEAT(JPT,J)
    ZSEAS(J)=XSEAS(JPT,J)
    ZSEAU(J)=XSEAU(JPT,J)
    ZSEAV(J)=XSEAV(JPT,J)    
    ZSEAE(J)=XSEAE(JPT,J)
    !
    ZSEAU_REL(J)=XSEAU_REL(JPT,J)
    ZSEAV_REL(J)=XSEAV_REL(JPT,J)
    ZSEAT_REL(J)=XSEAT_REL(JPT,J)
    ZSEAS_REL(J)=XSEAS_REL(JPT,J)
    !
    IF (J>=BINF .AND. ZSEAE(J)>=(ZEMIN*SQRT(2.))) ZSEAHMO=ZSEAHMO-XDZ1(J)
  ENDDO

 !precalculation of DRHO
  DO J=BINF-1,BSUP
    ZU(J)=0.
    ZV(J)=0.
    ZT(J)=0.
    ZS(J)=0.
    ZE(J)=0.
    ADVT(J)=0.
    ADVS(J)=0.
    ADVU(J)=0.
    ADVV(J)=0.
    ADVE(J)=0.
    ZZDRHO(J)=(ZSEAT(J)-ZT1)*(ZT2+ZT3*(ZSEAT(J)-ZT1)) + ZS2*(ZSEAS(J)-ZS1)
    ZUDTREL(J)=0.
    ZVDTREL(J)=0.
    ZTDTREL(J)=0.
    ZSDTREL(J)=0.
    ZDTFSOL(J)=0.
  ENDDO
  ZDTFNSOL=0.
!
! Control print
!IF (LREL_CUR)  WRITE(*,*) "WARNING :: Damping on current will be done"
!IF (LREL_TS)  WRITE(*,*) "WARNING :: Relaxation on T, S ocean will be done "
!IF (LDIAPYCNAL)  WRITE(*,*) "WARNING :: diapycnal mixing has been activated"
!IF (LFLX_CORR) WRITE(*,*) "WARNING :: ocean fluxes correctin has been activated"
!
!------------------------------------------------------------------------------
!
!       2.     Oceanic vertical mixing 
!              -----------------------
!!
!!       2.a    Diapycnal mixing
!!              ----------------
!!
  DO J=BSUP-1,BINF,-1
    ZKMES(J) =0.
    ZKMEWM(J)=0.
    ZKMEWS(J)=0.
    ZKMED(J) =0.
    IF ((ZTEST==0.).AND.(ZSEAE(J)>=ZSEUIL)) THEN
      IKHML=J
      ZTEST=1.
    ENDIF
  ENDDO
  DO J=IKHML,BSUP-1
    ZDRHODZ=(ZZDRHO(J)-ZZDRHO(J+1))/XDZ1(J)
    ZDU=ZSEAU(J+1)-ZSEAU(J)
    ZDV=ZSEAV(J+1)-ZSEAV(J)
!! Modif PP - HG : flag diapycnal
    IF (LDIAPYCNAL) THEN
      IF((ZDU*ZDU+ZDV*ZDV).LE.1.E-7) THEN 
        ZRICH = 0.8
      ELSE
        ZRICH = -ZALG*ZDRHODZ/(ZDU**2+ZDV**2)/XK4(J)
      ENDIF    
!coefficient de mélange aux ondes internes
      ZKMEWM(J)=1.E-3
      ZKMEWS(J)=1.E-4
!coefficient de mélange du au cisaillement
      IF(ZRICH>7.) THEN
        ZKMES(J) = 0.
      ELSEIF(ZRICH>=0.) THEN
        ZKMES(J) = 5.E-3*(1.-(ZRICH/0.7)*(ZRICH/0.7))**3
      ELSE
        ZKMES(J) = 5.E-3
      ENDIF
    ENDIF
    !plm si ldiapycnal=F zkmes non modofie et zrich ne sert a rien !
    !plm ELSE
    !plm  ZRICH=-ZALG*ZDRHODZ / (ZDU**2+ZDV**2) / XK4(J)   
  ENDDO
!
!       2.b    Mixing length and coefficient 
!              -----------------------------
!
  DO J=BINF,BSUP-1

    ZEE=ZSEAE(J)**2/ZALG
    ZXROD=(ZZDRHO(J)+ ZZDRHO(J+1))*0.5   

    ZXLME=0.
    ZPOT=0.
    JIN=J
    DO JJ=J+1,BSUP
      ZAUX=ZPOT + XDZ2(JJ)*(ZZDRHO(JJ)-ZXROD)
      IF (ZAUX<ZEE) THEN
        JIN=JJ
        ZPOT=ZAUX
        ZXLME=ZXLME+XDZ2(JJ)
      ENDIF
    ENDDO
    IF (JIN==J) THEN
      ZXLME=2.*(XZ2(J)-XZHOC(J+1))/(ZZDRHO(J+1)-ZXROD)/ZALG
      ZXLME=SQRT(ZXLME)*ZSEAE(J)
    ELSE
      IF (JIN/=BSUP) THEN
        ZXDL=(ZEE-ZPOT)/(ZZDRHO(JIN+1)-ZXROD)
        ZXLME=ZXLME+ZXDL
      ENDIF
    ENDIF

    ZXLPE=0.
    ZPOT=0.
    JIN=J
    DO JJ=J-1,NOCKMIN,-1
      ZAUX=ZPOT + XDZ2(JJ+1)*(ZXROD-ZZDRHO(JJ+1))
      IF (ZAUX<ZEE) THEN
        JIN=JJ
        ZPOT=ZAUX
        ZXLPE=ZXLPE+XDZ2(JJ+1)
      ENDIF
    ENDDO
    IF (JIN==J) THEN
      ZXLPE=2.*(XZ2(J)-XZHOC(J))/(ZZDRHO(J)-ZXROD)/ZALG
      ZXLPE=SQRT(ZXLPE)*ZSEAE(J)
    ELSE 
      IF (JIN/=NOCKMIN) THEN
        ZXDL=- (ZEE-ZPOT)/(ZZDRHO(JIN)-ZXROD)
        ZXLPE=ZXLPE+ZXDL
      ENDIF
    ENDIF

    ZLE(J)=SQRT(ZXLME*ZXLPE)
    ZKMEL(J)=XCKL*ZLE(J)*ZSEAE(J)
  
  ENDDO

  ZLE(BSUP)=ZLE(BSUP-1)
  ZKMEL(BSUP)=XCKL*ZLE(BSUP)*ZSEAE(BSUP)

  !first coef at all levels: because needed at j+1 further
  ZKMELM(BINF)=ZKMEL(BINF)
  DO J=BINF+1,BSUP
    ZKMELM(J)= (ZKMEL(J)+ZKMEL(J-1))/2.
  ENDDO

!--------------------------------------------------------------------------------------------------------
!
!!       2.c    Numerical resolution of evolution equations
!!              -------------------------------------------
!
  DO J=BINF,BSUP
    IF (LREL_CUR) THEN
      ZUDTREL(J) =  - (ZSEAU(J)-ZSEAU_REL(J))  / XTAU_REL 
      ZVDTREL(J) =  - (ZSEAV(J)-ZSEAV_REL(J))  / XTAU_REL 
    ENDIF
    ! flux solaire
    ZDTFSOL(J) = XRAY(J)*ZFSOL/XDZ2(J) 
  ENDDO
!
  IF (LFLX_CORR) THEN
    ! NO relaxation
    ! Barnier correction on surface fluxes          
    ! flux non solaire corrige
    ZTDTREL(BINF) =  XQCORR *(ZSEAT_REL(BINF)-ZSEAT(BINF)) / (XRHOSW*XCPSW)
    ZFNSOL = ZFNSOL + ZTDTREL(BINF)
  ELSEIF (LREL_TS) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RELAXATION IS MADE INSTEAD OF FLUX CORRECTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    DO J=BINF,BSUP
      ! flux non solaire
      ZTDTREL(J) =  - (ZSEAT(J)-ZSEAT_REL(J)) / XTAU_REL
      ZSDTREL(J) =  - (ZSEAS(J)-ZSEAS_REL(J)) / XTAU_REL
    ENDDO
  ENDIF

  ! flux non solaire
  ZDTFNSOL = ZFNSOL/XDZ2(BINF) 
!
!simple coef a, b, c at NOCKMIN+1 level
  ZC(BINF)= XOCEAN_TSTEP*XK2(BINF)*(ZKMEL(BINF)+ZKMES(BINF)+ZKMEWS(BINF)+ZKMED(BINF))              
  ZB(BINF)= 1.-ZC(BINF)
!v coef a, b, c at NOCKMIN+1 level
  ZCI(BINF)=XOCEAN_TSTEP*XK2(BINF)*(ZKMEL(BINF)+ZKMES(BINF)+ZKMEWM(BINF)+ZKMED(BINF))
  ZBI(BINF)= 1.-ZCI(BINF)
!c coef a, b, c at NOCKMIN+1 level
  ZCU(BINF)=ZCI(BINF)*(1.,0.)
  ZBU(BINF)=ZBI(BINF)*(1.,0.)+(0.,1.)*ZF*XGAMA * XOCEAN_TSTEP
!omt (temperature), oms (salinity), omu from ZB and ZBC at NOCKMIN+1 level
  ZOMT(BINF) = 1./ZB(BINF)
  ZOMS(BINF) = 1./ZB(BINF)
  ZOMU(BINF) = 1./ZBU(BINF)
!coef yt, ys and yu at NOCKMIN+1 level
  ZYT(BINF) = ZSEAT(BINF) + XOCEAN_TSTEP * (ZDTFNSOL + ZDTFSOL(BINF) + ADVT(BINF))
  ZYS(BINF) = ZSEAS(BINF) + XOCEAN_TSTEP * ( ZSEAS(BINF)*ZSFTEAU      / XDZ2(BINF) + ADVS(BINF))
! Modif. PP. relaxation towards ref. profile
! If LREL_TS, relaxation on temp, salinity  
  IF (LREL_TS) THEN
    ZYT(BINF) = ZYT(BINF) + XOCEAN_TSTEP * ZTDTREL(BINF)  
    ZYS(BINF) = ZYS(BINF) + XOCEAN_TSTEP * ZSDTREL(BINF) 
  ENDIF
  ZUC(BINF,1)=ZSEAU(BINF)*(1.,0.)+ZSEAV(BINF)*(0.,1.)
  ZYU(BINF) = ZUC(BINF,1) + XOCEAN_TSTEP * ( ZUC(BINF,1) * ZF*(1-XGAMA) * (0.,-1.) - &
         (ZSFU*(1.,0.) + ZSFV*(0.,1.))/XDZ2(BINF)/XRHOSW + ADVU(BINF)*(1.,0.) + ADVV(BINF)*(0.,1.) )
! damping on current if LREL_CUR=T  in explicit scheme         
  IF (LREL_CUR) ZYU(BINF) = ZYU(BINF) + XOCEAN_TSTEP*(ZUDTREL(BINF)*(1.,0.)+ ZVDTREL(BINF)*(0.,1.))
!wt, ws, wu from OMT, OMS, OMU, ZY
  ZWT(BINF) = ZOMT(BINF)*ZYT(BINF)
  ZWS(BINF) = ZOMS(BINF)*ZYS(BINF)
  ZWU(BINF) = ZOMU(BINF)*ZYU(BINF)
!*Evolution of the square root of TKE computation of density gradient
  ZDRHODZ=(ZZDRHO(BINF)-ZZDRHO(BINF+1))/XDZ1(BINF)
  ZDU=ZSEAU(BINF+1)-ZSEAU(BINF)
  ZDV=ZSEAV(BINF+1)-ZSEAV(BINF)
  ZPTH(BINF) = XCKL*ZLE(BINF) * ZALG*ZDRHODZ
  ZPDY(BINF) = XCKL*ZLE(BINF) * XK4(BINF)*(ZDU**2+ZDV**2)
!coef a2, b2, c2 at NOCKMIN and BSUP levels
  ZC2(BINF)=XOCEAN_TSTEP*XK3(BINF)*ZKMELM(BINF+1)
  ZB2(BINF)=1. - ZC2(BINF) + XOCEAN_TSTEP*ZSEAE(BINF) * 1./ZLE(BINF)/XZCE                            
!OME at NOCKMIN from B2 
  ZOME(BINF) = 1./ZB2(BINF)
!Y3 coef at NOCKMIN
  ZYE(BINF)=ZSEAE(BINF) + XOCEAN_TSTEP * &
             (0.5 *ZSEAE(BINF)**2/ZLE(BINF)/XZCE &
             + ZEWS/XDZ1(BINF) + ADVE(BINF)) + ZPTH(BINF) + ZPDY(BINF) 
!WE from OME and Y3
  ZWE(BINF) = ZOME(BINF)*ZYE(BINF)

!!
!loop on levels
  DO J=BINF+1,BSUP-1
    ZC(J)     = XOCEAN_TSTEP * XK2(J) * (ZKMEL(J)   + ZKMES(J)   + ZKMEWS(J)   + ZKMED(J)) 
    ZA(J)     = XOCEAN_TSTEP * XK1(J) * (ZKMEL(J-1) + ZKMES(J-1) + ZKMEWS(J-1) + ZKMED(J-1)) 
    ZB(J)     = 1. - ZA(J) - ZC(J)

    ZCI(J)    = XOCEAN_TSTEP * XK2(J) * (ZKMEL(J)   + ZKMES(J)   + ZKMEWM(J)   +ZKMED(J)) 
    ZAI(J)    = XOCEAN_TSTEP * XK1(J) * (ZKMEL(J-1) + ZKMES(J-1) + ZKMEWM(J-1) + ZKMED(J-1)) 
    ZBI(J)    = 1.-ZAI(J)-ZCI(J)

    ZCU(J)    = ZCI(J)*(1.,0.)
    ZAU(J)    = ZAI(J)*(1.,0.)    
    ZBU(J)    = ZBI(J)*(1.,0.)+(0.,1.) * ZF*XGAMA * XOCEAN_TSTEP
    
    ZOMT(J)   = 1./(ZB(J)  - ZA(J)  * ZOMT(J-1) * ZC(J-1)) 
    ZOMU(J)   = 1./(ZBU(J) - ZAU(J) * ZOMU(J-1) * ZCU(J-1)) 
    ZOMS(J)   = 1./(ZB(J)  - ZA(J)  * ZOMS(J-1) * ZC(J-1)) 

    ZYT(J)    = ZSEAT(J) + XOCEAN_TSTEP * (ZDTFSOL(J) + ADVT(J))
    ZYS(J)    = ZSEAS(J) + XOCEAN_TSTEP * (                         ADVS(J))
    IF (LREL_TS) THEN
      ZYT(J) = ZYT(J) + XOCEAN_TSTEP * ZTDTREL(J)
      ZYS(J) = ZYS(J) + XOCEAN_TSTEP * ZSDTREL(J)
    ENDIF
    ZUC(J,1)  = ZSEAU(J)*(1.,0.) + ZSEAV(J)*(0.,1.)
    ZYU(J)    = ZUC(J,1) + XOCEAN_TSTEP * (ZUC(J,1)*ZF*(1.-XGAMA)*(0.,-1.) + ADVU(J)*(1.,0.) + ADVV(J)*(0.,1.))
    ! damping on current if LREL_CUR=T explicit
    ! Pareil V,UREF da,s relawxation
    IF (LREL_CUR) ZYU(J) = ZYU(J) + XOCEAN_TSTEP * (ZUDTREL(J)*(1.,0.)+ZVDTREL(J)*(0.,1.))  
    ZWT(J)    = ZOMT(J) * (ZYT(J)- ZA(J) *ZWT(J-1))
    ZWS(J)    = ZOMS(J) * (ZYS(J)- ZA(J) *ZWS(J-1))
    ZWU(J)    = ZOMU(J) * (ZYU(J)- ZAU(J)*ZWU(J-1))

    ZDRHODZ   = (ZZDRHO(J)-ZZDRHO(J+1))/XDZ1(J)

    ZDU       = ZSEAU(J+1)-ZSEAU(J)
    ZDV       = ZSEAV(J+1)-ZSEAV(J)

    ZPTH(J)   = XCKL*ZLE(J)*ZALG*ZDRHODZ
    ZPDY(J)   = XCKL*ZLE(J)*XK4(J)*(ZDU**2+ZDV**2)

    ZC2(J)    = XOCEAN_TSTEP*XK3(J)*ZKMELM(J+1)
    ZA2(J)    = XOCEAN_TSTEP*XK2(J)*ZKMELM(J)
    ZB2(J)    = 1. - ZA2(J) - ZC2(J) + XOCEAN_TSTEP * ZSEAE(J) * 1./ZLE(J)/XZCE

    ZOME(J)   = 1. / (ZB2(J) - ZA2(J)*ZOME(J-1)*ZC2(J-1)) 
    ZYE(J)    = ZSEAE(J) + XOCEAN_TSTEP * (0.5 * ZSEAE(J)**2/ZLE(J)/XZCE + ADVE(J)) + ZPTH(J) + ZPDY(J)
    ZWE(J)    = ZOME(J) * (ZYE(J) - ZA2(J)*ZWE(J-1))
  ENDDO
!
!what leaves at BSUP level
  ZA(BSUP)  = 0.
  ZB(BSUP)  = 1.-ZA(BSUP)
  ZAI(BSUP)= 0.
  ZBI(BSUP)  = 1.-ZAI(BSUP)
  ZAU(BSUP)  = ZAI(BSUP)*(1.,0.)
  ZBU(BSUP)  = ZBI(BSUP)  *(1.,0.)+XOCEAN_TSTEP * ZF*XGAMA*(0.,1.)
  ZOMT(BSUP)   = 1./(ZB(BSUP)-ZA(BSUP)*ZOMT(BSUP-1)*ZC(BSUP-1))
  ZOMU(BSUP)=1./(ZBU(BSUP)-ZAU(BSUP)*ZOMU(BSUP-1)*ZCU(BSUP-1))
  ZOMS(BSUP)=1./(ZB(BSUP)-ZA(BSUP)*ZOMS(BSUP-1)*ZC(BSUP-1))
  ZYT(BSUP)=ZSEAT(BSUP) + XOCEAN_TSTEP * ( ZDTFSOL(BSUP) + ADVT(BSUP) ) 
  ZYS(BSUP)=ZSEAS(BSUP) + XOCEAN_TSTEP * ( ADVS(BSUP) )   
  IF (LREL_TS) THEN
    ZYT(BSUP) = ZYT(BSUP) + XOCEAN_TSTEP * ZTDTREL(BSUP)  
    ZYS(BSUP) = ZYS(BSUP) + XOCEAN_TSTEP * ZSDTREL(BSUP)
  ENDIF
  ZUC(BSUP,1)=ZSEAU(BSUP)*(1.,0.)+ZSEAV(BSUP)*(0.,1.)
  ZYU(BSUP)=ZUC(BSUP,1) + XOCEAN_TSTEP * ( ZUC(BSUP,1) * ZF*(1.-XGAMA)*(0.,-1.) + ADVU(BSUP)*(1.,0.) + ADVV(BSUP)*(0.,1.) )
  ! damping on current if LREL_CUR=T explicit
  ! Pareil V,UREF da,s relawxation
  IF (LREL_CUR) ZYU(BSUP) = ZYU(BSUP) + XOCEAN_TSTEP * (ZUDTREL(BSUP)*(1.,0.)+ZVDTREL(BSUP)*(0.,1.))
  ZWT(BSUP)=ZOMT(BSUP)*(ZYT(BSUP)-ZA(BSUP)*ZWT(BSUP-1))
  ZWS(BSUP)=ZOMS(BSUP)*(ZYS(BSUP)-ZA(BSUP)*ZWS(BSUP-1))
  ZWU(BSUP)=ZOMU(BSUP)*(ZYU(BSUP)-ZAU(BSUP)*ZWU(BSUP-1))
  ZPTH(BSUP)=ZPTH(BSUP-1)
  ZPDY(BSUP)=ZPDY(BSUP-1)
  ZA2(BSUP)= 0.
  ZB2(BSUP)=1. - ZA2(BSUP) + XOCEAN_TSTEP * ZSEAE(BSUP)/ZLE(BSUP)/XZCE
  ZOME(BSUP)=1./(ZB2(BSUP)-ZA2(BSUP)*ZOME(BSUP-1)*ZC2(BSUP-1))
  ZYE(BSUP)    = ZSEAE(BSUP) + XOCEAN_TSTEP * (0.5 * ZSEAE(BSUP)**2/ZLE(BSUP)/XZCE + ADVE(BSUP)) + ZPTH(BSUP)+ZPDY(BSUP) 
  ZWE(BSUP)=ZOME(BSUP)*(ZYE(BSUP)-ZA(BSUP)*ZWE(BSUP-1))
!
!---------------------------------------------------------------------------------------------------
!
  ZUC(BSUP,2)=ZWU(BSUP)
  ZT(BSUP)=ZWT(BSUP)
  ZS(BSUP)=ZWS(BSUP)
  ZE(BSUP)=ZWE(BSUP)
  DO J=BSUP-1,BINF,-1
    ZUC(J,2)=ZWU(J)-ZOMU(J)*ZCU(J)*ZUC(J+1,2)
    ZT(J)=ZWT(J)-ZOMT(J)*ZC(J)*ZT(J+1)
    ZS(J)=ZWS(J)-ZOMS(J)*ZC(J)*ZS(J+1)
    ZE(J)=ZWE(J)-ZOME(J)*ZC2(J)*ZE(J+1)
  ENDDO
!
!
  DO J=BINF,BSUP
    ZU(J)  = REAL(ZUC(J,2))
    ZV(J)  = AIMAG(ZUC(J,2))
    ZE(J)  = MAX(ZEMIN,ZE(J))
  ! Transformation to preserve E <EMAX; secure if mixt crash
    ZE(J)  = MIN(ZE(J),ZEMAX)
  !bilan TKE
    ZTENDE(J)=(ZE(J)*ZE(J)-ZSEAE(J)**2)/XOCEAN_TSTEP
    ZDIFFV(J)=ZTENDE(J) - ZSEAE(J)*(ZPDY(J) + ZPTH(J))
    !
    ZSEAU(J)  = ZU(J)
    ZSEAV(J)  = ZV(J)
    ZSEAT(J)  = ZT(J)
    ZSEAS(J)  = ZS(J)
    ZSEAE(J)  = ZE(J)
  ENDDO
!
!------------------------------------------------------------------------------
!!       3.     New oceanic profiles
!!              --------------------
!!
  IF (LPROGSST) XSEATEND(JPT)=(ZT(BINF)-ZSEAT(BINF))/XOCEAN_TSTEP
  ZSEAU(NOCKMIN)  = ZU(BINF)
  ZSEAV(NOCKMIN)  = ZV(BINF)
  ZSEAT(NOCKMIN)  = ZT(BINF)
  ZSEAS(NOCKMIN)  = ZS(BINF)
  ZSEAE(NOCKMIN)  = ZE(BINF)

  !bathymetrie
  DO J=BINF,BSUP
    IF (XSEABATH(JPT,J)==0.) THEN
      ZSEAU(J)  = ZSEAU(J-1)
      ZSEAV(J)  = ZSEAV(J-1)
      ZSEAT(J)  = ZSEAT(J-1) 
      ZSEAS(J)  = ZSEAS(J-1)
      ZSEAE(J)  = ZSEAE(J-1)
    ENDIF
  ENDDO
!
!SST diagnosticed with 1D oceanic model
  ZSEATEMP=ZSEAT(BINF)+XTT
!
  XSEAHMO(JPT)=ZSEAHMO
  PSEATEMP(JPT)=ZSEATEMP

  DO J=BINF,BSUP
    XLE(JPT,J)=ZLE(J)
    XLK(JPT,J)=ZLE(J)
    XKMEL(JPT,J)=ZKMEL(J)
    XKMELM(JPT,J)=ZKMELM(J)
  ENDDO

  DO J=BINF-1,BSUP
    XSEAT(JPT,J)=ZSEAT(J)
    XSEAS(JPT,J)=ZSEAS(J)
    XSEAU(JPT,J)=ZSEAU(J)
    XSEAV(JPT,J)=ZSEAV(J)
    XSEAE(JPT,J)=ZSEAE(J)
    XDTFSOL(JPT,J) = ZDTFSOL(J)
  ENDDO

  XDTFNSOL(JPT) = ZDTFNSOL

ENDDO
!  
IF (LHOOK) CALL DR_HOOK('MIXTL_N',1,ZHOOK_HANDLE)
!
!!-------------------------------------------------------------------------------
!!
!ENDDO !end of iterations on sea surfex grid points
!!------------------------------------------------------------------------------
!
END SUBROUTINE MIXTL_n
