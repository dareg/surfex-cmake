!option! -O nomove
!****---------------------------------------------------------------------------
!****   CACSTS : INITIALISE LES CHAMPS DE SURFACE
!****   ------
!****  Auteurs :   CB 01/91, BU 05/92, VC 05/93, DG 03/94, PA 09/95, DG 05/96
!****      Mod : E. Bazile 01/97 Soustraction du biais moyen de temperature
!****                            et/ou d'humidite pour les increments utilises
!****                            pour l'analyse de l'eau du sol
!****      Mod : D. Giard  03/99 ACSOL -> ACSOLW
!****            E. Bazile , F. Bouyssel : Remplacement du logique LLLACW par
!****                une fonction continue ZDACW (LSOLV).
!****      Mod : F. Taillefer 09/02 : mise a jour constantes surface selon SST
!****      Mod : F. Bouyssel 02/04 : Seuil utilisant l'angle zenithal solaire
!****      Mod : E. Bazile 01/2007 : Parametre pour la correction PSNS et WPI
!        M.Hamrud      01-Jul-2006  Revised surface fields
!        A.Trojakova   27-Jun-2007 bugfixing ZV10M (surface pointers)
!        F. Bouyssel    27-Mar-2011  Use of REPS2 instead of REPS3 for ZNEI
!****---------------------------------------------------------------------------
!
SUBROUTINE OI_CACSTS(KNBPT,PT2INC,PH2INC,PWGINC,PWS_O,                      &
                     KDAT,KSSSSS,                                           &
                     PTP,PWP,PTL,PSNS,PTS,PWS,                              &
                     PTCLS,PHCLS,PUCLS,PVCLS,PSSTC,PWPINC1,PWPINC2,PWPINC3, &
                     PT2MBIAS,PH2MBIAS,                                     &
                     PRRCL,PRRSL,PRRCN,PRRSN,PATMNEB,PEVAP,PEVAPTR,         &
                     PITM,PVEG,PALBF,PEMISF,PZ0F,                           &
                     PIVEG,PARG,PD2,PSAB,PLAI,PRSMIN,PZ0H,                  &
                     PTSC,PTPC,PWSC,PWPC,PSNC,PGELAT,PGELAM,PGEMU)  
!
!****---------------------------------------------------------------------------
!**  BUT : INITIALISE LES CHAMPS DE SURFACE PROGNOSTIQUES
!**  ---
!**  SEQUENCE D'APPEL :
!**  ----------------
!**        CALL CACSTS(....)

!**  ARGUMENTS D'ENTREE :
!**  ------------------        
!**                               
!**        - EXPLICITE - 
!**                      KNBPT  : nombre reel de points traites
!**                      PT2INC : increment d'analyse de T2m
!**                      PH2INC : increment d'analyse de Hu2m
!**                      PSP_SB,PSP_SG,PSP_RR,PSD_VF,PSD_VV,PSD_VX,PSP_CI,PSP_X2    : 
!**                      buffer des champs pdg de l'analyse
!**                      PGELAM, PGELAT, PGEMU : coordonnees geographiques

!**  ARGUMENTS DE SORTIE : 
!**  -------------------
!**  EXTERNES : CAVEGI (FCTVEG) - ACSOLW - TSL
!**  --------

!**  ALGORITHME : - INITIALISE LA TEMPERATURE DE SURFACE.
!**  ----------   - INITIALISE LA TEMPERATURE PROFONDE.
!**               - INITIALISE LE RESERVOIR DE SURFACE.
!**               - INITIALISE LE RESERVOIR PROFOND.
!**               - CORRIGE LA QUANTITE DE NEIGE.
!***-----------------------------------------------------------------
!
USE MODD_CSTS,       ONLY : XG, XTT, XRHOLW, XDAY
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_ASSIM
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_OI_CAVEGI
USE MODI_OI_ACSOLW
USE MODI_OI_JACOBIANS
USE MODI_OI_TSL
USE MODI_OI_FCTVEG
USE MODI_OI_KALMAN_GAIN
!
IMPLICIT NONE
!
INTEGER,INTENT(IN)    :: KNBPT, KDAT, KSSSSS
!
REAL   ,INTENT(IN)    :: PT2INC(KNBPT) 
REAL   ,INTENT(IN)    :: PH2INC(KNBPT) 
REAL   ,INTENT(IN)    :: PWGINC(KNBPT)
REAL   ,INTENT(IN)    :: PWS_O(KNBPT)
REAL   ,INTENT(INOUT) :: PTP(KNBPT)
REAL   ,INTENT(INOUT) :: PWP(KNBPT)
REAL   ,INTENT(INOUT) :: PTL(KNBPT)
REAL   ,INTENT(INOUT) :: PSNS(KNBPT)
REAL   ,INTENT(INOUT) :: PTS(KNBPT)
REAL   ,INTENT(INOUT) :: PWS(KNBPT)
REAL   ,INTENT(INOUT) :: PTCLS(KNBPT) 
REAL   ,INTENT(INOUT) :: PHCLS(KNBPT) 
REAL   ,INTENT(INOUT) :: PUCLS(KNBPT)
REAL   ,INTENT(INOUT) :: PVCLS(KNBPT)
REAL   ,INTENT(INOUT) :: PSSTC(KNBPT)
REAL   ,INTENT(INOUT) :: PWPINC1(KNBPT)
REAL   ,INTENT(INOUT) :: PWPINC2(KNBPT)
REAL   ,INTENT(INOUT) :: PWPINC3(KNBPT)
REAL   ,INTENT(INOUT) :: PT2MBIAS(KNBPT)
REAL   ,INTENT(INOUT) :: PH2MBIAS(KNBPT)
REAL   ,INTENT(IN)    :: PRRCL(KNBPT)
REAL   ,INTENT(IN)    :: PRRSL(KNBPT)
REAL   ,INTENT(IN)    :: PRRCN(KNBPT)
REAL   ,INTENT(IN)    :: PRRSN(KNBPT)
REAL   ,INTENT(IN)    :: PATMNEB(KNBPT)
REAL   ,INTENT(IN)    :: PEVAP(KNBPT)
REAL   ,INTENT(IN)    :: PEVAPTR(KNBPT)
REAL   ,INTENT(IN)    :: PITM(KNBPT) 
REAL   ,INTENT(IN) :: PVEG(KNBPT) 
REAL   ,INTENT(INOUT) :: PALBF(KNBPT)
REAL   ,INTENT(INOUT) :: PEMISF(KNBPT)
REAL   ,INTENT(INOUT) :: PZ0F(KNBPT)
REAL   ,INTENT(INOUT) :: PIVEG(KNBPT)
REAL   ,INTENT(INOUT) :: PARG(KNBPT)
REAL   ,INTENT(INOUT) :: PD2(KNBPT)
REAL   ,INTENT(INOUT) :: PSAB(KNBPT) 
REAL   ,INTENT(INOUT) :: PLAI(KNBPT)
REAL   ,INTENT(INOUT) :: PRSMIN(KNBPT)
REAL   ,INTENT(INOUT) :: PZ0H(KNBPT)
REAL   ,INTENT(IN)    :: PTSC(KNBPT)
REAL   ,INTENT(IN)    :: PTPC(KNBPT)
REAL   ,INTENT(IN)    :: PWSC(KNBPT)
REAL   ,INTENT(IN)    :: PWPC(KNBPT)
REAL   ,INTENT(IN)    :: PSNC(KNBPT) 
REAL   ,INTENT(IN)    :: PGELAT(KNBPT) 
REAL   ,INTENT(IN)    :: PGELAM(KNBPT) 
REAL   ,INTENT(IN)    :: PGEMU(KNBPT)
!
REAL, DIMENSION(24) :: ZVGAT1,ZVGAT2,ZVGAT3,ZVGBT1,ZVGBT2,ZVGBT3,ZVGCT1,ZVGCT2
REAL, DIMENSION(24) :: ZVGAH1,ZVGAH2,ZVGAH3,ZVGBH1,ZVGBH2,ZVGBH3,ZVGCH1,ZVGCH2
REAL, DIMENSION(24) :: ZSIGT2MP,ZSIGHP2
!
REAL, DIMENSION(KNBPT) :: ZIVEG
REAL, DIMENSION(KNBPT) :: ZWFC, ZWPMX, ZWSAT, ZWSMX, ZWWILT
REAL, DIMENSION(KNBPT) :: ZDWG_DWG, ZDWG_DW2
!
REAL :: ZECHGU, ZNEI, ZCLI, ZPD, ZCLIMCA
REAL :: ZTSC, ZTPC, ZWSC, ZWPC, ZSNC
REAL :: ZV10M, ZPRECIP, ZWPI, ZDACW, ZDACW2, ZMU0, ZMU0M
!
REAL :: ZVGST,ZVGSH,ZVGPT1,ZVGPH1,ZVGPT2,ZVGPH2,ZG1,ZG2,ZG3,ZG4
!
REAL :: ZZT, ZZH, ZLAISRS, ZTEFF, ZHEFF
REAL :: ZCWPH, ZCWPT, ZT2D, ZH2D
REAL :: ZWSD, ZWPD, ZWPDX
!
REAL :: ZWSA, ZWSMIN, ZWPA, ZWPMIN
REAL :: ZGEL, ZSNA, ZMSN, ZK1, ZK2
!
INTEGER :: IH, JROF
LOGICAL :: GSGOBS
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!--------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('OI_CACSTS',0,ZHOOK_HANDLE)
!
ZECHGU = REAL(NECHGU) * 3600.
!
!**  1. Initialisation des polynomes bruts et des champs de reference.    
!
 CALL OI_CAVEGI(ZVGAT1,ZVGAT2,ZVGAT3,ZVGBT1,ZVGBT2,ZVGBT3,ZVGCT1,ZVGCT2, &
                ZVGAH1,ZVGAH2,ZVGAH3,ZVGBH1,ZVGBH2,ZVGBH3,ZVGCH1,ZVGCH2, &
                ZSIGT2MP,ZSIGHP2,GSGOBS) 
!
!
!*   1.2  Initialisation des variables intermediaires
!
DO JROF = 1,KNBPT
  ZIVEG(JROF) = ANINT(PIVEG(JROF))
ENDDO
!
 CALL OI_ACSOLW (1, KNBPT,                        &
                 PARG, PD2, PWS, ZIVEG, PSAB,     &
                 LLDHMT,                          &
                 ZWFC, ZWPMX, ZWSAT, ZWSMX, ZWWILT)  
!
! Analytical Jacobians for WG assimilation
!
 CALL OI_JACOBIANS (KNBPT,PWS_O,PSAB,PARG,PD2,PWP,ZDWG_DWG,ZDWG_DW2) 
!
!**---------------------------------------------------------------------
!**  - 2 - Calcul des champs analyses.

DO JROF = 1,KNBPT

! analyse de surface ou rappel vers la climatologie, sur terre
!    IF (PITM(JROF) > 0.5.AND.(RCLIMCA >= 0.0).AND.PWS(JROF)/=XUNDEF) THEN
  IF ( PWS(JROF)/=XUNDEF ) THEN
!
! stockage des champs prevus
    ZNEI = MAX(0.0,PSNS(JROF)/(PSNS(JROF)+XWCRIN))
! mise a jour des champs climatologiques    
    ZCLI = XRCLIMCA /(1.0+XRCLIMN*ZNEI)

    IF ( .NOT. LCLIM ) THEN
      ZTSC = PTS (JROF)
      ZTPC = PTP (JROF)
      ZWSC = PWS (JROF)
      ZWPC = PWP (JROF)
      ZSNC = PSNS(JROF)
    ELSE
      ZTSC = PTSC(JROF)
      ZTPC = PTPC(JROF)
      ZWSC = PWSC(JROF) * ZWSMX(JROF)
      ZWPC = PWPC(JROF) * ZWPMX(JROF)
      ZSNC = PSNC(JROF)
    ENDIF
    
!-----------------------------------------------------------------------------------
!
!*   2.1  Analyse de temperature
!
! report de l'increment de temperature a 2m sur Ts et Tp avec amortissement
!
    IF ( NNEIGT<=0 .OR. ZNEI<XREPS2 ) THEN
      ZPD = 1.0
    ELSEIF ( XSNEIGT<XREPS3 ) THEN
      ZPD = 0.0
    ELSE
      ZPD = (1.0-MIN(ZNEI,XSNEIGT)/XSNEIGT)**NNEIGT
    ENDIF

    PTS(JROF) =  PTS(JROF)  + PT2INC(JROF)*ZPD
    PTP(JROF) =  PTP(JROF)  + PT2INC(JROF)*ZPD/(XSODELX(1)/XSODELX(0))

! Rappel de Ts
    ZCLIMCA = XRCLIMTS * ZCLI
    PTS(JROF) = (1.0-ZCLIMCA)*PTS(JROF) + ZCLIMCA*ZTSC

! Rappel de Tp
    ZCLIMCA = XRCLIMTP * XRCLIMCA
    PTP(JROF) = (1.0-ZCLIMCA)*PTP(JROF )+ ZCLIMCA*ZTPC

!-----------------------------------------------------------------------------------
!
!*   2.O  Initialisations pour l'analyse de surface
!
! conditions locales d'analyse effective des champs de surface
! calcul du temps solaire local utile
!
    CALL OI_TSL(KDAT,KSSSSS,PGELAT(JROF),PGELAM(JROF),ZMU0,ZMU0M,IH)
!
    ZV10M = SQRT(PUCLS(JROF)**2+PVCLS(JROF)**2)
!
    ZPRECIP = MAX(0.,PRRCL(JROF))+ MAX(0.,PRRSL(JROF)) &
            + MAX(0.,PRRCN(JROF))+ MAX(0.,PRRSN(JROF))  
!
    IF (LFGEL) THEN
      ZWPI = PTL(JROF)
    ELSE
      ZWPI = 0.0
    ENDIF
!
! Surface water forcing to the superficial reservoir 
!
    ZDACW =  MIN(1.0,MAX(0.0,ABS(REAL(NINT(ZIVEG(JROF))-NTVGLA))) )  &
           * MIN(1.0,MAX(0.0,REAL(IH)))                              &
           * MIN(1.0,MAX(0.0,REAL(NIDJ)/REAL(NMINDJ)))               &
           * MIN(1.0,MAX(0.0,1.0-ZV10M/(XV10MX+XREPS3)))             &
           * MIN(1.0,MAX(0.0,1.0-ZPRECIP/(XSPRECIP+XREPS3)))         &
           * MIN(1.0,MAX(0.0,1.0-ZWPI/XSICE))  
!
! coefficients : dependance par rapport a l'angle zenithal solaire
!
    IF ( XSMU0>XREPS3 ) THEN
      ZPD = 0.5 * (1.0+TANH(XSMU0*(ZMU0M-0.5)))
    ELSE
      ZPD = 1.0
    ENDIF
    ZDACW = ZDACW * ZPD
!
! coefficients : dependance par rapport a l'evaporation de surface
!
!*    Seuil d'evaporation min. pour analyse de W (SEVAP en mm/jour)
    IF ( XSEVAP>XREPS3 ) THEN
      ZPD = MIN(1.0,MAX(0.0,PEVAP(JROF)/(-XSEVAP/XDAY)))
    ELSE
      ZPD = 1.0
    ENDIF
    ZDACW = ZDACW * ZPD
!
! coefficients : dependance par rapport a la nebulosite
!
    IF ( XANEBUL>XREPS3 ) THEN
      ZPD = 1.0 - XANEBUL*(PATMNEB(JROF)/ZECHGU)**NNEBUL
    ELSE
      ZPD = 1.0
    ENDIF
    ZDACW = ZDACW * ZPD
!
! coefficients : dependance par rapport a la couverture neigeuse

    IF ( NNEIGW<=0 .OR. ZNEI<XREPS2 ) THEN
      ZPD = 1.0
    ELSEIF ( XSNEIGW<XREPS3 ) THEN
      ZPD = 0.0
    ELSE
      ZPD = ( 1.0 - MIN(ZNEI,XSNEIGW)/XSNEIGW)**NNEIGW
    ENDIF
    ZDACW = ZDACW * ZPD
!
    ZDACW2 = MIN(1.0,MAX(0.0,1.0-(ZPRECIP+ABS(PEVAP(JROF)))/(XSPRECIP2+XREPS3)))  &
           * MIN(1.0,MAX(0.0,1.0-ZWPI/XSICE))
!
    ZDACW2 = ZDACW2 * ZPD

! increments de temperature et d'humidite relative a 2m utiles
!

!*   2.2  Analyse d'humidite par interpolation optimale pour ISBA


! coefficients : dependance principale par rapport a la vegetation
!
!  fctveg.h 
!****---------------Calcul de ZWSD et ZWPD------------------------------------
!
    CALL OI_FCTVEG(IH,PVEG(JROF),                                            &
                   ZVGAT1,ZVGAT2,ZVGAT3,ZVGBT1,ZVGBT2,ZVGBT3,ZVGCT1,ZVGCT2,  &
                   ZVGAH1,ZVGAH2,ZVGAH3,ZVGBH1,ZVGBH2,ZVGBH3,ZVGCH1,ZVGCH2,  &
                   ZSIGT2MP,ZSIGHP2,                                         &
                   ZG1,ZG2,ZG3,ZG4,                                          &
                   ZVGST,ZVGSH,ZVGPT1,ZVGPH1,ZVGPT2,ZVGPH2)  
!
! coefficients : dependance par rapport aux erreurs d'observation
! nb - in our case GSGOBS=.F.
!
    IF ( GSGOBS ) THEN
      ZZT = ZG1 / ZG2
      ZZH = ZG3 / ZG4 
    ELSE
      ZZT = 1.0
      ZZH = 1.0
    ENDIF
    
! coefficients : dependance par rapport a la texture
!
    ZPD = (ZWFC(JROF)-ZWWILT(JROF))/XADWR

! calcul des increments bruts pour ws=Ws/ds/ro, wp=Wp/dp/ro
! coefficients finaux
!
    ZZT = ZZT * ZPD * ZDACW
    ZZH = ZZH * ZPD * ZDACW
!
    ZVGST = ZVGST * ZZT
    ZVGSH = ZVGSH * ZZH

    ZLAISRS = PLAI(JROF)/MAX(1.0,PRSMIN(JROF))    
    ZCWPT   = ( ZVGPT1 + ZLAISRS*ZVGPT2 ) * ZZT
    ZCWPH   = ( ZVGPH1 + ZLAISRS*ZVGPH2 ) * ZZH
!
! limitation eventuelle des increments de T2m et H2m
! limitation de la valeur absolue des increments
!
    ZT2D = PT2INC(JROF)
    ZH2D = PH2INC(JROF)
    IF (XSIGT2MO < 0.0) ZT2D=MAX(XSIGT2MO,MIN(-XSIGT2MO,ZT2D))
    IF (XSIGH2MO < 0.0) ZH2D=MAX(XSIGH2MO,MIN(-XSIGH2MO,ZH2D))

! retrait du biais moyen
! soustraction du biais moyen si SCOEF(T/H) <> 1
!
    PT2MBIAS(JROF) = PT2MBIAS(JROF)*(1.0-XSCOEFT) + ZT2D*XSCOEFT
    PH2MBIAS(JROF) = PH2MBIAS(JROF)*(1.0-XSCOEFH) + ZH2D*XSCOEFH

! si le biais courant est inferieur au biais moyen on le met a zero
!                IF (ABS(ZT2D).LT.ABS(PSP_CI(JROF,YSP_CI%YCI(12)%MP0)) ZTEFF = 0.
!                IF (ABS(ZH2D).LT.ABS(PSP_CI(JROF,YSP_CI%YCI(13)%MP0)) ZHEFF = 0.
! si le biais courant est inferieur au biais effectif on le garde
!
    IF ( XSCOEFT/= 0.0 .OR. XSCOEFH/=0.0 ) THEN
      ZTEFF = ZT2D - PT2MBIAS(JROF)
      ZHEFF = ZH2D - PH2MBIAS(JROF)
      IF (ABS(ZT2D) < ABS(ZTEFF)) ZTEFF = ZT2D
      IF (ABS(ZH2D) < ABS(ZHEFF)) ZHEFF = ZH2D
      ZT2D = ZTEFF
      ZH2D = ZHEFF
    ENDIF

! increments bruts
!
    IF ( LOBS2M .AND. (.NOT.LOBSWG .OR. PWGINC(JROF)==0.0) ) THEN
      ZWSD = XRSCALDW * (ZVGST*ZT2D + ZVGSH*ZH2D)
      ZWPD = XRSCALDW * (ZCWPT*ZT2D + ZCWPH*ZH2D)
    ELSEIF ( LOBSWG ) THEN
      CALL OI_KALMAN_GAIN(ZDWG_DWG(JROF),ZDWG_DW2(JROF),PD2(JROF),ZK1,ZK2)
      ZWSD = ZK1*ZDACW2*PWGINC(JROF)
      ZWPD = ZK2*ZDACW2*PWGINC(JROF)
    ELSE
      ZWSD = 0.0
      ZWPD = 0.0
    ENDIF

! limitations sur les corrections
! pas d'analyse de ws si pas d'evaporation sur sol nu
!
    IF (PEVAP(JROF)-PEVAPTR(JROF)>= 0.0 .AND. .NOT.LOBSWG) ZWSD = 0.0

!===============================================================
! Lower limit for Wp set to  Wwilt instead of veg*Wwilt
!===============================================================
!        ZZVEG = 1.0
!=========================WS et WP================================

! analyse de wp limitee pour assurer veg*wwilt <= wp <= SWFC*wfc
    IF ( LIMVEG ) CALL GET_ZW(ZWPD,PWP(JROF),PD2(JROF))

    ! lissage des increments d'analyse de wp

    IF ( LISSEW ) THEN
      ZWPDX = ZWPD
      IF ( NLISSEW >= 3 ) THEN
        ZWPD = 0.25*(PWPINC3(JROF)+PWPINC2(JROF)+PWPINC1(JROF)+ZWPDX)    
      ELSE
        ZWPD = 0.0
      ENDIF
      IF ( NLISSEW >= 2 ) PWPINC3(JROF) = PWPINC2(JROF)
      IF ( NLISSEW >= 1 ) PWPINC2(JROF) = PWPINC1(JROF)
      PWPINC1(JROF) = ZWPDX
    ENDIF

! analyse de ws limitee pour assurer veg*wwilt <= ws <= SWFC*wfc
    IF ( LIMVEG ) CALL GET_ZW(ZWSD,PWS(JROF),XRD1)

! report des increments sur Ws, Wp
!
    ZWSA = PWS(JROF) + ZWSD * XRD1 * XRHOLW
    ZWSMIN = XREPS1 * XRD1 * XRHOLW
    PWS(JROF) = MAX(ZWSMIN,MIN(ZWSMX(JROF),ZWSA)) 

! contenu en eau total minimum
!
    ZWPA = PWP(JROF) + ZWPD * PD2(JROF) * XRHOLW
    ZWPMIN = MAX(PWS(JROF), XREPS1 * PD2(JROF) * XRHOLW)
    PWP(JROF) = MAX(ZWPMIN,MIN(ZWPMX(JROF),ZWPA))

! Rappel de Ws
    ZCLIMCA = XRCLIMWS * ZCLI
    ZCLIMCA = ZCLIMCA*PVEG(JROF) + MIN(1.0,XRCLIMV*ZCLIMCA)*(1.0-PVEG(JROF))
    PWS(JROF) = (1.0-ZCLIMCA)*PWS(JROF) + ZCLIMCA*ZWSC

! Rappel de Wp
    ZCLIMCA = XRCLIMWP * ZCLI
    ZCLIMCA = ZCLIMCA*PVEG(JROF) + MIN(1.0,XRCLIMV*ZCLIMCA)*(1.0-PVEG(JROF))
    IF ( LFGEL ) THEN
      ZGEL = ZWPI / MAX(PWP(JROF)+ZWPI,XREPS3)
      ZWPC = ZWPC * (1.0 - MAX(0.0,MIN(1.0,ZGEL)))
      ZWPC = MAX(ZWPMIN,ZWPC)
    ENDIF
    PWP(JROF) = (1.0-ZCLIMCA)*PWP(JROF) + ZCLIMCA*ZWPC

! rappel de Sn avec correction eventuelle pour fonte
!
    ZSNA  = (1.0-XRCLIMCA)*PSNS(JROF) + XRCLIMCA*ZSNC
    ZMSN  = MAX(0.0, XRSNSA/21600.*ZECHGU*(PTCLS(JROF)-XTT))**XRSNSB
    PSNS(JROF) = MAX (ZSNA-ZMSN,0.0)

    IF (LFGEL) THEN
      ZMSN = MAX(0.0, XRWPIA/21600. * ZECHGU*(PTCLS(JROF)-XTT))**XRWPIB
      PTL(JROF)=MAX (ZWPI-ZMSN ,0.0)
      PWP(JROF)=PWP(JROF)-PTL(JROF)+ZWPI
    ENDIF

!*   2.5  Rappel de SST, sur mer

!  ELSEIF ( PITM(JROF) <= 0.5 .AND. RCLISST /= 0. .AND. LCLIM ) THEN  
!    PTS(JROF) = (1.0-RCLISST)*PTS(JROF) + RCLISST *PSSTC(JROF)    
!    PTP(JROF) = PTS(JROF)
!   PWS(JROF) = XUNDEF
!   PWP(JROF) = XUNDEF
!   PTL(JROF) = 0.0
  ENDIF

!*   2.6  Mise a jour des constantes de surface sur mer,
!*        en fonction de la banquise

  IF ( PITM(JROF) <= 0.5 ) THEN
    IF ( PTS(JROF) <= XTMERGL ) THEN
      PALBF (JROF) = XSALBB
      PEMISF(JROF) = XSEMIB
      PZ0F  (JROF) = XSZZ0B*XG
      PZ0H  (JROF) = XRZHZ0G * XSZZ0B*XG
    ELSE
      PALBF (JROF) = XSALBM
      PEMISF(JROF) = XSEMIM
    ENDIF
  ENDIF

ENDDO

IF (LHOOK) CALL DR_HOOK('OI_CACSTS',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE GET_ZW(PWD,PIW,PD)
!
REAL, INTENT(INOUT) :: PWD
REAL, INTENT(IN) :: PIW, PD
REAL :: ZWR, ZWD1, ZWD2
!
ZWR = PIW / (PD*XRHOLW)
!
IF ( ZWR > ZWFC(JROF)*XSWFC ) THEN
  IF ( LHUMID ) THEN
    PWD = MIN(0.0,PWD)
  ELSE
    PWD = 0.0
  ENDIF
ELSEIF ( ZWR < ZWWILT(JROF)*PVEG(JROF)) THEN
  IF (LHUMID) THEN
    PWD = MAX(0.0,PWD)
  ELSE
    PWD = 0.0
  ENDIF
ELSE
  ZWD1 = ZWWILT(JROF)*PVEG(JROF) - ZWR
  ZWD2 = ZWFC(JROF)*XSWFC   - ZWR
  PWD = MAX(ZWD1,MIN(ZWD2,PWD))
ENDIF

END SUBROUTINE
!
!**---------------------------------------------------------------------
END SUBROUTINE OI_CACSTS
