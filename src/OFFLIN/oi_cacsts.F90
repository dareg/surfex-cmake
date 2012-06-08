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
!        F. Bouyssel    27-Mar-2011  Use of REPS2 instead of REPS3 for ZNEIG
!****---------------------------------------------------------------------------
!
SUBROUTINE OI_CACSTS(KNBPT,PT2INC,PH2INC,PWGINC,PWS_O,                        &
                       IDAT,NSSSSS,                                           &
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
INTEGER,INTENT(IN)    :: KNBPT, IDAT, NSSSSS
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
REAL   ,INTENT(INOUT) :: PVEG(KNBPT) 
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
REAL  :: VGAT1(24),VGAT2(24),VGAT3(24)
REAL  :: VGBT1(24),VGBT2(24),VGBT3(24)
REAL  :: VGCT1(24),VGCT2(24)
REAL  :: VGAH1(24),VGAH2(24),VGAH3(24)
REAL  :: VGBH1(24),VGBH2(24),VGBH3(24)
REAL  :: VGCH1(24),VGCH2(24)
REAL  :: SIGT2MP(24),SIGHP2(24)
!
REAL  :: VGST,VGSH,VGPT1,VGPH1,VGPT2,VGPH2,G1,G2,G3,G4
!
REAL :: ZITS(KNBPT), ZITP(KNBPT), ZIWS(KNBPT), ZIWP(KNBPT)
REAL :: ZIVEG(KNBPT), ZWFC(KNBPT), ZWWILT(KNBPT)
REAL :: ZWSMX(KNBPT), ZWPMX(KNBPT), ZWPI(KNBPT), ZWSAT(KNBPT)
REAL :: ZISN(KNBPT)
REAL :: ZDWG_DWG(KNBPT), ZDWG_DW2(KNBPT)
!
REAL    :: ZCLIM, ZCLIMCA, ZCOEF, ZCWPH, ZCWPT, ZCWSH,     &
            ZCWST, ZDW, ZECHGU, ZEVAP, ZGEL,                &
            ZH2D, ZHEFF, ZHSA, ZHSP, ZLAISRS, ZMSN, ZNEIG,  &
            ZPDN, ZPDS, ZPDT, ZPRECIP, ZSNA, ZSNC, ZT2D,    &
            ZTEFF, ZTINER, ZTPC, ZTSC, ZV10M, ZVEG, ZWP,    &
            ZWPA, ZWPC, ZWPD, ZWPD1, ZWPD2, ZWPDX, ZWPMIN,  &
            ZWPR, ZWSA, ZWSC, ZWSD, ZWSD1, ZWSD2, ZWSMSPI,  &
            ZWSR, ZZH, ZZT, ZDACW, ZMU0, ZMU0M, ZPDM, ZPDV, &
            ZK1,  ZK2, ZDACW2, ZPR_EVA, ZWSMIN, ZZVEG
!
INTEGER :: IH, JROF
LOGICAL :: LSGOBS
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!--------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('OI_CACSTS',0,ZHOOK_HANDLE)
ZECHGU = REAL(NECHGU) * 3600.
!
!**  1. Initialisation des polynomes bruts et des champs de reference.    
!
 CALL OI_CAVEGI(VGAT1,VGAT2,VGAT3,VGBT1,VGBT2,VGBT3,VGCT1,VGCT2,       &
                 VGAH1,VGAH2,VGAH3,VGBH1,VGBH2,VGBH3,VGCH1,VGCH2,      &
                 SIGT2MP,SIGHP2,LSGOBS)  
!
!
!*    Seuil d'evaporation min. pour analyse de W (SEVAP en mm/jour)

ZEVAP  = -SEVAP / XDAY

!*   1.2  Initialisation des variables intermediaires
!

DO JROF = 1 , KNBPT
  ZIVEG(JROF) = ANINT(PIVEG(JROF))
  IF (LFGEL) THEN
    ZWPI(JROF) = PTL(JROF)
  ELSE
    ZWPI (JROF) = 0.0
  ENDIF
ENDDO
!
 CALL OI_ACSOLW (1,KNBPT,                          &
                  PARG,PD2,PWS,ZIVEG,PSAB,         &
                  LLDHMT,                          &
                  ZWFC, ZWPMX, ZWSAT, ZWSMX, ZWWILT)  
!
! Analytical Jacobians for WG assimilation
!
 CALL OI_JACOBIANS (KNBPT,PWS_O,PSAB,PARG,PD2,PWP,ZDWG_DWG,ZDWG_DW2) 
!
!**---------------------------------------------------------------------
!**  - 2 - Calcul des champs analyses.

DO JROF = 1 , KNBPT

! stockage des champs prevus

    ZITS(JROF) = PTS(JROF)
    ZITP(JROF) = PTP(JROF)
    ZIWS(JROF) = PWS(JROF)
    ZIWP(JROF) = PWP(JROF)
    ZVEG       = PVEG(JROF)
    ZISN(JROF) = PSNS(JROF)
    ZNEIG      = MAX(0.0,PSNS(JROF)/(PSNS(JROF)+WCRIN))
!
! analyse de surface ou rappel vers la climatologie, sur terre

!    IF (PITM(JROF) > 0.5.AND.(RCLIMCA >= 0.0).AND.PWS(JROF)/=XUNDEF) THEN
    IF (PWS(JROF)/=XUNDEF) THEN

! analyse de surface

!*   2.O  Initialisations pour l'analyse de surface

! conditions locales d'analyse effective des champs de surface
! calcul du temps solaire local utile
!
        ZV10M=SQRT(PUCLS(JROF)**2+PVCLS(JROF)**2)
        ZPRECIP = MAX(0.,PRRCL(JROF))+ MAX(0.,PRRSL(JROF)) &
                 + MAX(0.,PRRCN(JROF))+ MAX(0.,PRRSN(JROF))  
!
! Surface water forcing to the superficial reservoir 
!
        ZPR_EVA = ZPRECIP + ABS(PEVAP(JROF))  
!
        CALL OI_TSL(IDAT,NSSSSS,PGELAT(JROF),PGELAM(JROF),ZMU0,ZMU0M,IH)
        ZMU0 = ZMU0M
!
        ZDACW = MIN(1.0,MAX(0.0,ABS(REAL(NINT(ZIVEG(JROF))-NTVGLA))))   &
               * MIN(1.0,MAX(0.0,REAL(IH)))                              &
               * MIN(1.0,MAX(0.0,REAL(IDJ)/REAL(MINDJ)))                 &
               * MIN(1.0,MAX(0.0,1.0-ZV10M/(V10MX+REPS3)))               &
               * MIN(1.0,MAX(0.0,1.0-ZPRECIP/(SPRECIP+REPS3)))           &
               * MIN(1.0,MAX(0.0,1.0-ZWPI(JROF)/SICE))  
!
        ZDACW2 = MIN(1.0,MAX(0.0,1.0-ZPR_EVA/(SPRECIP2+REPS3)))         &
                * MIN(1.0,MAX(0.0,1.0-ZWPI(JROF)/SICE))  
!
! coefficients : dependance par rapport a l'angle zenithal solaire
!
        IF ( SMU0 > REPS3 ) THEN
          ZPDM=0.5*(1.0+TANH(SMU0*(ZMU0-0.5)))
        ELSE
          ZPDM=1.0
        ENDIF
        ZDACW = ZDACW * ZPDM
!
! coefficients : dependance par rapport a l'evaporation de surface
!
        IF ( SEVAP > REPS3 ) THEN
          ZPDV=MIN(1.0,MAX(0.0,PEVAP(JROF)/ZEVAP))
        ELSE
          ZPDV=1.0
        ENDIF
        ZDACW = ZDACW * ZPDV
!
! coefficients : dependance par rapport a la nebulosite
!
        IF ( ANEBUL > REPS3 ) THEN
          ZPDN=1.0-ANEBUL*(PATMNEB(JROF)/ZECHGU)**NNEBUL
        ELSE
          ZPDN=1.0
        ENDIF
        ZDACW = ZDACW * ZPDN
!
! increments de temperature et d'humidite relative a 2m utiles
!
        ZT2D = PT2INC(JROF)
        ZH2D = PH2INC(JROF)
!
!*   2.1  Analyse de temperature
!
! report de l'increment de temperature a 2m sur Ts et Tp avec amortissement
!
!
        ZTINER = SODELX(1)/SODELX(0)
        IF (NNEIGT <= 0.OR. ZNEIG < REPS2) THEN
          ZPDT= 1.0
        ELSEIF (SNEIGT < REPS3) THEN
          ZPDT= 0.0
        ELSE
          ZPDT= (1.0- MIN(ZNEIG,SNEIGT)/SNEIGT)**NNEIGT
        ENDIF

        PTS(JROF) =  PTS(JROF)  + ZT2D*ZPDT
        PTP(JROF) =  PTP(JROF)  + ZT2D*ZPDT/ZTINER

!*   2.2  Analyse d'humidite par interpolation optimale pour ISBA


! coefficients : dependance principale par rapport a la vegetation
!
!  fctveg.h 
!****-------------------------------------------------------------------
!
        CALL OI_FCTVEG(IH,ZVEG,                                               &
                        VGAT1,VGAT2,VGAT3,VGBT1,VGBT2,VGBT3,VGCT1,VGCT2,      &
                        VGAH1,VGAH2,VGAH3,VGBH1,VGBH2,VGBH3,VGCH1,VGCH2,      &
                        SIGT2MP,SIGHP2,                                       &
                        G1,G2,G3,G4,                                          &
                        VGST,VGSH,VGPT1,VGPH1,VGPT2,VGPH2)  
!
        ZLAISRS = PLAI(JROF)/MAX(1.0,PRSMIN(JROF))
        ZCWST   = VGST
        ZCWSH   = VGSH
        ZCWPT   = VGPT1 + ZLAISRS*VGPT2
        ZCWPH   = VGPH1 + ZLAISRS*VGPH2

! coefficients : dependance par rapport a la texture

        ZDW = (ZWFC(JROF)-ZWWILT(JROF))/ADWR

! coefficients : dependance par rapport aux erreurs d'observation
! nb - in our case LSGOBS=.F.

        IF ( LSGOBS ) THEN
          ZZT = G1 / G2
          ZZH = G3 / G4 
        ELSE
          ZZT = 1.0
          ZZH = 1.0
        ENDIF

! coefficients : dependance par rapport a la couverture neigeuse

        IF (NNEIGW <= 0.OR. ZNEIG < REPS2) THEN
          ZPDS= 1.0
        ELSEIF (SNEIGW < REPS3) THEN
          ZPDS= 0.0
        ELSE
          ZPDS= (1.0- MIN(ZNEIG,SNEIGW)/SNEIGW)**NNEIGW
        ENDIF
        ZDACW = ZDACW * ZPDS
!
        ZDACW2 = ZDACW2 * ZPDS

! calcul des increments bruts pour ws=Ws/ds/ro, wp=Wp/dp/ro
! coefficients finaux
 
        ZCWST = ZCWST * ZDW * ZZT * ZDACW
        ZCWSH = ZCWSH * ZDW * ZZH * ZDACW
        ZCWPT = ZCWPT * ZDW * ZZT * ZDACW
        ZCWPH = ZCWPH * ZDW * ZZH * ZDACW

! limitation eventuelle des increments de T2m et H2m
! limitation de la valeur absolue des increments

        IF (SIGT2MO < 0.0) ZT2D=MAX(SIGT2MO,MIN(-SIGT2MO,ZT2D))
        IF (SIGH2MO < 0.0) ZH2D=MAX(SIGH2MO,MIN(-SIGH2MO,ZH2D))

! retrait du biais moyen
! soustraction du biais moyen si SCOEF(T/H) <> 1

        PT2MBIAS(JROF)= PT2MBIAS(JROF)*(1.0-SCOEFT)+ZT2D*SCOEFT    
        PH2MBIAS(JROF)= PH2MBIAS(JROF)*(1.0-SCOEFH)+ZH2D*SCOEFH    
        ZTEFF = ZT2D - PT2MBIAS(JROF)
        ZHEFF = ZH2D - PH2MBIAS(JROF)

! si le biais courant est inferieur au biais moyen on le met a zero
!                IF (ABS(ZT2D).LT.ABS(PSP_CI(JROF,YSP_CI%YCI(12)%MP0)) ZTEFF = 0.
!                IF (ABS(ZH2D).LT.ABS(PSP_CI(JROF,YSP_CI%YCI(13)%MP0)) ZHEFF = 0.
! si le biais courant est inferieur au biais effectif on le garde

        IF ( (SCOEFT /= 0.0) .OR. (SCOEFH /= 0.0) ) THEN
          IF (ABS(ZT2D) < ABS(ZTEFF)) ZTEFF = ZT2D
          IF (ABS(ZH2D) < ABS(ZHEFF)) ZHEFF = ZH2D
          ZT2D = ZTEFF
          ZH2D = ZHEFF
        ENDIF

! increments bruts

        IF (LOBSWG) THEN
          CALL OI_KALMAN_GAIN(ZDWG_DWG(JROF),ZDWG_DW2(JROF),PD2(JROF),ZK1,ZK2)
          ZWSD = ZK1*ZDACW2*PWGINC(JROF)
          ZWPD = ZK2*ZDACW2*PWGINC(JROF)
          IF (LOBS2M) THEN
            IF (PWGINC(JROF) == 0.0) THEN
              ZWSD = RSCALDW*(ZCWST*ZT2D + ZCWSH*ZH2D)
              ZWPD = RSCALDW*(ZCWPT*ZT2D + ZCWPH*ZH2D)                    
            ENDIF
          ENDIF 
        ELSEIF (LOBS2M) THEN
          ZWSD = RSCALDW*(ZCWST*ZT2D + ZCWSH*ZH2D)
          ZWPD = RSCALDW*(ZCWPT*ZT2D + ZCWPH*ZH2D)        
        ELSE
          ZWSD = 0.0
          ZWPD = 0.0
        ENDIF

! limitations sur les corrections
! pas d'analyse de ws si pas d'evaporation sur sol nu

        IF (PEVAP(JROF)-PEVAPTR(JROF) >= 0.0 .AND. .NOT.LOBSWG)THEN
          ZWSD = 0.0
        ENDIF

        ZZVEG = ZVEG
!===============================================================
! Lower limit for Wp set to  Wwilt instead of veg*Wwilt
!===============================================================
!        ZZVEG = 1.0
!===============================================================

! analyse de wp limitee pour assurer veg*wwilt <= wp <= SWFC*wfc

        IF ( LIMVEG ) THEN

          ZWPR = ZIWP(JROF)/(PD2(JROF)*XRHOLW)
          IF ( ZWPR > ZWFC(JROF)*SWFC ) THEN
            IF ( LHUMID ) THEN
              ZWPD = MIN(0.0,ZWPD)
            ELSE
              ZWPD = 0.0
            ENDIF
          ELSEIF ( ZWPR < ZWWILT(JROF)*ZZVEG ) THEN
            IF ( LHUMID ) THEN
              ZWPD = MAX(0.0,ZWPD)
            ELSE
              ZWPD = 0.0
            ENDIF
          ELSE
            ZWPD1 = ZWWILT(JROF)*ZZVEG -ZWPR
            ZWPD2 = ZWFC(JROF)*SWFC   -ZWPR
            ZWPD = MAX(ZWPD1,MIN(ZWPD2,ZWPD))
          ENDIF

! analyse de ws limitee pour assurer veg*wwilt <= ws <= SWFC*wfc

          ZWSR = ZIWS(JROF)/(RD1*XRHOLW)
          IF ( ZWSR > ZWFC(JROF)*SWFC ) THEN
            IF ( LHUMID ) THEN
              ZWSD = MIN(0.0,ZWSD)
            ELSE
              ZWSD = 0.0
            ENDIF
          ELSEIF ( ZWSR < ZWWILT(JROF)*ZVEG) THEN
            IF (LHUMID) THEN
              ZWSD = MAX(0.0,ZWSD)
            ELSE
              ZWSD = 0.0
            ENDIF
          ELSE
            ZWSD1 = ZWWILT(JROF)*ZVEG -ZWSR
            ZWSD2 = ZWFC(JROF)*SWFC   -ZWSR
            ZWSD = MAX(ZWSD1,MIN(ZWSD2,ZWSD))
          ENDIF
        ENDIF

! lissage des increments d'analyse de wp

        IF ( LISSEW ) THEN
          ZWPDX = ZWPD
          IF ( NLISSEW >= 3 ) THEN
            ZWPD =.25*(PWPINC3(JROF)+PWPINC2(JROF)+PWPINC1(JROF)+ZWPDX)    
          ELSE
            ZWPD = 0.0
          ENDIF
          IF ( NLISSEW >= 2 ) THEN
            PWPINC3(JROF)=PWPINC2(JROF)
          ENDIF
          IF ( NLISSEW >= 1 ) THEN
            PWPINC2(JROF)=PWPINC1(JROF)
          ENDIF
          PWPINC1(JROF)=ZWPDX
        ENDIF

! report des increments sur Ws, Wp

        ZWSA = PWS(JROF) + ZWSD*RD1*XRHOLW
        ZWPA = PWP(JROF) + ZWPD*PD2(JROF)*XRHOLW
        ZWSMIN = REPS1*RD1*XRHOLW      
        PWS(JROF) = MAX(ZWSMIN,MIN(ZWSMX(JROF),ZWSA))        

! contenu en eau total minimum

        ZWPMIN = MAX(PWS(JROF),REPS1*PD2(JROF)*XRHOLW)
        PWP(JROF) = MAX(ZWPMIN,MIN(ZWPMX(JROF),ZWPA))

!*   2.4  Rappel vers la climatologie


! mise a jour des champs climatologiques

        IF ( .NOT. LCLIM ) THEN
          ZTSC = ZITS(JROF)
          ZTPC = ZITP(JROF)
          ZWSC = ZIWS(JROF)
          ZWPC = ZIWP(JROF)
          ZSNC = PSNS(JROF)
        ELSE
          ZTSC = PTSC(JROF)
          ZTPC = PTPC(JROF)
          ZWSC = PWSC(JROF) * ZWSMX(JROF)
          ZWPC = PWPC(JROF) * ZWPMX(JROF)
          ZSNC = PSNC(JROF)
        ENDIF

        ZCLIM = RCLIMCA /(1.0+RCLIMN*ZNEIG)

! Rappel de Ts
        ZCLIMCA = RCLIMTS * ZCLIM
        PTS(JROF) = (1.0-ZCLIMCA)*PTS(JROF)+    ZCLIMCA  *ZTSC
! Rappel de Tp
        ZCLIMCA = RCLIMTP * RCLIMCA
        PTP(JROF) = (1.0-ZCLIMCA)*PTP(JROF)+    ZCLIMCA  *ZTPC
! Rappel de Ws
        ZCLIMCA = RCLIMWS * ZCLIM
        ZCLIMCA = ZCLIMCA* ZVEG + MIN(1.0,RCLIMV*ZCLIMCA)* (1.0-ZVEG)
        PWS(JROF) = (1.0-ZCLIMCA)*PWS(JROF)+    ZCLIMCA  *ZWSC
! Rappel de Wp
        ZCLIMCA = RCLIMWP * ZCLIM
        ZCLIMCA = ZCLIMCA* ZVEG + MIN(1.0,RCLIMV*ZCLIMCA)* (1.0-ZVEG)
        IF ( LFGEL ) THEN
          ZWP = PWP(JROF)
          ZGEL = ZWPI(JROF) / MAX(ZWP+ZWPI(JROF),REPS3)
          ZWPC = ZWPC * (1.0 - MAX(0.0,MIN(1.0,ZGEL)))
          ZWPC = MAX(ZWPMIN,ZWPC)
        ENDIF
        PWP(JROF) = (1.0-ZCLIMCA)*PWP(JROF) + ZCLIMCA*ZWPC

! rappel de Sn avec correction eventuelle pour fonte

        ZSNA = (1.0-RCLIMCA)*PSNS(JROF) + RCLIMCA*ZSNC
        ZCOEF= RSNSA/21600. * ZECHGU
        ZMSN = MAX (0.0, ZCOEF*(PTCLS(JROF)-XTT))**RSNSB
        PSNS(JROF) = MAX (ZSNA-ZMSN ,0.0)

        IF (LFGEL) THEN
          ZCOEF= RWPIA/21600. * ZECHGU
          ZMSN = MAX (0.0, ZCOEF*(PTCLS(JROF)-XTT))**RWPIB
          PTL(JROF)=MAX (ZWPI(JROF)-ZMSN ,0.0)
          PWP(JROF)=PWP(JROF)-PTL(JROF)+ZWPI(JROF)
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
    IF ( PTS(JROF) <= TMERGL ) THEN
      PALBF(JROF)   = SALBB
      PEMISF(JROF)  = SEMIB
      PZ0F(JROF)    = SZZ0B*XG
      PZ0H(JROF)    = RZHZ0G*SZZ0B*XG
    ELSE
      PALBF(JROF)   = SALBM
      PEMISF(JROF)  = SEMIM
    ENDIF
  ENDIF

ENDDO
IF (LHOOK) CALL DR_HOOK('OI_CACSTS',1,ZHOOK_HANDLE)
!
!**---------------------------------------------------------------------
END SUBROUTINE OI_CACSTS
