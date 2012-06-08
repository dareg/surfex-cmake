C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACSIM_MT (FA,  KREP, KRANG, PCHAME, PCHAMS, 
     S                   KPULAS, KSTRON )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C     traitement des champs en coefficients spectraux, preparatoire
C     au codage GRIB.
C              ( Coefficients Spectraux, Integration Methodique ! )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANG  (Entree) ==> Rang de l'unite logique;
C    ( Tableau ) PCHAME (Entree) ==> Champ en coef. spectraux en entree;
C    ( Tableau ) PCHAMS (Sortie) ==> Champ en sortie, partie a coder;
C                KPULAS (Sortie) ==> Puissance de laplacien utilisee.
C                KSTRON (Entree) ==> Niveau de sous-troncature non
C                                    compactee.
C*
C       En mode multi-taches, il doit y avoir verrouillage du fichier
C     concerne avant l'appel au sous-programme.
C
C     Modifications
C     -------------
C
C  Juillet 1998, J. Clochard, SCEM/TTI/DAO:
C
C    -Reinitialisation de tableaux utilises pour le calcul iteratif
C     au changement de sens de balayage.
C    -Plus de "IF" pour le calcul d'extrema dans le cas ALADIN.
C    -Diagnostic plus precis en mode "mise au point".
C
C  Octobre 1998, J. Clochard, SCEM/TTI/DAO:
C
C    -Ajout de l'argument d'appel KSTRON pour compatibilite avec
C     evaluation dynamique (eventuelle) de la sous-troncature en
C     fonction de la troncature et du nombre de bits par valeur
C     compactee.
C
C  Avril   2004, D. Paradis,  DSI/DEV:
C
C    -Initialisations des tableaux XLAPxDx et FLAP1Dx faites
C     en debut de routine par appel a FAIXLA et FAIFLA.
C
C  April  2009, F. Vana and NEC:
C
C    - OpenMP directives
C
C  March 2010: J. Masek - fix of precomputed optimal Laplacian power
C              F. Vana  - simplification of IFC_SMAX,IFC_SMIN for
C                              better performance
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANG, KPULAS, KSTRON
C
      REAL (KIND=JPDBLR) PCHAME (*), PCHAMS (*)
C
      INTEGER IDIMNC, IRANGC, ITRONC, IPUFLA, IDMOPL, JN, JM, J
      INTEGER IMLIM, IOFF, IM, IMOD, INDLAP, INDZ, ILONG, IDECAL, IMINI
      INTEGER IMAXI, ILCHAM, INBITS, IMTRONC, IMODPL, JIND
      INTEGER IMEILL, JSENS, INDICE, IPUISS, IPOSEX, JMODPL
      INTEGER IPLUS, IMOINS, IPUISX, IPUIS2, IRAPOR, IPUISR, INIMES
      INTEGER INUMER, IDEB, IFIN, IXLOPA
      INTEGER IPULAS (0:1)
C
      REAL (KIND=JPDBLR) ZMIN, ZMAX, ZERRXI, ZERRXF, ZBIGVA
      REAL (KIND=JPDBLR) ZMINI (FA%JPXTRO,0:2),ZMAXI (FA%JPXTRO,0:2)
      REAL (KIND=JPDBLR) Z(4*FA%JPXTRO*FA%JPXTRO,2), ZECART (2,0:1)
C
      LOGICAL LLARPE,LLMLAM
C
      INTEGER ISMIN_1, ISMAX_1
      EXTERNAL ISMIN_1, ISMAX_1
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACSIM_MT',0,ZHOOK_HANDLE)
      IDIMNC=0
      ZBIGVA=HUGE(ZBIGVA)
C
      IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
C
C Si ce n'est pas encore fait, initialisation des tableaux XLAP... et FA%FLAP1D.
C
      IF (FA%LIXLAP) THEN
        CALL FAIXLA_MT (FA)
        FA%LIXLAP = .FALSE.
      ENDIF
      IF (FA%LIFLAP(KRANG)) THEN
        CALL FAIFLA_MT(FA, KRANG)
        FA%LIFLAP(KRANG) = .FALSE.
      ENDIF
C
      IRANGC=FA%NUCADR(KRANG)
      ITRONC=FA%MTRONC(IRANGC)
      IXLOPA=FA%NXLOPA(IRANGC)
      LLMLAM=FA%LIMLAM(IRANGC)
C
      IF (LLMLAM) IMTRONC=FA%NOZPAR(2,IRANGC)
      IF (ITRONC.LE.KSTRON) THEN
        KREP=-88
        GOTO 1001
      ELSEIF (LLMLAM.AND.IMTRONC.LE.KSTRON) THEN
        KREP=-88
        GOTO 1001
      ELSEIF (LLMLAM.AND.(IMTRONC.GT.3*ITRONC.OR.
     S    ITRONC.GT.3*IMTRONC)) THEN
C Il s'agit d'un garde-fou, modifiable (ne pas oublier FARCIS et FAPULA)
        KREP=-114
        GOTO 1001
      ELSE       
        KREP=0
      ENDIF
C
      IPUFLA=FA%NPUFLA(KRANG)
      IMODPL=FA%NMFDPL(KRANG)
C
      IF (LLMLAM) THEN
         ILCHAM=FA%NSFLAM(IRANGC)
         IDIMNC=4*(1+ITRONC+IMTRONC+(KSTRON*(KSTRON-1))/2)
CDP      IDIMNC=FA%NOZPAR(5,IRANGC)+4*KSTRON-1
      ELSE
         ILCHAM=(1+ITRONC)**2
         IDIMNC=(1+KSTRON)**2
      ENDIF
C**
C     2.  -  DETERMINATION DE LA "MEILLEURE" PUISSANCE DE LAPLACIEN
C            POUR LA PARTIE DU CHAMP QUI SERA COMPACTEE EN "GRIB".
C-----------------------------------------------------------------------
C
      IF (IMODPL.EQ.0) THEN
C
C           On elimine le cas ou aucune modulation de la puissance
C         de laplacien n'est possible.
C
        KPULAS=IPUFLA
        GOTO 300
      ENDIF
C*
C     2.1 -  AMORCAGE DU PROCESSUS ITERATIF: CALCUL DES EXTREMA DU CHAMP
C            MULTIPLIE PAR LA PUISSANCE DE LAPLACIEN NOMINALE DU FICHIER
C            ( le traitement est decoupe nombre d'onde "n" par "n" )
C-----------------------------------------------------------------------
C
C       Calcul des extrema du champ d'entree (partie a compacter),
C     pour chaque nombre d'onde "n".
C
      IF (LLMLAM) THEN
        ZMIN=ZBIGVA
        ZMAX=-ZBIGVA
!$OMP PARALLEL DO IF(FA%LOPENMP)
!$OMP&PRIVATE(JN,IMLIM,IDEB,IFIN,JIND,IOFF,IM,IMOD,INDLAP,INDZ)
!$OMP&REDUCTION(MAX:ZMAX) REDUCTION(MIN:ZMIN)
        DO 2110 JN=1,ITRONC
         IMLIM=KSTRON-JN
         IDEB=MAX(FA%NOZPAR(2*JN+3,IRANGC)+4*(1+IMLIM),
     S          FA%NOZPAR(2*JN+3,IRANGC)+4)
         IFIN=FA%NOZPAR(2*JN+4,IRANGC)
        DO 2110 JIND=IDEB,IFIN
          IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
          IM=IOFF/4
          IMOD=MOD(IOFF,4)
C
          INDLAP=((JN-1)*FA%JPXTRO)+IM
          INDZ=IMOD*FA%JPXTRO*FA%JPXTRO+INDLAP
          Z(INDZ,1)=PCHAME(JIND)*FA%FLAP1DA(INDLAP,KRANG)
          ZMAX=MAX (ZMAX,Z(INDZ,1))
          ZMIN=MIN (ZMIN,Z(INDZ,1))
C
 2110   CONTINUE
!$OMP END PARALLEL DO
      ELSE
        DO 211 JN=KSTRON+1,ITRONC
        ILONG=2*JN+1
        IDECAL=JN**2
        IMAXI=ISMAX_1 (ILONG,PCHAME(IDECAL+1))
        ZMAXI(JN,0)=PCHAME(IDECAL+IMAXI)
        IMINI=ISMIN_1 (ILONG,PCHAME(IDECAL+1))
        ZMINI(JN,0)=PCHAME(IDECAL+IMINI)
  211   CONTINUE
C
C
C
        DO 212 JN=KSTRON+1,ITRONC
        ZMAXI(JN,1)=ZMAXI(JN,0)*FA%FLAP1D(JN,KRANG)
        ZMINI(JN,1)=ZMINI(JN,0)*FA%FLAP1D(JN,KRANG)
  212   CONTINUE
C
C
        IMAXI=KSTRON+ISMAX_1 (ITRONC-KSTRON,ZMAXI(KSTRON+1,1))
        IMINI=KSTRON+ISMIN_1 (ITRONC-KSTRON,ZMINI(KSTRON+1,1))
        ZMIN=ZMINI(IMINI,1)
        ZMAX=ZMAXI(IMAXI,1)
      ENDIF
C
      INBITS=FA%NBFCSP(KRANG)
      LLARPE=FA%NFGRIB(KRANG).EQ.2
C
      IF (ZMAX.LE.ZMIN) THEN
C
C           On elimine le cas trivial du champ constant,
C         eventuellement apres transformation...
C
        KPULAS=IPUFLA
        GOTO 300
      ENDIF
C
C        Calcul de l'erreur de compactage initiale.
C
      CALL FAXION_MT (FA, PCHAME,IPUFLA,IDIMNC,ILCHAM,ZMIN,
     S             ZMAX,INBITS,LLARPE,ZERRXI,LLMLAM,FA%NOZPAR(1,IRANGC),
     S             KSTRON,ITRONC,IXLOPA)
      IMEILL=0
      ZECART(2,IMEILL)=ZERRXI
C*
C     2.3 -  BOUCLE SUR LES DEGRES DE MODULATION POSSIBLES,
C            PAR INCREMENTS DE PUISSANCE VALANT +1 (ESSAYE EN PREMIER)
C            PUIS (-1).
C-----------------------------------------------------------------------
C
      DO 239 JSENS=1,-1,-2
      INDICE=(1-JSENS)/2
      IPUISS=IPUFLA
      ZECART(1,INDICE)=ZERRXI
      IPOSEX=2
C
      IF (JSENS.EQ.-1) THEN
C
C       Compte-tenu du caractere "incremental" du calcul des extrema
C       pour des puissances successives, on doit reinitialiser lors du
C       changement de sens de balayage ZMAXI et ZMINI pour le cas ARPEGE
C       et Z pour le cas ALADIN.
C
        IF (LLMLAM) THEN
C
          ZMIN=ZBIGVA
          ZMAX=-ZBIGVA
!$OMP PARALLEL DO IF(FA%LOPENMP)
!$OMP&PRIVATE(JN,IMLIM,IDEB,IFIN,JIND,IOFF,IM,IMOD,INDLAP,INDZ)
          DO 2311 JN=1,ITRONC
          IMLIM=KSTRON-JN
          IDEB=MAX(FA%NOZPAR(2*JN+3,IRANGC)+4*(1+IMLIM),
     S           FA%NOZPAR(2*JN+3,IRANGC)+4)
          IFIN=FA%NOZPAR(2*JN+4,IRANGC)
          DO 2311 JIND=IDEB,IFIN
          IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
          IM=IOFF/4
          IMOD=MOD(IOFF,4)
C
             INDLAP=((JN-1)*FA%JPXTRO)+IM
             INDZ=IMOD*FA%JPXTRO*FA%JPXTRO+INDLAP
             Z(INDZ,1)=PCHAME(JIND)*FA%FLAP1DA(INDLAP,KRANG)
C
 2311     CONTINUE
!$OMP END PARALLEL DO
C
        ELSE
C
          DO 2312 JN=KSTRON+1,ITRONC
          ZMAXI(JN,1)=ZMAXI(JN,0)*FA%FLAP1D(JN,KRANG)
          ZMINI(JN,1)=ZMINI(JN,0)*FA%FLAP1D(JN,KRANG)
 2312     CONTINUE
C
        ENDIF
C
      ENDIF
C
      DO 238 JMODPL=1,IMODPL
      IPUISS=IPUISS+JSENS
C
      IF (LLMLAM) THEN
        ZMIN=ZBIGVA
        ZMAX=-ZBIGVA
!$OMP PARALLEL DO IF(FA%LOPENMP)
!$OMP&PRIVATE(JN,IMLIM,IDEB,IFIN,JIND,IOFF,IM,IMOD,INDLAP,INDZ)
!$OMP&REDUCTION(MAX:ZMAX) REDUCTION(MIN:ZMIN)
        DO 2310 JN=1,ITRONC
        IMLIM=KSTRON-JN
        IDEB=MAX(FA%NOZPAR(2*JN+3,IRANGC)+4*(1+IMLIM),
     S         FA%NOZPAR(2*JN+3,IRANGC)+4)
        IFIN=FA%NOZPAR(2*JN+4,IRANGC)
!ocl novrec
        DO 2313 JIND=IDEB,IFIN
        IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
        IM=IOFF/4
        IMOD=MOD(IOFF,4)
C
           INDLAP=((JN-1)*FA%JPXTRO)+IM
           INDZ=IMOD*FA%JPXTRO*FA%JPXTRO+INDLAP
           Z(INDZ,IPOSEX)=Z(INDZ,3-IPOSEX)*
     S     FA%XLAP1DA(INDLAP,INDICE)
           ZMAX=MAX (ZMAX,Z(INDZ,IPOSEX))
           ZMIN=MIN (ZMIN,Z(INDZ,IPOSEX))
C
 2313   CONTINUE
 2310   CONTINUE
!$OMP END PARALLEL DO
      ELSE 
        DO 231 JN=KSTRON+1,ITRONC
        ZMAXI(JN,IPOSEX)=ZMAXI(JN,3-IPOSEX)*FA%XLAP1D(JN,INDICE)
        ZMINI(JN,IPOSEX)=ZMINI(JN,3-IPOSEX)*FA%XLAP1D(JN,INDICE)
  231   CONTINUE
C
        IMAXI=KSTRON+ISMAX_1 (ITRONC-KSTRON,ZMAXI(KSTRON+1,IPOSEX))
        IMINI=KSTRON+ISMIN_1 (ITRONC-KSTRON,ZMINI(KSTRON+1,IPOSEX))
        ZMIN=ZMINI(IMINI,IPOSEX)
        ZMAX=ZMAXI(IMAXI,IPOSEX)
      ENDIF
C
      IF (ZMAX.LE.ZMIN) THEN
C
C           On elimine le cas du champ constant...
C
        KPULAS=IPUISS
        GOTO 240
      ENDIF
C
C        Calcul de la nouvelle erreur de compactage.
C
      CALL FAXION_MT (FA, PCHAME,IPUISS,IDIMNC,ILCHAM,ZMIN,ZMAX,INBITS,
     S             LLARPE,ZECART(IPOSEX,INDICE),LLMLAM,
     S             FA%NOZPAR(1,IRANGC),KSTRON,ITRONC,IXLOPA)
C
      IF (ZECART(IPOSEX,INDICE).GE.ZECART(3-IPOSEX,INDICE)) THEN
C
C        Ecart pas meilleur que celui calcule precedemment, on s'arrete.
C
        IPULAS(INDICE)=IPUISS-JSENS
        GOTO 239
      ENDIF
C
      IPOSEX=3-IPOSEX
  238 CONTINUE
C
C        On a epuise les degres de modulation possibles... on plafonne.
C                    (pour un sens de balayage)
C
      IPULAS(INDICE)=IPUISS
  239 CONTINUE
C
C        Choix du meilleur resultat obtenu dans les 2 sens de balayage.
C
      IPLUS=1+MOD (IPULAS(0)-IPUFLA,2)
      IMOINS=1+MOD (IPUFLA-IPULAS(1),2)
C
      IF (ZECART(IPLUS,0).LE.ZECART(IMOINS,1)) THEN
        IMEILL=0
      ELSE
        IMEILL=1
      ENDIF
C
      KPULAS=IPULAS(IMEILL)
C
  240 CONTINUE
C*
C     2.4 -  DIAGNOSTICS EVENTUELS, EN MODE MISE AU POINT SEULEMENT.
C-----------------------------------------------------------------------
C
      IF (FA%LFAMOP) THEN
        ZERRXF=MIN (ZECART(1,IMEILL),ZECART(2,IMEILL))
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         'FACSIM - Erreur Initiale (P=',IPUFLA,') ',ZERRXI,
     S         ', Finale (P=',KPULAS,') ', ZERRXF
      ENDIF
C**
C     3.  -  TRANSFORMATION DE LA PARTIE A COMPACTER DU CHAMP.
C-----------------------------------------------------------------------
C
  300 CONTINUE
C
C        On fait des multiplications plutot que des divisions,
C     et on essaie d'eviter l'exponentiation.
C
      IF (KPULAS.EQ.0) THEN
C
        IF (LLMLAM) THEN
!$OMP PARALLEL DO PRIVATE(JN,JIND) IF(FA%LOPENMP)
          DO 3010 JN=0,ITRONC
          DO 3010 JIND=FA%NOZPAR(2*JN+3,IRANGC),FA%NOZPAR(2*JN+4,IRANGC)
          PCHAMS(JIND)=PCHAME(JIND)
 3010     CONTINUE
!$OMP END PARALLEL DO
        ELSE
          DO 301 J=IDIMNC+1,ILCHAM
          PCHAMS(J)=PCHAME(J)
  301     CONTINUE
        ENDIF
C
      ELSE
        IPUISX=IABS (KPULAS)
C
        IF (KPULAS.GT.0) THEN
          INDICE=0
        ELSE
          INDICE=1
        ENDIF
C
        IF (IPUISX.LE.FA%JPUILA) THEN
C
          IF (LLMLAM) THEN
!$OMP PARALLEL DO PRIVATE(JN,JIND,IOFF,IM,INDLAP) IF(FA%LOPENMP)
            DO 3020 JN=1,ITRONC
            DO 3020 JIND=FA%NOZPAR(2*JN+3,IRANGC)+4,
     S                   FA%NOZPAR(2*JN+4,IRANGC)
            IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
            IM=IOFF/4
            INDLAP=((JN-1)*FA%JPXTRO)+IM
            PCHAMS(JIND)=PCHAME(JIND)*FA%XLAP2DA(INDLAP,IPUISX,INDICE)
 3020       CONTINUE
!$OMP END PARALLEL DO
          ELSE
            DO 302 J=IDIMNC+1,ILCHAM
            PCHAMS(J)=PCHAME(J)*FA%XLAP2D(J,IPUISX,INDICE)
  302       CONTINUE
          ENDIF
C
        ELSEIF (IPUISX.LE.2*FA%JPUILA) THEN
          IPUIS2=IPUISX/2
C
          IF (IPUISX.EQ.2*IPUIS2) THEN
C
            IF (LLMLAM) THEN
!$OMP PARALLEL DO PRIVATE(JN,JIND,IOFF,IM,INDLAP) IF(FA%LOPENMP)
              DO 3030 JN=1,ITRONC
              DO 3030 JIND=FA%NOZPAR(2*JN+3,IRANGC)+4,
     S                     FA%NOZPAR(2*JN+4,IRANGC)
              IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
              IM=IOFF/4
              INDLAP=((JN-1)*FA%JPXTRO)+IM
              PCHAMS(JIND)=PCHAME(JIND)*
     S                     FA%XLAP2DA(INDLAP,IPUIS2,INDICE)**2
 3030         CONTINUE
!$OMP END PARALLEL DO
            ELSE
              DO 303 J=IDIMNC+1,ILCHAM
              PCHAMS(J)=PCHAME(J)*FA%XLAP2D(J,IPUIS2,INDICE)**2
  303         CONTINUE
            ENDIF
C
          ELSE
C
            IF (LLMLAM) THEN
!$OMP PARALLEL DO PRIVATE(JN,JIND,IOFF,IM,INDLAP) IF(FA%LOPENMP)
              DO 3040 JN=1,ITRONC
              DO 3040 JIND=FA%NOZPAR(2*JN+3,IRANGC)+4,
     S                     FA%NOZPAR(2*JN+4,IRANGC)     
              IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
              IM=IOFF/4
              INDLAP=((JN-1)*FA%JPXTRO)+IM
              PCHAMS(JIND)=PCHAME(JIND)*
     S                FA%XLAP2DA(INDLAP,FA%JPUILA,INDICE)
     S                *FA%XLAP2DA(INDLAP,IPUISX-FA%JPUILA,INDICE)
 3040         CONTINUE
!$OMP END PARALLEL DO

            ELSE
              DO 304 J=IDIMNC+1,ILCHAM
              PCHAMS(J)=PCHAME(J)*FA%XLAP2D(J,FA%JPUILA,INDICE)
     S                *FA%XLAP2D(J,IPUISX-FA%JPUILA,INDICE)
  304         CONTINUE
            ENDIF
          ENDIF
C
        ELSE
          IRAPOR=1+(IPUISX-1)/FA%JPUILA
          IPUISR=IPUISX/IRAPOR
C
          IF (IPUISX.EQ.IRAPOR*IPUISR) THEN
C 
            IF (LLMLAM) THEN
!$OMP PARALLEL DO PRIVATE(JN,JIND,IOFF,IM,INDLAP) IF(FA%LOPENMP)
              DO 3050 JN=1,ITRONC
              DO 3050 JIND=FA%NOZPAR(2*JN+3,IRANGC)+4,
     S                     FA%NOZPAR(2*JN+4,IRANGC)
              IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
              IM=IOFF/4
              INDLAP=((JN-1)*FA%JPXTRO)+IM
              PCHAMS(JIND)=PCHAME(JIND)*
     S            FA%XLAP2DA(INDLAP,IPUISR,INDICE)**IRAPOR
 3050         CONTINUE
!$OMP END PARALLEL DO
            ELSE
              DO 305 J=IDIMNC+1,ILCHAM
              PCHAMS(J)=PCHAME(J)*FA%XLAP2D(J,IPUISR,INDICE)**IRAPOR
  305         CONTINUE
            ENDIF
C 
          ELSE
C
            IF (LLMLAM) THEN
!$OMP PARALLEL DO PRIVATE(JN,JIND,IOFF,IM,INDLAP) IF(FA%LOPENMP)
              DO 3060 JN=1,ITRONC
              DO 3060 JIND=FA%NOZPAR(2*JN+3,IRANGC)+4,
     S                     FA%NOZPAR(2*JN+4,IRANGC)
              IOFF=JIND-FA%NOZPAR(2*JN+3,IRANGC)
              IM=IOFF/4
              INDLAP=((JN-1)*FA%JPXTRO)+IM
              PCHAMS(JIND)=PCHAME(JIND)*
     S          FA%XLAP2DA(INDLAP,FA%JPUILA,INDICE)**(IRAPOR-1)*
     S          FA%XLAP2DA(INDLAP,IPUISX-FA%JPUILA*(IRAPOR-1),INDICE)
 3060         CONTINUE 
!$OMP END PARALLEL DO
            ELSE
              DO 306 J=IDIMNC+1,ILCHAM
              PCHAMS(J)=PCHAME(J)*
     S                  FA%XLAP2D(J,FA%JPUILA,INDICE)**(IRAPOR-1)
     S                *FA%XLAP2D(J,IPUISX-FA%JPUILA*(IRAPOR-1),INDICE)
  306         CONTINUE
            ENDIF
C
          ENDIF
C
        ENDIF
C
      ENDIF
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE EVENTUELLE,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      LLFATA=LLMOER (KREP,KRANG)
C
      IF (FA%LFAMOP.OR.LLFATA) THEN
        INIMES=2
        CLNSPR='FACSIM'
        INUMER=FA%JPNIIL
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4,
     S         '', PCHAME(1)='',G12.5,'', PCHAMS('',I3,'')='',G12.5,
     S         '', KPULAS='',I3)')
     S     KREP,KRANG,PCHAME(1),IDIMNC+1,PCHAMS(IDIMNC+1),KPULAS
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S               CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FACSIM_MT',1,ZHOOK_HANDLE)
      END

