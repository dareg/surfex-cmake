C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAXION_MT (FA,  PCHAME, KPUISS, KDIMNC, KLCHAM, PMIN, 
     S                      PMAX, KNBITS, LDARPE, PECART, LDMLAM, 
     S                      KNOZPA, KSTROF, KTRONC, KXLOPA )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C     calcul de l'erreur totale par le codage GRIB
C     sur un champ en coefficients spectraux, en mode "complex packing".
C**
C    Arguments :
C    ( Tableau ) PCHAME (Entree) ==> Champ en coef. spectraux en entree;
C                KPUISS (Entree) ==> Puissance de laplacien a utiliser;
C                KDIMNC (Entree) ==> Taille de la zone non codee;
C                KLCHAM (Entree) ==> Longueur totale du champ en c.spx.;
C                PMIN   (Entree) ==> Minimum precalcule du champ module;
C                PMAX   (Entree) ==> Maximum      "      "    "     "  ;
C                KNBITS (Entree) ==> Nombre de bits par valeur a coder;
C                LDARPE (Entree) ==> Vrai si GRIB non standard;
C                PECART (Sortie) ==> Erreur relative maximum commise,
C                                    en valeur absolue;
C                LDMLAM (Entree) ==> Vrai si fichier aladin
C    (Tableau)   KNOZPA (Entree) ==> Nombre d'onde zonal maxi par parallele;
C                                    (du pole nord vers l'equateur seulement)
C                KSTROF (Entree) ==> Sous-troncature non compactee EFFECTIVE (coef. spectraux)
C                KTRONC (Entree) ==> Troncature (NS pour Aladin)
C                KXLOPA (Entree) ==> Nombre maximum de longitudes par cercle de latitudes
C
C     Modifications
C     -------------
C
C  Juillet 1998, J. Clochard, SCEM/TTI/DAO:
C
C    -Le critere a minimiser est mis sous la forme de fonctions
C    -Pour le cas LDARPE=.TRUE., utilisation des memes gardes-fous
C     que dans CODEGA.
C
C  Septem. 2000, D. Paradis, SCEM/TTI/DEV:
C
C    -Ajustement de la puissance Dolby pour des nb d'onde (JN)
C     uniquement inferieurs a MIN(FA%MTRONC(),(FA%NXLOPA()-1)/3): on
C     s'affranchit de l'influence des grands nb d'onde sur des
C     grilles lineaires pour le choix du Dolby.
C
C  Mars 2003, R. El Khatib CNRM/GMAP:
C
C    -Optimisation : ISAMAX remplace par MAXLOC
C
C  March 2010: J. Masek - fix of precomputed optimal Laplacian power
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KPUISS, KDIMNC, KLCHAM, KNBITS, KSTROF
      INTEGER KTRONC, KXLOPA
C
      INTEGER KNOZPA (FA%JPXIND)
      INTEGER ILOC(1)
C
      REAL (KIND=JPDBLR) PMIN, PMAX, PECART
      REAL (KIND=JPDBLR) PCHAME(KLCHAM)
      REAL (KIND=JPDBLR) ZERR  (KLCHAM)
      REAL (KIND=JPDBLR) ZECART_LOC
C
      LOGICAL LDARPE, LDMLAM
C
      INTEGER JN, JIND, ISCALE, J, INDICE, IPUISX, IOFF, IM
      INTEGER INDLAP, IRAPOR, IPUISR, IMLIM, IEXP, IMANT, IPUIS2
      INTEGER IZZZZZ, IDEB, IFIN, ITRDOL, ILCHADO
C
      REAL (KIND=JPDBLR) ZREFER, ZDIFFR, ZCOMPA, ZAUXI1, ZAUXI2, ZS
      REAL (KIND=JPDBLR) ZYPR, ZYP0, ZCRITR, ZCRIT0, ZCOEFF
#include "facom_mt.h"
C
C
C
C       Fonctions intermediaires utilisees pour alleger l'ecriture.
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      ZYPR ( IZZZZZ, ZCOEFF ) = ZREFER + ZAUXI2 * ANINT ( 
     S                  ZAUXI1 * ( PCHAME(IZZZZZ)/ZCOEFF - ZREFER ) )
C
      ZYP0 ( IZZZZZ ) = ZREFER + ZAUXI2 * ANINT (
     S                  ZAUXI1 * ( PCHAME(IZZZZZ) - ZREFER ) )
C
C       Fonction "critere" a minimiser lorsque la modulation de la
C       puissance de laplacien est possible.
C
      ZCRITR ( IZZZZZ, ZCOEFF ) = ( ZYPR(IZZZZZ,ZCOEFF)*ZCOEFF - 
     S                      PCHAME(IZZZZZ) )**2
C
C       Meme chose, mais pour la puissance 0.
C
      ZCRIT0 ( IZZZZZ ) = (ZYP0(IZZZZZ) - PCHAME(IZZZZZ) )**2
      IF (LHOOK) CALL DR_HOOK('FAXION_MT',0,ZHOOK_HANDLE)
C**
C     0.  -  INITIALISATIONS FAITES AU PREMIER APPEL DU SOUS-PROGRAMME.
C-----------------------------------------------------------------------
C
      IF (FA%FAXION_LLPREA) THEN
        FA%FAXION_ZEPSIL=TINY(FA%FAXION_ZEPSIL)
        FA%FAXION_ISCALX=99
        FA%FAXION_LLPREA=.FALSE.
      ENDIF
C**
C     1.  -  CALCUL DIRECT, EN ESSAYANT DE MINIMISER LES OPERATIONS.
C            ( a l'aide des puissances precalculees dans FA%XLAP2D )
C-----------------------------------------------------------------------
C
      IF (LDARPE) THEN
C
        ZDIFFR=PMAX-PMIN
C
        IF ( ZDIFFR .LE. FA%FAXION_ZEPSIL ) THEN
          ZAUXI1=MIN ( ABS (PMIN), ABS (PMAX) )
          IF ( ZAUXI1 .LE. FA%FAXION_ZEPSIL ) ZAUXI1=0.
          PMAX=SIGN (ZAUXI1,PMAX)
          PMIN=PMAX
          ZAUXI1=0.
          ZAUXI2=0.
        ELSE
          ZCOMPA=2**KNBITS-1
          ZAUXI1=ZCOMPA/ZDIFFR
          ZAUXI2=ZDIFFR/ZCOMPA
        ENDIF
C
        ZREFER=PMIN
C
      ELSE
C
        CALL CONFI (PMIN,IEXP,IMANT,ZREFER)
        ZS = (PMAX-ZREFER)/(2**(KNBITS+1)-1)
        ZAUXI1=1.
        ZAUXI2=2.
        IF (ZS.NE.0.) ZS = LOG(ZS)/LOG(ZAUXI2) + ZAUXI2
        ISCALE = MIN(INT(ZS),INT(ZS+SIGN(ZAUXI1,ZS)))
        ISCALE = MAX(-FA%FAXION_ISCALX,MIN(FA%FAXION_ISCALX,ISCALE))
        ZAUXI1 = ZAUXI2**(-ISCALE)
        ZAUXI2 = ZAUXI2**ISCALE
C
      ENDIF
C
C  Calcul de la borne superieure de la plage des nb d'onde
C  servant pour l'ajustement du Dolby (borne inferieure=KSTROF)
C  valable uniquement dans le cas Arpege pour ne pas prendre
C  en compte les coeff spectraux relatifs aux nb d'onde non
C  quadratiques (grille lineaire).
C
      IF (LDARPE.AND..NOT.LDMLAM) THEN
        ITRDOL = MIN ( KTRONC , (KXLOPA-1)/3 )
        ILCHADO = (ITRDOL+1)**2
      ENDIF
C
      IF (KPUISS.EQ.0) THEN
C
        IF (LDMLAM) THEN 
!$OMP PARALLEL DO PRIVATE(JN,JIND) IF(FA%LOPENMP)
          DO 1010 JN=0,KTRONC
          DO 1010 JIND=KNOZPA(2*JN+3),KNOZPA(2*JN+4)
          ZERR(JIND)=ZCRIT0(JIND)
 1010     CONTINUE
!$OMP END PARALLEL DO
        ELSE
          DO 101 J=KDIMNC+1,ILCHADO
          ZERR(J)=ZCRIT0(J)
  101     CONTINUE
        ENDIF
C
      ELSE
        IPUISX=IABS (KPUISS)
C
        IF (KPUISS.GT.0) THEN
          INDICE=1
        ELSE
          INDICE=0
        ENDIF
C
        IF (IPUISX.LE.FA%JPUILA) THEN
C
C        Dans ce cas precis, on pourrait aussi remplacer la division
C     par ZCOEFF par une multiplication par FA%XLAP2D(J,IPUISX,1-INDICE).
C     Il y aurait gain en calcul, mais reference memoire supplementaire.
C
          IF (LDMLAM) THEN
#if !defined(RS6K)
!$OMP PARALLEL DO PRIVATE(JN,JIND,IOFF,IM,INDLAP,ZCOEFF)
!$OMP&  IF(FA%LOPENMP)
#endif
            DO 1020 JN=1,KTRONC
            DO 1020 JIND=KNOZPA(2*JN+3)+4,KNOZPA(2*JN+4)
            IOFF=JIND-KNOZPA(2*JN+3)
            IM=IOFF/4 
            INDLAP=((JN-1)*FA%JPXTRO)+IM
            ZERR(JIND)=ZCRITR(JIND,FA%XLAP2DA(INDLAP,IPUISX,INDICE))
 1020       CONTINUE
#if !defined(RS6K)
!$OMP END PARALLEL DO
#endif
          ELSE
            DO 102 J=KDIMNC+1,ILCHADO
            ZERR(J)=ZCRITR(J,FA%XLAP2D(J,IPUISX,INDICE))
  102       CONTINUE
          ENDIF
C
        ELSEIF (IPUISX.LE.2*FA%JPUILA) THEN
          IPUIS2=IPUISX/2
C
          IF (IPUISX.EQ.2*IPUIS2) THEN
C
            IF (LDMLAM) THEN
#if !defined(RS6K)
!$OMP PARALLEL DO PRIVATE(JN,JIND,IOFF,IM,INDLAP,ZCOEFF)
!$OMP&  IF(FA%LOPENMP)
#endif
              DO 1030 JN=1,KTRONC
              DO 1030 JIND=KNOZPA(2*JN+3)+4,KNOZPA(2*JN+4)
              IOFF=JIND-KNOZPA(2*JN+3)
              IM=IOFF/4 
              INDLAP=((JN-1)*FA%JPXTRO)+IM
              ZERR(JIND)=ZCRITR(JIND,
     S                   FA%XLAP2DA(INDLAP,IPUIS2,INDICE)**2)
 1030         CONTINUE
#if !defined(RS6K)
!$OMP END PARALLEL DO
#endif
            ELSE
              DO 103 J=KDIMNC+1,ILCHADO
              ZERR(J)=ZCRITR(J,FA%XLAP2D(J,IPUIS2,INDICE)**2)
  103         CONTINUE
            ENDIF
C
          ELSE
C
          IF (LDMLAM) THEN
#if !defined(RS6K)
!$OMP PARALLEL DO PRIVATE(JN,JIND,IOFF,IM,INDLAP,ZCOEFF)
!$OMP&  IF(FA%LOPENMP)
#endif
             DO 1040 JN=1,KTRONC
             DO 1040 JIND=KNOZPA(2*JN+3)+4,KNOZPA(2*JN+4)
              IOFF=JIND-KNOZPA(2*JN+3)
              IM=IOFF/4 
              INDLAP=((JN-1)*FA%JPXTRO)+IM
              ZERR(JIND)=ZCRITR(JIND,FA%XLAP2DA(INDLAP,FA%JPUILA,INDICE)
     S             *FA%XLAP2DA(INDLAP,IPUISX-FA%JPUILA,INDICE))
 1040         CONTINUE
#if !defined(RS6K)
!$OMP END PARALLEL DO
#endif
          ELSE
              DO 104 J=KDIMNC+1,ILCHADO
              ZERR(J)=ZCRITR(J,FA%XLAP2D(J,FA%JPUILA,INDICE)
     S             *FA%XLAP2D(J,IPUISX-FA%JPUILA,INDICE))
  104       CONTINUE
          ENDIF
C
          ENDIF
C
        ELSE
          IRAPOR=1+(IPUISX-1)/FA%JPUILA
          IPUISR=IPUISX/IRAPOR
C
          IF (IPUISX.EQ.IRAPOR*IPUISR) THEN
C
            IF (LDMLAM) THEN
#if !defined(RS6K)
!$OMP PARALLEL DO PRIVATE(JN,JIND,IOFF,IM,INDLAP,ZCOEFF)
!$OMP&  IF(FA%LOPENMP)
#endif
              DO 1050 JN=1,KTRONC
              DO 1050 JIND=KNOZPA(2*JN+3)+4,KNOZPA(2*JN+4)
              IOFF=JIND-KNOZPA(2*JN+3)
              IM=IOFF/4 
              INDLAP=((JN-1)*FA%JPXTRO)+IM
              ZERR(JIND)=ZCRITR(JIND,
     S               FA%XLAP2DA(INDLAP,IPUISR,INDICE)**IRAPOR)
 1050         CONTINUE
#if !defined(RS6K)
!$OMP END PARALLEL DO
#endif
            ELSE
              DO 105 J=KDIMNC+1,ILCHADO
              ZERR(J)=ZCRITR(J,FA%XLAP2D(J,IPUISR,INDICE)**IRAPOR)
  105         CONTINUE
            ENDIF
C
          ELSE
C
           IF (LDMLAM) THEN
#if !defined(RS6K)
!$OMP PARALLEL DO PRIVATE(JN,JIND,IOFF,IM,INDLAP,ZCOEFF)
!$OMP&  IF(FA%LOPENMP)
#endif
              DO 1060 JN=1,KTRONC
              DO 1060 JIND=KNOZPA(2*JN+3)+4,KNOZPA(2*JN+4)
              IOFF=JIND-KNOZPA(2*JN+3)
              IM=IOFF/4 
              INDLAP=((JN-1)*FA%JPXTRO)+IM
              ZERR(JIND)=ZCRITR(JIND,
     S           FA%XLAP2DA(INDLAP,FA%JPUILA,INDICE)**(IRAPOR-1)
     S          *FA%XLAP2DA(INDLAP,IPUISX-FA%JPUILA*(IRAPOR-1),INDICE))
 1060         CONTINUE
#if !defined(RS6K)
!$OMP END PARALLEL DO
#endif
           ELSE
              DO 106 J=KDIMNC+1,ILCHADO
              ZERR(J)=ZCRITR(J,
     S             FA%XLAP2D(J,FA%JPUILA,INDICE)**(IRAPOR-1)
     S            *FA%XLAP2D(J,IPUISX-FA%JPUILA*(IRAPOR-1),INDICE))
  106         CONTINUE
            ENDIF
C
          ENDIF
C
        ENDIF
C
      ENDIF
C
      PECART=0.
      IF (LDMLAM) THEN
!$OMP PARALLEL IF(FA%LOPENMP)
!$OMP& PRIVATE(JN,IMLIM,IDEB,IFIN,JIND,ZECART_LOC)
!$OMP& REDUCTION(+:PECART)
!$OMP DO
        DO 1070 JN=1,KTRONC
        IMLIM=KSTROF-JN
        IDEB=MAX(KNOZPA(2*JN+3)+4,KNOZPA(2*JN+3)+4*(1+IMLIM))
        IFIN=KNOZPA(2*JN+4)
        ZECART_LOC=0.
        DO 1071 JIND=IDEB,IFIN
          ZECART_LOC=ZECART_LOC+ZERR(JIND)
 1071   CONTINUE
        PECART=PECART+ZECART_LOC
C
 1070   CONTINUE
!$OMP END DO 
!$OMP END PARALLEL
      ELSE
        DO J=KDIMNC+1,ILCHADO
          PECART=PECART+ZERR(J)
        ENDDO
      ENDIF
C
      IF (FA%LFAMOP) THEN
        IF (LDARPE.AND..NOT.LDMLAM) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FAXION: KPUISS=', KPUISS,', KDIMNC=',KDIMNC,
     S           ', KLCHAM=', KLCHAM, ', ITRDOL=', ITRDOL, ', ILCHADO=',
     S           ILCHADO, ', PMIN=', PMIN, ', PMAX=', PMAX
        ELSE
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FAXION: KPUISS=', KPUISS,', KDIMNC=',KDIMNC,
     S           ', KLCHAM=', KLCHAM,
     S           ', PMIN=', PMIN, ', PMAX=', PMAX
        ENDIF
        WRITE (UNIT=FA%NULOUT,FMT=*)'FAXION: PECART=', PECART
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAXION_MT',1,ZHOOK_HANDLE)
      END


