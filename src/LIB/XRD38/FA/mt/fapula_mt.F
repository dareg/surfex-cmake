C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAPULA_MT (FA,  KREP, KRANG, PSPEC, KPULAP )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C     calcul de la PUissance de LAplacien "optimale" pour aplatir
C     le spectre hors de la sous-troncature.
C
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANG  (Entree) ==> Rang de l'unite logique;
C    ( Tableau ) PSPEC  (Entree) ==> Champ en coef. spectraux en entree;
C                KPULAP (Sortie) ==> Puissance de laplacien calculee.
C*
C
C     Creation
C     --------
C  Avril 1999, D. Paradis, SCEM/TTI/DEV:
C     Ce ss-programme est derive de la routine CALCOP du logiciel GRIBEX
C     du CEPMMT. Il retourne la puissance de laplacien (exprimee en un
C     nombre entier de milliemes) qui aplatit le mieux le spectre des
C     coefficients spectraux (hors sous-troncature).
C     Cette puissance de laplacien, P, sera l'exposant du facteur
C     .(n*(n+1))**P   pour les champs ARPEGE ou n est le nombre d'onde total
C     .(n**2+m**2)**P pour les champs ALADIN ou n est le nombre d'onde meridien
C                                            et m est le nombre d'onde zonal
C     Ce facteur sera applique aux coefficients spectraux hors sous-troncature
C     pour en reduire l'amplitude avant le codage GRIB.
C
C
C     Method.
C     -------
C
C     For ARPEGE case, consider F(n,m) = (n*(n+1))**(-P) * G(n,m),
C     where F(n,m) is the original spectral field and n the total wavenumber.
C     For ALADIN case (LAM), consider F(n,m) = (n**2+m**2)**(-P) * G(n,m),
C     where F(n,m) is the original spectral field, n the meridian
C     wavenumber and m the zonal wavenumber.
C     The aim is to minimize G in a 1 norm with respect to P.  This can only
C     partially be achieved.
C
C     1) For ARPEGE case, what we do is to compute H(n), where
C         H(n) = max(F(n,m))    with respect to m.
C
C        We then perform a least square fit for the equation
C            log(H(n)) = beta0+beta1*log(n*(n+1)).
C
C        To ensure a better fit for the lower end of the spectrum, we
C        apply an (arbitrary) weighting function before fitting
C             W(n) = 1.0 / (n-ISTRON+1)
C
C     2) For ALADIN case, what we do is to compute H(in), where
C         H(in) = max(F(n,m))      where n and m verify: k = n**2 + m**2
C           and "in" is the rank of k among all the reached values by k:
C           we are only interested in k values which verify
C           k = n**2 + m**2, with (n,m) belonging to spectrum without
C           the non-compacted sub-truncation area. We call inmax the number
C           of these k values. "in" belongs to interval 1,inmax.
C
C        We then perform a least square fit for the equation
C            log(H(in)) = beta0+beta1*log(k).
C
C        To ensure a better fit for the lower end of the spectrum, we
C        apply an (arbitrary) weighting function before fitting
C           W(in) = 1.0 / in
C
C     Reference.
C     ----------
C
C     Seber, G.A.F. (1979). Linear Regression Analyses.
C     John Wiley and Sons
C
C     ECMWF Research Department documentation of the IFS
C
C
C     Modifications
C     -------------
C
#include "precision.h"
C
C
C     Subroutine arguments
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANG, KPULAP
C
      REAL (KIND=JPDBLR) PSPEC(*)
C
C     Local arguments
C
      INTEGER IRANGC, ITRONC, IMTRONC, INMAX, INMIN
      INTEGER IKMAX, IDEB, IFIN, JN, JM, JK, IK, IN
      INTEGER IMLIM, JIND, IOFF, INUMER, INIMES, ISTRON
      INTEGER, DIMENSION(:), ALLOCATABLE    :: ITAB1, ITAB2
C
      REAL (KIND=JPDBLR) ZEPS, ZZ, ZXMW, ZYMW, ZWSUM, ZX, ZY
      REAL (KIND=JPDBLR) ZK, ZP, ZBETA1, ZSUM1, ZSUM2
      REAL (KIND=JPDBLR), PARAMETER :: NINE = 9.999
      REAL (KIND=JPDBLR), DIMENSION(:), ALLOCATABLE    :: ZW, ZNORM
C
      LOGICAL LLMLAM
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAPULA_MT',0,ZHOOK_HANDLE)
      IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
C
      ISTRON=FA%NSTROF(KRANG)
      IRANGC=FA%NUCADR(KRANG)
      ITRONC=FA%MTRONC(IRANGC)
      LLMLAM=FA%LIMLAM(IRANGC)
C
      IF (LLMLAM) IMTRONC=FA%NOZPAR(2,IRANGC)
      IF (ITRONC.LE.ISTRON) THEN
        KREP=-88
        GOTO 1001
      ELSEIF (LLMLAM.AND.IMTRONC.LE.ISTRON) THEN
        KREP=-88
        GOTO 1001
      ELSEIF (LLMLAM.AND.(IMTRONC.GT.3*ITRONC.OR.
     S    ITRONC.GT.3*IMTRONC)) THEN
C Il s'agit d'un garde-fou, modifiable (ne pas oublier FARCIS et FACSIM)
        KREP=-114
        GOTO 1001
      ELSE       
        KREP=0
      ENDIF
C
      ZEPS  = 1.0D-15
C
      IF (LLMLAM) THEN
        IKMAX = ITRONC*ITRONC + IMTRONC*IMTRONC
C
C  Pour le cas ALADIN, il faut creer les tableaux ITAB1 et ITAB2 :
C  ITAB1(JK) = 0 sauf pour JK tel qu'il existe (n,m) appartenant au spectre
C    hors de la sous-troncature non-compactee verifiant JK = n**2 + m**2:
C    on a alors ITAB1(JK) = 1
C  ITAB2 regroupe toutes les valeurs de JK precedentes.
C
C  ITAB1 est necessaire a la creation de ITAB2 et ITAB2 permet de boucler
C  directement sur les valeurs de JK utiles, et de dimensionner au plus
C  juste les tableaux ZW et ZNORM.
C
        ALLOCATE ( ITAB1(IKMAX) )
        ITAB1 = 0
CDP#ifdef SCALAR
        ALLOCATE ( ITAB2((ITRONC-1)*(IMTRONC-1)) )
CDP#else
CDP        ALLOCATE ( ITAB2(0:(ITRONC-1)*(IMTRONC-1)) )
CDP        ITAB2 = 0
CDP#endif
        DO JM = 1,IMTRONC
          IMLIM=ISTRON-JM
          IDEB=MAX(FA%NOMPAR(2*JM+3,IRANGC)+4*(1+IMLIM),
     S            FA%NOMPAR(2*JM+3,IRANGC)+4)
          IFIN=FA%NOMPAR(2*JM+4,IRANGC)
          DO JIND=IDEB,IFIN
            IOFF=JIND-FA%NOMPAR(2*JM+3,IRANGC)
            JN=IOFF/4
            JK=(JN**2 + JM**2)
C  On conserve les valeurs de JK atteintes
            ITAB1(JK) = 1
          ENDDO
        ENDDO
C  On construit le tableau ITAB2 des valeurs de JK atteintes
        IK = 0
        DO JK = 1,IKMAX
CDP#ifdef SCALAR
          IF (ITAB1(JK).GT.0) THEN
            IK = IK+1
            ITAB2(IK) = JK
            ITAB1(JK) = IK
          ENDIF
CDP#else
CDP       IK = IK + ITAB1(JK)
CDP       ITAB2(IK) = ITAB2(IK) + JK*ITAB1(JK)
CDP       ITAB1(JK) = IK
CDP#endif
        ENDDO
        INMIN = 1
        INMAX = IK
C
C  CAS ARPEGE
C
      ELSE
        INMIN = ISTRON + 1
        INMAX = ITRONC
      ENDIF
C
C**
C     2.  -  CALCUL DES NORMES ET DES POIDS
C-----------------------------------------------------------------------
C
      ALLOCATE ( ZW(INMAX), ZNORM(INMAX) )
      ZNORM = 0.0D0
C
      IF (LLMLAM) THEN
C
C  CAS ALADIN
C
        DO JM = 1,IMTRONC
          IMLIM=ISTRON-JM
          IDEB=MAX(FA%NOMPAR(2*JM+3,IRANGC)+4*(1+IMLIM),
     S            FA%NOMPAR(2*JM+3,IRANGC)+4)
          IFIN=FA%NOMPAR(2*JM+4,IRANGC)
          DO JIND=IDEB,IFIN
            IOFF=JIND-FA%NOMPAR(2*JM+3,IRANGC)
            JN=IOFF/4
            JK=(JN**2 + JM**2)
            ZNORM(ITAB1(JK)) = MAX(ZNORM(ITAB1(JK)), ABS(PSPEC(JIND)))
          ENDDO
        ENDDO
        DEALLOCATE ( ITAB1 )
C
C  CAS ARPEGE
C
      ELSE
C
        DO JN = INMIN , INMAX
          IDEB = JN**2+1
          IFIN = (1+JN)**2
          DO JIND = IDEB , IFIN
            ZNORM(JN) = MAX(ZNORM(JN), ABS(PSPEC(JIND)))
          ENDDO
        ENDDO
C
      ENDIF
C
      ZZ = REAL(INMAX-INMIN+1,JPDBLR)
      DO IN = INMIN, INMAX
        ZW(IN) = ZZ / REAL(IN-INMIN+1,JPDBLR)
C
C     Ensure the norms have a value which is not too small in case of
C     problems with math functions (e.g. LOG).
C
        IF(ZNORM(IN).LT.ZEPS) THEN
          ZNORM(IN) = ZEPS
          ZW(IN) = 100.D0*ZEPS
        ENDIF
      ENDDO
C*
C     3.  -  AJUSTEMENT PAR LES MOINDRES CARRES POUR DETERMINER LA PENTE
C            (ZBETA1).
C-----------------------------------------------------------------------
C
C
 300  CONTINUE
C
      ZXMW  = 0.0D0
      ZYMW  = 0.0D0
      ZWSUM = 0.0D0
C
C     Sum the weighted X and Ys.
      DO JN = INMIN, INMAX
        IF (LLMLAM) THEN
          ZX =  LOG(REAL(ITAB2(JN),JPDBLR))
        ELSE
          ZX =  LOG(REAL(JN*(JN+1),JPDBLR))
        ENDIF
        ZY    = LOG(ZNORM(JN))
        ZXMW  = ZXMW+ZX*ZW(JN)
        ZYMW  = ZYMW+ZY*ZW(JN)
        ZWSUM = ZWSUM+ZW(JN)
      ENDDO
C
C     Form mean weighted X and Y.
      ZXMW  = ZXMW / ZWSUM
      ZYMW  = ZYMW / ZWSUM
      ZSUM1 = 0.0D0
      ZSUM2 = 0.0D0
C
C     Perform a least square fit for the equation
      DO JN = INMIN, INMAX
        IF (LLMLAM) THEN
          ZX =  LOG(REAL(ITAB2(JN),JPDBLR))
        ELSE
          ZX =  LOG(REAL(JN*(JN+1),JPDBLR))
        ENDIF
        ZY    = LOG(ZNORM(JN))
        ZSUM1 = ZSUM1+ZW(JN)*(ZY-ZYMW)*(ZX-ZXMW)
        ZSUM2 = ZSUM2+ZW(JN)*(ZX-ZXMW)**2
      ENDDO
C
C
      DEALLOCATE ( ZW, ZNORM )
      IF (LLMLAM) DEALLOCATE ( ITAB2 )
C
C     On peut en tirer la pente recherchee:
      ZBETA1 = ZSUM1 / ZSUM2
C
C     Et finalement, la puissance de laplacien, arrondie en un
C     nombre entier de milliemes.
      ZP = -ZBETA1
      ZP = MAX(-NINE, MIN(NINE, ZP))
      KPULAP = NINT( ZP * 1000.0 )
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
        CLNSPR='FAPULA'
        INUMER=FA%JPNIIL
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4,
     S         '', ISTRON='',I4,'', KPULAP='',I6)')
     S     KREP,KRANG,ISTRON,KPULAP
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S               CLNSPR,CLNSPR,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAPULA_MT',1,ZHOOK_HANDLE)
      END

