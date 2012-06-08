C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACDEC_MT (FA, KREP, PA, PMIN, KNBIT, KDEC)
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!****
!      Sous-programme de calcul du FACteur d'echelle DECimal optimal
!     pour coder un champ d'amplitude PA donnee, sur KNBIT bits avec
!     une precision maximale, compte tenu d'un intervalle fixe
!     de facteurs decimaux.
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                                    0: OK; -1: pb rencontre
!                PA     (Entree) ==> Amplitude du champ a compacter;
!                PMIN   (Entree) ==> Minimum du champ a compacter;
!                KNBIT  (Entree) ==> Nb de bits servant au compactage;
!                KDEC   (Sortie) ==> Facteur d'echelle decimal optimal;
!
!     Modifications:
!      R. El Khatib 12-Aug-2009 Bugfix for compatibility with Gribex
!
      IMPLICIT NONE
#include "precision.h"
      TYPE(FA_COM) :: FA
      REAL (KIND=JPDBLR) PA, PMIN
      INTEGER KREP, KNBIT, KDEC
!
      REAL (KIND=JPDBLR) XNBINT, XTINYR4, XHUGER4
      INTEGER IDECMIN, IDECMAX, INBINT
      INTEGER JDEC, IE, INUTIL, INUMAX, IEMAX
      INTEGER INU0, IE0
#include "facom_mt.h"
!
!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACDEC_MT',0,ZHOOK_HANDLE)
      IF (KNBIT.LE.0 .OR. KNBIT.GT.64) THEN
        KREP=-1
        WRITE (UNIT=FA%NULOUT,FMT=*)'****'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '**** FACDEC: ERROR, bits number out of range 1-64'
        WRITE (UNIT=FA%NULOUT,FMT=*)'****         KNBIT = ',KNBIT
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '**** ! Optimal decimal scale factor is not computed !'
        WRITE (UNIT=FA%NULOUT,FMT=*)'****'
        IF (LHOOK) CALL DR_HOOK('FACDEC_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      IF (ABS(PA) .LT. TINY(PA)) THEN
        KREP = 0
        KDEC = 0
        IF (FA%LFAMOP) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)'////'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '//// FACDEC: WARNING, range of the field is null :',
     S                        PA
          WRITE (UNIT=FA%NULOUT,FMT=*)'////'
        ENDIF
        IF (LHOOK) CALL DR_HOOK('FACDEC_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      KREP    = 0
      INUMAX  = 0
      IDECMIN = -15
      IDECMAX =   5
      XTINYR4 = TINY(1._4)
      XHUGER4 = HUGE(1._4)
      INBINT  = 2**KNBIT -1
      XNBINT  = FLOAT(INBINT)
! Cas du facteur decimal nul (reference a calculer dans tous les cas)
      JDEC    = 0
      CALL FACTEC_MT(FA, KREP, PA, KNBIT, JDEC, IE0, INU0)
      
!
!**
!     2.  -  BOUCLE TESTANT LES FACTEURS DECIMAUX
!-----------------------------------------------------------------------
!
!!OCL SCALAR
      DO JDEC = IDECMIN, IDECMAX
!
! 3 tests sur la pertinence de ce facteur decimal
!
! 0/ PA * 10**JDEC > 1.E-12
!    Du fait d'un test exagere dans gribex 
        IF (PA * 10._JPDBLR**JDEC .LE. 1.E-12) CYCLE
! 1/ PMIN * 10**JDEC > TINY(real*4) pour decodage eventuel
!    de la valeur de reference en arithmetique 32 bits IEEE.
        IF (ABS(PMIN) .GT. TINY(PMIN)) THEN
          IF (LOG10(ABS(PMIN))+JDEC .LE. LOG10( XTINYR4 ) ) CYCLE
        ELSE
!  instruction "bidon" pour que l'optimisation ne fasse pas
!  l'evaluation de LOG10(ABS(PMIN)) par anticipation (Pb sur VPP5000)
          PMIN = 0._JPDBLR
        ENDIF
! 2/ PA * 10**JDEC inclus dans le domaine de validite des real*8
        IF (ABS(LOG10(ABS(PA))+JDEC) .GE. RANGE(PA)) CYCLE
!
        CALL FACTEC_MT(FA, KREP, PA, KNBIT, JDEC, IE, INUTIL)
        IF (KREP.NE.0) CYCLE
! 3/ PMIN*10**JDEC + (2**KNBIT-1)*2**IE < HUGE(real*4) pour decodage
!    eventuel du max du champ en arithmetique 32 bits IEEE.
        IF (PMIN*10._JPDBLR**JDEC + XNBINT*2._JPDBLR**IE .GE. 
     S                          XHUGER4) CYCLE
! 4/ Le facteur d'echelle binaire doit pouvoir etre code sur 1 octet,
!    il doit donc etre compris entre -126 et 127:
        IF (IE.LT.-126 .OR. IE.GT.127) CYCLE
!
        IF (INUTIL.GT.INUMAX) THEN
          INUMAX = INUTIL
          IEMAX  = IE
          KDEC   = JDEC
        ENDIF
      ENDDO
!
      IF (INUMAX.EQ.0) THEN
        KREP=-1
        WRITE (UNIT=FA%NULOUT,FMT=*)'****'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '**** FACDEC: all the decimal factors comprised between'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '**** ',IDECMIN,' and ',IDECMAX,' are rejected !!'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '**** Range and min of the field are :', PA, PMIN
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '****'
        WRITE (UNIT=FA%NULOUT,FMT=*)'****'
        IF (LHOOK) CALL DR_HOOK('FACDEC_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      IF (FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         'FACDEC: champ d''amplitude ',PA,' ,de minimum ',PMIN
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '        => fact decimal opt de',KDEC,
     S          ' ,pour 1 fact binaire de ',IEMAX
        WRITE (UNIT=FA%NULOUT,FMT='(1X,A,I3,A,I9,A,I9,A,F5.1,A,E11.4)')
     S         ' Eff des',
     S         KNBIT,' bits = ',INUMAX,' sur ',INBINT,' soit: ',
     S         FLOAT(INUMAX)*100./XNBINT,' % et une precision de ',
     S         10.**(-KDEC)*2.**(IEMAX-1)
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '        a comparer, si le fact decimal = 0, avec'
        WRITE (UNIT=FA%NULOUT,FMT='(1X,A,I9,A,I9,A,F5.1,A,E11.4)')
     S         ' une efficacite de ',
     S         INU0,' sur ',INBINT,' soit: ',FLOAT(INU0)*100./XNBINT,
     S         ' % et une precision de ',2.**(IE0-1)
      ENDIF
!
      IF (LHOOK) CALL DR_HOOK('FACDEC_MT',1,ZHOOK_HANDLE)
      END

