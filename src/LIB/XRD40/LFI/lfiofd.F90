! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe LFI

SUBROUTINE LFIOFD_MT64               &
&                     (LFI, KFACMD )
USE LFIMOD, ONLY : LFICOM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!        Ce sous-programme permet d'obtenir la valeur courante
!     du FACTEUR MULTIPLICATIF par DEFAUT du logiciel LFI.
!**
!        ARGUMENT : KFACMD (Sortie) ==> Facteur multiplicatif par defaut
!
!        Le facteur multiplicatif est, pour une unite logique LFI, le
!     rapport entre la longueur d'enregistrement du fichier et LFI%JPRECL.
!
!
TYPE(LFICOM) :: LFI
INTEGER (KIND=JPLIKB) KFACMD, INIMES, IREP, INUMER
CHARACTER(LEN=LFI%JPLSPX) CLNSPR
CHARACTER(LEN=LFI%JPLMES) CLMESS
CHARACTER(LEN=LFI%JPLFTX) CLACTI

!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('LFIOFD_MT64',0,ZHOOK_HANDLE)
CLACTI=''
IF (LFI%LFIOFD_LLPREA) THEN
  CALL LFIINI_MT64                 &
&                 (LFI, 2_JPLIKB )
  LFI%LFIOFD_LLPREA=.FALSE.
ENDIF
!
IF (LFI%LMULTI) CALL LFIVER_MT64                       &
&                               (LFI, LFI%VERGLA,'ON')
KFACMD=LFI%MFACTU(0)
IF (LFI%LMULTI) CALL LFIVER_MT64                        &
&                               (LFI, LFI%VERGLA,'OFF')
INIMES=LFI%NIMESG
!
IF (INIMES.EQ.2) THEN
  IREP=0
  INUMER=LFI%JPNIL
  CLNSPR='LFIOFD'
  WRITE (UNIT=CLMESS,FMT='(''KFACMD='',I3)') KFACMD
  CALL LFIEMS_MT64                                  &
&                 (LFI, INUMER,INIMES,IREP,.FALSE., &
&                  CLMESS,CLNSPR,CLACTI)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('LFIOFD_MT64',1,ZHOOK_HANDLE)
END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE LFIOFD64           &
&           (KFACMD)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KFACMD                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIOFD_MT64            &
&           (LFI, KFACMD)

END SUBROUTINE

SUBROUTINE LFIOFD             &
&           (KFACMD)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KFACMD                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIOFD_MT             &
&           (LFI, KFACMD)

END SUBROUTINE

SUBROUTINE LFIOFD_MT             &
&           (LFI, KFACMD)
USE LFIMOD, ONLY : LFICOM
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
TYPE (LFICOM)          LFI                                    ! INOUT
INTEGER (KIND=JPLIKM)  KFACMD                                 !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  IFACMD                                 !   OUT
! Convert arguments


CALL LFIOFD_MT64            &
&           (LFI, IFACMD)

KFACMD     = INT (    IFACMD, JPLIKM)

END SUBROUTINE

!INTF KFACMD          OUT 
