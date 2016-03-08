!
!     ##########################
      SUBROUTINE DISPATCH_WG (IP, IMX, IR, KI, PWG, PWGI, PDG)
!     ##########################
!
!!
!!    PURPOSE
!!    -------
!      from  AVERAGE_DIAG_MISC_ISBA_n
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      
!!    REFERENCE
!!    ---------
!!     
!!    AUTHOR
!!    ------
!!
!!       ELYAZIDI/HEYMES/RISTOR * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original  02/2011 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE MODD_SURF_PAR,  ONLY : XUNDEF, NUNDEF
USE MODD_ISBA_PAR,      ONLY : XWGMIN
USE MODD_COUPLING_TOPD, ONLY :  XATOP
!
USE YOMHOOK   ,     ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,     ONLY : JPRB
!
!
IMPLICIT NONE
!
!*      0.1    declarations of argumentsXPATCH
!
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: IMX
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
 INTEGER, INTENT(IN)               :: KI
 REAL, DIMENSION(:,:), INTENT(IN) :: PWG
 REAL, DIMENSION(:,:), INTENT(IN) :: PWGI
 REAL, DIMENSION(:,:), INTENT(IN) :: PDG
!      
!*      0.2    declarations of local variables
 INTEGER                         :: JJ, JLAYER, JPATCH ! loop indexes
 INTEGER                         :: IDEPTH 
 INTEGER                         :: INI, INP
 REAL                            :: ZWORK,ZTMP, ZWORK2
 REAL, DIMENSION(SIZE(IP%XPATCH,1)) :: ZSUMPATCH
 REAL, DIMENSION(SIZE(IP%XPATCH,1),SIZE(IP%XPATCH,2)) :: ZFRAC_PATCH2
 REAL, DIMENSION(SIZE(IP%XPATCH,1),SIZE(IP%XPATCH,2)) :: ZFRAC_PATCH3
 REAL, DIMENSION(SIZE(PWG,1),SIZE(PWG,2)) :: ZWG_CTL
 !
 REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DISPATCH_WG',0,ZHOOK_HANDLE)
!
INI=SIZE(IP%XPATCH,1)
INP=SIZE(IP%XPATCH,2)
!
 !DO JPATCH=1,INP
 !
 ! write(*,*) 'In dispatch XPATCH (1)',JPATCH,XPATCH(1,JPATCH),XWG(1,2,JPATCH)
 !ENDDO
!write(*,*) 'In dispatch wg ,KI,INI,INP ',KI,INI,INP
IF (INP/=1)THEN
 DO JPATCH=1,INP
   DO JJ=1,INI
    IF ((IP%XPATCH(JJ,JPATCH)/=XUNDEF).AND.(IP%XPATCH(JJ,JPATCH)/=0.)&
        .AND.(XATOP(JJ)==1.)) THEN
     WHERE (IR%XWG(JJ,:,JPATCH)/=XUNDEF)
      IR%XWG(JJ,:,JPATCH) = PWG(JJ,:) 
      IR%XWGI(JJ,:,JPATCH)= PWGI(JJ,:) 
      IMX%XDG(JJ,:,JPATCH) = PDG (JJ,:)
     ENDWHERE
    ENDIF
   ENDDO
 ENDDO

ELSE
 IR%XWG(:,:,1) = PWG(:,:) 
 IR%XWGI(:,:,1)= PWGI(:,:) 
 IMX%XDG(:,:,1) = PDG (:,:)
ENDIF
!
WHERE (IR%XWG(:,:,:)<XWGMIN) 
  IR%XWG(:,:,:)=XWGMIN
ENDWHERE
 !     
IF (LHOOK) CALL DR_HOOK('DISPATCH_WG',1,ZHOOK_HANDLE)

END SUBROUTINE DISPATCH_WG


