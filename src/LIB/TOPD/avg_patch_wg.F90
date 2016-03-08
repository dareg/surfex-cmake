!
!     ##########################
      SUBROUTINE AVG_PATCH_WG (IMX, IP, IR, KI, PWG, PWGI, PDG)
!     ##########################
!
!!
!!    PURPOSE
!!    -------
!      from  AVERAGE_DIAG_MISC_ISBA_n
!!     ONLY for 3L cases!!
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
USE YOMHOOK   ,     ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,     ONLY : JPRB
!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments

!
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: IMX
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
 INTEGER, INTENT(IN)               :: KI
 REAL, DIMENSION(:,:), INTENT(OUT) :: PWG
 REAL, DIMENSION(:,:), INTENT(OUT) :: PWGI
 REAL, DIMENSION(:,:), INTENT(OUT) :: PDG
!      
!*      0.2    declarations of local variables
 INTEGER                         :: JJ, JLAYER, JPATCH ! loop indexes
 INTEGER                         :: IDEPTH 
 INTEGER                         :: INI, INP
 REAL                            :: ZWORK 
REAL, DIMENSION(SIZE(IP%XPATCH,1)) :: ZSUMPATCH
 !
 REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AVG_PATCH_WG',0,ZHOOK_HANDLE)
!
INI=SIZE(IP%XPATCH,1)
INP=SIZE(IP%XPATCH,2)

ZSUMPATCH(:) = 0.0
DO JPATCH=1,INP
   DO JJ=1,INI
      ZSUMPATCH(JJ) = ZSUMPATCH(JJ) + IP%XPATCH(JJ,JPATCH)
   END DO
END DO

PWG(:,:) =0.0
PWGI(:,:)=0.0
PDG(:,:) =0.0
!
! 
IF (INP/=1)THEN
  DO JPATCH=1,INP
     DO JJ=1,INI     
        IF(ZSUMPATCH(JJ) > 0.)THEN
!
          ZWORK=MAX(0.0,IMX%XDG(JJ,3,JPATCH)-IMX%XDG(JJ,2,JPATCH))
          PWG(JJ,1)  = PWG(JJ,1)  + IP%XPATCH(JJ,JPATCH) * IR%XWG(JJ,1,JPATCH)  * IMX%XDG (JJ,1,JPATCH) 
          PWG(JJ,2)  = PWG(JJ,2)  + IP%XPATCH(JJ,JPATCH) * IR%XWG(JJ,2,JPATCH)  * IMX%XDG (JJ,2,JPATCH) 
          PWG(JJ,3)  = PWG(JJ,3)  + IP%XPATCH(JJ,JPATCH) * IR%XWG(JJ,3,JPATCH)  * ZWORK
          PWGI(JJ,1) = PWGI(JJ,1) + IP%XPATCH(JJ,JPATCH) * IR%XWGI(JJ,1,JPATCH) * IMX%XDG (JJ,1,JPATCH) 
          PWGI(JJ,2) = PWGI(JJ,2) + IP%XPATCH(JJ,JPATCH) * IR%XWGI(JJ,2,JPATCH) * IMX%XDG (JJ,2,JPATCH) 
          PWGI(JJ,3) = PWGI(JJ,3) + IP%XPATCH(JJ,JPATCH) * IR%XWGI(JJ,3,JPATCH) * ZWORK
          ! 
          PDG(JJ,1) = PDG(JJ,1) + IP%XPATCH(JJ,JPATCH) * IMX%XDG(JJ,1,JPATCH)
          PDG(JJ,2) = PDG(JJ,2) + IP%XPATCH(JJ,JPATCH) * IMX%XDG (JJ,2,JPATCH)
          PDG(JJ,3) = PDG(JJ,3) + IP%XPATCH(JJ,JPATCH) * IMX%XDG (JJ,3,JPATCH)
          !
!          
       ENDIF
     ENDDO
  ENDDO     
!     
 WHERE (PDG(:,1)>0.0)
    PWG(:,1)  = PWG(:,1)  / PDG(:,1)
    PWGI(:,1)  = PWGI(:,1)  / PDG(:,1)
 ENDWHERE
 WHERE (PDG(:,2)>0.0)
    PWG(:,2)  = PWG(:,2)  / PDG(:,2)
    PWGI(:,2)  = PWGI(:,2)  / PDG(:,2)
 ENDWHERE
 WHERE (PDG(:,3)-PDG(:,2)>0.0)
    PWG(:,3)  = PWG(:,3)  / (PDG(:,3)-PDG(:,2))
    PWGI(:,3)  = PWGI(:,3)  / (PDG(:,3)-PDG(:,2))
 ENDWHERE
ELSE
    PWG(:,:)  = IR%XWG(:,:,1)
    PWGI(:,:) = IR%XWGI(:,:,1)
    PDG (:,:) = IMX%XDG (:,:,1)

ENDIF 
!

IF (LHOOK) CALL DR_HOOK('AVG_PATCH_WG',1,ZHOOK_HANDLE)

END SUBROUTINE AVG_PATCH_WG


