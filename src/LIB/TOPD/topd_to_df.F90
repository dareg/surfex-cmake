!
!     ##########################
      SUBROUTINE TOPD_TO_DF(KI,PWG)
!     ##########################
!
!!
!!    PURPOSE
!!    -------
!     This routines updates the soil water content of ISBA DIF afeter TOPODYN
!     lateral distribution  
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
!
USE MODD_ISBA_n,        ONLY : NGROUND_LAYER, NPATCH,&
                               XWG, XDG, XDZG, XDG2, &
                               XWSAT, NWG_LAYER, &
                               NSIZE_NATURE_P
USE MODD_SURF_PAR,      ONLY : XUNDEF, NUNDEF
USE MODD_COUPLING_TOPD, ONLY : XTOTBV_IN_MESH, XFRAC_D3
USE MODD_ISBA_PAR,      ONLY : XWGMIN
!
USE YOMHOOK   ,         ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,         ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments

 INTEGER, INTENT(IN)              :: KI
 REAL, DIMENSION(:,:), INTENT(IN) :: PWG
!      
!*      0.2    declarations of local variables
REAL                              :: ZWORK          ! numbers of layers in root and deep zones
INTEGER                           :: IDEPTH
INTEGER                           :: JI, JLAYER, JPATCH ! loop indexes
REAL(KIND=JPRB)                   :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TOPD_TO_DF',0,ZHOOK_HANDLE)
!
DO JPATCH=1,NPATCH

 IF (NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE

 DO JLAYER = 1,NGROUND_LAYER

  DO JI=1,KI

  IDEPTH=NWG_LAYER(JI,JPATCH)

     IF(JLAYER<=IDEPTH.AND.IDEPTH/=NUNDEF.AND.(XTOTBV_IN_MESH(JI)/=0.0).AND.(XTOTBV_IN_MESH(JI)/=XUNDEF)) THEN

     ! root layers
      IF (XDZG(JI,JLAYER,JPATCH)/=XUNDEF.AND.XDG2(JI,JPATCH)/=XUNDEF.AND.XDG(JI,JLAYER,JPATCH)/=XUNDEF)&! 
        ZWORK=MIN(XDZG(JI,JLAYER,JPATCH),MAX(0.0,XDG2(JI,JPATCH)-XDG(JI,JLAYER,JPATCH)+XDZG(JI,JLAYER,JPATCH)))

      IF ((PWG(JI,2)/=XUNDEF).AND.(ZWORK>0.).AND.(ZWORK/=XUNDEF))& 
        XWG(JI,JLAYER,JPATCH)=MIN(MAX(PWG(JI,2),XWGMIN),XWSAT(JI,JLAYER)) 

     ! deep layers
     IF ((XFRAC_D3(JI)/=0.0).AND.(XFRAC_D3(JI)/=XUNDEF)) THEN     

      IF (XDZG(JI,JLAYER,JPATCH)/=XUNDEF.AND.XDG2(JI,JPATCH)/=XUNDEF.AND.XDG(JI,JLAYER,JPATCH)/=XUNDEF) &
        ZWORK=MIN(XDZG(JI,JLAYER,JPATCH),MAX(0.0,XDG(JI,JLAYER,JPATCH)-XDG2(JI,JPATCH)))

      IF ((PWG(JI,3)/=XUNDEF).AND.(ZWORK>0.).AND.(ZWORK/=XUNDEF)) &! 
        XWG(JI,JLAYER,JPATCH)=MIN(MAX(PWG(JI,3),XWGMIN),XWSAT(JI,JLAYER))

     ENDIF

    ENDIF

  ENDDO
 ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('TOPD_TO_DF',1,ZHOOK_HANDLE)

END SUBROUTINE TOPD_TO_DF


