!
!     ##########################
      SUBROUTINE DG_DFTO3L(KI,PDG)
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
!
USE MODD_ISBA_n,    ONLY : NGROUND_LAYER, XROOTFRAC, &
                           XDG, XSOILWGHT, XPATCH, &
                           NSIZE_NATURE_P, XDZG,   &
                           NWG_LAYER, XDG2, NPATCH
USE MODD_SURF_PAR,  ONLY : XUNDEF, NUNDEF
USE MODD_PACK_ISBA, ONLY : NK_WG_LAYER
USE YOMHOOK   ,     ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,     ONLY : JPRB
!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments

 INTEGER, INTENT(IN)               :: KI
 REAL, DIMENSION(:,:), INTENT(OUT) :: PDG
!      
!*      0.2    declarations of local variables
 INTEGER                         :: JJ, JLAYER, JPATCH ! loop indexes
 INTEGER                         :: IDEPTH 
 INTEGER                         :: INI, INP
 REAL                            :: ZWORK 
 !
 REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DG_DFTO3L',0,ZHOOK_HANDLE)
INI=SIZE(XPATCH,1)
INP=SIZE(XPATCH,2)
!
PDG(:,:)=0.0
!
    DO JPATCH=1,INP
!  
      IF (NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
      DO JLAYER = 1,NGROUND_LAYER
        DO JJ=1,INI
          IDEPTH=NWG_LAYER(JJ,JPATCH)
          IF(JLAYER<=IDEPTH.AND.IDEPTH/=NUNDEF.AND.XPATCH(JJ,JPATCH)/=XUNDEF)THEN
            !
            PDG      (JJ,1) = PDG      (JJ,1) + XDG(JJ,1,JPATCH) * XPATCH(JJ,JPATCH) 
            ! ISBA-FR-DG2 comparable soil wetness index, liquid water and ice contents
            ZWORK=MIN(XDZG(JJ,JLAYER,JPATCH),MAX(0.0,XDG2(JJ,JPATCH)-XDG(JJ,JLAYER,JPATCH)+XDZG(JJ,JLAYER,JPATCH)))
            PDG      (JJ,2) = PDG      (JJ,2) + ZWORK * XPATCH(JJ,JPATCH) 
            !
            ! ISBA-FR-DG3 comparable soil wetness index, liquid water and ice contents
            ZWORK=MIN(XDZG(JJ,JLAYER,JPATCH),MAX(0.0,XDG(JJ,JLAYER,JPATCH)-XDG2(JJ,JPATCH)))
            PDG      (JJ,3) = PDG      (JJ,3) + ZWORK * XPATCH(JJ,JPATCH) 
            !
          ENDIF
        ENDDO
      ENDDO
!
    ENDDO
    ! 
     PDG (:,3) =  PDG (:,2) + PDG (:,3)
     WHERE (PDG(:,:)==0.0)
             PDG(:,:)=XUNDEF
     ENDWHERE
!
IF (LHOOK) CALL DR_HOOK('DG_DFTO3L',1,ZHOOK_HANDLE)

END SUBROUTINE DG_DFTO3L


