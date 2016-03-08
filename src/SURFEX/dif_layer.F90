!#############################################################
SUBROUTINE DIF_LAYER(KLU, IO, IP, MX   )  
!#############################################################
!
!!****  *DIF_LAYER_n* - routine to initialize dif numbers of layers
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!    S. Faroux
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2012!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_SGH_PAR,        ONLY : XHORT_DEPTH
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, INTENT(IN) :: KLU
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: MX
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(KLU) :: ZWORK
INTEGER, DIMENSION(KLU,IO%NPATCH) :: IWORK
INTEGER :: JLAYER, JPATCH, JILU, IDEPTH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('DIF_LAYER',0,ZHOOK_HANDLE)
!
DO JLAYER = 1, IO%NGROUND_LAYER
  !DO JPATCH=1,SIZE(IP%XPATCH,2)
  !  DO JILU=1,SIZE(IP%XPATCH,1)
  !    IF (IP%XPATCH(JILU,JPATCH)/=0.) THEN
  !      IF (MX%XROOTFRAC(JILU,JLAYER,JPATCH)<0..OR.MX%XROOTFRAC(JILU,JLAYER,JPATCH)>1.) then
  !              print*,JILU,JLAYER,JPATCH,MX%XROOTFRAC(JILU,JLAYER,JPATCH)
  !      endif
  !    endif
  !  enddo
  !enddo
  IF (ANY((MX%XROOTFRAC(:,JLAYER,:)<0. .OR. MX%XROOTFRAC(:,JLAYER,:)>1.) .AND. IP%XPATCH(:,:).NE.0.)) &
    CALL ABOR1_SFX('DIF_LAYER: WITH CISBA=DIF ROOTFRAC MUST BE DEFINED')
ENDDO
!
IP%XDZG     (:,:,:) = XUNDEF
IP%XDZDIF   (:,:,:) = XUNDEF
IP%XSOILWGHT(:,:,:) = 0.0
!
DO JPATCH=1,IO%NPATCH
!
  IF (IP%NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
!
!*   soil layers thicknesses
  IP%XDZG(:,1,JPATCH) = MX%XDG(:,1,JPATCH)
  DO JLAYER=2,IO%NGROUND_LAYER
    DO JILU=1,KLU
      IP%XDZG(JILU,JLAYER,JPATCH) = MX%XDG(JILU,JLAYER,JPATCH) - MX%XDG(JILU,JLAYER-1,JPATCH)
    ENDDO
  ENDDO
!
!*   distance between consecuative layer mid-points
  DO JLAYER=1,IO%NGROUND_LAYER
    DO JILU=1,KLU
      IF(JLAYER<IO%NGROUND_LAYER)THEN
        IP%XDZDIF(JILU,JLAYER,JPATCH)=0.5*(IP%XDZG(JILU,JLAYER,JPATCH)+IP%XDZG(JILU,JLAYER+1,JPATCH))
      ELSE
        IP%XDZDIF(JILU,JLAYER,JPATCH)=0.5*IP%XDZG(JILU,JLAYER,JPATCH) 
      ENDIF
    ENDDO
  ENDDO 
! 
ENDDO
!
! Horton runoff parameter
!
IWORK(:,:) = MX%NWG_LAYER(:,:)
!
DO JPATCH=1,IO%NPATCH
!  
  IF( IP%NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
!
  DO JILU=1,KLU
    IDEPTH = MX%NWG_LAYER(JILU,JPATCH)
    IF (IDEPTH==NUNDEF) IDEPTH = IO%NGROUND_LAYER
    DO JLAYER=1,IDEPTH-1
      IF(MX%XDG(JILU,JLAYER,JPATCH)<XHORT_DEPTH) IWORK(JILU,JPATCH)=JLAYER+1
    ENDDO
  ENDDO
!
END DO
!  
IO%NLAYER_HORT=MAXVAL(IWORK(:,:),IWORK(:,:)/=NUNDEF)
!  
! Dunne runoff parameter
!
IWORK(:,:)=MX%NWG_LAYER(:,:)
!
DO JPATCH=1,IO%NPATCH
!  
  IF (IP%NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
!
  DO JILU=1,KLU
    IF(IP%XPATCH(JILU,JPATCH)>0.0)THEN 
      IDEPTH = MX%NWG_LAYER(JILU,JPATCH)    
      IF(MX%XDROOT(JILU,JPATCH)>0.0.AND.MX%XDROOT(JILU,JPATCH)/=XUNDEF)THEN
        IP%XRUNOFFD(JILU,JPATCH) = MX%XDG(JILU,1,JPATCH)
        DO JLAYER=1,IDEPTH-1
          IF(MX%XROOTFRAC(JILU,JLAYER,JPATCH)<0.90)THEN
            IP%XRUNOFFD(JILU,JPATCH) = MX%XDG(JILU,JLAYER+1,JPATCH)
          ENDIF
        ENDDO
      ELSE
        IP%XRUNOFFD(JILU,JPATCH) = MIN(0.6,MX%XDG2(JILU,JPATCH))
      ENDIF
    ENDIF
  ENDDO
!
  ZWORK(:) = 0.0
  DO JLAYER=1,IO%NGROUND_LAYER
    DO JILU=1,KLU
      IF(IP%XPATCH(JILU,JPATCH)>0.0)THEN
        IDEPTH=MX%NWG_LAYER(JILU,JPATCH)
        IF(JLAYER<=IDEPTH)THEN
          ZWORK    (JILU              ) = ZWORK(JILU) + IP%XDZG(JILU,JLAYER,JPATCH)  
          IP%XSOILWGHT(JILU,JLAYER,JPATCH) = MIN(IP%XDZG(JILU,JLAYER,JPATCH), &
                                          MAX(0.0,IP%XRUNOFFD(JILU,JPATCH)-ZWORK(JILU)+IP%XDZG(JILU,JLAYER,JPATCH)))
        ENDIF
        IF(MX%XDG(JILU,JLAYER,JPATCH)<IP%XRUNOFFD(JILU,JPATCH))THEN
          IWORK(JILU,JPATCH)=JLAYER+1
        ENDIF
      ENDIF
    ENDDO
  ENDDO
!  
END DO
!
IO%NLAYER_DUN=MAXVAL(IWORK(:,:),IWORK(:,:)/=NUNDEF)
!
IF (LHOOK) CALL DR_HOOK('DIF_LAYER',1,ZHOOK_HANDLE)
!
END SUBROUTINE DIF_LAYER
