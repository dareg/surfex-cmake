!     #########
    SUBROUTINE ALLOCATE_PHYSIO (IO, IM, KLU, KVEGTYPE )
!   ##########################################################################
!
!!****  *ALLOCATE_PHYSIO* - 
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
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    xx/xxxx
!!      Modified 10/2014 P. Samuelsson  MEB
!
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_t
!
USE MODD_TYPE_DATE_SURF
!
USE MODD_AGRI,        ONLY : LAGRIP
!
USE MODD_TREEDRAG,       ONLY : LTREEDRAG
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PARAM_t), INTENT(INOUT) :: IM
!
INTEGER, INTENT(IN) :: KLU
INTEGER, INTENT(IN) :: KVEGTYPE
!
INTEGER               :: ISIZE_LMEB_PATCH  ! Number of patches with MEB=true
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! Mask and number of grid elements containing patches/tiles:
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_PHYSIO',0,ZHOOK_HANDLE)
!
ISIZE_LMEB_PATCH=COUNT(IO%LMEB_PATCH(:))
!
ALLOCATE(IM%X%XVEGTYPE                (KLU,KVEGTYPE            ))
!
ALLOCATE(IM%X%XDG                     (KLU,IO%NGROUND_LAYER,IO%NPATCH)) 
ALLOCATE(IM%X%XD_ICE                  (KLU,IO%NPATCH              )) 
!
ALLOCATE(IM%T%XLAI                    (KLU,IO%NPATCH              )) 
ALLOCATE(IM%T%XVEG                    (KLU,IO%NPATCH              )) 
ALLOCATE(IM%T%XZ0                     (KLU,IO%NPATCH              )) 
ALLOCATE(IM%T%XEMIS                   (KLU,IO%NPATCH              )) 
!
ALLOCATE(IM%T%XRSMIN                  (KLU,IO%NPATCH              )) 
ALLOCATE(IM%T%XGAMMA                  (KLU,IO%NPATCH              )) 
ALLOCATE(IM%T%XWRMAX_CF               (KLU,IO%NPATCH              )) 
ALLOCATE(IM%T%XRGL                    (KLU,IO%NPATCH              )) 
ALLOCATE(IM%T%XCV                     (KLU,IO%NPATCH              )) 
ALLOCATE(IM%T%XALBNIR_VEG             (KLU,IO%NPATCH              )) 
ALLOCATE(IM%T%XALBVIS_VEG             (KLU,IO%NPATCH              )) 
ALLOCATE(IM%T%XALBUV_VEG              (KLU,IO%NPATCH              )) 
!
ALLOCATE(IM%X%XZ0_O_Z0H               (KLU,IO%NPATCH              )) 
!
IF (ISIZE_LMEB_PATCH>0 .OR. IO%CPHOTO/='NON') THEN
  ALLOCATE(IM%T%XBSLAI                  (KLU,IO%NPATCH              )) 
ELSE
  ALLOCATE(IM%T%XBSLAI     (0,0))  
ENDIF
! - vegetation: Ags parameters ('AGS', 'LAI', 'AST', 'LST', 'NIT' options)
!
IF (IO%CPHOTO/='NON'.OR.LTREEDRAG) THEN
  ALLOCATE(IM%X%XH_TREE                 (KLU,IO%NPATCH              ))
ELSE
  ALLOCATE(IM%X%XH_TREE                 (0,0                     ))
ENDIF
!
IF (IO%CPHOTO/='NON') THEN
  ALLOCATE(IM%X%XRE25                   (KLU,IO%NPATCH              )) 
  ALLOCATE(IM%X%XDMAX                   (KLU,IO%NPATCH              ))  
  ALLOCATE(IM%T%XLAIMIN                 (KLU,IO%NPATCH              )) 
  ALLOCATE(IM%T%XSEFOLD                 (KLU,IO%NPATCH              )) 
  ALLOCATE(IM%T%XGMES                   (KLU,IO%NPATCH              )) 
  ALLOCATE(IM%T%XGC                     (KLU,IO%NPATCH              )) 
  IF (IO%CPHOTO/='AGS' .AND. IO%CPHOTO/='LAI') THEN
    ALLOCATE(IM%T%XF2I                    (KLU,IO%NPATCH              ))
    ALLOCATE(IM%T%LSTRESS                 (KLU,IO%NPATCH              )) 
    IF (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='NCB') THEN
      ALLOCATE(IM%T%XCE_NITRO               (KLU,IO%NPATCH              )) 
      ALLOCATE(IM%T%XCF_NITRO               (KLU,IO%NPATCH              )) 
      ALLOCATE(IM%T%XCNA_NITRO              (KLU,IO%NPATCH              ))  
    ELSE
      ALLOCATE(IM%T%XCE_NITRO    (0,0))
      ALLOCATE(IM%T%XCF_NITRO    (0,0))
      ALLOCATE(IM%T%XCNA_NITRO   (0,0))
 
    ENDIF
  ELSE
    ALLOCATE(IM%T%XF2I   (0,0))
    ALLOCATE(IM%T%LSTRESS(0,0))
    ALLOCATE(IM%T%XCE_NITRO    (0,0))
    ALLOCATE(IM%T%XCF_NITRO    (0,0))
    ALLOCATE(IM%T%XCNA_NITRO   (0,0))
  ENDIF
ELSE
  ALLOCATE(IM%X%XRE25      (0,0))
  ALLOCATE(IM%X%XDMAX      (0,0))  
  ALLOCATE(IM%T%XLAIMIN    (0,0))
  ALLOCATE(IM%T%XSEFOLD    (0,0))  
  ALLOCATE(IM%T%XGMES      (0,0))
  ALLOCATE(IM%T%XGC        (0,0))
  ALLOCATE(IM%T%XF2I   (0,0))
  ALLOCATE(IM%T%LSTRESS(0,0))
  ALLOCATE(IM%T%XCE_NITRO    (0,0))
  ALLOCATE(IM%T%XCF_NITRO    (0,0))
  ALLOCATE(IM%T%XCNA_NITRO   (0,0))
ENDIF  
!
! - Irrigation, seeding and reaping
!
IF (LAGRIP .AND. (IO%CPHOTO == 'LAI' .OR. IO%CPHOTO == 'LST' .OR. IO%CPHOTO == 'NIT' .OR. IO%CPHOTO == 'NCB'))  THEN
  ALLOCATE(IM%I%TSEED                  (KLU,IO%NPATCH              )) 
  ALLOCATE(IM%I%TREAP                  (KLU,IO%NPATCH              )) 
  ALLOCATE(IM%I%XWATSUP                 (KLU,IO%NPATCH              )) 
  ALLOCATE(IM%I%XIRRIG                  (KLU,IO%NPATCH              ))
ELSE
  ALLOCATE(IM%I%TSEED     (0,0))
  ALLOCATE(IM%I%TREAP     (0,0))
  ALLOCATE(IM%I%XWATSUP    (0,0))
  ALLOCATE(IM%I%XIRRIG     (0,0))        
ENDIF
!
! - ISBA-DF scheme
!
IF(IO%CISBA=='DIF')THEN
  ALLOCATE(IM%X%XROOTFRAC  (KLU,IO%NGROUND_LAYER,IO%NPATCH))
  ALLOCATE(IM%X%NWG_LAYER  (KLU,IO%NPATCH))
  ALLOCATE(IM%X%XDROOT     (KLU,IO%NPATCH))
  ALLOCATE(IM%X%XDG2       (KLU,IO%NPATCH))
ELSE  
  ALLOCATE(IM%X%XROOTFRAC  (0,0,0))
  ALLOCATE(IM%X%NWG_LAYER  (0,0)  )
  ALLOCATE(IM%X%XDROOT     (0,0)  )        
  ALLOCATE(IM%X%XDG2       (0,0)  )        
ENDIF
!
ALLOCATE(IM%M%XGNDLITTER (KLU,IO%NPATCH))
ALLOCATE(IM%M%XRGLGV     (KLU,IO%NPATCH))
ALLOCATE(IM%M%XGAMMAGV   (KLU,IO%NPATCH))
ALLOCATE(IM%M%XRSMINGV   (KLU,IO%NPATCH))
ALLOCATE(IM%M%XROOTFRACGV(KLU,IO%NGROUND_LAYER,IO%NPATCH))
ALLOCATE(IM%M%XWRMAX_CFGV(KLU,IO%NPATCH))
ALLOCATE(IM%M%XLAIGV     (KLU,IO%NPATCH))
ALLOCATE(IM%M%XZ0LITTER  (KLU,IO%NPATCH))
ALLOCATE(IM%M%XH_VEG     (KLU,IO%NPATCH))
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_PHYSIO',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOCATE_PHYSIO
