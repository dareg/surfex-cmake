!     #########
    SUBROUTINE ALLOCATE_TEB_GARDEN(HINIT,KLU,KVEGTYPE,KGROUND_LAYER,KPATCH,KLW, &
                                     KDIMTAB)  
!   ##########################################################################
!
USE MODD_TEB_GARDEN_n
USE MODD_DIAG_TEB_GARDEN_n
USE MODD_AGRI_GARDEN_n
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
CHARACTER(LEN=3),INTENT(IN)  :: HINIT     ! choice of fields to initialize
INTEGER, INTENT(IN) :: KLU
INTEGER, INTENT(IN) :: KVEGTYPE
INTEGER, INTENT(IN) :: KGROUND_LAYER
INTEGER, INTENT(IN) :: KPATCH
INTEGER, INTENT(IN) :: KLW
INTEGER, INTENT(IN) :: KDIMTAB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! Mask and number of grid elements containing patches/tiles:
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GARDEN',0,ZHOOK_HANDLE)
ALLOCATE(XPATCH                  (KLU,KPATCH              ))
ALLOCATE(XVEGTYPE                (KLU,KVEGTYPE            ))
ALLOCATE(XVEGTYPE_PATCH          (KLU,KVEGTYPE,KPATCH     ))
ALLOCATE(NSIZE_NATURE_P          (KPATCH                  ))
ALLOCATE(NR_NATURE_P             (KLU,KPATCH              ))
!
!-------------------------------------------------------------------------------
!
! Averaged Surface radiative parameters:
!
ALLOCATE(XALBNIR_DRY             (KLU                     )) 
ALLOCATE(XALBVIS_DRY             (KLU                     )) 
ALLOCATE(XALBUV_DRY              (KLU                     )) 
ALLOCATE(XALBNIR_WET             (KLU                     )) 
ALLOCATE(XALBVIS_WET             (KLU                     )) 
ALLOCATE(XALBUV_WET              (KLU                     )) 
ALLOCATE(XALBNIR_SOIL            (KLU,KPATCH              )) 
ALLOCATE(XALBVIS_SOIL            (KLU,KPATCH              )) 
ALLOCATE(XALBUV_SOIL             (KLU,KPATCH              )) 
ALLOCATE(XEMIST                  (KLU                     )) 
!
! Subgrid orography parameters
!
!ALLOCATE(XAOSIP                  (KLU                     )) 
!ALLOCATE(XAOSIM                  (KLU                     )) 
!ALLOCATE(XAOSJP                  (KLU                     )) 
!ALLOCATE(XAOSJM                  (KLU                     )) 
!ALLOCATE(XHO2IP                  (KLU                     )) 
!ALLOCATE(XHO2IM                  (KLU                     )) 
!ALLOCATE(XHO2JP                  (KLU                     )) 
!ALLOCATE(XHO2JM                  (KLU                     )) 
ALLOCATE(XZ0EFFIP                (KLU,KPATCH              )) 
ALLOCATE(XZ0EFFIM                (KLU,KPATCH              )) 
ALLOCATE(XZ0EFFJP                (KLU,KPATCH              )) 
ALLOCATE(XZ0EFFJM                (KLU,KPATCH              )) 
ALLOCATE(XZ0EFFJPDIR             (KLU                     )) 
ALLOCATE(XZ0REL                  (KLU                     )) 
!ALLOCATE(XSSO_SLOPE              (KLU                     )) 
!ALLOCATE(XSSO_STDEV              (KLU                     )) 
!
!-------------------------------------------------------------------------------
!
! Input Parameters, per patch:
!
! - vegetation + bare soil:
!
ALLOCATE(XZ0_O_Z0H               (KLU,KPATCH              )) 
ALLOCATE(XALBNIR                 (KLU,KPATCH              )) 
ALLOCATE(XALBVIS                 (KLU,KPATCH              )) 
ALLOCATE(XALBUV                  (KLU,KPATCH              )) 
ALLOCATE(XEMIS                   (KLU,KPATCH              )) 
ALLOCATE(XZ0                     (KLU,KPATCH              )) 
!
! - vegetation:
!
ALLOCATE(XALBNIR_VEG             (KLU,KPATCH              )) 
ALLOCATE(XALBVIS_VEG             (KLU,KPATCH              )) 
ALLOCATE(XALBUV_VEG              (KLU,KPATCH              )) 
!
! - vegetation: default option (Jarvis) and general parameters:
!
ALLOCATE(XVEG                    (KLU,KPATCH              )) 
ALLOCATE(XWRMAX_CF               (KLU,KPATCH              )) 
ALLOCATE(XGAMMA                  (KLU,KPATCH              )) 
ALLOCATE(XCV                     (KLU,KPATCH              )) 
ALLOCATE(XRGL                    (KLU,KPATCH              )) 
ALLOCATE(XRSMIN                  (KLU,KPATCH              )) 
ALLOCATE(XROOTFRAC               (KLU,KGROUND_LAYER,KPATCH)) 
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags parameters ('AGS', 'LAI', 'AST', 'LST', 'NIT' options)
!
ALLOCATE(XABC                    (3                       )) 
ALLOCATE(XPOI                    (3                       )) 
ALLOCATE(XBSLAI                  (KLU,KPATCH              )) 
ALLOCATE(XLAIMIN                 (KLU,KPATCH              )) 
ALLOCATE(XSEFOLD                 (KLU,KPATCH              )) 
ALLOCATE(XH_TREE                 (KLU,KPATCH              )) 
ALLOCATE(XANF                    (KLU,KPATCH              )) 
ALLOCATE(XANMAX                  (KLU,KPATCH              )) 
ALLOCATE(XFZERO                  (KLU,KPATCH              )) 
ALLOCATE(XEPSO                   (KLU,KPATCH              )) 
ALLOCATE(XGAMM                   (KLU,KPATCH              )) 
ALLOCATE(XQDGAMM                 (KLU,KPATCH              )) 
ALLOCATE(XGMES                   (KLU,KPATCH              )) 
ALLOCATE(XRE25                   (KLU,KPATCH              )) 
ALLOCATE(XQDGMES                 (KLU,KPATCH              )) 
ALLOCATE(XT1GMES                 (KLU,KPATCH              )) 
ALLOCATE(XT2GMES                 (KLU,KPATCH              )) 
ALLOCATE(XAMAX                   (KLU,KPATCH              )) 
ALLOCATE(XQDAMAX                 (KLU,KPATCH              )) 
ALLOCATE(XT1AMAX                 (KLU,KPATCH              )) 
ALLOCATE(XT2AMAX                 (KLU,KPATCH              )) 
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Stress parameters ('AST', 'LST', 'NIT' options)
!
ALLOCATE(LSTRESS                 (KLU,KPATCH              )) 
ALLOCATE(XF2I                    (KLU,KPATCH              )) 
ALLOCATE(XGC                     (KLU,KPATCH              )) 
ALLOCATE(XAH                     (KLU,KPATCH              )) 
ALLOCATE(XBH                     (KLU,KPATCH              )) 
ALLOCATE(XDMAX                   (KLU,KPATCH              )) 
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Nitrogen-model parameters ('NIT' option)
!
ALLOCATE(XCE_NITRO               (KLU,KPATCH              )) 
ALLOCATE(XCF_NITRO               (KLU,KPATCH              )) 
ALLOCATE(XCNA_NITRO              (KLU,KPATCH              )) 
ALLOCATE(XBSLAI_NITRO            (KLU,KPATCH              )) 
!
!-------------------------------------------------------------------------------
!
! - soil: primary parameters
!
!ALLOCATE(XSAND                   (KLU,KGROUND_LAYER       )) 
!ALLOCATE(XCLAY                   (KLU,KGROUND_LAYER       )) 
!ALLOCATE(XRUNOFFB                (KLU                     )) 
!ALLOCATE(XWDRAIN                 (KLU                     )) 
ALLOCATE(XTAUICE                 (KLU                     )) 
ALLOCATE(XGAMMAT                 (KLU                     )) 
ALLOCATE(XDG                     (KLU,KGROUND_LAYER,KPATCH)) 
ALLOCATE(XRUNOFFD                (KLU,KPATCH              )) 
!
!-------------------------------------------------------------------------------
!
! - soil: Secondary parameters: hydrology
!
ALLOCATE(XC1SAT                  (KLU,KPATCH              )) 
ALLOCATE(XC2REF                  (KLU,KPATCH              )) 
ALLOCATE(XC3                     (KLU,2,KPATCH            )) 
ALLOCATE(XC4B                    (KLU                     )) 
ALLOCATE(XC4REF                  (KLU,KPATCH              )) 
ALLOCATE(XACOEF                  (KLU                     )) 
ALLOCATE(XPCOEF                  (KLU                     )) 
ALLOCATE(XWFC                    (KLU,KGROUND_LAYER       )) 
ALLOCATE(XWWILT                  (KLU,KGROUND_LAYER       )) 
ALLOCATE(XWSAT                   (KLU,KGROUND_LAYER       )) 
ALLOCATE(XBCOEF                  (KLU,KGROUND_LAYER       )) 
ALLOCATE(XCONDSAT                (KLU,KGROUND_LAYER,KPATCH))
ALLOCATE(XCONDSAT_EXP            (KLU,KGROUND_LAYER,KPATCH))
ALLOCATE(XEXP_DIF                (KLU,KGROUND_LAYER,KPATCH))
ALLOCATE(XMPOTSAT                (KLU,KGROUND_LAYER       )) 
!
!-------------------------------------------------------------------------------
!
! - soil: Secondary parameters: thermal 
!
ALLOCATE(XCGSAT                  (KLU                     )) 
ALLOCATE(XHCAPSOIL               (KLU,KGROUND_LAYER       )) 
ALLOCATE(XCONDDRY                (KLU,KGROUND_LAYER       )) 
ALLOCATE(XCONDSLD                (KLU,KGROUND_LAYER       )) 
ALLOCATE(XTDEEP                  (KLU                     )) 
ALLOCATE(XGAMMAT                 (KLU                     )) 
!
!
! Prognostic variables:
!
! - Soil and vegetation heat and water:
!
ALLOCATE(XPCPS                   (KLU,KPATCH              )) 
ALLOCATE(XPLVTT                  (KLU,KPATCH              )) 
ALLOCATE(XPLSTT                  (KLU,KPATCH              )) 
!
! - Vegetation: Ags Prognostic (YPHOTO = ('LAI', 'LST', or 'NIT') or prescribed (YPHOTO='NON', 'AGS' or 'LST')
!
ALLOCATE(XLAI                    (KLU,KPATCH              )) 
!
! - SGH scheme
!                                   
ALLOCATE(XD_ICE                  (KLU,KPATCH              )) 
ALLOCATE(XKSAT_ICE               (KLU,KPATCH              )) 
!
! - Irrigation, seeding and reaping
!
ALLOCATE(TSEED                   (KLU,KPATCH              )) 
ALLOCATE(TREAP                   (KLU,KPATCH              )) 
ALLOCATE(XWATSUP                 (KLU,KPATCH              )) 
ALLOCATE(XIRRIG                  (KLU,KPATCH              ))
!
!-------------------------------------------------------------------------------
!
IF (HINIT=='ALL') THEN
!
! Prognostic variables:
!
! - Soil and vegetation heat and water:
!
ALLOCATE(XWR                     (KLU,KPATCH              )) 
ALLOCATE(XTG                     (KLU,KGROUND_LAYER,KPATCH)) 
ALLOCATE(XWG                     (KLU,KGROUND_LAYER,KPATCH)) 
ALLOCATE(XWGI                    (KLU,KGROUND_LAYER,KPATCH)) 
ALLOCATE(XRESA                   (KLU,KPATCH              )) 
!
! - Vegetation: Ags Prognostic (YPHOTO = 'LAI', 'LST', 'AGS' or 'LST')
!
ALLOCATE(XAN                     (KLU,KPATCH              )) 
ALLOCATE(XANDAY                  (KLU,KPATCH              )) 
ALLOCATE(XANFM                   (KLU,KPATCH              )) 
ALLOCATE(XLE                     (KLU,KPATCH              ))
!
! - Vegetation (Ags 'NIT' 'NCB' option):
!
ALLOCATE(XBIOMASS                (KLU,NNBIOMASS,KPATCH    ))
ALLOCATE(XRESP_BIOMASS           (KLU,NNBIOMASS,KPATCH    ))
!
! - Vegetation (ISBA-CC, YRESPSL = 'CNT'):
!
ALLOCATE(XLITTER                 (KLU,NNLITTER,NNLITTLEVS,KPATCH))
ALLOCATE(XSOILCARB               (KLU,NNSOILCARB,KPATCH))
ALLOCATE(XLIGNIN_STRUC           (KLU,NNLITTLEVS,KPATCH))
!
!-------------------------------------------------------------------------------
!
! - Snow and flood fractions and total albedo at time t:
!
ALLOCATE(XPSN                    (KLU,NPATCH              ))
ALLOCATE(XPSNG                   (KLU,NPATCH              ))
ALLOCATE(XPSNV                   (KLU,NPATCH              ))
ALLOCATE(XPSNV_A                 (KLU,NPATCH              ))
!
ALLOCATE(LIRRIGATE               (KLU,KPATCH              ))
ALLOCATE(LIRRIDAY                (KLU,KPATCH              ))
ALLOCATE(NIRRINUM                (KLU,KPATCH              ))
ALLOCATE(XTHRESHOLDSPT           (KLU,KPATCH              ))
!
ENDIF
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GARDEN',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOCATE_TEB_GARDEN
