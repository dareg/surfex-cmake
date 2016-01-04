!     #########
      SUBROUTINE WRITESURF_PGD_ISBA_n (DGU, &
                                        DTI, DTZ, IG, I, U, &
                                       HPROGRAM)
!     ################################################
!
!!****  *WRITESURF_PGD_ISBA_n* - writes ISBA physiographic fields
!!                        
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003 
!!      P. Le Moigne 12/2004 : add type of photosynthesis 
!!      B. Decharme  06/2009 : add topographic index statistics
!!      A.L. Gibelin 04/2009 : dimension NBIOMASS for ISBA-A-gs
!!      B. Decharme  07/2011 : delete argument HWRITE
!!      B. Decharme  07/2012 : files of data for permafrost area and for SOC top and sub soil
!!                   11/2013 : same for groundwater distribution
!!                   11/2014 : Write XSOILGRID as a series of real 
!!      P. Samuelsson 10/2014 : MEB
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
!
!
#ifdef SFX_OL
USE MODN_IO_OFFLINE, ONLY : LWR_VEGTYPE
#endif
!
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_t
USE MODD_ISBA_GRID_n, ONLY : ISBA_GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
!
USE MODE_WRITE_SURF_COV, ONLY : WRITE_SURF_COV
!
USE MODI_WRITE_SURF
USE MODI_WRITE_GRID
USE MODI_WRITESURF_PGD_ISBA_PAR_n
USE MODI_WRITESURF_PGD_TSZ0_PAR_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
!
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
!
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: DTZ
TYPE(ISBA_GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
CHARACTER(LEN=4 ) :: YLVL
!
INTEGER :: JJ, JLAYER
INTEGER :: ISIZE_LMEB_PATCH  ! Number of patches with MEB=true
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!-------------------------------------------------------------------------------
!
!
!* soil scheme option
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_ISBA_N',0,ZHOOK_HANDLE)
YRECFM='ISBA'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%CISBA,IRESP,HCOMMENT=YCOMMENT)
!
!* Pedo-transfert function
!
YRECFM='PEDOTF'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%CPEDOTF,IRESP,HCOMMENT=YCOMMENT)
!
!* type of photosynthesis
!
YRECFM='PHOTO'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%CPHOTO,IRESP,HCOMMENT=YCOMMENT)
!
!* new radiative transfert
!
YRECFM='TR_ML'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%LTR_ML,IRESP,HCOMMENT=YCOMMENT)
!
!* threshold to remove little fractions of patches
!
YRECFM='RM_PATCH'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%XRM_PATCH,IRESP,HCOMMENT=YCOMMENT)

!* number of soil layers
!
YRECFM='GROUND_LAYER'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%NGROUND_LAYER,IRESP,HCOMMENT=YCOMMENT)
!
!* Reference grid for DIF
!
IF(I%O%CISBA=='DIF') THEN
  DO JLAYER=1,I%O%NGROUND_LAYER
     WRITE(YLVL,'(I4)') JLAYER     
     YRECFM='SOILGRID'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     YCOMMENT='Depth of ISBA soilgrid layer '//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%XSOILGRID(JLAYER),IRESP,HCOMMENT=YCOMMENT)
  END DO 
ENDIF
!
!* number of biomass pools
!
YRECFM='NBIOMASS'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%NNBIOMASS,IRESP,HCOMMENT=YCOMMENT)
!
!* number of tiles
!
YRECFM='PATCH_NUMBER'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%NPATCH,IRESP,HCOMMENT=YCOMMENT)
!
!* flag indicating if fields are computed from ecoclimap or not
!
YRECFM='ECOCLIMAP'
YCOMMENT=YRECFM
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%LECOCLIMAP,IRESP,HCOMMENT=YCOMMENT)
!
!* logical vector indicating for which patches MEB should be applied
!
YRECFM='MEB_PATCH'
YCOMMENT='(LOGICAL LIST)'
CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%LMEB_PATCH(:),IRESP,HCOMMENT=YCOMMENT,HDIR='-')
!
ISIZE_LMEB_PATCH = COUNT(I%O%LMEB_PATCH(:))
!
IF (ISIZE_LMEB_PATCH>0)THEN
!
!* flag indicating if forcing is from observed measurements or not
!
   YRECFM='FORC_MEASURE'
   YCOMMENT=YRECFM
   CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%LFORC_MEASURE,IRESP,HCOMMENT=YCOMMENT)
!
!* flag indicating if litter layer is used or not
!
   YRECFM='MEB_LITTER'
   YCOMMENT=YRECFM
   CALL WRITE_SURF(DGU, U, &
                HPROGRAM,YRECFM,I%O%LMEB_LITTER,IRESP,HCOMMENT=YCOMMENT)
!
ENDIF
!
!*       2.     Physiographic data fields:
!               -------------------------
!
!* cover classes
!
YRECFM='COVER_LIST'
YCOMMENT='(LOGICAL LIST)'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%LCOVER(:),IRESP,HCOMMENT=YCOMMENT,HDIR='-')
!
#ifdef SFX_OL
IF (LWR_VEGTYPE) THEN
  YRECFM='VEGTYPE'
  YCOMMENT='(X_Y_VEGTYPE)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,I%M%X%XVEGTYPE,IRESP,HCOMMENT=YCOMMENT)
ENDIF
#endif
!
!* orography
!
YRECFM='ZS'
YCOMMENT='ZS'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XZS(:),IRESP,HCOMMENT=YCOMMENT)
!
!* latitude, longitude
!
 CALL WRITE_GRID(DGU, U, &
                 HPROGRAM,IG%CGRID,IG%XGRID_PAR,IG%XLAT,IG%XLON,IG%XMESH_SIZE,IRESP,I%P%XZ0EFFJPDIR)
!
!
!* clay fraction
!
!
YRECFM='CLAY'
YCOMMENT='X_Y_CLAY'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XCLAY(:,1),IRESP,HCOMMENT=YCOMMENT)
!
!* sand fraction
!
YRECFM='SAND'
YCOMMENT='X_Y_SAND'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XSAND(:,1),IRESP,HCOMMENT=YCOMMENT)
!
!* soil organic carbon
!
YRECFM='SOCP'
YCOMMENT=''
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%LSOCP,IRESP,HCOMMENT=YCOMMENT)
!
IF(I%O%LSOCP)THEN
  !        
  YCOMMENT='X_Y_SOC'
  YRECFM='SOC_TOP'
  CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XSOC(:,1),IRESP,HCOMMENT=YCOMMENT)
  YRECFM='SOC_SUB'
  CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XSOC(:,2),IRESP,HCOMMENT=YCOMMENT)
  !
ENDIF
!
!* permafrost distribution
!
YRECFM='PERMAFROST'
YCOMMENT=''
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%LPERM,IRESP,HCOMMENT=YCOMMENT)
!
IF(I%O%LPERM)THEN
  YCOMMENT='X_Y_PERM'
  YRECFM='PERM'
  CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XPERM(:),IRESP,HCOMMENT=YCOMMENT)
ENDIF
!
!* groundwater distribution
!
YRECFM='GWKEY'
YCOMMENT=''
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%LGW,IRESP,HCOMMENT=YCOMMENT)
!
IF(I%O%LGW)THEN
  YCOMMENT='X_Y_GWFRAC'
  YRECFM='GWFRAC'
  CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XGW(:),IRESP,HCOMMENT=YCOMMENT)
ENDIF
!
!SOILNOX
!
YRECFM='NO'
YCOMMENT=''
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%LNOF,IRESP,HCOMMENT=YCOMMENT)
!
IF (I%O%LNOF) THEN
  !
  YRECFM='PH'
  YCOMMENT='X_Y_PH'
  CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XPH(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='FERT'
  YCOMMENT='X_Y_FERT'
  CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XFERT(:),IRESP,HCOMMENT=YCOMMENT)
  !
ENDIF
!
!* subgrid-scale orography parameters to compute dynamical roughness length
!
YRECFM='AOSIP'
YCOMMENT='X_Y_AOSIP'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XAOSIP,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='AOSIM'
YCOMMENT='X_Y_AOSIM'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XAOSIM,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='AOSJP'
YCOMMENT='X_Y_AOSJP'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XAOSJP,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='AOSJM'
YCOMMENT='X_Y_AOSJM'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XAOSJM,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='HO2IP'
YCOMMENT='X_Y_HO2IP'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XHO2IP,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='HO2IM'
YCOMMENT='X_Y_HO2IM'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XHO2IM,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='HO2JP'
YCOMMENT='X_Y_HO2JP'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XHO2JP,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='HO2JM'
YCOMMENT='X_Y_HO2JM'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XHO2JM,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='SSO_SLOPE'
YCOMMENT='X_Y_SSO_SLOPE (-)'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XSSO_SLOPE,IRESP,HCOMMENT=YCOMMENT)
!
!* orographic runoff coefficient
!
YRECFM='RUNOFFB'
YCOMMENT='X_Y_RUNOFFB'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XRUNOFFB,IRESP,HCOMMENT=YCOMMENT)
!
!* subgrid drainage coefficient
!
YRECFM='WDRAIN'
YCOMMENT='X_Y_WDRAIN'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XWDRAIN,IRESP,HCOMMENT=YCOMMENT)
!
!* topographic index statistics
!
YRECFM='CTI'
YCOMMENT=''
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%O%LCTI,IRESP,HCOMMENT=YCOMMENT)
!
IF(I%O%LCTI)THEN
!
YRECFM='TI_MIN'
YCOMMENT='X_Y_TI_MIN'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XTI_MIN,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='TI_MAX'
YCOMMENT='X_Y_TI_MAX'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XTI_MAX,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='TI_MEAN'
YCOMMENT='X_Y_TI_MEAN'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XTI_MEAN,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='TI_STD'
YCOMMENT='X_Y_TI_STD'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XTI_STD,IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='TI_SKEW'
YCOMMENT='X_Y_TI_SKEW'
 CALL WRITE_SURF(DGU, U, &
                 HPROGRAM,YRECFM,I%P%XTI_SKEW,IRESP,HCOMMENT=YCOMMENT)
!
ENDIF
!
!-------------------------------------------------------------------------------
 CALL WRITESURF_PGD_ISBA_PAR_n(DGU, U, &
                               DTI, &
                               HPROGRAM)
IF (U%CNATURE=='TSZ0') CALL WRITESURF_PGD_TSZ0_PAR_n(DGU, U, &
                                                     DTZ, &
                                                     HPROGRAM)
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_PGD_ISBA_n
