!     #########
      SUBROUTINE READ_PGD_ISBA_n (CHI, DTCO, DTI, DTZ, DGU, GB, IG, I, &
                                  UG, U, USS, SV, &
                                  HPROGRAM,OLAND_USE)
!     #########################################
!
!!****  *READ_PGD_ISBA_n* - routine to initialise ISBA physiographic variables 
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
!!      P. Le Moigne  12/2004 : add type of photosynthesis
!!      B. Decharme      2008 : add XWDRAIN
!!      B. Decharme   06/2009 : add topographic index statistics
!!      A.L. Gibelin 04/2009 : dimension NBIOMASS for ISBA-A-gs
!!      B. Decharme  07/2012  : files of data for permafrost area and for SOC top and sub soil
!!                   11/2013  : same for groundwater distribution
!!                   11/2014  : Read XSOILGRID as a series of real 
!!      P. Samuelsson 10/2014 : MEB
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_ISBA_GRID_n, ONLY : ISBA_GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
USE MODD_SV_n, ONLY : SV_t
!
USE MODD_TYPE_DATE_SURF
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_ISBA_PAR,    ONLY : XOPTIMGRID
!
USE MODE_READ_SURF_COV, ONLY : READ_SURF_COV
!
USE MODI_INIT_IO_SURF_n
USE MODI_END_IO_SURF_n
!
USE MODI_READ_SURF
USE MODI_PACK_INIT
USE MODI_PACK_SSO
USE MODI_READ_LCOVER
USE MODI_READ_PGD_ISBA_PAR_n
USE MODI_READ_PGD_TSZ0_PAR_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_TYPE_DIM_n
USE MODI_READ_LECOCLIMAP
!
USE MODI_ABOR1_SFX
!
USE MODI_GET_LUOUT
USE MODI_PACK_SAME_RANK
USE MODI_GET_SURF_MASK_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: DTZ
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(ISBA_GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
TYPE(SV_t), INTENT(INOUT) :: SV
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! calling program
LOGICAL,           INTENT(IN)  :: OLAND_USE ! 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER, DIMENSION(:), POINTER :: IMASK  ! mask for packing from complete field to nature field
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK
!
CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=4 ) :: YLVL
!
INTEGER :: ISIZE_LMEB_PATCH  ! Number of patches with MEB=true
INTEGER :: ILU    ! expected physical size of full surface array
INTEGER :: ILUOUT ! output listing logical unit
INTEGER :: IRESP  ! Error code after redding
INTEGER :: JLAYER ! loop counter on layers
INTEGER :: IVERSION, IBUGFIX   ! surface version
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_PGD_ISBA_N',0,ZHOOK_HANDLE)
YRECFM='SIZE_NATURE'
 CALL GET_TYPE_DIM_n(DTCO, U, &
                     'NATURE',IG%NDIM)
!
YRECFM='VERSION'
 CALL READ_SURF(&
                HPROGRAM,YRECFM,IVERSION,IRESP)
!
YRECFM='BUG'
 CALL READ_SURF(&
                HPROGRAM,YRECFM,IBUGFIX,IRESP)
!
!*       2.     Dimension initializations:
!               -------------------------
!
!* soil scheme
!
YRECFM='ISBA'
 CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%CISBA,IRESP)
!
IF (IVERSION>=7) THEN
  !
  !* Pedo-transfert function
  !
  YRECFM='PEDOTF'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%CPEDOTF,IRESP)
  !
ELSE
  I%O%CPEDOTF = 'CH78'
ENDIF
!
!* type of photosynthesis
!
YRECFM='PHOTO'
 CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%CPHOTO,IRESP)
!
!* new radiative transfert
!
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=2) THEN
  !
  YRECFM='TR_ML'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%LTR_ML,IRESP)
  !
ELSE 
  I%O%LTR_ML = .FALSE.
ENDIF
!
!* threshold to remove little fractions of patches
!
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
  !
  YRECFM='RM_PATCH'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%XRM_PATCH,IRESP)
  !
ELSE 
  I%O%XRM_PATCH = 0.0
ENDIF
!
!* number of soil layers
!
YRECFM='GROUND_LAYER'
 CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%NGROUND_LAYER,IRESP)
!
!* Reference grid for DIF
!
IF(I%O%CISBA=='DIF') THEN
  ALLOCATE(I%O%XSOILGRID(I%O%NGROUND_LAYER))
  I%O%XSOILGRID=XUNDEF
  IF (IVERSION>=8) THEN
     DO JLAYER=1,I%O%NGROUND_LAYER
        WRITE(YLVL,'(I4)') JLAYER
        YRECFM='SOILGRID'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
        CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%XSOILGRID(JLAYER),IRESP)
     ENDDO    
  ELSEIF (IVERSION==7 .AND. IBUGFIX>=2) THEN
    YRECFM='SOILGRID'
    CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%XSOILGRID,IRESP,HDIR='-')
  ELSE
    I%O%XSOILGRID(1:I%O%NGROUND_LAYER)=XOPTIMGRID(1:I%O%NGROUND_LAYER)
  ENDIF
ELSE
  ALLOCATE(I%O%XSOILGRID(0))
ENDIF
!
!* number of biomass pools
!
IF (IVERSION>=6) THEN
  YRECFM='NBIOMASS'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%NNBIOMASS,IRESP)
ELSE
  SELECT CASE (I%O%CPHOTO)
    CASE ('AGS','LAI','AST','LST')
      I%O%NNBIOMASS = 1
    CASE ('NIT')
      I%O%NNBIOMASS = 3
    CASE ('NCB')
      I%O%NNBIOMASS = 6
  END SELECT
ENDIF
!
!* number of tiles
!
YRECFM='PATCH_NUMBER'
 CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%NPATCH,IRESP)
!
!* logical vector indicating for which patches MEB should be applied
!
ALLOCATE(I%O%LMEB_PATCH(I%O%NPATCH))
!
IF (IVERSION>=8) THEN
!
   YRECFM='MEB_PATCH'
   CALL READ_SURF(HPROGRAM,YRECFM,I%O%LMEB_PATCH(:),IRESP,HDIR='-')
!
   ISIZE_LMEB_PATCH = COUNT(I%O%LMEB_PATCH(:))
!
   IF (ISIZE_LMEB_PATCH>0)THEN
      YRECFM='FORC_MEASURE'
      CALL READ_SURF(HPROGRAM,YRECFM,I%O%LFORC_MEASURE,IRESP)
      YRECFM='MEB_LITTER'
      CALL READ_SURF(HPROGRAM,YRECFM,I%O%LMEB_LITTER,IRESP)
   ELSE
      I%O%LFORC_MEASURE=.FALSE.
      I%O%LMEB_LITTER  =.FALSE.           
   ENDIF
!
ELSE
   I%O%LMEB_PATCH(:)=.FALSE.
   I%O%LFORC_MEASURE=.FALSE.
   I%O%LMEB_LITTER  =.FALSE.
ENDIF
!
!
!*       3.     Physiographic data fields:
!               -------------------------
!
!
!*       3.1    Cover classes :
!               -------------
!
ALLOCATE(I%P%LCOVER(JPCOVER))
ALLOCATE(I%P%XZS(IG%NDIM))
ALLOCATE(IG%XLAT       (IG%NDIM))
ALLOCATE(IG%XLON       (IG%NDIM))
ALLOCATE(IG%XMESH_SIZE (IG%NDIM))
ALLOCATE(I%P%XZ0EFFJPDIR(IG%NDIM))
!
CALL PACK_INIT(DTCO,U,UG,HPROGRAM,'NATURE',IG%CGRID,IG%XGRID_PAR,     &
               I%P%LCOVER,I%P%XCOVER,I%P%XZS,IG%XLAT,IG%XLON,IG%XMESH_SIZE, I%P%XZ0EFFJPDIR )
!
!* clay fraction : attention, seul un niveau est present dans le fichier
!* on rempli tout les niveaux de  XCLAY avec les valeurs du fichiers
!
ALLOCATE(I%P%XCLAY(IG%NDIM,I%O%NGROUND_LAYER))
YRECFM='CLAY'
 CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XCLAY(:,1),IRESP)
DO JLAYER=2,I%O%NGROUND_LAYER
  I%P%XCLAY(:,JLAYER)=I%P%XCLAY(:,1)
END DO
!
!* sand fraction
!
ALLOCATE(I%P%XSAND(IG%NDIM,I%O%NGROUND_LAYER))
YRECFM='SAND'
 CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XSAND(:,1),IRESP)
DO JLAYER=2,I%O%NGROUND_LAYER
  I%P%XSAND(:,JLAYER)=I%P%XSAND(:,1)
END DO
!
!* Soil organic carbon profile
!
IF (IVERSION>7 .OR. (IVERSION==7 .AND. IBUGFIX>=3)) THEN
   YRECFM='SOCP'
   CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%LSOCP,IRESP)
ELSE
   I%O%LSOCP=.FALSE.
ENDIF
!
IF(I%O%LSOCP)THEN
!  
  ALLOCATE(I%P%XSOC (IG%NDIM,I%O%NGROUND_LAYER))
!
  YRECFM='SOC_TOP'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XSOC(:,1),IRESP)
  YRECFM='SOC_SUB'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XSOC(:,2),IRESP)
!
  DO JLAYER=2,I%O%NGROUND_LAYER
    I%P%XSOC (:,JLAYER)=I%P%XSOC (:,2)
  END DO
!
ELSE
!  
  ALLOCATE(I%P%XSOC (0,1))
!
ENDIF
!
!* permafrost distribution
!
IF (IVERSION>7 .OR. (IVERSION==7 .AND. IBUGFIX>=3)) THEN
   YRECFM='PERMAFROST'
   CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%LPERM,IRESP)
ELSE
   I%O%LPERM=.FALSE.
ENDIF
!
IF(I%O%LPERM)THEN
!  
  ALLOCATE(I%P%XPERM (IG%NDIM))
!
  YRECFM='PERM'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XPERM(:),IRESP)
!
ELSE
!  
  ALLOCATE(I%P%XPERM (0))
!
ENDIF
!
!* groundwater distribution
!
IF (IVERSION>=8) THEN
   YRECFM='GWKEY'
   CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%LGW,IRESP)
ELSE
   I%O%LGW=.FALSE.
ENDIF
!
IF(I%O%LGW)THEN
!  
  ALLOCATE(I%P%XGW (IG%NDIM))
!
  YRECFM='GWFRAC'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XGW(:),IRESP)
!
ELSE
!  
  ALLOCATE(I%P%XGW (0))
!
ENDIF
!
IF (IVERSION>7 .OR. (IVERSION==7 .AND. IBUGFIX>=3)) THEN
   YRECFM='NO'
   CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%LNOF,IRESP)
ELSE
   I%O%LNOF = .FALSE.
ENDIF
!
!SOILNOX
!
IF (CHI%LCH_NO_FLUX) THEN
  !
  IF (I%O%LNOF) THEN
    !
    ALLOCATE(I%P%XPH(IG%NDIM))
    YRECFM='PH'
    CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XPH(:),IRESP)
    !
    ALLOCATE(I%P%XFERT(IG%NDIM))
    YRECFM='FERT'
    CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XFERT(:),IRESP)
    !
  ELSE
    CALL ABOR1_SFX("READ_PGD_ISBAn: WITH LCH_NO_FLUX=T, PH AND FERT FIELDS ARE GIVEN AT PGD STEP")
  ENDIF
  !
ELSE
  ALLOCATE(I%P%XPH (0))
  ALLOCATE(I%P%XFERT(0))
END IF
!
!* subgrid-scale orography parameters to compute dynamical roughness length
!
ALLOCATE(I%P%XAOSIP(IG%NDIM))
ALLOCATE(I%P%XAOSIM(IG%NDIM))
ALLOCATE(I%P%XAOSJP(IG%NDIM))
ALLOCATE(I%P%XAOSJM(IG%NDIM))
ALLOCATE(I%P%XHO2IP(IG%NDIM))
ALLOCATE(I%P%XHO2IM(IG%NDIM))
ALLOCATE(I%P%XHO2JP(IG%NDIM))
ALLOCATE(I%P%XHO2JM(IG%NDIM))
ALLOCATE(I%P%XSSO_SLOPE(IG%NDIM))
ALLOCATE(I%P%XSSO_STDEV(IG%NDIM))
!
 CALL PACK_SSO(USS,HPROGRAM,U%NR_NATURE, I%P%XAOSIP, I%P%XAOSIM, I%P%XAOSJP, I%P%XAOSJM,  &
               I%P%XHO2IP, I%P%XHO2IM, I%P%XHO2JP, I%P%XHO2JM, I%P%XSSO_SLOPE,I%P%XSSO_STDEV)
!
!* orographic runoff coefficient
!
ALLOCATE(I%P%XRUNOFFB(IG%NDIM))
YRECFM='RUNOFFB'
 CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XRUNOFFB,IRESP)
!
!* subgrid drainage coefficient
!
ALLOCATE(I%P%XWDRAIN(IG%NDIM))
IF (IVERSION<=3) THEN
  I%P%XWDRAIN = 0.
ELSE
  YRECFM='WDRAIN'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XWDRAIN,IRESP)
ENDIF
!
!* topographic index statistics
!
IF(I%O%CRUNOFF=='SGH ' .AND. IVERSION>=5) THEN 
!
  YRECFM='CTI'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%O%LCTI,IRESP)        
!
  IF (.NOT.I%O%LCTI) CALL ABOR1_SFX("READ_PGD_ISBA_n:WITH CRUNOFF=SGH, CTI MAPS MUST BE GIVEN TO PGD")
  !
  ALLOCATE(I%P%XTI_MIN(IG%NDIM))
  ALLOCATE(I%P%XTI_MAX(IG%NDIM))
  ALLOCATE(I%P%XTI_MEAN(IG%NDIM))
  ALLOCATE(I%P%XTI_STD(IG%NDIM))
  ALLOCATE(I%P%XTI_SKEW(IG%NDIM))
!
  YRECFM='TI_MIN'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XTI_MIN,IRESP)
!
  YRECFM='TI_MAX'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XTI_MAX,IRESP)
!
  YRECFM='TI_MEAN'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XTI_MEAN,IRESP)
!
  YRECFM='TI_STD'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XTI_STD,IRESP)
!
  YRECFM='TI_SKEW'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,I%P%XTI_SKEW,IRESP)
!
ELSE
!
  ALLOCATE(I%P%XTI_MIN(0))
  ALLOCATE(I%P%XTI_MAX(0))
  ALLOCATE(I%P%XTI_MEAN(0))
  ALLOCATE(I%P%XTI_STD(0))
  ALLOCATE(I%P%XTI_SKEW(0))
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!* biogenic chemical emissions
!
IF (CHI%LCH_BIO_FLUX) THEN
  ALLOCATE(ZWORK(U%NSIZE_FULL,1))
  !
  CALL END_IO_SURF_n(HPROGRAM)
CALL INIT_IO_SURF_n(DTCO, DGU, U, &
                      HPROGRAM,'FULL  ','SURF  ','READ ')
  !
  CALL GET_LUOUT(HPROGRAM,ILUOUT)
  ALLOCATE(IMASK(IG%NDIM))
  ILU=0
  CALL GET_SURF_MASK_n(DTCO, U, &
                       'NATURE',IG%NDIM,IMASK,ILU,ILUOUT)
  ALLOCATE(GB%XISOPOT(IG%NDIM))
  ALLOCATE(GB%XMONOPOT(IG%NDIM))
  !
  ZWORK(:,:) = 0.  
  YRECFM='E_ISOPOT'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,ZWORK,IRESP)
  CALL PACK_SAME_RANK(IMASK,ZWORK(:,1),GB%XISOPOT(:))
  !
  ZWORK(:,:) = 0.  
  YRECFM='E_MONOPOT'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,ZWORK,IRESP)
  CALL PACK_SAME_RANK(IMASK,ZWORK(:,1),GB%XMONOPOT(:))
  !
  CALL END_IO_SURF_n(HPROGRAM)
CALL INIT_IO_SURF_n(DTCO, DGU, U, &
                      HPROGRAM,'NATURE','ISBA  ','READ ')
  !
  DEALLOCATE(ZWORK)
ELSE
  ALLOCATE(GB%XISOPOT (0))
  ALLOCATE(GB%XMONOPOT(0))
END IF
!
!-------------------------------------------------------------------------------
!
!*       4.     Physiographic data fields not to be computed by ecoclimap
!               ---------------------------------------------------------
!
 CALL READ_LECOCLIMAP(&
                      HPROGRAM,I%O%LECOCLIMAP)
!
 CALL READ_PGD_ISBA_PAR_n(DTCO, U, &
                          DTI, IG, I%O, &
                          HPROGRAM,IG%NDIM,OLAND_USE)
IF (U%CNATURE == 'TSZ0') CALL READ_PGD_TSZ0_PAR_n(&
                                                  DTZ, &
                                                  HPROGRAM)
!
IF (LHOOK) CALL DR_HOOK('READ_PGD_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_PGD_ISBA_n
