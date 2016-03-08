!     #########
      SUBROUTINE READ_PGD_TEB_GREENROOF_PAR_n (DGR, GRO, GRP, KDIM, HPROGRAM)
!     ################################################
!
!!****  *READ_PGD_TEB_GREENROOF_PAR_n* - reads ISBA physiographic fields
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
!!      C. de Munck  02/2012 : added parameterisation for sedum species under NVT_TROG 
!-------------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
!
USE MODD_CSTS,                 ONLY : XDAY
USE MODD_SURF_PAR,             ONLY : XUNDEF
USE MODD_DATA_COVER_PAR,       ONLY : NVEGTYPE, NVT_GRAS, NVT_TROG
!paramètres ci-dessus à initialiser pour les GR (sauf XPAR_OM_GR, XPAR_SAND_GR, XPAR_CLAY_GR qui sont lues) 
USE MODD_PREP_TEB_GREENROOF,   ONLY : NGRID_LEVEL, XGRID_SOIL
!
USE MODI_READ_SURF
USE MODI_VEG_FROM_LAI
USE MODI_Z0V_FROM_LAI
USE MODI_EMIS_FROM_VEG
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
!
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DGR
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GRO
TYPE(ISBA_PGD_t), INTENT(INOUT) :: GRP
INTEGER, INTENT(IN) :: KDIM
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER                               :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12)                     :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100)                    :: YCOMMENT       ! Comment string
INTEGER                               :: JI             ! loop index
INTEGER                               :: JTIME          ! loop index
INTEGER                               :: JLAYER         ! loop index
!
LOGICAL :: GAGRI_TO_GRASS
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.    Reading of PGD file
!              --------------------
!
IF (LHOOK) CALL DR_HOOK('READ_PGD_TEB_GREENROOF_PAR_N',0,ZHOOK_HANDLE)
!
GAGRI_TO_GRASS=.FALSE.
!
YRECFM='GR_NTIME'
 CALL READ_SURF(HPROGRAM,YRECFM,DGR%NTIME,IRESP)
!
YRECFM='GR_LAYER'
 CALL READ_SURF( HPROGRAM,YRECFM,GRO%NGROUND_LAYER,IRESP)
!
! Read type of green roof
YRECFM='D_TYPE_GR'
 CALL READ_SURF(HPROGRAM,YRECFM,GRO%CTYP_COV,IRESP)
!
! Read green roof OM fraction
ALLOCATE(GRP%XSOC    (KDIM,GRO%NGROUND_LAYER))
DO JLAYER=1,GRO%NGROUND_LAYER
  !WRITE(YRECFM,FMT='(A8,I1.1)') 'D_OM_GR0',JLAYER
  WRITE(YRECFM,FMT='(A7,I2.2)') 'D_OM_GR',JLAYER
  CALL READ_SURF(HPROGRAM,YRECFM,GRP%XSOC(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
! Read green roof SAND fraction
ALLOCATE(GRP%XSAND   (KDIM,GRO%NGROUND_LAYER))
DO JLAYER=1,GRO%NGROUND_LAYER
  !WRITE(YRECFM,FMT='(A10,I1.1)') 'D_SAND_GR0',JLAYER
  WRITE(YRECFM,FMT='(A9,I2.2)') 'D_SAND_GR',JLAYER
  CALL READ_SURF(HPROGRAM,YRECFM,GRP%XSAND(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
! Read green roof CLAY fraction
ALLOCATE(GRP%XCLAY   (KDIM,GRO%NGROUND_LAYER))
DO JLAYER=1,GRO%NGROUND_LAYER
  !WRITE(YRECFM,FMT='(A10,I1.1)') 'D_CLAY_GR0',JLAYER
  WRITE(YRECFM,FMT='(A9,I2.2)') 'D_CLAY_GR',JLAYER
  CALL READ_SURF(HPROGRAM,YRECFM,GRP%XCLAY(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
! Read green roof LAI
ALLOCATE(DGR%XPAR_LAI    (KDIM,DGR%NTIME,1))
DO JTIME=1,DGR%NTIME
  WRITE(YRECFM,FMT='(A8,I2.2)') 'D_LAI_GR',JTIME
  CALL READ_SURF(HPROGRAM,YRECFM,DGR%XPAR_LAI(:,JTIME,1),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!
!-------------------------------------------------------------------------------
!
!*       2.    Definition of ISBA parameters
!              -----------------------------
!
ALLOCATE(DGR%XPAR_VEG        (KDIM,DGR%NTIME,1))
ALLOCATE(DGR%XPAR_RSMIN      (KDIM,1))
ALLOCATE(DGR%XPAR_GAMMA      (KDIM,1))
ALLOCATE(DGR%XPAR_WRMAX_CF   (KDIM,1))
ALLOCATE(DGR%XPAR_RGL        (KDIM,1))
ALLOCATE(DGR%XPAR_CV         (KDIM,1))
ALLOCATE(DGR%XPAR_DG         (KDIM,GRO%NGROUND_LAYER,1))
ALLOCATE(DGR%XPAR_ROOTFRAC   (KDIM,GRO%NGROUND_LAYER,1))
ALLOCATE(DGR%XPAR_DICE       (KDIM,1))
ALLOCATE(DGR%XPAR_Z0         (KDIM,DGR%NTIME,1))
ALLOCATE(DGR%XPAR_Z0_O_Z0H   (KDIM,1))
ALLOCATE(DGR%XPAR_ALBNIR_VEG (KDIM,1))
ALLOCATE(DGR%XPAR_ALBVIS_VEG (KDIM,1))
ALLOCATE(DGR%XPAR_ALBUV_VEG  (KDIM,1))
ALLOCATE(DGR%XPAR_ALBNIR_SOIL(KDIM,1))
ALLOCATE(DGR%XPAR_ALBVIS_SOIL(KDIM,1))
ALLOCATE(DGR%XPAR_ALBUV_SOIL (KDIM,1))
ALLOCATE(DGR%XPAR_EMIS       (KDIM,DGR%NTIME,1))
ALLOCATE(DGR%XPAR_VEGTYPE    (KDIM,NVEGTYPE))
ALLOCATE(DGR%XPAR_GMES       (KDIM,1))
ALLOCATE(DGR%XPAR_RE25       (KDIM,1))
ALLOCATE(DGR%XPAR_BSLAI      (KDIM,1))
ALLOCATE(DGR%XPAR_LAIMIN     (KDIM,1))
ALLOCATE(DGR%XPAR_SEFOLD     (KDIM,1))
ALLOCATE(DGR%XPAR_GC         (KDIM,1))
ALLOCATE(DGR%XPAR_DMAX       (KDIM,1))
ALLOCATE(DGR%XPAR_F2I        (KDIM,1))
ALLOCATE(DGR%LPAR_STRESS    (KDIM,1))
ALLOCATE(DGR%XPAR_H_TREE     (KDIM,1))
ALLOCATE(DGR%XPAR_CE_NITRO   (KDIM,1))
ALLOCATE(DGR%XPAR_CF_NITRO   (KDIM,1))
ALLOCATE(DGR%XPAR_CNA_NITRO  (KDIM,1))
!
DGR%XPAR_VEG          (:,:,:) = XUNDEF
DGR%XPAR_RSMIN          (:,:) = XUNDEF
DGR%XPAR_GAMMA          (:,:) = XUNDEF
DGR%XPAR_WRMAX_CF       (:,:) = XUNDEF
DGR%XPAR_RGL            (:,:) = XUNDEF
DGR%XPAR_CV             (:,:) = XUNDEF
DGR%XPAR_DG           (:,:,:) = XUNDEF
DGR%XPAR_DICE           (:,:) = XUNDEF
DGR%XPAR_ROOTFRAC     (:,:,:) = XUNDEF
DGR%XPAR_Z0           (:,:,:) = XUNDEF
DGR%XPAR_Z0_O_Z0H       (:,:) = XUNDEF
DGR%XPAR_ALBNIR_VEG     (:,:) = XUNDEF
DGR%XPAR_ALBVIS_VEG     (:,:) = XUNDEF
DGR%XPAR_ALBUV_VEG      (:,:) = XUNDEF
DGR%XPAR_ALBNIR_SOIL    (:,:) = XUNDEF
DGR%XPAR_ALBVIS_SOIL    (:,:) = XUNDEF
DGR%XPAR_ALBUV_SOIL     (:,:) = XUNDEF
DGR%XPAR_EMIS         (:,:,:) = XUNDEF
DGR%XPAR_VEGTYPE      (:,:) = XUNDEF
DGR%XPAR_GMES           (:,:) = XUNDEF
DGR%XPAR_RE25           (:,:) = XUNDEF
DGR%XPAR_BSLAI          (:,:) = XUNDEF
DGR%XPAR_LAIMIN         (:,:) = XUNDEF
DGR%XPAR_SEFOLD         (:,:) = XUNDEF
DGR%XPAR_GC             (:,:) = XUNDEF
DGR%XPAR_DMAX           (:,:) = XUNDEF
DGR%XPAR_F2I            (:,:) = XUNDEF
DGR%LPAR_STRESS        (:,:) = .FALSE.
DGR%XPAR_H_TREE         (:,:) = XUNDEF
DGR%XPAR_CE_NITRO       (:,:) = XUNDEF
DGR%XPAR_CF_NITRO       (:,:) = XUNDEF
DGR%XPAR_CNA_NITRO      (:,:) = XUNDEF
!
!---------------------------------------------------------------------------
! Vegtypes adapted to greenroofs:
!--------------------------------
! NPATCH = 1 
! 2D cases : all greenroofs have same vegetation (defined by CTYP_GR)
! (CTYP_GR == 'GRASS') <=> NVT_GRAS (10)
!  ** OR **
! (CTYP_GR == 'SEDUM') <=> NVT_TROG (11)
! NB1: => no aggregation of vegetype parameters needed 
! NB2: Functions existing for gardens are used for initial greenroofs
!      This will need to be refined specifically for greenroofs
!
DGR%XPAR_VEGTYPE(:,:) = 0.
IF (GRO%CTYP_COV == 'GRASS') DGR%XPAR_VEGTYPE(:, NVT_GRAS) = 1.
IF (GRO%CTYP_COV == 'SEDUM') DGR%XPAR_VEGTYPE(:, NVT_TROG) = 1.
!--------------------------------------------------------------------------
!
! Critical normilized soil water content for stress parameterisation
DGR%XPAR_F2I(:,:) = 0.3
!
! Ratio between roughness length for momentum and heat
DGR%XPAR_Z0_O_Z0H(:,:) = 10.
!
! Defensive/offensive strategy (1/0)
DGR%LPAR_STRESS(:,:) = .FALSE. 
!
DO JI=1,KDIM
! 
! Vegetation albedo: near-IR, visible, and UV albedo
! * Will need to be adapted to greenroof GRASS and SEDUM species *
! * vérifier si/où l'abedo ds l'UV est utilisé *
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_ALBNIR_VEG(JI,:)= 0.3
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_ALBNIR_VEG(JI,:)= 0.154 ! mesures ONERA/Doya (2011)

 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_ALBVIS_VEG(JI,:)= 0.10
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_ALBVIS_VEG(JI,:)= 0.154 ! mesures ONERA/Doya (2011)

 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_ALBUV_VEG(JI,:) = 0.0800
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_ALBUV_VEG(JI,:) = 0.1250
!
! Min stomatal resistance  
 !IF(XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  XPAR_RSMIN(JI)= 40 (dans isba & garden)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_RSMIN(JI,:)= 120  ! for GRASS
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_RSMIN(JI,:)= 150. ! for SEDUM
 !IF(XPAR_VEGTYPE(JI,NVT_TROG)>0. )  XPAR_RSMIN(JI)= 120.
! 
! Gamma parameter 
! (* Check if values needs to be refined for GRASS and SEDUM *)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_GAMMA(JI,:)= 0.
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_GAMMA(JI,:)= 0.
!
! Wrmax_cf 
! (* Check if needs to be refined for GRASS and SEDUM greenroofs *)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_WRMAX_CF(JI,:)= 0.2
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_WRMAX_CF(JI,:)= 0.2
!
! Rgl 
! (* Check if needs to be refined for GRASS and SEDUM greenroofs *)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_RGL(JI,:)= 100.
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_RGL(JI,:)= 100.
!
! Cv 
! (* Check if needs to be refined for GRASS and SEDUM greenroofs *)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_CV(JI,:)= 2.E-5
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_CV(JI,:)= 2.E-5
!
!! Mesophyll conductance (m s-1) 
! (* Check if needs to be refined for GRASS and SEDUM greenroofs *)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_GMES(JI,:)= 0.020
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_GMES(JI,:)= 0.020
 !IF(XPAR_VEGTYPE(JI,NVT_TROG)>0. )  XPAR_GMES(JI)= 0.003
!
! Ecosystem Respiration (kg/kg.m.s-1)
! (* Check if needs to be refined for GRASS and SEDUM greenroofs *)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0.  )  DGR%XPAR_RE25(JI,:)= 3.0E-7
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG  )>0.)  DGR%XPAR_RE25(JI,:)= 3.0E-7
!
! Cuticular conductance (m s-1)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_GC(JI,:)= 0.00025
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_GC(JI,:)= 0.00025        
!
! Ratio d(biomass)/d(lai) (kg/m2)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_BSLAI(JI,:)= 0.36
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_BSLAI(JI,:)= 0.06
!
! Maximum air saturation deficit tolerate by vegetation (kg/kg)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_DMAX(JI,:)= 0.1
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_DMAX(JI,:)= 0.1
!
! e-folding time for senescence (days)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_SEFOLD(JI,:)=  90.* XDAY
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_SEFOLD(JI,:)=  60.* XDAY
!
! Minimum LAI (m2/m2)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_LAIMIN (JI,:) = 0.3
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_LAIMIN (JI,:) = 0.3
!
! Leaf aera ratio sensitivity to nitrogen concentration
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_CE_NITRO(JI,:)= 5.56
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_CE_NITRO(JI,:)= 3.79
!
! Lethal minimum value of leaf area ratio
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_CF_NITRO(JI,:)=  6.73
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DGR%XPAR_CF_NITRO(JI,:)=  9.84
!
! Nitrogen concentration of active biomass
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DGR%XPAR_CNA_NITRO(JI,:)= 1.9
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG )>0.)  DGR%XPAR_CNA_NITRO(JI,:)= 1.3
!
! Depth of greenroof ground layers
 DGR%XPAR_DG(JI, 1,:) = XGRID_SOIL(NGRID_LEVEL - 5)
 DGR%XPAR_DG(JI, 2,:) = XGRID_SOIL(NGRID_LEVEL - 4)
 DGR%XPAR_DG(JI, 3,:) = XGRID_SOIL(NGRID_LEVEL - 3)
 DGR%XPAR_DG(JI, 4,:) = XGRID_SOIL(NGRID_LEVEL - 2)
 DGR%XPAR_DG(JI, 5,:) = XGRID_SOIL(NGRID_LEVEL - 1)
 DGR%XPAR_DG(JI, 6,:) = XGRID_SOIL(NGRID_LEVEL - 0)
!
! Root fractions
 DGR%XPAR_ROOTFRAC(JI, 1,:)  = 0.04
 DGR%XPAR_ROOTFRAC(JI, 2,:)  = 0.36
 DGR%XPAR_ROOTFRAC(JI, 3,:)  = 0.68
 DGR%XPAR_ROOTFRAC(JI, 4,:)  = 1.
 DGR%XPAR_ROOTFRAC(JI, 5,:)  = 1.
 DGR%XPAR_ROOTFRAC(JI, 6,:)  = 1.
!
! Depth of the soil column for the calculation of the frozen soil fraction (m)
 DGR%XPAR_DICE(JI,:) = DGR%XPAR_DG(JI,1,:) 
!
DO JTIME=1,DGR%NTIME
! Leaf Area Index

! Fraction of vegetation on greenroof
!* Will need to be refined for greenroofs *)
  !XPAR_VEG (JI,1,JTIME) = VEG_FROM_LAI (XPAR_LAI_GR(JI,JTIME),   &
  !                                       XPAR_VEGTYPE(JI,:),GAGRI_TO_GRASS)  
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )   DGR%XPAR_VEG (JI,JTIME,:) = 0.9
 !IF(XPAR_VEGTYPE(JI,NVT_TROG)>0. )   XPAR_VEG (JI,JTIME) = 1.0
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )   DGR%XPAR_VEG (JI,JTIME,:) = 0.95

! Roughness length for momentum
!* Will need to be refined for greenroofs *)
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )   DGR%XPAR_Z0 (JI,JTIME,:) = 0.01
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )   DGR%XPAR_Z0 (JI,JTIME,:) = 0.01
 !                                        
! Emissivity
!* Will need to be refined for greenroofs *)
  !XPAR_EMIS (JI,1,JTIME) = EMIS_FROM_VEG (XPAR_VEG    (JI,1,JTIME),&
  !                                         XPAR_VEGTYPE(JI,:))  
 IF(DGR%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )   DGR%XPAR_EMIS (JI,JTIME,:) = 0.95 
 IF(DGR%XPAR_VEGTYPE(JI,NVT_TROG)>0. )   DGR%XPAR_EMIS (JI,JTIME,:) = 0.83 ! Feng. et al. (2010)

END DO
!
ENDDO
!
IF (LHOOK) CALL DR_HOOK('READ_PGD_TEB_GREENROOF_PAR_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_PGD_TEB_GREENROOF_PAR_n
