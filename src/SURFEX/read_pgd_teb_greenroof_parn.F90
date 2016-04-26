!     #########
      SUBROUTINE READ_PGD_TEB_GREENROOF_PAR_n (DTI, IO, S, K, KDIM, HPROGRAM)
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
USE MODD_ISBA_n, ONLY : ISBA_S_t, ISBA_K_t
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
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_K_t), INTENT(INOUT) :: K
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
 CALL READ_SURF(HPROGRAM,YRECFM,DTI%NTIME,IRESP)
!
! Read type of green roof
YRECFM='D_TYPE_GR'
 CALL READ_SURF(HPROGRAM,YRECFM,IO%CTYP_COV,IRESP)
!
! Read green roof OM fraction
DO JLAYER=1,IO%NGROUND_LAYER
  !WRITE(YRECFM,FMT='(A8,I1.1)') 'D_OM_GR0',JLAYER
  WRITE(YRECFM,FMT='(A7,I2.2)') 'D_OM_GR',JLAYER
  CALL READ_SURF(HPROGRAM,YRECFM,S%XSOC(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
! Read green roof SAND fraction
DO JLAYER=1,IO%NGROUND_LAYER
  !WRITE(YRECFM,FMT='(A10,I1.1)') 'D_SAND_GR0',JLAYER
  WRITE(YRECFM,FMT='(A9,I2.2)') 'D_SAND_GR',JLAYER
  CALL READ_SURF(HPROGRAM,YRECFM,K%XSAND(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
! Read green roof CLAY fraction
DO JLAYER=1,IO%NGROUND_LAYER
  !WRITE(YRECFM,FMT='(A10,I1.1)') 'D_CLAY_GR0',JLAYER
  WRITE(YRECFM,FMT='(A9,I2.2)') 'D_CLAY_GR',JLAYER
  CALL READ_SURF(HPROGRAM,YRECFM,K%XCLAY(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
! Read green roof LAI
ALLOCATE(DTI%XPAR_LAI    (KDIM,DTI%NTIME,1))
DO JTIME=1,DTI%NTIME
  WRITE(YRECFM,FMT='(A8,I2.2)') 'D_LAI_GR',JTIME
  CALL READ_SURF(HPROGRAM,YRECFM,DTI%XPAR_LAI(:,JTIME,1),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!
!-------------------------------------------------------------------------------
!
!*       2.    Definition of ISBA parameters
!              -----------------------------
!
ALLOCATE(DTI%XPAR_VEG        (KDIM,DTI%NTIME,1))
ALLOCATE(DTI%XPAR_RSMIN      (KDIM,1))
ALLOCATE(DTI%XPAR_GAMMA      (KDIM,1))
ALLOCATE(DTI%XPAR_WRMAX_CF   (KDIM,1))
ALLOCATE(DTI%XPAR_RGL        (KDIM,1))
ALLOCATE(DTI%XPAR_CV         (KDIM,1))
ALLOCATE(DTI%XPAR_DG         (KDIM,IO%NGROUND_LAYER,1))
ALLOCATE(DTI%XPAR_ROOTFRAC   (KDIM,IO%NGROUND_LAYER,1))
ALLOCATE(DTI%XPAR_DICE       (KDIM,1))
ALLOCATE(DTI%XPAR_Z0         (KDIM,DTI%NTIME,1))
ALLOCATE(DTI%XPAR_Z0_O_Z0H   (KDIM,1))
ALLOCATE(DTI%XPAR_ALBNIR_VEG (KDIM,1))
ALLOCATE(DTI%XPAR_ALBVIS_VEG (KDIM,1))
ALLOCATE(DTI%XPAR_ALBUV_VEG  (KDIM,1))
ALLOCATE(DTI%XPAR_ALBNIR_SOIL(KDIM,1))
ALLOCATE(DTI%XPAR_ALBVIS_SOIL(KDIM,1))
ALLOCATE(DTI%XPAR_ALBUV_SOIL (KDIM,1))
ALLOCATE(DTI%XPAR_EMIS       (KDIM,DTI%NTIME,1))
ALLOCATE(DTI%XPAR_VEGTYPE    (KDIM,NVEGTYPE))
ALLOCATE(DTI%XPAR_GMES       (KDIM,1))
ALLOCATE(DTI%XPAR_RE25       (KDIM,1))
ALLOCATE(DTI%XPAR_BSLAI      (KDIM,1))
ALLOCATE(DTI%XPAR_LAIMIN     (KDIM,1))
ALLOCATE(DTI%XPAR_SEFOLD     (KDIM,1))
ALLOCATE(DTI%XPAR_GC         (KDIM,1))
ALLOCATE(DTI%XPAR_DMAX       (KDIM,1))
ALLOCATE(DTI%XPAR_F2I        (KDIM,1))
ALLOCATE(DTI%LPAR_STRESS    (KDIM,1))
ALLOCATE(DTI%XPAR_H_TREE     (KDIM,1))
ALLOCATE(DTI%XPAR_CE_NITRO   (KDIM,1))
ALLOCATE(DTI%XPAR_CF_NITRO   (KDIM,1))
ALLOCATE(DTI%XPAR_CNA_NITRO  (KDIM,1))
!
DTI%XPAR_VEG          (:,:,:) = XUNDEF
DTI%XPAR_RSMIN          (:,:) = XUNDEF
DTI%XPAR_GAMMA          (:,:) = XUNDEF
DTI%XPAR_WRMAX_CF       (:,:) = XUNDEF
DTI%XPAR_RGL            (:,:) = XUNDEF
DTI%XPAR_CV             (:,:) = XUNDEF
DTI%XPAR_DG           (:,:,:) = XUNDEF
DTI%XPAR_DICE           (:,:) = XUNDEF
DTI%XPAR_ROOTFRAC     (:,:,:) = XUNDEF
DTI%XPAR_Z0           (:,:,:) = XUNDEF
DTI%XPAR_Z0_O_Z0H       (:,:) = XUNDEF
DTI%XPAR_ALBNIR_VEG     (:,:) = XUNDEF
DTI%XPAR_ALBVIS_VEG     (:,:) = XUNDEF
DTI%XPAR_ALBUV_VEG      (:,:) = XUNDEF
DTI%XPAR_ALBNIR_SOIL    (:,:) = XUNDEF
DTI%XPAR_ALBVIS_SOIL    (:,:) = XUNDEF
DTI%XPAR_ALBUV_SOIL     (:,:) = XUNDEF
DTI%XPAR_EMIS         (:,:,:) = XUNDEF
DTI%XPAR_VEGTYPE      (:,:) = XUNDEF
DTI%XPAR_GMES           (:,:) = XUNDEF
DTI%XPAR_RE25           (:,:) = XUNDEF
DTI%XPAR_BSLAI          (:,:) = XUNDEF
DTI%XPAR_LAIMIN         (:,:) = XUNDEF
DTI%XPAR_SEFOLD         (:,:) = XUNDEF
DTI%XPAR_GC             (:,:) = XUNDEF
DTI%XPAR_DMAX           (:,:) = XUNDEF
DTI%XPAR_F2I            (:,:) = XUNDEF
DTI%LPAR_STRESS        (:,:) = .FALSE.
DTI%XPAR_H_TREE         (:,:) = XUNDEF
DTI%XPAR_CE_NITRO       (:,:) = XUNDEF
DTI%XPAR_CF_NITRO       (:,:) = XUNDEF
DTI%XPAR_CNA_NITRO      (:,:) = XUNDEF
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
DTI%XPAR_VEGTYPE(:,:) = 0.
IF (IO%CTYP_COV == 'GRASS') DTI%XPAR_VEGTYPE(:, NVT_GRAS) = 1.
IF (IO%CTYP_COV == 'SEDUM') DTI%XPAR_VEGTYPE(:, NVT_TROG) = 1.
!--------------------------------------------------------------------------
!
! Critical normilized soil water content for stress parameterisation
DTI%XPAR_F2I(:,:) = 0.3
!
! Ratio between roughness length for momentum and heat
DTI%XPAR_Z0_O_Z0H(:,:) = 10.
!
! Defensive/offensive strategy (1/0)
DTI%LPAR_STRESS(:,:) = .FALSE. 
!
DO JI=1,KDIM
! 
! Vegetation albedo: near-IR, visible, and UV albedo
! * Will need to be adapted to greenroof GRASS and SEDUM species *
! * vérifier si/où l'abedo ds l'UV est utilisé *
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_ALBNIR_VEG(JI,:)= 0.3
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_ALBNIR_VEG(JI,:)= 0.154 ! mesures ONERA/Doya (2011)

 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_ALBVIS_VEG(JI,:)= 0.10
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_ALBVIS_VEG(JI,:)= 0.154 ! mesures ONERA/Doya (2011)

 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_ALBUV_VEG(JI,:) = 0.0800
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_ALBUV_VEG(JI,:) = 0.1250
!
! Min stomatal resistance  
 !IF(XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  XPAR_RSMIN(JI)= 40 (dans isba & garden)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_RSMIN(JI,:)= 120  ! for GRASS
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_RSMIN(JI,:)= 150. ! for SEDUM
 !IF(XPAR_VEGTYPE(JI,NVT_TROG)>0. )  XPAR_RSMIN(JI)= 120.
! 
! Gamma parameter 
! (* Check if values needs to be refined for GRASS and SEDUM *)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_GAMMA(JI,:)= 0.
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_GAMMA(JI,:)= 0.
!
! Wrmax_cf 
! (* Check if needs to be refined for GRASS and SEDUM greenroofs *)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_WRMAX_CF(JI,:)= 0.2
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_WRMAX_CF(JI,:)= 0.2
!
! Rgl 
! (* Check if needs to be refined for GRASS and SEDUM greenroofs *)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_RGL(JI,:)= 100.
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_RGL(JI,:)= 100.
!
! Cv 
! (* Check if needs to be refined for GRASS and SEDUM greenroofs *)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_CV(JI,:)= 2.E-5
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_CV(JI,:)= 2.E-5
!
!! Mesophyll conductance (m s-1) 
! (* Check if needs to be refined for GRASS and SEDUM greenroofs *)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_GMES(JI,:)= 0.020
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_GMES(JI,:)= 0.020
 !IF(XPAR_VEGTYPE(JI,NVT_TROG)>0. )  XPAR_GMES(JI)= 0.003
!
! Ecosystem Respiration (kg/kg.m.s-1)
! (* Check if needs to be refined for GRASS and SEDUM greenroofs *)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0.  )  DTI%XPAR_RE25(JI,:)= 3.0E-7
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG  )>0.)  DTI%XPAR_RE25(JI,:)= 3.0E-7
!
! Cuticular conductance (m s-1)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_GC(JI,:)= 0.00025
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_GC(JI,:)= 0.00025        
!
! Ratio d(biomass)/d(lai) (kg/m2)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_BSLAI(JI,:)= 0.36
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_BSLAI(JI,:)= 0.06
!
! Maximum air saturation deficit tolerate by vegetation (kg/kg)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_DMAX(JI,:)= 0.1
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_DMAX(JI,:)= 0.1
!
! e-folding time for senescence (days)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_SEFOLD(JI,:)=  90.* XDAY
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_SEFOLD(JI,:)=  60.* XDAY
!
! Minimum LAI (m2/m2)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_LAIMIN (JI,:) = 0.3
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_LAIMIN (JI,:) = 0.3
!
! Leaf aera ratio sensitivity to nitrogen concentration
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_CE_NITRO(JI,:)= 5.56
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_CE_NITRO(JI,:)= 3.79
!
! Lethal minimum value of leaf area ratio
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_CF_NITRO(JI,:)=  6.73
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )  DTI%XPAR_CF_NITRO(JI,:)=  9.84
!
! Nitrogen concentration of active biomass
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )  DTI%XPAR_CNA_NITRO(JI,:)= 1.9
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG )>0.)  DTI%XPAR_CNA_NITRO(JI,:)= 1.3
!
! Depth of greenroof ground layers
 DTI%XPAR_DG(JI, 1,:) = XGRID_SOIL(NGRID_LEVEL - 5)
 DTI%XPAR_DG(JI, 2,:) = XGRID_SOIL(NGRID_LEVEL - 4)
 DTI%XPAR_DG(JI, 3,:) = XGRID_SOIL(NGRID_LEVEL - 3)
 DTI%XPAR_DG(JI, 4,:) = XGRID_SOIL(NGRID_LEVEL - 2)
 DTI%XPAR_DG(JI, 5,:) = XGRID_SOIL(NGRID_LEVEL - 1)
 DTI%XPAR_DG(JI, 6,:) = XGRID_SOIL(NGRID_LEVEL - 0)
!
! Root fractions
 DTI%XPAR_ROOTFRAC(JI, 1,:)  = 0.04
 DTI%XPAR_ROOTFRAC(JI, 2,:)  = 0.36
 DTI%XPAR_ROOTFRAC(JI, 3,:)  = 0.68
 DTI%XPAR_ROOTFRAC(JI, 4,:)  = 1.
 DTI%XPAR_ROOTFRAC(JI, 5,:)  = 1.
 DTI%XPAR_ROOTFRAC(JI, 6,:)  = 1.
!
! Depth of the soil column for the calculation of the frozen soil fraction (m)
 DTI%XPAR_DICE(JI,:) = DTI%XPAR_DG(JI,1,:) 
!
DO JTIME=1,DTI%NTIME
! Leaf Area Index

! Fraction of vegetation on greenroof
!* Will need to be refined for greenroofs *)
  !XPAR_VEG (JI,1,JTIME) = VEG_FROM_LAI (XPAR_LAI_GR(JI,JTIME),   &
  !                                       XPAR_VEGTYPE(JI,:),GAGRI_TO_GRASS)  
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )   DTI%XPAR_VEG (JI,JTIME,:) = 0.9
 !IF(XPAR_VEGTYPE(JI,NVT_TROG)>0. )   XPAR_VEG (JI,JTIME) = 1.0
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )   DTI%XPAR_VEG (JI,JTIME,:) = 0.95

! Roughness length for momentum
!* Will need to be refined for greenroofs *)
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )   DTI%XPAR_Z0 (JI,JTIME,:) = 0.01
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )   DTI%XPAR_Z0 (JI,JTIME,:) = 0.01
 !                                        
! Emissivity
!* Will need to be refined for greenroofs *)
  !XPAR_EMIS (JI,1,JTIME) = EMIS_FROM_VEG (XPAR_VEG    (JI,1,JTIME),&
  !                                         XPAR_VEGTYPE(JI,:))  
 IF(DTI%XPAR_VEGTYPE(JI,NVT_GRAS)>0. )   DTI%XPAR_EMIS (JI,JTIME,:) = 0.95 
 IF(DTI%XPAR_VEGTYPE(JI,NVT_TROG)>0. )   DTI%XPAR_EMIS (JI,JTIME,:) = 0.83 ! Feng. et al. (2010)

END DO
!
ENDDO
!
IF (LHOOK) CALL DR_HOOK('READ_PGD_TEB_GREENROOF_PAR_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_PGD_TEB_GREENROOF_PAR_n
