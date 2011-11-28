!     #########
      SUBROUTINE CONVERT_COVER_ISBA   (HISBA,KDECADE,PCOVER,HPHOTO,          &
                                         HSFTYPE,PVEG,                       &
                                         PLAI,PRSMIN,PGAMMA,PWRMAX_CF,       &
                                         PRGL,PCV,PDG,PD_ICE,PZ0,PZ0_O_Z0H,  &
                                         PALBNIR_VEG,PALBVIS_VEG,PALBUV_VEG, &
                                         PEMIS_ECO,                          &
                                         PVEGTYPE,PROOTFRAC,                 &
                                         PGMES,PBSLAI,PLAIMIN,PSEFOLD,PGC,   &
                                         PDMAX, PF2I, OSTRESS, PH_TREE,PRE25,&
                                         PCE_NITRO, PCF_NITRO, PCNA_NITRO,   &
                                         TPSEED, TPREAP, PWATSUP, PIRRIG     )  
!     ##############################################################
!
!!**** *CONVERT_COVER* convert surface cover classes into secondary 
!!                     physiographic variables for ISBA
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
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
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original   01/2004
!!    
!!      P Le Moigne 09/2005 AGS modifs of L. Jarlan
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_COVER,     ONLY : XDATA_LAI, XDATA_H_TREE, XDATA_VEGTYPE,   &
                                  XDATA_VEG, XDATA_Z0, XDATA_Z0_O_Z0H,    &
                                  XDATA_EMIS_ECO, XDATA_GAMMA, XDATA_CV,  &
                                  XDATA_RGL, XDATA_RSMIN, XDATA_DG,       &
                                  XDATA_ALBNIR_VEG, XDATA_ALBVIS_VEG,     &
                                  XDATA_ALBUV_VEG, XDATA_DICE,            &
                                  XDATA_GMES, XDATA_BSLAI, XDATA_LAIMIN,  &
                                  XDATA_SEFOLD, XDATA_GC, XDATA_WRMAX_CF, &
                                  XDATA_CUM_ROOT_FRAC, XDATA_STRESS,      &
                                  XDATA_DMAX, XDATA_F2I, XDATA_RE25,      &
                                  XDATA_CE_NITRO, XDATA_CF_NITRO,         &
                                  XDATA_CNA_NITRO,                        &
                                  XDATA_GMES_ST, XDATA_BSLAI_ST,          &
                                  XDATA_SEFOLD_ST, XDATA_GC_ST,           &
                                  XDATA_DMAX_ST, XDATA_WATSUP,            &
                                  TDATA_SEED, TDATA_REAP,XDATA_IRRIG,     &
                                  XDATA_ROOT_DEPTH, XDATA_GROUND_DEPTH,   &
                                  XDATA_ROOT_EXTINCTION, XDATA_ROOT_LIN  
!
USE MODD_DATA_COVER_n,   ONLY : XDATA_NATURE, XDATA_GARDEN
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE, JPCOVER
USE MODD_TYPE_DATE_SURF
!
USE MODI_AV_PGD
USE MODI_INI_DATA_SOIL
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=*),       INTENT(IN)    :: HISBA   ! type of soil (Force-Restore OR Diffusion)
INTEGER,                INTENT(IN)    :: KDECADE
REAL, DIMENSION(:,:),   INTENT(IN)    :: PCOVER
CHARACTER(LEN=*),       INTENT(IN)    :: HPHOTO  ! type of photosynthesis
CHARACTER(LEN=*),       INTENT(IN)    :: HSFTYPE ! nature / garden
!
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PVEG
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PLAI
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PRSMIN
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PGAMMA
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PWRMAX_CF
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PRGL
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PCV
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)   :: PDG
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PD_ICE
REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(OUT)   :: PROOTFRAC
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PZ0
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PZ0_O_Z0H
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PALBNIR_VEG
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PALBVIS_VEG
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PALBUV_VEG
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PEMIS_ECO
!
REAL, DIMENSION(:,:), OPTIONAL, INTENT(OUT)   :: PVEGTYPE
!
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PGMES
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PRE25
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PBSLAI
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PLAIMIN
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PSEFOLD
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PGC
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PDMAX
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PF2I
LOGICAL, DIMENSION(:,:),OPTIONAL, INTENT(OUT)   :: OSTRESS
!
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PH_TREE
!
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PCE_NITRO
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PCF_NITRO
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PCNA_NITRO
!
TYPE(DATE_TIME), DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: TPSEED
TYPE(DATE_TIME), DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: TPREAP
!
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PWATSUP
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PIRRIG
!
!*    0.2    Declaration of local variables
!            ------------------------------
! calculation of veg from lai in the pixel
REAL, DIMENSION (:,:), ALLOCATABLE           :: ZWORK      ! work array
REAL, DIMENSION (JPCOVER, NVEGTYPE)          :: ZDATA_DG   !        "
REAL, DIMENSION (JPCOVER, NVEGTYPE)          :: ZDATA_ROOT !        "
!
CHARACTER(LEN=3)  :: YTREE, YNAT, YLAI, YVEG
!
INTEGER :: JLAYER ! loop counter on surface layers
INTEGER :: JVEG   ! loop counter on vegetation types
INTEGER :: JPATCH ! loop counter on patches
!
INTEGER :: IGROUND! number of ground layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*    2.      SECONDARY VARIABLES
!             -------------------
!
IF (LHOOK) CALL DR_HOOK('CONVERT_COVER_ISBA',0,ZHOOK_HANDLE)
IF (HSFTYPE=='NAT') THEN
  YNAT='NAT'
  YTREE='TRE'
  YLAI='LAI'
  YVEG='VEG'
ELSEIF (HSFTYPE=='GRD') THEN
  YNAT='GRD'
  YTREE='GRT'
  YLAI='GRL'
  YVEG='GRV'
ENDIF
!
!*    2.1     fields on natural surfaces only, taking into account patches/ 
!             -------------------------------
!
!
IF (PRESENT(PH_TREE)) &
  CALL AV_PGD (PH_TREE ,PCOVER ,XDATA_H_TREE (:,:) ,YTREE,'ARI')  
!
DO JVEG=1,NVEGTYPE
  IF (PRESENT(PVEGTYPE)) &
    CALL AV_PGD (PVEGTYPE(:,JVEG),PCOVER ,XDATA_VEGTYPE(:,JVEG),YNAT,'ARI')  
END DO
!
! vegetation fraction
! -------------------
!
IF (PRESENT(PVEG)) &
  CALL AV_PGD (PVEG ,PCOVER ,XDATA_VEG (:,KDECADE,:),YNAT,'ARI')  
!
! Leaf Aera Index
! ---------------
!
IF (PRESENT(PLAI)) &
  CALL AV_PGD (PLAI ,PCOVER ,XDATA_LAI (:,KDECADE,:),YVEG,'ARI',KDECADE=KDECADE)  
!
! roughness length
! ----------------
!
IF (PRESENT(PZ0)) &
  CALL AV_PGD (PZ0 ,PCOVER ,XDATA_Z0 (:,KDECADE,:),YNAT,'CDN')  
!
IF (PRESENT(PZ0_O_Z0H)) &
  CALL AV_PGD (PZ0_O_Z0H ,PCOVER ,XDATA_Z0_O_Z0H (:,:),YNAT,'ARI')  
!
!
!emis-eco
!--------
!
IF (PRESENT(PEMIS_ECO)) &
  CALL AV_PGD (PEMIS_ECO ,PCOVER ,XDATA_EMIS_ECO (:,KDECADE,:),YNAT,'ARI')  
! 
!---------------------------------------------------------------------------------
! 
!* 1/Rsmin
!
IF (PRESENT(PRSMIN)) THEN
  IF (SIZE(PRSMIN)>0) &
    CALL AV_PGD (PRSMIN,PCOVER ,XDATA_RSMIN,YLAI,'INV',KDECADE=KDECADE)  
END IF
!
!* other vegetation parameters
!
IF (PRESENT(PGAMMA)) &
  CALL AV_PGD (PGAMMA     ,PCOVER ,XDATA_GAMMA   (:,:),YVEG,'ARI',KDECADE=KDECADE)  
IF (PRESENT(PWRMAX_CF)) &
  CALL AV_PGD (PWRMAX_CF  ,PCOVER ,XDATA_WRMAX_CF(:,:),YVEG,'ARI',KDECADE=KDECADE)  
!
!
IF (PRESENT(PRGL)) &
  CALL AV_PGD (PRGL       ,PCOVER ,XDATA_RGL   (:,:),YVEG,'ARI',KDECADE=KDECADE)  
IF (PRESENT(PCV)) &
  CALL AV_PGD (PCV        ,PCOVER ,XDATA_CV    (:,:),YVEG,'INV',KDECADE=KDECADE)  
!
!---------------------------------------------------------------------------------
!
!* soil layers
!  -----------
!
IF (PRESENT(PDG)) THEN
  !
  !* definition of soil layers and cumulative root fraction
  !
  IGROUND = SIZE(PDG,2)
  !
  IF (.NOT.ALLOCATED(XDATA_DG)) &
    ALLOCATE(XDATA_DG(SIZE(XDATA_ROOT_DEPTH,1),IGROUND,SIZE(XDATA_ROOT_DEPTH,2)))
  IF (.NOT.ALLOCATED(XDATA_CUM_ROOT_FRAC)) &
    ALLOCATE(XDATA_CUM_ROOT_FRAC(SIZE(XDATA_ROOT_DEPTH,1),IGROUND,SIZE(XDATA_ROOT_DEPTH,2)))

  CALL INI_DATA_SOIL(HISBA,PSURF=XDATA_NATURE,PSURF2=XDATA_GARDEN, &
    PROOT_DEPTH=XDATA_ROOT_DEPTH,PGROUND_DEPTH=XDATA_GROUND_DEPTH,&
    PROOT_EXTINCTION=XDATA_ROOT_EXTINCTION,PROOT_LIN=XDATA_ROOT_LIN,&
    PDG_OUT=XDATA_DG,PCUM_ROOT_FRAC=XDATA_CUM_ROOT_FRAC)
  !
  !* soil layers
  !
  DO JLAYER=1,IGROUND
    CALL AV_PGD (PDG(:,JLAYER,:),PCOVER,XDATA_DG(:,JLAYER,:),YNAT,'ARI')
  ENDDO

END IF

!
!* cumulative root fraction
!
IF (PRESENT(PROOTFRAC)) THEN
  ALLOCATE(ZWORK(SIZE(PROOTFRAC,1),SIZE(PROOTFRAC,3)))
  DO JLAYER=1,IGROUND
    DO JVEG=1,NVEGTYPE
      ZDATA_ROOT(:,JVEG) = XDATA_CUM_ROOT_FRAC(:,JLAYER,JVEG)
    END DO
    CALL AV_PGD (ZWORK(:,:),PCOVER ,ZDATA_ROOT(:,:),YNAT,'ARI')
    DO JPATCH=1,SIZE(PDG,3)
      PROOTFRAC(:,JLAYER,JPATCH) = ZWORK(:,JPATCH)
    END DO
  END DO
  DEALLOCATE(ZWORK)
END IF
!
!* soil ice for runoff
!  -------------------
!
IF (PRESENT(PD_ICE)) &
 CALL AV_PGD (PD_ICE,PCOVER ,XDATA_DICE(:,:),YNAT,'ARI')
!
!---------------------------------------------------------------------------------
!
IF (PRESENT(PALBNIR_VEG)) &
  CALL AV_PGD (PALBNIR_VEG,PCOVER ,XDATA_ALBNIR_VEG(:,:),YVEG,'ARI',KDECADE=KDECADE)  
IF (PRESENT(PALBVIS_VEG)) &
  CALL AV_PGD (PALBVIS_VEG,PCOVER ,XDATA_ALBVIS_VEG(:,:),YVEG,'ARI',KDECADE=KDECADE)  
IF (PRESENT(PALBUV_VEG)) &
  CALL AV_PGD (PALBUV_VEG, PCOVER ,XDATA_ALBUV_VEG (:,:),YVEG,'ARI',KDECADE=KDECADE)  
! parameters for "stress option"
IF (HPHOTO == 'AST' .OR. HPHOTO == 'LST' .OR. HPHOTO == 'NIT' .OR.  HPHOTO == 'NCB')  THEN

   IF (PRESENT(PGMES)) THEN
     IF (SIZE(PGMES)>0) &
       CALL AV_PGD (PGMES  ,PCOVER ,XDATA_GMES_ST  (:,:),YVEG,'ARI',KDECADE=KDECADE)  
   END IF

   IF (PRESENT(PBSLAI)) THEN
     IF (SIZE(PBSLAI)>0) &
       CALL AV_PGD (PBSLAI ,PCOVER ,XDATA_BSLAI_ST (:,:),YVEG,'ARI',KDECADE=KDECADE)  
   END IF

   IF (PRESENT(PSEFOLD)) THEN
     IF (SIZE(PSEFOLD)>0) &
       CALL AV_PGD (PSEFOLD,PCOVER ,XDATA_SEFOLD_ST(:,:),YVEG,'ARI',KDECADE=KDECADE)  
   END IF

   IF (PRESENT(PGC)) THEN
     IF (SIZE(PGC)>0) &
       CALL AV_PGD (PGC    ,PCOVER ,XDATA_GC_ST    (:,:),YVEG,'ARI',KDECADE=KDECADE)  
   END IF

   IF (PRESENT(PDMAX)) THEN
     IF (SIZE(PDMAX)>0) &
       CALL AV_PGD (PDMAX  ,PCOVER ,XDATA_DMAX_ST  (:,:),YTREE,'ARI',KDECADE=KDECADE)  
   END IF

ELSE

   IF (PRESENT(PGMES)) THEN
     IF (SIZE(PGMES)>0) &
       CALL AV_PGD (PGMES  ,PCOVER ,XDATA_GMES  (:,:),YVEG,'ARI',KDECADE=KDECADE)  
   END IF

   IF (PRESENT(PBSLAI)) THEN
     IF (SIZE(PBSLAI)>0) &
       CALL AV_PGD (PBSLAI ,PCOVER ,XDATA_BSLAI (:,:),YVEG,'ARI',KDECADE=KDECADE)  
   END IF

   IF (PRESENT(PSEFOLD)) THEN
     IF (SIZE(PSEFOLD)>0) &
       CALL AV_PGD (PSEFOLD,PCOVER ,XDATA_SEFOLD(:,:),YVEG,'ARI',KDECADE=KDECADE)  
   END IF

   IF (PRESENT(PGC)) THEN
     IF (SIZE(PGC)>0) &
       CALL AV_PGD (PGC    ,PCOVER ,XDATA_GC    (:,:),YVEG,'ARI',KDECADE=KDECADE)  
   END IF

   IF (PRESENT(PDMAX)) THEN
     IF (SIZE(PDMAX)>0) &
       CALL AV_PGD (PDMAX  ,PCOVER ,XDATA_DMAX  (:,:),YTREE,'ARI',KDECADE=KDECADE)  
   END IF

ENDIF

IF (PRESENT(PRE25)) THEN
  IF (SIZE(PRE25)>0) &
    CALL AV_PGD (PRE25  ,PCOVER ,XDATA_RE25  (:,:),YNAT,'ARI')  
END IF

IF (PRESENT(PLAIMIN)) THEN
  IF (SIZE(PLAIMIN)>0) &
    CALL AV_PGD (PLAIMIN,PCOVER ,XDATA_LAIMIN(:,:),YVEG,'ARI',KDECADE=KDECADE)  
END IF

IF (PRESENT(PCE_NITRO)) THEN
  IF (SIZE(PCE_NITRO)>0) &
    CALL AV_PGD (PCE_NITRO  ,PCOVER ,XDATA_CE_NITRO  (:,:),YVEG,'ARI',KDECADE=KDECADE)  
END IF

IF (PRESENT(PCF_NITRO)) THEN
  IF (SIZE(PCF_NITRO)>0) &
    CALL AV_PGD (PCF_NITRO  ,PCOVER ,XDATA_CF_NITRO  (:,:),YVEG,'ARI',KDECADE=KDECADE)  
END IF

IF (PRESENT(PCNA_NITRO)) THEN
  IF (SIZE(PCNA_NITRO)>0) &
    CALL AV_PGD (PCNA_NITRO  ,PCOVER ,XDATA_CNA_NITRO(:,:),YVEG,'ARI',KDECADE=KDECADE)  
END IF

IF (PRESENT(PF2I)) THEN
  IF (SIZE(PF2I)>0) &
    CALL AV_PGD (PF2I   ,PCOVER ,XDATA_F2I   (:,:),YVEG,'ARI',KDECADE=KDECADE)  
END IF
!
IF (PRESENT(OSTRESS)) THEN
  IF (SIZE(OSTRESS)>0) THEN
    ALLOCATE(ZWORK(SIZE(OSTRESS,1),SIZE(OSTRESS,2)))
    CALL AV_PGD (ZWORK,PCOVER ,XDATA_STRESS(:,:),YVEG,'ARI',KDECADE=KDECADE)
    WHERE (ZWORK<0.5) 
      OSTRESS = .FALSE.
    ELSEWHERE
      OSTRESS = .TRUE.
    END WHERE
    DEALLOCATE(ZWORK)
  END IF
END IF
!
IF (HPHOTO == 'LAI' .OR. HPHOTO == 'LST' .OR. HPHOTO == 'NIT')  THEN
   !
   ! date of seeding
   ! ---------------
   !
   IF (PRESENT(TPSEED)) THEN
     IF (SIZE(TPSEED)>0) &
       CALL AV_PGD (TPSEED ,PCOVER ,TDATA_SEED(:,:),YVEG,'MAJ',KDECADE=KDECADE)  
   END IF
   !
   ! date of reaping
   ! ---------------
   !
   IF (PRESENT(TPREAP)) THEN
     IF (SIZE(TPREAP)>0) &
       CALL AV_PGD (TPREAP ,PCOVER ,TDATA_REAP(:,:),YVEG,'MAJ',KDECADE=KDECADE)  
   END IF
   !
   ! fraction of irrigated surface
   ! ---------------------------
   !
   IF (PRESENT(PIRRIG)) THEN
     IF (SIZE(PIRRIG)>0) &
       CALL AV_PGD (PIRRIG ,PCOVER ,XDATA_IRRIG(:,:),YVEG,'ARI',KDECADE=KDECADE)  
   END IF
   !
   ! water supply for irrigation
   ! ---------------------------
   !
   IF (PRESENT(PWATSUP)) THEN
     IF (SIZE(PWATSUP)>0) &
       CALL AV_PGD (PWATSUP ,PCOVER ,XDATA_WATSUP(:,:),YVEG,'ARI',KDECADE=KDECADE)  
   END IF
!
END IF
IF (LHOOK) CALL DR_HOOK('CONVERT_COVER_ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CONVERT_COVER_ISBA
