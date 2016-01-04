!     #########
      SUBROUTINE INIT_FROM_DATA_GRDN_n (DTGD, &
                                        KDECADE, HPHOTO,                               &
                                               PVEG,                                  &
                                               PLAI,PRSMIN,PGAMMA,PWRMAX_CF,          &
                                               PRGL,PCV,PDG,PD_ICE,PZ0,PZ0_O_Z0H,     &
                                               PALBNIR_VEG,PALBVIS_VEG,PALBUV_VEG,    &
                                               PEMIS,                                 &
                                               PVEGTYPE,PROOTFRAC,                    &
                                               PGMES,PBSLAI,PLAIMIN,PSEFOLD,PGC,      &
                                               PDMAX, PF2I, OSTRESS,                  &
                                               PH_TREE, PRE25,                        &
                                               PCE_NITRO, PCF_NITRO, PCNA_NITRO,      &
                                               PALBNIR_SOIL,PALBVIS_SOIL,PALBUV_SOIL  )  
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
!     
!!      P Le Moigne 09/2005 AGS modifs of L. Jarlan
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!

!
!
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
!
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTGD
!
INTEGER,                INTENT(IN)    :: KDECADE
 CHARACTER(LEN=*),       INTENT(IN)    :: HPHOTO  ! type of photosynthesis
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
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PEMIS
!
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PVEGTYPE
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
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PALBNIR_SOIL
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PALBVIS_SOIL
REAL, DIMENSION(:,:),   OPTIONAL, INTENT(OUT)   :: PALBUV_SOIL
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: ITIME
INTEGER :: ILUOUT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*    1.      TIME INITIALIZATION
!             -------------------
!
! data every month
IF (LHOOK) CALL DR_HOOK('INIT_FROM_DATA_GRDN_N',0,ZHOOK_HANDLE)
IF (DTGD%NTIME==12) THEN
  ITIME = (KDECADE+2)/3    
ELSEIF (DTGD%NTIME==1) THEN
  ITIME = 1
ENDIF
!
!*    2.      SECONDARY VARIABLES
!             -------------------
!
!*    2.1     fields on natural surfaces only, taking into account patches/ 
!             -------------------------------
!
!
IF (PRESENT(PH_TREE)) THEN
  IF (SIZE(PH_TREE)>0) PH_TREE = DTGD%XPAR_H_TREE
ENDIF
!
IF (PRESENT(PVEGTYPE)) PVEGTYPE = DTGD%XPAR_VEGTYPE
!
! vegetation fraction
! -------------------
!
IF (PRESENT(PVEG)) PVEG =  DTGD%XPAR_VEG (:,ITIME,:)
!
! Leaf Aera Index
! ---------------
!
IF (PRESENT(PLAI)) PLAI = DTGD%XPAR_LAI (:,ITIME,:)
!
! roughness length
! ----------------
!
IF (PRESENT(PZ0)) PZ0 =  DTGD%XPAR_Z0 (:,ITIME,:)
!
IF (PRESENT(PZ0_O_Z0H)) PZ0_O_Z0H = DTGD%XPAR_Z0_O_Z0H
!
!
!emis-eco
!--------
!
IF (PRESENT(PEMIS)) PEMIS =  DTGD%XPAR_EMIS (:,ITIME,:)
! 
!---------------------------------------------------------------------------------
! 
!* 1/Rsmin
!
IF (PRESENT(PRSMIN)) THEN
  IF (SIZE(PRSMIN)>0) PRSMIN = DTGD%XPAR_RSMIN
END IF
!
!* other vegetation parameters
!
IF (PRESENT(PGAMMA)) PGAMMA = DTGD%XPAR_GAMMA
IF (PRESENT(PWRMAX_CF)) PWRMAX_CF = DTGD%XPAR_WRMAX_CF
!
!
IF (PRESENT(PRGL)) PRGL = DTGD%XPAR_RGL
IF (PRESENT(PCV)) PCV = DTGD%XPAR_CV
!
!---------------------------------------------------------------------------------
!
!* soil layers
!  -----------
!
IF (PRESENT(PDG)) PDG = DTGD%XPAR_DG
!
!* cumulative root fraction
!
IF (PRESENT(PROOTFRAC)) THEN
  IF (SIZE(PROOTFRAC)>0) PROOTFRAC = DTGD%XPAR_ROOTFRAC
ENDIF
!
!* soil ice for runoff
!
IF (PRESENT(PD_ICE)) PD_ICE = DTGD%XPAR_DICE
!
!---------------------------------------------------------------------------------
IF (PRESENT(PALBNIR_VEG)) PALBNIR_VEG = DTGD%XPAR_ALBNIR_VEG
IF (PRESENT(PALBVIS_VEG)) PALBVIS_VEG = DTGD%XPAR_ALBVIS_VEG
IF (PRESENT(PALBUV_VEG)) PALBUV_VEG = DTGD%XPAR_ALBUV_VEG

IF (PRESENT(PALBNIR_SOIL)) PALBNIR_SOIL = DTGD%XPAR_ALBNIR_SOIL
IF (PRESENT(PALBVIS_SOIL)) PALBVIS_SOIL = DTGD%XPAR_ALBVIS_SOIL
IF (PRESENT(PALBUV_SOIL)) PALBUV_SOIL = DTGD%XPAR_ALBUV_SOIL

IF (PRESENT(PGMES)) THEN
  IF (SIZE(PGMES)>0) PGMES = DTGD%XPAR_GMES
END IF

IF (PRESENT(PBSLAI)) THEN
  IF (SIZE(PBSLAI)>0) PBSLAI = DTGD%XPAR_BSLAI
END IF

IF (PRESENT(PSEFOLD)) THEN
  IF (SIZE(PSEFOLD)>0) PSEFOLD = DTGD%XPAR_SEFOLD
END IF

IF (PRESENT(PGC)) THEN
  IF (SIZE(PGC)>0) PGC = DTGD%XPAR_GC
END IF

IF (PRESENT(PDMAX)) THEN
  IF (SIZE(PDMAX)>0) PDMAX = DTGD%XPAR_DMAX
END IF

IF (PRESENT(PRE25)) THEN
  IF (SIZE(PRE25)>0) PRE25 = DTGD%XPAR_RE25
END IF

IF (PRESENT(PLAIMIN)) THEN
  IF (SIZE(PLAIMIN)>0) PLAIMIN = DTGD%XPAR_LAIMIN
END IF

IF (PRESENT(PCE_NITRO)) THEN
  IF (SIZE(PCE_NITRO)>0) PCE_NITRO = DTGD%XPAR_CE_NITRO
END IF

IF (PRESENT(PCF_NITRO)) THEN
  IF (SIZE(PCF_NITRO)>0) PCF_NITRO = DTGD%XPAR_CF_NITRO
END IF

IF (PRESENT(PCNA_NITRO)) THEN
  IF (SIZE(PCNA_NITRO)>0) PCNA_NITRO = DTGD%XPAR_CNA_NITRO
END IF

IF (PRESENT(PF2I)) THEN
  IF (SIZE(PF2I)>0) PF2I = DTGD%XPAR_F2I
END IF
!
IF (PRESENT(OSTRESS)) THEN
  IF (SIZE(OSTRESS)>0) OSTRESS = DTGD%LPAR_STRESS
END IF
IF (LHOOK) CALL DR_HOOK('INIT_FROM_DATA_GRDN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_FROM_DATA_GRDN_n
