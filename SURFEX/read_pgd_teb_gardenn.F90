!     #########
      SUBROUTINE READ_PGD_TEB_GARDEN_n(HPROGRAM)
!     #########################################
!
!!****  *READ_PGD_TEB_GARDEN_n* - routine to initialise ISBA physiographic variables 
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003 
!!      P. Le Moigne  12/2004 : add type of photosynthesis
!!      B. Decharme      2008 : add XWDRAIN
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TEB_GARDEN_n,    ONLY : NPATCH, CISBA, CDIF, CPEDOTF,                &
                                   CPHOTO, CRUNOFF, XCLAY, XSAND,             &
                                   NGROUND_LAYER, NNBIOMASS,                  &
                                   XAOSIP, XAOSIM, XAOSJP, XAOSJM,            &
                                   XHO2IP, XHO2IM, XHO2JP, XHO2JM,            &
                                   XSSO_SLOPE, XSSO_STDEV, XRUNOFFB,          &
                                   XZ0EFFJPDIR, XWDRAIN, LPAR_GARDEN
USE MODD_GR_BIOG_GARDEN_n,ONLY : XISOPOT, XMONOPOT
USE MODD_CH_TEB_n,        ONLY : LCH_BIO_FLUX
USE MODD_TEB_GRID_n,      ONLY : NDIM
USE MODD_SURF_PAR,        ONLY : NVERSION
!
USE MODI_READ_PGD_TEB_GARDEN_PAR_n
USE MODI_READ_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_TYPE_DIM_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! calling program
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! Error code after redding
!
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
!
!
CHARACTER(LEN=3)  :: HISBA          ! ISBA version
CHARACTER(LEN=3)  :: HPHOTO         ! ISBA photosynthesis version
INTEGER           :: KGROUND_LAYER  ! number of ground layers (ISBA)
INTEGER           :: KPATCH         ! number of patches (ISBA)          
INTEGER           :: JLAYER         ! loop counter on layers
INTEGER           :: IVERSION       ! surface version
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_PGD_TEB_GARDEN_N',0,ZHOOK_HANDLE)
YRECFM='SIZE_TOWN'
CALL GET_TYPE_DIM_n('TOWN  ',NDIM)
!
!*       2.     Definition of version
!               ---------------------
!
IVERSION = NVERSION
!
!
!*       2.     Initialisation of ISBA options
!               ------------------------------
!
!
YRECFM='TWN_ISBA'
CALL READ_SURF(HPROGRAM,YRECFM,CISBA,IRESP)
!
YRECFM='TWN_THETAPSI'
CALL READ_SURF(HPROGRAM,YRECFM,CDIF,IRESP)
!
YRECFM='TWN_PEDOTF'
CALL READ_SURF(HPROGRAM,YRECFM,CPEDOTF,IRESP)
!
YRECFM='TWN_PHOTO'
CALL READ_SURF(HPROGRAM,YRECFM,CPHOTO,IRESP)
!
YRECFM='TWN_LAYER'
CALL READ_SURF(HPROGRAM,YRECFM,NGROUND_LAYER,IRESP)
!
!* number of biomass pools
!
YRECFM='TWN_NBIOMASS'
CALL READ_SURF(HPROGRAM,YRECFM,NNBIOMASS,IRESP)
!
NPATCH = 1
!
!
!*       3.     Physiographic data fields:
!               -------------------------
!
ALLOCATE(XAOSIP    (NDIM)              )
ALLOCATE(XAOSIM    (NDIM)              )
ALLOCATE(XAOSJP    (NDIM)              )
ALLOCATE(XAOSJM    (NDIM)              )
ALLOCATE(XHO2IP    (NDIM)              )
ALLOCATE(XHO2IM    (NDIM)              )
ALLOCATE(XHO2JP    (NDIM)              )
ALLOCATE(XHO2JM    (NDIM)              )
ALLOCATE(XSSO_SLOPE(NDIM)              )
ALLOCATE(XSSO_STDEV(NDIM)              )
!
XAOSIP     (:) = 0.
XAOSIM     (:) = 0.
XAOSJP     (:) = 0.
XAOSJM     (:) = 0.
XHO2IP     (:) = 0.
XHO2IM     (:) = 0.
XHO2JP     (:) = 0.
XHO2JM     (:) = 0.
XSSO_SLOPE (:) = 0.
XSSO_STDEV (:) = 0.
!
!
!* clay fraction : attention, seul un niveau est present dans le fichier
!* on rempli tout les niveaux de  XCLAY avec les valeurs du fichiers
!
ALLOCATE(XCLAY(NDIM,NGROUND_LAYER))
YRECFM='TWN_CLAY'
CALL READ_SURF(HPROGRAM,YRECFM,XCLAY(:,1),IRESP)
DO JLAYER=2,NGROUND_LAYER
  XCLAY(:,JLAYER)=XCLAY(:,1)
END DO
!
!* sand fraction
!
ALLOCATE(XSAND(NDIM,NGROUND_LAYER))
YRECFM='TWN_SAND'
CALL READ_SURF(HPROGRAM,YRECFM,XSAND(:,1),IRESP)
DO JLAYER=2,NGROUND_LAYER
  XSAND(:,JLAYER)=XSAND(:,1)
END DO
!
!* orographic runoff coefficient
!
ALLOCATE(XRUNOFFB(NDIM))
YRECFM='TWN_RUNOFFB'
CALL READ_SURF(HPROGRAM,YRECFM,XRUNOFFB,IRESP)
!
!* subgrid drainage coefficient
!
ALLOCATE(XWDRAIN(NDIM))
IF (IVERSION<=3) THEN
  XWDRAIN = 0.
ELSE
  YRECFM='TWN_WDRAIN'
  CALL READ_SURF(HPROGRAM,YRECFM,XWDRAIN,IRESP)
ENDIF
!
!-------------------------------------------------------------------------------
!
!* biogenic chemical emissions
!
IF (LCH_BIO_FLUX) THEN
  ALLOCATE(XISOPOT(NDIM))
  YRECFM='EMIS_ISOPOT'
  CALL READ_SURF(HPROGRAM,YRECFM,XISOPOT,IRESP)
  !
  ALLOCATE(XMONOPOT(NDIM))
  YRECFM='EMIS_MONOPOT'
  CALL READ_SURF(HPROGRAM,YRECFM,XMONOPOT,IRESP)
ELSE
  ALLOCATE(XISOPOT (0))
  ALLOCATE(XMONOPOT(0))
END IF
!
!-------------------------------------------------------------------------------
!
!*       4.     Physiographic data fields not to be computed by ecoclimap
!               ---------------------------------------------------------
!
LPAR_GARDEN = .FALSE.
IF (IVERSION>=7) THEN
  YRECFM='PAR_GARDEN'
  CALL READ_SURF(HPROGRAM,YRECFM,LPAR_GARDEN,IRESP)
END IF
!
IF (LPAR_GARDEN) CALL READ_PGD_TEB_GARDEN_PAR_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('READ_PGD_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_PGD_TEB_GARDEN_n
