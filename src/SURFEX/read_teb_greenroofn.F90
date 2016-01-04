!     #########
      SUBROUTINE READ_TEB_GREENROOF_n (DTCO, U, GRM, &
                                       HPROGRAM,HPATCH)
!     ##################################
!
!!****  *READ_TEB_GREENROOF_n* - routine to initialise ISBA variables
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
!!    based on read_teb_greenroofn
!!
!!    AUTHOR
!!    ------
!!      C. de Munck & A. Lemonsu *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2011
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURFEX_n, ONLY : TEB_GREENROOF_MODEL_t
!
USE MODD_CO2V_PAR,          ONLY : XANFMINIT, XCONDCTMIN
!                                
USE MODD_SURF_PAR,          ONLY : XUNDEF
USE MODD_SNOW_PAR,          ONLY : XZ0SN
!
USE MODI_READ_SURF
!
USE MODI_READ_GR_SNOW
!
!
!
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
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(TEB_GREENROOF_MODEL_t), INTENT(INOUT) :: GRM

!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! calling program
 CHARACTER(LEN=3),  INTENT(IN)  :: HPATCH   ! current TEB patch identificator
!
!*       0.2   Declarations of local variables
!              -------------------------------
INTEGER           :: ILU                             ! 1D physical dimension
INTEGER           :: IRESP                           ! Error code after redding
INTEGER           :: IWORK                           ! Work integer
INTEGER           :: JLAYER, JNBIOMASS               ! loop counter on layers
 CHARACTER(LEN=30) :: YRECFM                          ! Name of the article to be read
 CHARACTER(LEN=4)  :: YLVL
REAL, DIMENSION(:),ALLOCATABLE  :: ZWORK             ! 2D array to write data in file
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_TEB_GREENROOF_N',0,ZHOOK_HANDLE)
YRECFM='SIZE_TOWN'
 CALL GET_TYPE_DIM_n(DTCO, U, &
                     'TOWN  ',ILU)
!
!
!*       2.     Prognostic fields:
!               -----------------
!
ALLOCATE(ZWORK(ILU))
!
!* soil temperatures
!
IWORK = GRM%TV%O%NGROUND_LAYER
!
DO JLAYER=1,IWORK
  WRITE(YLVL,'(I2)') JLAYER
  YRECFM=HPATCH//'GR_TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,ZWORK(:),IRESP)
  GRM%TV%R%CUR%XTG(:,JLAYER,1) = ZWORK
END DO
!
!
!* soil liquid water content
!
DO JLAYER=1,GRM%TV%O%NGROUND_LAYER
  WRITE(YLVL,'(I2)') JLAYER
  YRECFM=HPATCH//'GR_WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,ZWORK(:),IRESP)
  GRM%TV%R%CUR%XWG(:,JLAYER,1) = ZWORK
END DO
!
!* soil ice water content
!
DO JLAYER=1,GRM%TV%O%NGROUND_LAYER
  WRITE(YLVL,'(I2)') JLAYER
  YRECFM=HPATCH//'GR_WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,ZWORK(:),IRESP)
  GRM%TV%R%CUR%XWGI(:,JLAYER,1) = ZWORK
END DO
!
!* water intercepted on leaves
!
YRECFM=HPATCH//'GR_WR'
YRECFM=ADJUSTL(YRECFM)
 CALL READ_SURF(&
                 HPROGRAM,YRECFM,GRM%TV%R%CUR%XWR(:,1),IRESP)
!
!* Leaf Area Index
!
IF (GRM%TV%O%CPHOTO=='LAI' .OR. GRM%TV%O%CPHOTO=='LST' .OR. GRM%TV%O%CPHOTO=='NIT' .OR. GRM%TV%O%CPHOTO=='NCB') THEN
  YRECFM = HPATCH//'GR_LAI'
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(&
                 HPROGRAM,YRECFM,GRM%TV%M%T%CUR%XLAI(:,1),IRESP)
END IF
!
!* snow mantel
!
 CALL READ_GR_SNOW(&
                   HPROGRAM,'GR',HPATCH,ILU,1,GRM%TV%R%CUR%TSNOW  )! GROO:GreenROOf 
!
!-------------------------------------------------------------------------------
!
!*       4.  Semi-prognostic variables
!            -------------------------
!
!* aerodynamical resistance
!
YRECFM = HPATCH//'GR_RESA'
YRECFM=ADJUSTL(YRECFM)
GRM%TV%R%CUR%XRESA(:,1) = 100.
 CALL READ_SURF(&
                 HPROGRAM,YRECFM,GRM%TV%R%CUR%XRESA(:,1),IRESP)
!
GRM%TV%R%CUR%XLE(:,1) = XUNDEF
!
!* ISBA-AGS variables
!
IF (GRM%TV%O%CPHOTO/='NON') THEN
  GRM%TV%R%CUR%XAN(:,1)    = 0.
  GRM%TV%R%CUR%XANDAY(:,1) = 0.
  GRM%TV%R%CUR%XANFM(:,1)  = XANFMINIT
  GRM%TV%R%CUR%XLE(:,1)    = 0.
END IF
!
IF (GRM%TV%O%CPHOTO=='AGS' .OR. GRM%TV%O%CPHOTO=='AST') THEN
  GRM%TV%R%CUR%XBIOMASS(:,:,1)      = 0.
  GRM%TV%R%CUR%XRESP_BIOMASS(:,:,1) = 0.
ELSEIF (GRM%TV%O%CPHOTO=='LAI' .OR. GRM%TV%O%CPHOTO=='LST') THEN
  GRM%TV%R%CUR%XBIOMASS(:,1,1)      = GRM%TV%M%T%CUR%XBSLAI(:,1) * GRM%TV%M%T%CUR%XLAI(:,1)
  GRM%TV%R%CUR%XRESP_BIOMASS(:,:,1) = 0.
ELSEIF (GRM%TV%O%CPHOTO=='NIT') THEN
  GRM%TV%R%CUR%XBIOMASS(:,:,1) = 0.
  DO JNBIOMASS=1,GRM%TV%O%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    YRECFM=HPATCH//'GR_BIOMA'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YRECFM=ADJUSTL(YRECFM)
    CALL READ_SURF(&
                 HPROGRAM,YRECFM,GRM%TV%R%CUR%XBIOMASS(:,JNBIOMASS,1),IRESP)
  END DO

  GRM%TV%R%CUR%XRESP_BIOMASS(:,:,1) = 0.
  DO JNBIOMASS=2,GRM%TV%O%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    YRECFM=HPATCH//'GR_RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YRECFM=ADJUSTL(YRECFM)
    CALL READ_SURF(&
                 HPROGRAM,YRECFM,GRM%TV%R%CUR%XRESP_BIOMASS(:,JNBIOMASS,1),IRESP)
  END DO
ENDIF
!
!
DEALLOCATE(ZWORK)
IF (LHOOK) CALL DR_HOOK('READ_TEB_GREENROOF_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_TEB_GREENROOF_n
