!     #########
      SUBROUTINE READ_TEB_GREENROOF_n (DTCO, U, GRO, GRR, GRMT, HPROGRAM,HPATCH)
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
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_TIME_t
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
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GRO
TYPE(ISBA_PROG_t), INTENT(INOUT) :: GRR
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: GRMT
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
 CALL GET_TYPE_DIM_n(DTCO, U, 'TOWN  ',ILU)
!
!
!*       2.     Prognostic fields:
!               -----------------
!
ALLOCATE(ZWORK(ILU))
!
!* soil temperatures
!
IWORK = GRO%NGROUND_LAYER
!
DO JLAYER=1,IWORK
  WRITE(YLVL,'(I2)') JLAYER
  YRECFM=HPATCH//'GR_TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:),IRESP)
  GRR%XTG(:,JLAYER,1) = ZWORK
END DO
!
!
!* soil liquid water content
!
DO JLAYER=1,GRO%NGROUND_LAYER
  WRITE(YLVL,'(I2)') JLAYER
  YRECFM=HPATCH//'GR_WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:),IRESP)
  GRR%XWG(:,JLAYER,1) = ZWORK
END DO
!
!* soil ice water content
!
DO JLAYER=1,GRO%NGROUND_LAYER
  WRITE(YLVL,'(I2)') JLAYER
  YRECFM=HPATCH//'GR_WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:),IRESP)
  GRR%XWGI(:,JLAYER,1) = ZWORK
END DO
!
!* water intercepted on leaves
!
YRECFM=HPATCH//'GR_WR'
YRECFM=ADJUSTL(YRECFM)
 CALL READ_SURF(HPROGRAM,YRECFM,GRR%XWR(:,1),IRESP)
!
!* Leaf Area Index
!
IF (GRO%CPHOTO=='LAI' .OR. GRO%CPHOTO=='LST' .OR. GRO%CPHOTO=='NIT' .OR. GRO%CPHOTO=='NCB') THEN
  YRECFM = HPATCH//'GR_LAI'
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(HPROGRAM,YRECFM,GRMT%XLAI(:,1),IRESP)
END IF
!
!* snow mantel
!
 CALL READ_GR_SNOW(HPROGRAM,'GR',HPATCH,ILU,1,GRR%TSNOW  )! GROO:GreenROOf 
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
GRR%XRESA(:,1) = 100.
 CALL READ_SURF(HPROGRAM,YRECFM,GRR%XRESA(:,1),IRESP)
!
GRR%XLE(:,1) = XUNDEF
!
!* ISBA-AGS variables
!
IF (GRO%CPHOTO/='NON') THEN
  GRR%XAN(:,1)    = 0.
  GRR%XANDAY(:,1) = 0.
  GRR%XANFM(:,1)  = XANFMINIT
  GRR%XLE(:,1)    = 0.
END IF
!
IF (GRO%CPHOTO=='AGS' .OR. GRO%CPHOTO=='AST') THEN
  GRR%XBIOMASS(:,:,1)      = 0.
  GRR%XRESP_BIOMASS(:,:,1) = 0.
ELSEIF (GRO%CPHOTO=='LAI' .OR. GRO%CPHOTO=='LST') THEN
  GRR%XBIOMASS(:,1,1)      = GRMT%XBSLAI(:,1) * GRMT%XLAI(:,1)
  GRR%XRESP_BIOMASS(:,:,1) = 0.
ELSEIF (GRO%CPHOTO=='NIT') THEN
  GRR%XBIOMASS(:,:,1) = 0.
  DO JNBIOMASS=1,GRO%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    YRECFM=HPATCH//'GR_BIOMA'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YRECFM=ADJUSTL(YRECFM)
    CALL READ_SURF(HPROGRAM,YRECFM,GRR%XBIOMASS(:,JNBIOMASS,1),IRESP)
  END DO

  GRR%XRESP_BIOMASS(:,:,1) = 0.
  DO JNBIOMASS=2,GRO%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    YRECFM=HPATCH//'GR_RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YRECFM=ADJUSTL(YRECFM)
    CALL READ_SURF(HPROGRAM,YRECFM,GRR%XRESP_BIOMASS(:,JNBIOMASS,1),IRESP)
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
