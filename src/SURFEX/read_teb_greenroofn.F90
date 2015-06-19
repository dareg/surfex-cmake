!     #########
      SUBROUTINE READ_TEB_GREENROOF_n(HPROGRAM,HPATCH)
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
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_CO2V_PAR,          ONLY : XANFMINIT, XCONDCTMIN
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TGRP => TEB_GREENROOF_PGD
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TGRPE => TEB_GREENROOF_PGD_EVOL
USE MODD_TEB_GREENROOF_n, ONLY : TGR => TEB_GREENROOF
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
IWORK = TGRO%NLAYER_GR
!
DO JLAYER=1,IWORK
  WRITE(YLVL,'(I2)') JLAYER
  YRECFM=HPATCH//'GR_TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,ZWORK(:),IRESP)
  TGR%CUR%XTG(:,JLAYER) = ZWORK
END DO
!
!
!* soil liquid water content
!
DO JLAYER=1,TGRO%NLAYER_GR
  WRITE(YLVL,'(I2)') JLAYER
  YRECFM=HPATCH//'GR_WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,ZWORK(:),IRESP)
  TGR%CUR%XWG(:,JLAYER) = ZWORK
END DO
!
!* soil ice water content
!
DO JLAYER=1,TGRO%NLAYER_GR
  WRITE(YLVL,'(I2)') JLAYER
  YRECFM=HPATCH//'GR_WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,ZWORK(:),IRESP)
  TGR%CUR%XWGI(:,JLAYER) = ZWORK
END DO
!
!* water intercepted on leaves
!
YRECFM=HPATCH//'GR_WR'
YRECFM=ADJUSTL(YRECFM)
 CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,TGR%CUR%XWR(:),IRESP)
!
!* Leaf Area Index
!
IF (TVG%CPHOTO=='LAI' .OR. TVG%CPHOTO=='LST' .OR. TVG%CPHOTO=='NIT' .OR. TVG%CPHOTO=='NCB') THEN
  YRECFM = HPATCH//'GR_LAI'
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,TGRPE%CUR%XLAI(:),IRESP)
END IF
!
!* snow mantel
!
 CALL READ_GR_SNOW(HPROGRAM,'GR',HPATCH,ILU,1,TGR%CUR%TSNOW  )! GROO:GreenROOf 
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
TGR%CUR%XRESA(:) = 100.
 CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,TGR%CUR%XRESA(:),IRESP)
!
TGR%CUR%XLE(:) = XUNDEF
!
!* ISBA-AGS variables
!
IF (TVG%CPHOTO/='NON') THEN
  TGR%CUR%XAN(:)    = 0.
  TGR%CUR%XANDAY(:) = 0.
  TGR%CUR%XANFM(:)  = XANFMINIT
  TGR%CUR%XLE(:)    = 0.
END IF
!
IF (TVG%CPHOTO=='AGS' .OR. TVG%CPHOTO=='AST') THEN
  TGR%CUR%XBIOMASS(:,:)      = 0.
  TGR%CUR%XRESP_BIOMASS(:,:) = 0.
ELSEIF (TVG%CPHOTO=='LAI' .OR. TVG%CPHOTO=='LST') THEN
  TGR%CUR%XBIOMASS(:,1)      = TGRP%XBSLAI(:) * TGRPE%CUR%XLAI(:)
  TGR%CUR%XRESP_BIOMASS(:,:) = 0.
ELSEIF (TVG%CPHOTO=='NIT') THEN
  TGR%CUR%XBIOMASS(:,:) = 0.
  DO JNBIOMASS=1,TVG%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    YRECFM=HPATCH//'GR_BIOMA'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,TGR%CUR%XBIOMASS(:,JNBIOMASS),IRESP)
  END DO

  TGR%CUR%XRESP_BIOMASS(:,:) = 0.
  DO JNBIOMASS=2,TVG%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    YRECFM=HPATCH//'GR_RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL READ_SURF(IOB, &
                 HPROGRAM,YRECFM,TGR%CUR%XRESP_BIOMASS(:,JNBIOMASS),IRESP)
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
