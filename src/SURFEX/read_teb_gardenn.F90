!     #########
      SUBROUTINE READ_TEB_GARDEN_n(HPROGRAM)
!     ##################################
!
!!****  *READ_TEB_GARDEN_n* - routine to initialise ISBA variables
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!!
!!      READ_SURF for general reading : 08/2003 (S.Malardel)
!!      B. Decharme  2008    : Floodplains
!!      B. Decharme  01/2009 : Optional Arpege deep soil temperature read
!!      B. Decharme  09/2012 : suppress NWG_LAYER (parallelization problems)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_CO2V_PAR,       ONLY : XANFMINIT, XCONDCTMIN
USE MODD_TEB_GARDEN_n,   ONLY : NGROUND_LAYER, NPATCH,              &
                                CPHOTO, CRESPSL, NNBIOMASS,         &
                                XTG, XWG, XWGI, XWR, XLAI, TSNOW,   &
                                XRESA, XANFM, XANF, XAN, XLE, XANDAY,&
                                XBSLAI, XBIOMASS, XRESP_BIOMASS,    &
                                XLITTER, XSOILCARB, XLIGNIN_STRUC,  &
                                NNLITTER, NNLITTLEVS, NNSOILCARB  
!                                
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_SNOW_PAR,       ONLY : XZ0SN
!
USE MODI_READ_SURF
!
USE MODI_INIT_IO_SURF_n
USE MODI_SET_SURFEX_FILEIN
USE MODI_END_IO_SURF_n
USE MODI_TOWN_PRESENCE
USE MODI_ALLOCATE_GR_SNOW
USE MODI_READ_GR_SNOW
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
LOGICAL           :: GTOWN          ! town variables written in the file
INTEGER           :: ILU            ! 1D physical dimension
INTEGER           :: IRESP          ! Error code after redding
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=4)  :: YLVL
REAL, DIMENSION(:,:),ALLOCATABLE  :: ZWORK      ! 2D array to write data in file
!
INTEGER :: IWORK   ! Work integer
!
INTEGER :: JLAYER  ! loop counter on layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_TEB_GARDEN_N',0,ZHOOK_HANDLE)
YRECFM='SIZE_TOWN'
CALL GET_TYPE_DIM_n('TOWN  ',ILU)
!
!
!*       2.     Prognostic fields:
!               -----------------
!
ALLOCATE(ZWORK(ILU,NPATCH))
!* soil temperatures
!
IWORK=NGROUND_LAYER
ALLOCATE(XTG(ILU,IWORK,NPATCH))
!
DO JLAYER=1,IWORK
  WRITE(YLVL,'(I4)') JLAYER
  YRECFM='TWN_TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  XTG(:,JLAYER,:)=ZWORK
END DO
!
!
!* soil liquid water content
!
ALLOCATE(XWG(ILU,IWORK,NPATCH))
DO JLAYER=1,NGROUND_LAYER
  WRITE(YLVL,'(I4)') JLAYER
  YRECFM='TWN_WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  XWG(:,JLAYER,:)=ZWORK
END DO
!
!* soil ice water content
!
ALLOCATE(XWGI(ILU,IWORK,NPATCH))
DO JLAYER=1,NGROUND_LAYER
  WRITE(YLVL,'(I4)') JLAYER
  YRECFM='TWN_WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  CALL READ_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP)
  XWGI(:,JLAYER,:)=ZWORK
END DO
!
!* water intercepted on leaves
!
ALLOCATE(XWR(ILU,NPATCH))
YRECFM = 'TWN_WR'
CALL READ_SURF(HPROGRAM,YRECFM,XWR(:,:),IRESP)
!
!* Leaf Area Index
!
IF (CPHOTO=='LAI' .OR. CPHOTO=='LST' .OR. CPHOTO=='NIT' .OR. CPHOTO=='NCB') THEN
  XLAI(:,:) = 0.
END IF
!
!* snow mantel
!
CALL END_IO_SURF_n(HPROGRAM)
CALL SET_SURFEX_FILEIN(HPROGRAM,'PGD ')
CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
!
CALL TOWN_PRESENCE(HPROGRAM,GTOWN)
!
CALL END_IO_SURF_n(HPROGRAM)
CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP')
CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
!
IF (.NOT. GTOWN) THEN
  TSNOW%SCHEME='1-L'
  CALL ALLOCATE_GR_SNOW(TSNOW,ILU,NPATCH)
ELSE
  CALL READ_GR_SNOW(HPROGRAM,'GARD',ILU,NPATCH,TSNOW  )
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.  Semi-prognostic variables
!            -------------------------
!
!* aerodynamical resistance
!
ALLOCATE(XRESA(ILU,NPATCH))
YRECFM = 'TWN_RESA'
XRESA(:,:) = 100.
CALL READ_SURF(HPROGRAM,YRECFM,XRESA(:,:),IRESP)
!
ALLOCATE(XLE(ILU,NPATCH))
XLE(:,:) = XUNDEF
!
!* ISBA-AGS variables
!
IF (CPHOTO/='NON') THEN
  ALLOCATE(XAN   (ILU,NPATCH)) 
  ALLOCATE(XANDAY(ILU,NPATCH)) 
  ALLOCATE(XANFM (ILU,NPATCH))
  ALLOCATE(XANF  (ILU,NPATCH))  
  XAN(:,:)    = 0.
  XANDAY(:,:) = 0.
  XANFM(:,:)  = XANFMINIT
  XLE(:,:)    = 0.
ELSE
  ALLOCATE(XAN   (0,0)) 
  ALLOCATE(XANDAY(0,0)) 
  ALLOCATE(XANFM (0,0))
ENDIF
!
IF(CPHOTO/='NON') THEN
  ALLOCATE(XBIOMASS         (ILU,NNBIOMASS,NPATCH))
  ALLOCATE(XRESP_BIOMASS    (ILU,NNBIOMASS,NPATCH))
ELSE
  ALLOCATE(XBIOMASS         (0,0,0))
  ALLOCATE(XRESP_BIOMASS    (0,0,0))
END IF
!
IF (CPHOTO=='AGS' .OR. CPHOTO=='AST') THEN
  !
  XBIOMASS(:,:,:) = 0.
  XRESP_BIOMASS(:,:,:) = 0.
ELSEIF (CPHOTO=='LAI' .OR. CPHOTO=='LST') THEN
  !
  XBIOMASS(:,1,:) = XBSLAI(:,:) * XLAI(:,:)
  XRESP_BIOMASS(:,:,:) = 0.
ELSEIF (CPHOTO=='NIT') THEN
  !
  XBIOMASS(:,:,:) = 0.
  XRESP_BIOMASS(:,:,:) = 0.
ELSEIF (CPHOTO=='NCB') THEN
  !
  XBIOMASS(:,:,:) = 0.
  XRESP_BIOMASS(:,:,:) = 0.
  !
ENDIF
!
!
!*       6. Soil carbon
!
!
IF (CRESPSL=='CNT') THEN
  !
  ALLOCATE(XLITTER          (ILU,NNLITTER,NNLITTLEVS,NPATCH))
  ALLOCATE(XSOILCARB        (ILU,NNSOILCARB,NPATCH))
  ALLOCATE(XLIGNIN_STRUC    (ILU,NNLITTLEVS,NPATCH))  
  XLITTER(:,:,:,:) = 0.
  XSOILCARB(:,:,:) = 0.
  XLIGNIN_STRUC(:,:,:) = 0.
!
ENDIF
!
!
DEALLOCATE(ZWORK)
IF (LHOOK) CALL DR_HOOK('READ_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_TEB_GARDEN_n
