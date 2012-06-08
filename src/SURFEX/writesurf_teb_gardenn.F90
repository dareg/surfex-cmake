!     #########
      SUBROUTINE WRITESURF_TEB_GARDEN_n(HPROGRAM)
!     #####################################
!
!!****  *WRITESURF_TEB_GARDEN_n* - writes ISBA prognostic fields
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
!!      P. LeMoigne 12/2004 : correct dimensionning if more than 10 layers in
!!                            the soil (diffusion version)
!!      B. Decharme  2008    : Floodplains
!!      B. Decharme  01/2009 : Optional Arpege deep soil temperature write
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TEB_GARDEN_n, ONLY :   NGROUND_LAYER, CISBA, CPHOTO, CRESPSL,       &
                                NNBIOMASS, NNLITTER, NNSOILCARB, NNLITTLEVS, &
                                XTG, XWG, XWGI, XWR, XLAI, TSNOW,            &
                                XRESA, XAN, XANFM, XLE, XANDAY,              &
                                XRESP_BIOMASS, XBIOMASS, NWG_LAYER,          &
                                XLITTER, XSOILCARB, XLIGNIN_STRUC                                
!
USE MODD_SURF_PAR, ONLY : NUNDEF
!
USE MODI_WRITE_SURF
USE MODI_WRITESURF_GR_SNOW
USE MODD_DST_n
USE MODD_DST_SURF
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
CHARACTER(LEN=14) :: YFORM          ! Writing format
CHARACTER(LEN=4 ) :: YLVL
!
INTEGER :: JLAYER ! loop counter on soil layers
!
REAL, DIMENSION(:,:),ALLOCATABLE  :: ZWORK      ! 2D array to write data in file
!
INTEGER :: IWORK   ! Work integer
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
!*       2.     Prognostic fields:
!               -----------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_TEB_GARDEN_N',0,ZHOOK_HANDLE)
ALLOCATE(ZWORK(SIZE(XTG,1),SIZE(XTG,3)))
!* soil temperatures
!
IWORK=NGROUND_LAYER
!
DO JLAYER=1,IWORK
  WRITE(YLVL,'(I4)') JLAYER
  YRECFM='TWN_TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A11,I1.1,A4)'
  IF (JLAYER >= 10)  YFORM='(A11,I2.2,A4)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_TWN_TG',JLAYER,' (K)'
  ZWORK=XTG(:,JLAYER,:)
  CALL WRITE_SURF(HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
END DO
!
!
!* soil liquid water content
!
IF(CISBA=='DIF')THEN
  IWORK=MAXVAL(NWG_LAYER(:,:),NWG_LAYER(:,:)/=NUNDEF)
ELSE
  IWORK=NGROUND_LAYER
ENDIF
!
DO JLAYER=1,IWORK
  WRITE(YLVL,'(I4)') JLAYER
  YRECFM='TWN_WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A11,I1.1,A8)'
  IF (JLAYER >= 10)  YFORM='(A11,I2.2,A8)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_TWN_WG',JLAYER,' (m3/m3)'
  ZWORK=XWG(:,JLAYER,:)
  CALL WRITE_SURF(HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
END DO
!
!
!* soil ice water content
!
DO JLAYER=1,IWORK
  WRITE(YLVL,'(I4)') JLAYER
  YRECFM='TWN_WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YFORM='(A12,I1.1,A8)'
  IF (JLAYER >= 10)  YFORM='(A12,I2.2,A8)'
  WRITE(YCOMMENT,YFORM) 'X_Y_TWN_WGI',JLAYER,' (m3/m3)'
  ZWORK=XWGI(:,JLAYER,:)
  CALL WRITE_SURF(HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
END DO
!
DEALLOCATE(ZWORK)
!
!* water intercepted on leaves
!
YRECFM='TWN_WR'
YCOMMENT='X_Y_TWN_WR (kg/m2)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XWR(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!* Leaf Area Index
!
IF (CPHOTO/='NON' .AND. CPHOTO/='AGS' .AND. CPHOTO/='AST') THEN
  YRECFM='TWN_LAI'
  YCOMMENT='X_Y_TWN_LAI (m2/m2)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XLAI(:,:),IRESP,HCOMMENT=YCOMMENT)
END IF
!
!
!* aerodynamical resistance
!
YRECFM='TWN_RESA'
YCOMMENT='X_Y_TWN_RESA (s/m)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XRESA(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!* snow mantel
!
CALL WRITESURF_GR_SNOW(HPROGRAM,'GARD',TSNOW)
IF (LHOOK) CALL DR_HOOK('WRITESURF_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_TEB_GARDEN_n
