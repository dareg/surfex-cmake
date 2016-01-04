!     #########
      SUBROUTINE WRITESURF_TEB_GREENROOF_n (DGU, U, GRM, &
                                            HPROGRAM,HPATCH)
!     #####################################
!
!!****  *WRITESURF_TEB_GREENROOF_n* - writes ISBA prognostic fields
!!                        
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    Based on "writesurf_teb_gardenn" 
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
!!      A. Lemonsu & C. de Munck
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
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURFEX_n, ONLY : TEB_GREENROOF_MODEL_t
!
USE MODI_WRITE_SURF
USE MODI_WRITESURF_GR_SNOW
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
!
!
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(TEB_GREENROOF_MODEL_t), INTENT(INOUT) :: GRM
!
 CHARACTER(LEN=6),  INTENT(IN)     :: HPROGRAM ! program calling
 CHARACTER(LEN=3),  INTENT(IN)     :: HPATCH   ! current teb patch
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER                           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=30)                 :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100)                :: YCOMMENT       ! Comment string
 CHARACTER(LEN=14)                 :: YFORM          ! Writing format
 CHARACTER(LEN=4 )                 :: YLVL
!
INTEGER                           :: JLAYER         ! loop counter on soil layers
!
REAL, DIMENSION(:),  ALLOCATABLE  :: ZWORK          ! 2D array to write data in file
!
INTEGER                           :: IWORK          ! Work integer
!
INTEGER                           :: JSV, JNBIOMASS
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
!*       2.     Prognostic fields:
!               -----------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_TEB_GREENROOF_N',0,ZHOOK_HANDLE)
ALLOCATE(ZWORK(SIZE(GRM%TV%R%CUR%XTG,1)))
!
!
!* soil temperatures
!
IWORK=GRM%TV%O%NGROUND_LAYER
!
DO JLAYER=1,IWORK
  WRITE(YLVL,'(I2)') JLAYER
  YRECFM=HPATCH//'GR_TG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YRECFM=ADJUSTL(YRECFM)
  YFORM='(A13,I1.1,A4)'
  IF (JLAYER >= 10)  YFORM='(A13,I2.2,A4)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_TWN_TG_GR',JLAYER,' (K)'
  ZWORK=GRM%TV%R%CUR%XTG(:,JLAYER,1)
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
END DO
!
!* soil liquid water content
!
DO JLAYER=1,GRM%TV%O%NGROUND_LAYER
  WRITE(YLVL,'(I2)') JLAYER
  YRECFM=HPATCH//'GR_WG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YRECFM=ADJUSTL(YRECFM)
  YFORM='(A13,I1.1,A8)'
  IF (JLAYER >= 10)  YFORM='(A13,I2.2,A8)'
  WRITE(YCOMMENT,FMT=YFORM) 'X_Y_TWN_WG_GR',JLAYER,' (m3/m3)'
  ZWORK=GRM%TV%R%CUR%XWG(:,JLAYER,1)
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
END DO
!
!
!* soil ice water content
!
DO JLAYER=1,GRM%TV%O%NGROUND_LAYER
  WRITE(YLVL,'(I2)') JLAYER
  YRECFM=HPATCH//'GR_WGI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
  YRECFM=ADJUSTL(YRECFM)
  YFORM='(A14,I1.1,A8)'
  IF (JLAYER >= 10)  YFORM='(A14,I2.2,A8)'
  WRITE(YCOMMENT,YFORM) 'X_Y_GR_WGI',JLAYER,' (m3/m3)'
  ZWORK=GRM%TV%R%CUR%XWGI(:,JLAYER,1)
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
END DO
!
DEALLOCATE(ZWORK)
! 
!* water intercepted on leaves
!
YRECFM=HPATCH//'GR_WR'
YRECFM=ADJUSTL(YRECFM)
YCOMMENT='X_Y_TWN_WR_GR (kg/m2)'
 CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,GRM%TV%R%CUR%XWR(:,1),IRESP,HCOMMENT=YCOMMENT)
!
!* Leaf Area Index
!
IF (GRM%TV%O%CPHOTO/='NON' .AND. GRM%TV%O%CPHOTO/='AGS' .AND. GRM%TV%O%CPHOTO/='AST') THEN
  YRECFM=HPATCH//'GR_LAI'
  YRECFM=ADJUSTL(YRECFM)
  YCOMMENT='X_Y_GR_LAI (m2/m2)'
 CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,GRM%TV%M%T%CUR%XLAI(:,1),IRESP,HCOMMENT=YCOMMENT)
END IF
!
!
!* biomass
!
IF (GRM%TV%O%CPHOTO=='NIT') THEN
  DO JNBIOMASS=1,GRM%TV%O%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    YRECFM=HPATCH//'GR_BIOMA'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YRECFM=ADJUSTL(YRECFM)
    YFORM='(A11,I1.1,A8)'
    WRITE(YCOMMENT,FMT=YFORM) 'X_Y_BIOMASS',JNBIOMASS,' (kg/m2)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,GRM%TV%R%CUR%XBIOMASS(:,JNBIOMASS,1),IRESP,HCOMMENT=YCOMMENT)
  END DO
  !
  !
  DO JNBIOMASS=2,GRM%TV%O%NNBIOMASS
    WRITE(YLVL,'(I1)') JNBIOMASS
    YRECFM=HPATCH//'GR_RESPI'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YRECFM=ADJUSTL(YRECFM)
    YFORM='(A16,I1.1,A10)'
    WRITE(YCOMMENT,FMT=YFORM) 'X_Y_RESP_BIOMASS',JNBIOMASS,' (kg/m2/s)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,GRM%TV%R%CUR%XRESP_BIOMASS(:,JNBIOMASS,1),IRESP,HCOMMENT=YCOMMENT)
  END DO
END IF
!
!* aerodynamical resistance
!
!
YRECFM=HPATCH//'GR_RESA'
YRECFM=ADJUSTL(YRECFM)
YCOMMENT='X_Y_GR_RESA (s/m)'
 CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,GRM%TV%R%CUR%XRESA(:,1),IRESP,HCOMMENT=YCOMMENT)
!
!* snow mantel
!
YRECFM='GR'
 CALL WRITESURF_GR_SNOW(DGU, U, &
                        HPROGRAM,YRECFM,HPATCH,GRM%TV%R%CUR%TSNOW)
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_TEB_GREENROOF_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_TEB_GREENROOF_n
