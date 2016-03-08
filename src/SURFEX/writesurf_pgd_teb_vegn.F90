!     #########
      SUBROUTINE WRITESURF_PGD_TEB_VEG_n (HSELECT, KTIME, TGDO, TGDP, HPROGRAM)
!     ###############################################
!
!!****  *WRITE_PGD_TEB_VEG_n* - writes ISBA fields describing urban gardens
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
!!      A. Lemonsu & C. de Munck   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2011 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
!
USE MODI_WRITE_SURF
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
 CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSELECT
!
INTEGER, INTENT(IN) :: KTIME
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(ISBA_PGD_t), INTENT(INOUT) :: TGDP
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
CHARACTER(LEN=4 ) :: YLVL
!
INTEGER :: JJ, JLAYER
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_TEB_VEG_N',0,ZHOOK_HANDLE)
!
!* soil scheme option
!
YRECFM='GD_ISBA'
YCOMMENT=YRECFM
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,TGDO%CISBA,IRESP,HCOMMENT=YCOMMENT)
!
!* Reference grid for DIF
!
IF(TGDO%CISBA=='DIF') THEN
  DO JLAYER=1,TGDO%NGROUND_LAYER
     WRITE(YLVL,'(I4)') JLAYER     
     YRECFM='GD_SGRID'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     YCOMMENT='Depth of TEB Garden soilgrid layer '//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
     CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,TGDO%XSOILGRID(JLAYER),IRESP,HCOMMENT=YCOMMENT)
  END DO 
ENDIF
!
!* number of soil layers
!
YRECFM='GD_LAYER'
YCOMMENT=YRECFM
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,TGDO%NGROUND_LAYER,IRESP,HCOMMENT=YCOMMENT)
!
!* number of time data for vegetation characteristics (VEG, LAI, EMIS, Z0) 
!
YRECFM='GD_NTIME'
YCOMMENT=YRECFM
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,KTIME,IRESP,HCOMMENT=YCOMMENT)
!
! * clay fraction
!
YRECFM='GD_CLAY'
YCOMMENT='X_Y_GD_CLAY'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,TGDP%XCLAY(:,1),IRESP,HCOMMENT=YCOMMENT)
!        
! * sand fraction
!
YRECFM='GD_SAND'
YCOMMENT='X_Y_GD_SAND'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,TGDP%XSAND(:,1),IRESP,HCOMMENT=YCOMMENT)
!        
! * orographic runoff coefficient
!
YRECFM='GD_RUNOFFB'
YCOMMENT='X_Y_GD_RUNOFFB'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,TGDP%XRUNOFFB,IRESP,HCOMMENT=YCOMMENT)
!        
! * subgrid drainage coefficient
!
YRECFM='GD_WDRAIN'
YCOMMENT='X_Y_GD_WDRAIN'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,YRECFM,TGDP%XWDRAIN,IRESP,HCOMMENT=YCOMMENT)
!
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_TEB_VEG_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_PGD_TEB_VEG_n
