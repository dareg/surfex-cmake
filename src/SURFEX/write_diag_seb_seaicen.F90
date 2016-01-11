!     #########
      SUBROUTINE WRITE_DIAG_SEB_SEAICE_n (DTCO, DGU, U, DGS, DGSI, DGSIC, S, &
                                          HPROGRAM)
!     #################################
!
!!****  *WRITE_DIAG_SEB_SEAICE_n* - write the seaice diagnostic fields
!!
!!    PURPOSE
!!    -------
!!
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      S.Senesi                *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2014
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
!
USE MODD_SFX_OASIS,      ONLY : LCPL_SEAICE
!
!
!
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_t), INTENT(INOUT) :: DGU
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(DIAG_t), INTENT(INOUT) :: DGS
TYPE(DIAG_t), INTENT(INOUT) :: DGSI
TYPE(DIAG_t), INTENT(INOUT) :: DGSIC
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
CHARACTER(LEN=2)  :: YNUM
INTEGER           :: JSV, JSW
!
REAL(KIND=JPRB)   :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEB_SEAICE_N',0,ZHOOK_HANDLE)
!
CALL INIT_IO_SURF_n(DTCO, DGU, U, &
                    HPROGRAM,'SEA   ','SEAFLX','WRITE')
!
!
!*       8.2.     Richardson number :
!               -----------------
IF (DGS%N2M>=1) THEN
   !
   YRECFM='RI_SEAICE'
   YCOMMENT='X_Y_'//YRECFM
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XRI(:),IRESP,HCOMMENT=YCOMMENT)
   !
END IF
!
!*       8.3     Energy fluxes :
!               -------------
!
IF (DGS%LSURF_BUDGET) THEN

   YRECFM='RN_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XRN(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='H_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XH(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='LE_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XLE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='GFLX_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XGFLUX(:),IRESP,HCOMMENT=YCOMMENT)
   !
   IF (DGS%LRAD_BUDGET) THEN
      !
      YRECFM='SWU_SEAICE'
      YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
      !
      CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XSWU(:),IRESP,HCOMMENT=YCOMMENT)
      !
      YRECFM='LWU_SEAICE'
      YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
      !
      CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XLWU(:),IRESP,HCOMMENT=YCOMMENT)
      !
      DO JSW=1, SIZE(DGSI%XSWBU,2)
         YNUM=ACHAR(48+JSW)
         !
         YRECFM='SWU_SEAICE_'//YNUM
         YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
         !
         CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XSWBU(:,JSW),IRESP,HCOMMENT=YCOMMENT)
         !
      ENDDO
      !
   ENDIF
   !
   YRECFM='FMU_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (kg/ms2)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XFMU(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='FMV_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (kg/ms2)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XFMV(:),IRESP,HCOMMENT=YCOMMENT)
   !
END IF
!
IF (DGS%LSURF_BUDGETC) THEN
   !
   YRECFM='RNC_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSIC%XRN(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='HC_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSIC%XH(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='LEC_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSIC%XLE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='GFLXC_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSIC%XGFLUX(:),IRESP,HCOMMENT=YCOMMENT)
   IF (DGS%LRAD_BUDGET .OR. (DGS%LSURF_BUDGETC .AND. .NOT.DGU%LRESET_BUDGETC)) THEN
      !
      YRECFM='SWUC_SEAICE'
      YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
      !
      CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSIC%XSWU(:),IRESP,HCOMMENT=YCOMMENT)
      !
      YRECFM='LWUC_SEAICE'
      YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
      !
      CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSIC%XLWU(:),IRESP,HCOMMENT=YCOMMENT)
      !
   ENDIF
   !
   YRECFM='FMUC_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (kg/ms)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSIC%XFMU(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='FMVC_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (kg/ms)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSIC%XFMV(:),IRESP,HCOMMENT=YCOMMENT)
   !
END IF
!
!*       8.4     transfer coefficients
!               ---------------------
!
IF (DGS%LCOEF) THEN
   !
   YRECFM='CD_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (W/s2)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XCD(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='CH_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (W/s)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XCH(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='Z0_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (M)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XZ0(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='Z0H_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (M)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XZ0H(:),IRESP,HCOMMENT=YCOMMENT)
   !
END IF
!
!
!*       8.5     Surface humidity
!               ----------------
!
IF (DGS%LSURF_VARS) THEN
   YRECFM='QS_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (KG/KG)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XQS(:),IRESP,HCOMMENT=YCOMMENT)
   !
ENDIF
!

!
!*       8.6.     parameters at 2 and 10 meters :
!               -----------------------------
!
IF (DGS%N2M>=1) THEN
   !
   YRECFM='T2M_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (K)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XT2M(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='Q2M_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (KG/KG)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XQ2M(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='HU2M_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (-)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XHU2M(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='ZON10M_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (M/S)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XZON10M(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='MER10M_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (M/S)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XMER10M(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='W10M_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (M/S)'
   !
   CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGSI%XWIND10M(:),IRESP,HCOMMENT=YCOMMENT)
   !
END IF
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)

IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEB_SEAICE_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE WRITE_DIAG_SEB_SEAICE_n
