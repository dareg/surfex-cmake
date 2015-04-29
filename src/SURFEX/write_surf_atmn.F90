!     ####################################
      SUBROUTINE WRITE_SURF_ATM_n(HPROGRAM,HWRITE,OLAND_USE)
!     ####################################
!
!!****  *WRITE_SURF_ATM_n* - routine to write surface variables 
!!                           in their respective files or in file
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!!      Modified    06/2007, P.LeMoigne: do not write pgd fields in
!!                                       historical files
!!      Modified    03/2009, B.Decharme: keys for arrange cover
!!      Modified    04/2009, B.Decharme: write precipitation forcing into the restart file for ARPEGE/ALADIN run
!       Modified    06/2009, B.Decharme: flag to desactivate writing of horizontal grid 
!       Modified    08/2009, B.Decharme: BUDGETC for all tiles
!       Modified    07/2011, B.Decharme: delete write pgd fields
!       Modified    07/2011, B.Decharme: land_use key for writing semi-prognostic variables
!       Modified    05/2012, B.Decharme: supress LPROVAR_TO_DIAG to write prognostic fields if user want
!       Modified    05/2013, B.Decharme: WRITESURF_PRECIP becomes WRITESURF_CPL_GCM
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SSO_CANOPY_n, ONLY : SSCP => SSO_CANOPY
!
USE MODD_SURF_CONF,       ONLY : CPROGNAME
USE MODD_SURF_PAR,        ONLY : NVERSION, NBUGFIX
USE MODD_WRITE_SURF_ATM,  ONLY : LNOWRITE_CANOPY
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
USE MODD_SURF_ATM_SSO_n, ONLY : USS => SURF_ATM_SSO
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM  
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_WRITE_SEA_n
USE MODI_WRITE_INLAND_WATER_n
USE MODI_WRITE_NATURE_n
USE MODI_WRITE_TOWN_n
USE MODI_END_IO_SURF_n
USE MODI_WRITE_GRID
!
USE MODI_WRITESURF_ATM_CONF_n
USE MODI_WRITESURF_SSO_CANOPY_n
USE MODI_WRITESURF_CPL_GCM_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),    INTENT(IN)  :: HWRITE    ! 'PREP' : does not write SBL XUNDEF fields
!                                             ! 'ALL' : all fields are written
LOGICAL,             INTENT(IN)  :: OLAND_USE !
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
 CHARACTER(LEN=100) :: YCOMMENT
INTEGER            :: IRESP
LOGICAL            :: LSAVE_SELECT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_SURF_ATM_N',0,ZHOOK_HANDLE)
CPROGNAME = HPROGRAM
!
!*       1.     Configuration and cover fields:
!               ------------------------------
!
!
!         Initialisation for IO
!
 CALL INIT_IO_SURF_n(HPROGRAM,'FULL  ','SURF  ','WRITE')
!
LSAVE_SELECT=DGU%LSELECT
DGU%LSELECT     =.FALSE.
!
YCOMMENT='(-)'
 CALL WRITE_SURF(HPROGRAM,'VERSION',NVERSION,IRESP,YCOMMENT)
 CALL WRITE_SURF(HPROGRAM,'BUG    ',NBUGFIX ,IRESP,YCOMMENT)
 CALL WRITE_SURF(HPROGRAM,'STORAGETYPE',HWRITE,IRESP,YCOMMENT)
 CALL WRITE_SURF(HPROGRAM,'DIM_FULL  ',U%NDIM_FULL,IRESP,HCOMMENT=YCOMMENT)
!
YCOMMENT='s'
 CALL WRITE_SURF(HPROGRAM,'DTCUR',U%TTIME,IRESP,YCOMMENT)
!
DGU%LSELECT=LSAVE_SELECT
!
 CALL WRITE_GRID(HPROGRAM,UG%CGRID,UG%XGRID_PAR,UG%XLAT,UG%XLON,UG%XMESH_SIZE,IRESP)
!
 CALL WRITESURF_ATM_CONF_n(HPROGRAM)
!
IF (HWRITE/='PRE') CALL WRITESURF_SSO_CANOPY_n(SSCP, &
                                                         HPROGRAM,(USS%CROUGH=='BE04' .AND. .NOT. LNOWRITE_CANOPY))
!
 CALL WRITESURF_CPL_GCM_n(U, &
                          HPROGRAM)
!
YCOMMENT='flag for accumulated variables'
 CALL WRITE_SURF(HPROGRAM,'BUDC',DGU%LSURF_BUDGETC,IRESP,HCOMMENT=YCOMMENT)
!
IF (DGU%LSURF_BUDGETC) THEN
   YCOMMENT='time of beginning of accumulation'
   CALL WRITE_SURF(HPROGRAM,'TBUDC',DGU%TIME_BUDGETC,IRESP,HCOMMENT=YCOMMENT)   
END IF
!  
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
!
!*       3.     Sea
!               ---
!
IF (U%NDIM_SEA>0) CALL WRITE_SEA_n(HPROGRAM,HWRITE)
!
!
!*       4.     Inland water
!               ------------
!
IF (U%NDIM_WATER>0) CALL WRITE_INLAND_WATER_n(HPROGRAM,HWRITE)
!
!
!*       5.     Vegetation scheme
!               -----------------
!
IF (U%NDIM_NATURE>0) CALL WRITE_NATURE_n(HPROGRAM,HWRITE,OLAND_USE)
!
!
!*       6.     Urban scheme
!               ------------
!
IF (U%NDIM_TOWN>0) CALL WRITE_TOWN_n(HPROGRAM,HWRITE)
IF (LHOOK) CALL DR_HOOK('WRITE_SURF_ATM_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_SURF_ATM_n
