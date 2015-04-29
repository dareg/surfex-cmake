!     #########
      SUBROUTINE WRITE_DIAG_CH_AGGR_n (CHE, &
                                       HPROGRAM)
!     #################################
!
!!****  *WRITE_DIAG_CH_AGGR_n* - writes surface chemical emissions diagnostics
!!
!!    PURPOSE
!!    -------
!!
!!
!!**  METHOD
!!    ------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson & S. Queguiner  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2012
!!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_CH_EMIS_FIELD_n, ONLY : CH_EMIS_FIELD_t
!
USE MODD_CSTS,        ONLY : XAVOGADRO
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
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
TYPE(CH_EMIS_FIELD_t), INTENT(INOUT) :: CHE
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!

INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
!
INTEGER           :: JSPEC
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_CH_AGGR_n',0,ZHOOK_HANDLE)
 CALL INIT_IO_SURF_n(HPROGRAM,'FULL  ','SURF  ','WRITE')
!
!-------------------------------------------------------------------------------
!
!         Writes Emissions of all species
!
DO JSPEC=1,SIZE(CHE%TSEMISS)
  YRECFM = "E_"//TRIM(CHE%TSEMISS(JSPEC)%CNAME)
  YCOMMENT = "Emission data at time t (ppm*m/s)"
  CALL WRITE_SURF(HPROGRAM,YRECFM,CHE%TSEMISS(JSPEC)%XEMISDATA,IRESP,HCOMMENT=YCOMMENT)
END DO
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_CH_AGGR_n',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE WRITE_DIAG_CH_AGGR_n
