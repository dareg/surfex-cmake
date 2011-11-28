!     #########
      SUBROUTINE WRITESURF_TEB_n(HPROGRAM)
!     ####################################
!
!!****  *WRITE_TEB_n* - writes TEB fields
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_TEB_n,          ONLY : LGARDEN,                        &
                                  NROOF_LAYER, XT_ROOF, XWS_ROOF, &
                                  NROAD_LAYER, XT_ROAD, XWS_ROAD, &
                                  NWALL_LAYER, XT_WALL,           &
                                  XTI_BLD, XTI_ROAD,              &
                                  TSNOW_ROOF, TSNOW_ROAD,         &
                                  XT_CANYON, XQ_CANYON,           &
                                  TTIME  
!
USE MODI_WRITE_SURF
USE MODI_WRITESURF_GR_SNOW
USE MODI_WRITESURF_TEB_GARDEN_n
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
!
INTEGER :: JLAYER ! loop on surface layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!*       3.     Prognostic fields:
!               -----------------
!
!* roof temperatures
!

IF (LHOOK) CALL DR_HOOK('WRITESURF_TEB_N',0,ZHOOK_HANDLE)
DO JLAYER=1,NROOF_LAYER
  WRITE(YRECFM,'(A6,I1.1,A8)') 'T_ROOF',JLAYER,'       '
  WRITE(YCOMMENT,'(A10,I1.1,A4)') 'X_Y_T_ROOF',JLAYER,' (K)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XT_ROOF(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO

!
!* roof water content
!

YRECFM='WS_ROOF'
YCOMMENT='WS_ROOF (kg/m2)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XWS_ROOF(:),IRESP,HCOMMENT=YCOMMENT)
!
!* road temperatures
!

DO JLAYER=1,NROAD_LAYER
  WRITE(YRECFM,'(A6,I1.1,A8)') 'T_ROAD',JLAYER,'       '
  WRITE(YCOMMENT,'(A10,I1.1,A4)') 'X_Y_T_ROAD',JLAYER,' (K)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XT_ROAD(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!* road water content
!

YRECFM='WS_ROAD'
YCOMMENT='WS_ROAD (kg/m2)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XWS_ROAD(:),IRESP,HCOMMENT=YCOMMENT)
!
!* wall temperatures
!

DO JLAYER=1,NWALL_LAYER
  WRITE(YRECFM,'(A6,I1.1,A8)') 'T_WALL',JLAYER,'       '
  WRITE(YCOMMENT,'(A10,I1.1,A4)') 'X_Y_T_WALL',JLAYER,' (K)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XT_WALL(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!* internal building temperature
!
YRECFM='TI_BLD'
YCOMMENT='TI_BLD (K)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XTI_BLD(:),IRESP,HCOMMENT=YCOMMENT)
!
!* deep road temperature
!
YRECFM='TI_ROAD'
YCOMMENT='TI_ROAD (K)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XTI_ROAD(:),IRESP,HCOMMENT=YCOMMENT)
!
!* snow mantel

!
CALL WRITESURF_GR_SNOW(HPROGRAM,'ROOF',1,TSNOW_ROOF  )
!
CALL WRITESURF_GR_SNOW(HPROGRAM,'ROAD',1,TSNOW_ROAD  )
!
!-------------------------------------------------------------------------------
!
!*       4.     Semi-prognostic fields:
!               ----------------------
!
!* temperature of canyon air
!
YRECFM='T_CANYON'
YCOMMENT='T_CANYON (K)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XT_CANYON(:),IRESP,HCOMMENT=YCOMMENT)
!
!* humidity of canyon air
!
YRECFM='Q_CANYON'
YCOMMENT='Q_CANYON (kg/kg)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XQ_CANYON(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!*       5.  Time
!            ----
!
YRECFM='DTCUR'
YCOMMENT='s'
CALL WRITE_SURF(HPROGRAM,YRECFM,TTIME,IRESP,HCOMMENT=YCOMMENT)
!
!
!-------------------------------------------------------------------------------
!
!*       6.  §Urban green areas
!            ------------------
!

IF (LGARDEN) CALL WRITESURF_TEB_GARDEN_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITESURF_TEB_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_TEB_n
