!     #########
      SUBROUTINE WRITESURF_TEB_n(HPROGRAM,KPATCH,HWRITE)
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
!!      V. Masson   *Meteo France*
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
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_BEM_n, ONLY : B => BEM
USE MODD_TEB_PANEL_n, ONLY : TPN => TEB_PANEL
!
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
!
USE MODI_WRITE_SURF
USE MODI_WRITESURF_GR_SNOW
USE MODI_WRITESURF_TEB_GARDEN_n
USE MODI_WRITESURF_TEB_GREENROOF_n
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!RJ #ifdef SFX_MPI
!RJ INCLUDE "mpif.h"
!RJ #endif
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
INTEGER,           INTENT(IN)  :: KPATCH   ! current TEB patch
 CHARACTER(LEN=3),    INTENT(IN)  :: HWRITE    ! 'PREP' : does not write SBL XUNDEF fields
!                                             ! 'ALL' : all fields are written
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP           ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
 CHARACTER(LEN=3)  :: YPATCH         ! Patch identificator
 CHARACTER(LEN=7)  :: YDIR           ! Direction identificator
 CHARACTER(LEN=100):: YSTRING        ! Comment string
!
INTEGER :: JLAYER ! loop on surface layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_TEB_N',0,ZHOOK_HANDLE)
!
YPATCH='   '
IF (TOP%NTEB_PATCH>1) WRITE(YPATCH,FMT='(A,I1,A)') 'T',KPATCH,'_'
!
!
!*       2.     Option for road orientation:
!               ---------------------------
!
YCOMMENT='Option for Road orientation in TEB scheme'
 CALL WRITE_SURF(HPROGRAM,'ROAD_DIR',TOP%CROAD_DIR,IRESP,YCOMMENT)
YCOMMENT='Option for Wall representation in TEB scheme'
 CALL WRITE_SURF(HPROGRAM,'WALL_OPT',TOP%CWALL_OPT,IRESP,YCOMMENT)
!
!*       3.     Prognostic fields:
!               -----------------
!
!* roof temperatures
!

DO JLAYER=1,TOP%NROOF_LAYER
  WRITE(YRECFM,'(A3,A5,I1.1,A1)') YPATCH,'TROOF',JLAYER,' '
  WRITE(YCOMMENT,'(A9,I1.1,A4)') 'X_Y_TROOF',JLAYER,' (K)'
  YRECFM=ADJUSTL(YRECFM)
 CALL WRITE_SURF(HPROGRAM,YRECFM,T%XT_ROOF(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO

!
!* roof water content
!

YRECFM=YPATCH//'WS_ROOF'
YRECFM=ADJUSTL(YRECFM)
YCOMMENT='WS_ROOF (kg/m2)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,T%XWS_ROOF(:),IRESP,HCOMMENT=YCOMMENT)
!
!* road temperatures
!

DO JLAYER=1,TOP%NROAD_LAYER
  WRITE(YRECFM,'(A3,A5,I1.1,A1)') YPATCH,'TROAD',JLAYER,' '
  YRECFM=ADJUSTL(YRECFM)
  IF (TOP%CROAD_DIR=='UNIF' .OR. DTT%LDATA_ROAD_DIR) THEN
    YSTRING = 'X_Y_TROAD'
  ELSEIF (SIZE(T%XROAD_DIR)>0) THEN
    !* road direction is uniform spatially, one can then indicate it in the comment
    CALL ROAD_DIR(T%XROAD_DIR(1),YDIR)
    YSTRING=TRIM(YDIR)//' ROAD TEMP. LAYER '
  ELSE
    YSTRING='? ROAD TEMP. LAYER '
  ENDIF
  WRITE(YCOMMENT,'(A,I1.1,A4)') TRIM(YSTRING), JLAYER,' (K)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,T%XT_ROAD(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!* road water content
!

YRECFM=YPATCH//'WS_ROAD'
YRECFM=ADJUSTL(YRECFM)
YCOMMENT='WS_ROAD (kg/m2)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,T%XWS_ROAD(:),IRESP,HCOMMENT=YCOMMENT)
!
!* wall temperatures
!

DO JLAYER=1,TOP%NWALL_LAYER
 IF (TOP%CWALL_OPT=='UNIF') THEN
  WRITE(YRECFM,'(A3,A5,I1.1,A1)') YPATCH,'TWALL',JLAYER,' '
  YRECFM=ADJUSTL(YRECFM)
  WRITE(YCOMMENT,'(A9,I1.1,A4)') 'X_Y_TWALL',JLAYER,' (K)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,T%XT_WALL_A(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
 ELSE
  !* Wall A
  WRITE(YRECFM,'(A3,A6,I1.1)') YPATCH,'TWALLA',JLAYER
  YRECFM=ADJUSTL(YRECFM)
  IF (DTT%LDATA_ROAD_DIR) THEN
    YSTRING = 'X_Y_TWALL_A'
  ELSEIF (SIZE(T%XROAD_DIR)>0) THEN
    !* wall direction is uniform spatially, one can then indicate it in the comment
    CALL WALLA_DIR(T%XROAD_DIR(1),YDIR)
    YSTRING=TRIM(YDIR)//'-FACING WALL TEMP. LAYER '
  ELSE
    YSTRING='?-FACING WALL TEMP. LAYER '
  ENDIF
  WRITE(YCOMMENT,'(A,I1.1,A4)') TRIM(YSTRING), JLAYER,' (K)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,T%XT_WALL_A(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
  !
  !* Wall B
  WRITE(YRECFM,'(A3,A6,I1.1)') YPATCH,'TWALLB',JLAYER
  YRECFM=ADJUSTL(YRECFM)
  IF (DTT%LDATA_ROAD_DIR) THEN
    YSTRING = 'X_Y_TWALL_B'
  ELSEIF (SIZE(T%XROAD_DIR)>0) THEN
    !* wall direction is uniform spatially, one can then indicate it in the comment
    CALL WALLB_DIR(T%XROAD_DIR(1),YDIR)
    YSTRING=TRIM(YDIR)//'-FACING WALL TEMP. LAYER '
  ELSE
    YSTRING='?-FACING WALL TEMP. LAYER '
  ENDIF
  WRITE(YCOMMENT,'(A,I1.1,A4)') TRIM(YSTRING), JLAYER,' (K)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,T%XT_WALL_B(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
 END IF
END DO
!
!* internal building temperature
!
YRECFM=YPATCH//'TI_BLD'
YRECFM=ADJUSTL(YRECFM)
YCOMMENT='TI_BLD (K)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,B%XTI_BLD(:),IRESP,HCOMMENT=YCOMMENT)
!
!
!* outdoor window temperature
!
YRECFM=YPATCH//'T_WIN1'
YRECFM=ADJUSTL(YRECFM)
YCOMMENT='T_WIN1 (K)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,B%XT_WIN1(:),IRESP,HCOMMENT=YCOMMENT)
!
IF (TOP%CBEM=='BEM') THEN
!* internal building specific humidity
!
YRECFM=YPATCH//'QI_BLD'
YRECFM=ADJUSTL(YRECFM)
YCOMMENT='QI_BLD (kg/kg)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,B%XQI_BLD(:),IRESP,HCOMMENT=YCOMMENT)
!
  !
  !* indoor window temperature
  !
  YRECFM=YPATCH//'T_WIN2'
  YRECFM=ADJUSTL(YRECFM)
  YCOMMENT='T_WIN2 (K)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,B%XT_WIN2(:),IRESP,HCOMMENT=YCOMMENT)
  !
  !* floor temperatures
  !
  DO JLAYER=1,BOP%NFLOOR_LAYER
    WRITE(YRECFM,'(A3,A5,I1.1,A1)') YPATCH,'TFLOO',JLAYER,' '
    WRITE(YCOMMENT,'(A9,I1.1,A4)') 'X_Y_TFLOO',JLAYER,' (K)'
    YRECFM=ADJUSTL(YRECFM)
    CALL WRITE_SURF(HPROGRAM,YRECFM,B%XT_FLOOR(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
  END DO
  !
  !* internal th. mass temperature
  !
  DO JLAYER=1,BOP%NFLOOR_LAYER
    WRITE(YRECFM,'(A3,A5,I1.1,A1)') YPATCH,'TMASS',JLAYER,' '
    WRITE(YCOMMENT,'(A9,I1.1,A4)') 'X_Y_TMASS',JLAYER,' (K)'
    YRECFM=ADJUSTL(YRECFM)
    CALL WRITE_SURF(HPROGRAM,YRECFM,B%XT_MASS(:,JLAYER),IRESP,HCOMMENT=YCOMMENT)
  END DO        
  !
ENDIF
!
!* deep road temperature
!
YRECFM=YPATCH//'TI_ROAD'
YRECFM=ADJUSTL(YRECFM)
YCOMMENT='TI_ROAD (K)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,T%XTI_ROAD(:),IRESP,HCOMMENT=YCOMMENT)
!
!* snow mantel
!
YRECFM='RF'
 CALL WRITESURF_GR_SNOW(HPROGRAM,YRECFM,YPATCH,T%TSNOW_ROOF  )
!
YRECFM='RD'
 CALL WRITESURF_GR_SNOW(HPROGRAM,YRECFM,YPATCH,T%TSNOW_ROAD  )
!
!-------------------------------------------------------------------------------
!
!*       4.     Semi-prognostic fields:
!               ----------------------
!
!* temperature of canyon air
!
YRECFM=YPATCH//'TCANYON'
YRECFM=ADJUSTL(YRECFM)
YCOMMENT='T_CANYON (K)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,T%XT_CANYON(:),IRESP,HCOMMENT=YCOMMENT)
!
!* humidity of canyon air
!
YRECFM=YPATCH//'QCANYON'
YRECFM=ADJUSTL(YRECFM)
YCOMMENT='Q_CANYON (kg/kg)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,T%XQ_CANYON(:),IRESP,HCOMMENT=YCOMMENT)
!
!
!* Thermal solar panels present day production
!
IF (TOP%LSOLAR_PANEL) THEN
  YRECFM=YPATCH//'THER_PDAY'
  YRECFM=ADJUSTL(YRECFM)
  YCOMMENT='Thermal Solar Panels present day production (J/m2)'
  IF (.NOT. ASSOCIATED(TPN%XTHER_PRODC_DAY)) THEN
    ! for PREP cases
    ALLOCATE(TPN%XTHER_PRODC_DAY(SIZE(B%XTI_BLD)))
    TPN%XTHER_PRODC_DAY=0.
  END IF
  CALL WRITE_SURF(HPROGRAM,YRECFM,TPN%XTHER_PRODC_DAY(:),IRESP,HCOMMENT=YCOMMENT)
END IF
!-------------------------------------------------------------------------------
!
!*       5.  Time
!            ----
!
IF (KPATCH==1) THEN
  YRECFM='DTCUR'
  YCOMMENT='s'
  CALL WRITE_SURF(HPROGRAM,YRECFM,TOP%TTIME,IRESP,HCOMMENT=YCOMMENT)
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       6.  Urban green areas
!            ------------------
!
! Gardens
IF (TOP%LGARDEN) CALL WRITESURF_TEB_GARDEN_n(HPROGRAM,YPATCH)
!
! Grenn roofs
IF (TOP%LGREENROOF) CALL WRITESURF_TEB_GREENROOF_n(HPROGRAM,YPATCH)
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_TEB_N',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
CONTAINS
SUBROUTINE ROAD_DIR(PDIR,HDIR)
REAL,             INTENT(IN)  :: PDIR
 CHARACTER(LEN=7), INTENT(OUT) :: HDIR
REAL :: ZDIR
ZDIR=PDIR
IF (PDIR<0) ZDIR = PDIR +360.
IF (ZDIR>=  0.   .AND. ZDIR< 11.25) HDIR='N-S    '
IF (ZDIR>= 11.25 .AND. ZDIR< 33.75) HDIR='NNE-SSW'
IF (ZDIR>= 33.75 .AND. ZDIR< 56.25) HDIR='NE-SW'
IF (ZDIR>= 56.25 .AND. ZDIR< 78.75) HDIR='ENE-WSW'
IF (ZDIR>= 78.75 .AND. ZDIR<101.25) HDIR='E-W    '
IF (ZDIR>=101.25 .AND. ZDIR<123.75) HDIR='ESE-WNW'
IF (ZDIR>=123.75 .AND. ZDIR<146.25) HDIR='SE-NW  '
IF (ZDIR>=146.25 .AND. ZDIR<168.75) HDIR='SSE-NNW'
IF (ZDIR>=168.75 .AND. ZDIR<180.00) HDIR='N-S    '
END SUBROUTINE ROAD_DIR
SUBROUTINE WALLA_DIR(PDIR,HDIR)
REAL,             INTENT(IN)  :: PDIR
 CHARACTER(LEN=7), INTENT(OUT) :: HDIR
REAL :: ZDIR
ZDIR=PDIR
IF (PDIR<0) ZDIR = PDIR +360.
IF (ZDIR>=  0.   .AND. ZDIR< 11.25) HDIR='E      '
IF (ZDIR>= 11.25 .AND. ZDIR< 33.75) HDIR='ESE    '
IF (ZDIR>= 33.75 .AND. ZDIR< 56.25) HDIR='SE     ' 
IF (ZDIR>= 56.25 .AND. ZDIR< 78.75) HDIR='SSE    '
IF (ZDIR>= 78.75 .AND. ZDIR<101.25) HDIR='S      '
IF (ZDIR>=101.25 .AND. ZDIR<123.75) HDIR='SSW    '
IF (ZDIR>=123.75 .AND. ZDIR<146.25) HDIR='SW     '
IF (ZDIR>=146.25 .AND. ZDIR<168.75) HDIR='WSW    '
IF (ZDIR>=168.75 .AND. ZDIR<180.00) HDIR='W      '
END SUBROUTINE WALLA_DIR
SUBROUTINE WALLB_DIR(PDIR,HDIR)
REAL,             INTENT(IN)  :: PDIR
 CHARACTER(LEN=7), INTENT(OUT) :: HDIR
REAL :: ZDIR
ZDIR=PDIR
IF (PDIR<0) ZDIR = PDIR +360.
IF (ZDIR>=  0.   .AND. ZDIR< 11.25) HDIR='W      '
IF (ZDIR>= 11.25 .AND. ZDIR< 33.75) HDIR='WNW    '
IF (ZDIR>= 33.75 .AND. ZDIR< 56.25) HDIR='NW     ' 
IF (ZDIR>= 56.25 .AND. ZDIR< 78.75) HDIR='NNW    '
IF (ZDIR>= 78.75 .AND. ZDIR<101.25) HDIR='N      '
IF (ZDIR>=101.25 .AND. ZDIR<123.75) HDIR='NNE    '
IF (ZDIR>=123.75 .AND. ZDIR<146.25) HDIR='NE     '
IF (ZDIR>=146.25 .AND. ZDIR<168.75) HDIR='ENE    '
IF (ZDIR>=168.75 .AND. ZDIR<180.00) HDIR='E      '
END SUBROUTINE WALLB_DIR
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_TEB_n
