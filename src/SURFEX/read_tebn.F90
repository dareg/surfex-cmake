!     #########
      SUBROUTINE READ_TEB_n(HPROGRAM)
!     #########################################
!
!!****  *READ_TEB_n* - reads TEB fields
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
!
USE MODD_TEB_n,          ONLY : NROOF_LAYER, XT_ROOF, XWS_ROOF, &
                                  NROAD_LAYER, XT_ROAD, XWS_ROAD, &
                                  NWALL_LAYER, XT_WALL,           &
                                  XTI_BLD, XTI_ROAD,              &
                                  TSNOW_ROOF, TSNOW_ROAD,         &
                                  XT_CANYON, XQ_CANYON  
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
!
!*       0.2   Declarations of local variables
!     -------------------------------
!
LOGICAL           :: GTOWN          ! town variables written in the file
INTEGER           :: ILU            ! 1D physical dimension
INTEGER           :: IRESP          ! Error code after redding
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
INTEGER :: JLAYER  ! loop counter on layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_TEB_N',0,ZHOOK_HANDLE)
YRECFM='SIZE_TOWN'
CALL GET_TYPE_DIM_n('TOWN  ',ILU)
!
!*       2.     Prognostic fields:
!               -----------------
!
!* roof temperatures
!
ALLOCATE(XT_ROOF(ILU,NROOF_LAYER))
!
DO JLAYER=1,NROOF_LAYER
  WRITE(YRECFM,'(A6,I1.1,A9)') 'T_ROOF',JLAYER,'         '
CALL READ_SURF(HPROGRAM,YRECFM,XT_ROOF(:,JLAYER),IRESP)
END DO
!
!* roof water content
!
ALLOCATE(XWS_ROOF(ILU))
!
YRECFM='WS_ROOF'
CALL READ_SURF(HPROGRAM,YRECFM,XWS_ROOF(:),IRESP)
!
!* road temperatures
!
ALLOCATE(XT_ROAD(ILU,NROAD_LAYER))
!
DO JLAYER=1,NROAD_LAYER
  WRITE(YRECFM,'(A6,I1.1,A9)') 'T_ROAD',JLAYER,'         '
CALL READ_SURF(HPROGRAM,YRECFM,XT_ROAD(:,JLAYER),IRESP)
END DO
!
!* road water content
!
ALLOCATE(XWS_ROAD(ILU))
!
YRECFM='WS_ROAD'
CALL READ_SURF(HPROGRAM,YRECFM,XWS_ROAD(:),IRESP)
!
!* wall temperatures
!
ALLOCATE(XT_WALL(ILU,NWALL_LAYER))
!
DO JLAYER=1,NWALL_LAYER
  WRITE(YRECFM,'(A6,I1.1,A9)') 'T_WALL',JLAYER,'         '
CALL READ_SURF(HPROGRAM,YRECFM,XT_WALL(:,JLAYER),IRESP)
END DO
!
!* internal building temperature
!
ALLOCATE(XTI_BLD(ILU))
!
YRECFM='TI_BLD'
CALL READ_SURF(HPROGRAM,YRECFM,XTI_BLD(:),IRESP)
!
!* deep road temperature
!
ALLOCATE(XTI_ROAD(ILU))
!
YRECFM='TI_ROAD'
CALL READ_SURF(HPROGRAM,YRECFM,XTI_ROAD(:),IRESP)
!
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
  TSNOW_ROAD%SCHEME='1-L'
  CALL ALLOCATE_GR_SNOW(TSNOW_ROAD,ILU,1)
  TSNOW_ROOF%SCHEME='1-L'
  CALL ALLOCATE_GR_SNOW(TSNOW_ROOF,ILU,1)  
ELSE
  CALL READ_GR_SNOW(HPROGRAM,'ROAD',ILU,1,TSNOW_ROAD  )
  CALL READ_GR_SNOW(HPROGRAM,'ROOF',ILU,1,TSNOW_ROOF  )
END IF
!
!-------------------------------------------------------------------------------
!
!*       3.     Semi-prognostic fields:
!               ----------------------
!
!* temperature in canyon air
!
ALLOCATE(XT_CANYON(ILU))
XT_CANYON(:) = XT_ROAD(:,1)
!
YRECFM='T_CANYON'
CALL READ_SURF(HPROGRAM,YRECFM,XT_CANYON(:),IRESP)
!
!* water vapor in canyon air
!
ALLOCATE(XQ_CANYON(ILU))
XQ_CANYON(:) = 0.
!
YRECFM='Q_CANYON'
CALL READ_SURF(HPROGRAM,YRECFM,XQ_CANYON(:),IRESP)
IF (LHOOK) CALL DR_HOOK('READ_TEB_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_TEB_n
