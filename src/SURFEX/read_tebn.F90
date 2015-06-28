!     #########
      SUBROUTINE READ_TEB_n(HPROGRAM,KPATCH)
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
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
!
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_BEM_n, ONLY : B => BEM
USE MODD_TEB_PANEL_n, ONLY : TPN => TEB_PANEL
USE MODD_ASSIM, ONLY : LASSIM,XAT2M_TEB,NIPERT,NVAR
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
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_TYPE_DIM_n
USE MODD_SURF_PAR, ONLY : XUNDEF
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! calling program
INTEGER,           INTENT(IN)  :: KPATCH   ! current patch number
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
LOGICAL           :: GTOWN          ! town variables written in the file
INTEGER           :: ILU          ! 1D physical dimension
!
INTEGER           :: IRESP          ! Error code after redding
!
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=3)  :: YPATCH         ! suffix if more than 1 patch
!
INTEGER           :: IVERSION, IBUGFIX
LOGICAL           :: GOLD_NAME      ! name of temperatures in old versions of SURFEX
!
INTEGER :: JLAYER  ! loop counter on layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_TEB_N',0,ZHOOK_HANDLE)
YRECFM='SIZE_TOWN'
 CALL GET_TYPE_DIM_n(DTCO, U, &
                     'TOWN  ',ILU)
!
YPATCH='   '
IF (TOP%NTEB_PATCH>1) WRITE(YPATCH,FMT='(A1,I1,A1)') 'T',KPATCH,'_'
!  
 CALL READ_SURF(IOB, &
                HPROGRAM,'VERSION',IVERSION,IRESP)
 CALL READ_SURF(IOB, &
                HPROGRAM,'BUG',IBUGFIX,IRESP)
GOLD_NAME = (IVERSION<7 .OR. (IVERSION==7 .AND. IBUGFIX<=2))
!
!*       2.     Prognostic fields:
!               -----------------
!
!* roof temperatures
!
ALLOCATE(T%CUR%XT_ROOF(ILU,TOP%NROOF_LAYER))
!
DO JLAYER=1,TOP%NROOF_LAYER
  WRITE(YRECFM,'(A3,A5,I1.1)') YPATCH,'TROOF',JLAYER
  YRECFM=ADJUSTL(YRECFM)
  IF (GOLD_NAME) WRITE(YRECFM,'(A6,I1.1)') 'T_ROOF',JLAYER

 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,T%CUR%XT_ROOF(:,JLAYER),IRESP)
END DO
!
!* roof water content
!
ALLOCATE(T%CUR%XWS_ROOF(ILU))
!
YRECFM=YPATCH//'WS_ROOF'
YRECFM=ADJUSTL(YRECFM)
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,T%CUR%XWS_ROOF(:),IRESP)
!
!* road temperatures
!
ALLOCATE(T%CUR%XT_ROAD(ILU,TOP%NROAD_LAYER))
!
DO JLAYER=1,TOP%NROAD_LAYER
  WRITE(YRECFM,'(A3,A5,I1.1)') YPATCH,'TROAD',JLAYER
  YRECFM=ADJUSTL(YRECFM)
  IF (GOLD_NAME) WRITE(YRECFM,'(A6,I1.1)') 'T_ROAD',JLAYER
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,T%CUR%XT_ROAD(:,JLAYER),IRESP)
END DO
!
!* road water content
!
ALLOCATE(T%CUR%XWS_ROAD(ILU))
!
YRECFM=YPATCH//'WS_ROAD'
YRECFM=ADJUSTL(YRECFM)
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,T%CUR%XWS_ROAD(:),IRESP)
!
!* wall temperatures
!
ALLOCATE(T%CUR%XT_WALL_A(ILU,TOP%NWALL_LAYER))
ALLOCATE(T%CUR%XT_WALL_B(ILU,TOP%NWALL_LAYER))
!
DO JLAYER=1,TOP%NWALL_LAYER
  IF (TOP%CWALL_OPT=='UNIF' .OR. GOLD_NAME) THEN
    WRITE(YRECFM,'(A3,A5,I1.1)') YPATCH,'TWALL',JLAYER
    YRECFM=ADJUSTL(YRECFM)
    IF (GOLD_NAME) WRITE(YRECFM,'(A6,I1.1)') 'T_WALL',JLAYER
    CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,T%CUR%XT_WALL_A(:,JLAYER),IRESP)
    !
    T%CUR%XT_WALL_B = T%CUR%XT_WALL_A
  ELSE
    WRITE(YRECFM,'(A3,A6,I1.1)') YPATCH,'TWALLA',JLAYER
    YRECFM=ADJUSTL(YRECFM)
    CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,T%CUR%XT_WALL_A(:,JLAYER),IRESP)
    !
    WRITE(YRECFM,'(A3,A6,I1.1)') YPATCH,'TWALLB',JLAYER
    YRECFM=ADJUSTL(YRECFM)
    CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,T%CUR%XT_WALL_B(:,JLAYER),IRESP)
  END IF
END DO
!
!* internal building temperature
!
ALLOCATE(B%CUR%XTI_BLD(ILU))
!
YRECFM=YPATCH//'TI_BLD'
YRECFM=ADJUSTL(YRECFM)
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,B%CUR%XTI_BLD(:),IRESP)

!
!* outdoor window temperature
!
ALLOCATE(B%CUR%XT_WIN1(ILU))
!
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
   YRECFM=YPATCH//'T_WIN1'
   YRECFM=ADJUSTL(YRECFM)
   CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,B%CUR%XT_WIN1(:),IRESP)
ELSE
   B%CUR%XT_WIN1(:)=XUNDEF
ENDIF
!
!
!* internal building specific humidity
!
ALLOCATE(B%CUR%XQI_BLD(ILU))
!
IF (TOP%CBEM=='BEM' .AND. (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3)) THEN
   YRECFM=YPATCH//'QI_BLD'
   YRECFM=ADJUSTL(YRECFM)
   CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,B%CUR%XQI_BLD(:),IRESP)
ELSE
   B%CUR%XQI_BLD(:) = XUNDEF
ENDIF
!
IF (TOP%CBEM=='BEM') THEN
  !
  !* indoor window temperature
  !
  ALLOCATE(B%CUR%XT_WIN2(ILU))
  !
  YRECFM=YPATCH//'T_WIN2'
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,B%CUR%XT_WIN2(:),IRESP)        
  !
  !* floor temperatures
  !
  ALLOCATE(B%CUR%XT_FLOOR(ILU,BOP%NFLOOR_LAYER))
  !
  DO JLAYER=1,BOP%NFLOOR_LAYER
    WRITE(YRECFM,'(A3,A5,I1.1)') YPATCH,'TFLOO',JLAYER
    YRECFM=ADJUSTL(YRECFM)
    CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,B%CUR%XT_FLOOR(:,JLAYER),IRESP)
  END DO
  !
  !* mass temperatures
  !
  ALLOCATE(B%CUR%XT_MASS(ILU,BOP%NFLOOR_LAYER))
  !
  DO JLAYER=1,BOP%NFLOOR_LAYER
    WRITE(YRECFM,'(A3,A5,I1.1)') YPATCH,'TMASS',JLAYER
    YRECFM=ADJUSTL(YRECFM)
    CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,B%CUR%XT_MASS(:,JLAYER),IRESP)
  END DO
  !
ELSE 
  ALLOCATE(B%CUR%XT_WIN2(0))
  ALLOCATE(B%CUR%XT_FLOOR(0,0))
  ALLOCATE(B%CUR%XT_MASS(0,0))
ENDIF
!
!* deep road temperature
!
ALLOCATE(T%CUR%XTI_ROAD(ILU))
!
YRECFM=YPATCH//'TI_ROAD'
YRECFM=ADJUSTL(YRECFM)
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,T%CUR%XTI_ROAD(:),IRESP)
!
!
!* snow mantel
!
 CALL END_IO_SURF_n(HPROGRAM)
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PGD ')
 CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
!
 CALL TOWN_PRESENCE(IOB, &
                    HPROGRAM,GTOWN)
!
 CALL END_IO_SURF_n(HPROGRAM)
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP')
 CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
!
IF (.NOT. GTOWN) THEN
  T%CUR%TSNOW_ROAD%SCHEME='1-L'
  CALL ALLOCATE_GR_SNOW(T%CUR%TSNOW_ROAD,ILU,1)
  T%CUR%TSNOW_ROOF%SCHEME='1-L'
  CALL ALLOCATE_GR_SNOW(T%CUR%TSNOW_ROOF,ILU,1)  
ELSE
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
    CALL READ_GR_SNOW(IOB, &
                      HPROGRAM,'RD',YPATCH,ILU,1,T%CUR%TSNOW_ROAD  )
    CALL READ_GR_SNOW(IOB, &
                      HPROGRAM,'RF',YPATCH,ILU,1,T%CUR%TSNOW_ROOF  )
  ELSE
    CALL READ_GR_SNOW(IOB, &
                      HPROGRAM,'ROAD',YPATCH,ILU,1,T%CUR%TSNOW_ROAD  )
    CALL READ_GR_SNOW(IOB, &
                      HPROGRAM,'ROOF',YPATCH,ILU,1,T%CUR%TSNOW_ROOF  )
  ENDIF    
END IF
!
!-------------------------------------------------------------------------------
!
!*       3.     Semi-prognostic fields:
!               ----------------------
!
!* temperature in canyon air
!
ALLOCATE(T%CUR%XT_CANYON(ILU))
T%CUR%XT_CANYON(:) = T%CUR%XT_ROAD(:,1)
!
YRECFM=YPATCH//'TCANYON'
YRECFM=ADJUSTL(YRECFM)
IF (GOLD_NAME) YRECFM='T_CANYON'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,T%CUR%XT_CANYON(:),IRESP)
!
!* water vapor in canyon air
!
ALLOCATE(T%CUR%XQ_CANYON(ILU))
T%CUR%XQ_CANYON(:) = 0.
!
YRECFM=YPATCH//'QCANYON'
YRECFM=ADJUSTL(YRECFM)
IF (GOLD_NAME) YRECFM='Q_CANYON'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,T%CUR%XQ_CANYON(:),IRESP)
!
!* Thermal solar panels present day production
!
IF (TOP%LSOLAR_PANEL) THEN
  ALLOCATE(TPN%XTHER_PRODC_DAY(ILU))
  TPN%XTHER_PRODC_DAY(:) = 0.

  YRECFM=YPATCH//'THER_PDAY'
  YRECFM=ADJUSTL(YRECFM)
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TPN%XTHER_PRODC_DAY(:),IRESP)
END IF

IF ( LASSIM .AND. NIPERT/=NVAR+2 ) THEN
  ! Diagnostic fields for assimilation
  IF ( .NOT. ALLOCATED(XAT2M_TEB)) ALLOCATE(XAT2M_TEB(ILU))
  XAT2M_TEB=XUNDEF
  YRECFM='T2M'
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,XAT2M_TEB(:),IRESP)
ENDIF

!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READ_TEB_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_TEB_n
