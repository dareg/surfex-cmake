!     #########
      SUBROUTINE READ_GR_SNOW(HPROGRAM,HSURFTYPE,HPREFIX,     &
                              KLU,KPATCH,TPSNOW,HDIR,KVERSION,KBUGFIX)  
!     ##########################################################
!
!!****  *READ_GR_SNOW* - routine to read snow surface fields
!!
!!    PURPOSE
!!    -------
!       Initialize snow surface fields.
!
!!**  METHOD
!!    ------
!!    
!!    
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
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/01/99
!       F.solmon       06/00 adaptation for patch
!       V.Masson       01/03 new version of ISBA
!       B. Decharme    2008  If no WSNOW, WSNOW = XUNDEF
!-----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_TYPE_SNOW
!
USE MODI_READ_SURF
!
USE MODI_ALLOCATE_GR_SNOW
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_PREP_SNOW, ONLY : LSNOW_FRAC_TOT
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
 CHARACTER(LEN=6),   INTENT(IN)           :: HPROGRAM  ! calling program
 CHARACTER (LEN=*),  INTENT(IN)           :: HSURFTYPE ! generic name used for
                                                      ! snow characteristics
                                                      ! storage in file
 CHARACTER (LEN=3),  INTENT(IN)           :: HPREFIX   ! generic name for patch
!                                                     ! identification                      
INTEGER,            INTENT(IN)           :: KLU       ! horizontal size of snow var.
INTEGER,            INTENT(IN)           :: KPATCH    ! number of tiles
TYPE(SURF_SNOW)                          :: TPSNOW    ! snow characteristics
 CHARACTER (LEN=1),  INTENT(IN), OPTIONAL :: HDIR      ! type of reading
!                                                     ! HDIR = 'A' : entire field on All processors
!                                                     ! HDIR = 'H' : distribution on each processor
!
INTEGER,            INTENT(IN), OPTIONAL :: KVERSION
INTEGER,            INTENT(IN), OPTIONAL :: KBUGFIX
!
!*       0.2   declarations of local variables
!
 CHARACTER (LEN=7) :: YFMT0               ! format for writing
 CHARACTER (LEN=100) :: YFMT                ! format for writing
 CHARACTER(LEN=16)   :: YRECFM2 
 CHARACTER(LEN=12)   :: YRECFM              ! Name of the article to be read
 CHARACTER(LEN=4)    :: YNLAYER     !Format depending on the number of layers
 CHARACTER(LEN=1)    :: YDIR                ! type of reading
!
INTEGER             :: IRESP, JI, JP              ! Error code after redding
INTEGER             :: ISURFTYPE_LEN       ! 
INTEGER             :: JLAYER              ! loop counter
INTEGER             :: IVERSION, IBUGFIX
!
LOGICAL :: GVERSION
LOGICAL             :: GSNOW               ! snow written in the file
REAL(KIND=JPRB) :: ZHOOK_HANDLE, ZHOOK_HANDLE_OMP
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('READ_GR_SNOW_1',0,ZHOOK_HANDLE)
!
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
!-------------------------------------------------------------------------------
IF(PRESENT(KVERSION))THEN
  IVERSION=KVERSION
ELSE
  CALL READ_SURF(&
                 HPROGRAM,'VERSION',IVERSION,IRESP)
ENDIF
IF(PRESENT(KBUGFIX))THEN
  IBUGFIX=KBUGFIX
ELSE
  CALL READ_SURF(&
                 HPROGRAM,'BUG',IBUGFIX,IRESP)
ENDIF
!-------------------------------------------------------------------------------
!
GVERSION = (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3)
!
!*       1.    Type of snow scheme
!              -------------------
!
ISURFTYPE_LEN=LEN_TRIM(HSURFTYPE)
!
IF (IVERSION <=2 .OR. (IVERSION==3 .AND. IBUGFIX<=4)) THEN
  WRITE(YFMT,'(A5,I1,A4)')     '(A5,A',ISURFTYPE_LEN,',A5)'
  WRITE(YRECFM2,YFMT) 'SNOW_',HSURFTYPE,'_TYPE'
ELSE
  IF (IVERSION<7 .OR. IVERSION==7 .AND. IBUGFIX<3) THEN
    WRITE(YFMT,'(A5,I1,A4)')     '(A3,A',ISURFTYPE_LEN,',A5)'
    WRITE(YRECFM2,YFMT) 'SN_',HSURFTYPE,'_TYPE'
  ELSE
    WRITE(YFMT,'(A5,I1,A4)')     '(A3,A',ISURFTYPE_LEN,',A4)'
    WRITE(YRECFM2,YFMT) 'SN_',HSURFTYPE,'_TYP'
    YRECFM2=ADJUSTL(HPREFIX//YRECFM2)
  ENDIF
END IF
!
 CALL READ_SURF(HPROGRAM,YRECFM2,TPSNOW%SCHEME,IRESP)
!
!*       2.    Snow levels
!              -----------
!
!
IF (IVERSION <=2 .OR. (IVERSION==3 .AND. IBUGFIX<=4)) THEN
  WRITE(YFMT,'(A5,I1,A4)')     '(A5,A',ISURFTYPE_LEN,',A6)'
  WRITE(YRECFM2,YFMT) 'SNOW_',HSURFTYPE,'_LAYER'
ELSE
  WRITE(YFMT,'(A5,I1,A4)')     '(A3,A',ISURFTYPE_LEN,',A2)'
  WRITE(YRECFM2,YFMT) 'SN_',HSURFTYPE,'_N'
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) YRECFM2=ADJUSTL(HPREFIX//YRECFM2)
END IF
!
 CALL READ_SURF(HPROGRAM,YRECFM2,TPSNOW%NLAYER,IRESP)
!
!*       2.    Presence of snow fields in the file
!              -----------------------------------
!
IF (IVERSION >6 .OR. (IVERSION==6 .AND. IBUGFIX>=1)) THEN
  WRITE(YFMT,'(A5,I1,A1)')     '(A3,A',ISURFTYPE_LEN,')'
  WRITE(YRECFM,YFMT) 'SN_',HSURFTYPE
  IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) YRECFM=ADJUSTL(HPREFIX//YRECFM)
  CALL READ_SURF(HPROGRAM,YRECFM,GSNOW,IRESP)
ELSE
  IF (TPSNOW%NLAYER==0) THEN
    GSNOW = .FALSE.
    IF (TPSNOW%SCHEME=='D95' .OR. TPSNOW%SCHEME=='1-L' .OR. TPSNOW%SCHEME=='EBA') TPSNOW%NLAYER=1
    IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO'                          ) TPSNOW%NLAYER=3
  ELSE
    GSNOW = .TRUE.
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
!*       3.    Allocations
!              -----------
!
 CALL ALLOCATE_GR_SNOW(TPSNOW,KLU,KPATCH)
!
IF (.NOT. GSNOW) THEN
  IF (LHOOK) CALL DR_HOOK('READ_GR_SNOW_1',1,ZHOOK_HANDLE)
  RETURN
END IF
!-------------------------------------------------------------------------------
!
!*       4.    Additional key
!              ---------------
!
IF (IVERSION >= 7 .AND. HSURFTYPE=='VEG') CALL READ_SURF(HPROGRAM,'LSNOW_FRAC_T',LSNOW_FRAC_TOT,IRESP)
!
!-------------------------------------------------------------------------------
!
!*       5.    Snow reservoir
!              --------------
!
!
IF (TPSNOW%SCHEME=='1-L' .OR. TPSNOW%SCHEME=='D95' .OR. TPSNOW%SCHEME=='EBA' &
                         .OR. TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN 
  !
  WRITE(YFMT0,'(A5,I1,A1)') ',A1,A',ISURFTYPE_LEN
  !
  IF (GVERSION) THEN
    YFMT = '(A3'//YFMT0
  ELSE
    YFMT = '(A5'//YFMT0
  ENDIF
  CALL READ_LAYERS(GVERSION,TPSNOW%NLAYER,YDIR,HPREFIX,YFMT,"WSNOW",HSURFTYPE,TPSNOW%WSNOW)
  CALL READ_LAYERS(GVERSION,TPSNOW%NLAYER,YDIR,HPREFIX,YFMT,"RSNOW",HSURFTYPE,TPSNOW%RHO)
  !
  !*       7.    Snow temperature
  !              ----------------
  !
  IF (TPSNOW%SCHEME=='1-L') THEN
    !
    CALL READ_LAYERS(GVERSION,TPSNOW%NLAYER,YDIR,HPREFIX,YFMT,"TSNOW",HSURFTYPE,TPSNOW%T)
    !
  ENDIF
  !
  !*       8.    Heat content
  !              ------------
  !
  IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN
    !
    CALL READ_LAYERS(GVERSION,TPSNOW%NLAYER,YDIR,HPREFIX,YFMT,"HSNOW",HSURFTYPE,TPSNOW%HEAT)
    !
    IF (TPSNOW%SCHEME=='CRO') THEN
      !
      CALL READ_LAYERS(GVERSION,TPSNOW%NLAYER,YDIR,HPREFIX,YFMT,"SHIST",HSURFTYPE,TPSNOW%HIST)
      !
      !*       9.    Snow Gran1
      !              ------------
      !
      IF (GVERSION) THEN
        YFMT = "(A2"//YFMT0         
      ELSE
        YFMT = "(A5"//YFMT0
      ENDIF
      YFMT = YFMT//YNLAYER//')'
      CALL READ_LAYERS(GVERSION,TPSNOW%NLAYER,YDIR,HPREFIX,YFMT,"SGRAN",HSURFTYPE,TPSNOW%GRAN1,HREC2="1")
      CALL READ_LAYERS(GVERSION,TPSNOW%NLAYER,YDIR,HPREFIX,YFMT,"SGRAN",HSURFTYPE,TPSNOW%GRAN2,HREC2="2")
      !
      !*       12.    Age parameter
      !              -------------------
      !
      IF (GVERSION) THEN
        YFMT = "(A3"//YFMT0         
      ELSE
        YFMT = "(A4"//YFMT0
      ENDIF
      YFMT = YFMT//YNLAYER//')'
      CALL READ_LAYERS(GVERSION,TPSNOW%NLAYER,YDIR,HPREFIX,YFMT,"SAGE",HSURFTYPE,TPSNOW%AGE)
      !
    ELSEIF(TPSNOW%SCHEME=='3-L'.AND.IVERSION<8)THEN
      !
      DO JLAYER = 1,TPSNOW%NLAYER
        WHERE (TPSNOW%WSNOW(:,1,:) >= 0.0) 
          TPSNOW%AGE(:,JLAYER,:) = 0.0
        ELSEWHERE
          TPSNOW%AGE(:,JLAYER,:) = XUNDEF
        ENDWHERE
      ENDDO
      !
    END IF    
    !
  ENDIF
  !
  YFMT = TRIM(YFMT)//')'
  WRITE(YRECFM,YFMT) 'ASNOW','_',HSURFTYPE
  IF (GVERSION) YRECFM=ADJUSTL(HPREFIX//YRECFM) 
  CALL READ_SURF(HPROGRAM,YRECFM,TPSNOW%ALB(:,:),IRESP,HDIR=YDIR)
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('READ_GR_SNOW_1',1,ZHOOK_HANDLE)
!
IF (TPSNOW%SCHEME=='1-L' .OR. TPSNOW%SCHEME=='D95' .OR. TPSNOW%SCHEME=='EBA' &
                         .OR. TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN 
  !
!$OMP PARALLEL PRIVATE(ZHOOK_HANDLE_OMP)
IF (LHOOK) CALL DR_HOOK('READ_GR_SNOW_2',0,ZHOOK_HANDLE_OMP)
!$OMP DO PRIVATE(JI,JP,JLAYER)
  DO JI = 1,SIZE(TPSNOW%WSNOW,1)
    DO JP = 1,SIZE(TPSNOW%WSNOW,3)
      !
      IF (TPSNOW%WSNOW(JI,1,JP) == 0.0 ) THEN
        !
        TPSNOW%ALB(JI,JP) = XUNDEF
        !
        DO JLAYER = 1,TPSNOW%NLAYER
          !
          TPSNOW%RHO(JI,JLAYER,JP)=XUNDEF
          IF (TPSNOW%SCHEME=='1-L') THEN
            TPSNOW%T(JI,JLAYER,JP) = XUNDEF
          ELSEIF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN
            TPSNOW%HEAT(JI,JLAYER,JP) = XUNDEF
            IF (TPSNOW%SCHEME=='CRO') THEN
              TPSNOW%HIST (JI,JLAYER,JP) = XUNDEF
              TPSNOW%GRAN1(JI,JLAYER,JP) = XUNDEF
              TPSNOW%GRAN2(JI,JLAYER,JP) = XUNDEF
              TPSNOW%AGE  (JI,JLAYER,JP) = XUNDEF
            ENDIF
          ENDIF
          !
        ENDDO
      ENDIF
    ENDDO
  ENDDO
!$OMP ENDDO
IF (LHOOK) CALL DR_HOOK('READ_GR_SNOW_2',1,ZHOOK_HANDLE_OMP)
!$OMP END PARALLEL
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
SUBROUTINE READ_LAYERS(OVERSION,KNL,HDIRIN,HPREF,HFMT,HREC,HSURF,PTAB,HREC2)
!
IMPLICIT NONE
!
LOGICAL, INTENT(IN) :: OVERSION
INTEGER, INTENT(IN) :: KNL
 CHARACTER(LEN=*), INTENT(IN) :: HDIRIN
 CHARACTER(LEN=*), INTENT(IN) :: HPREF
 CHARACTER(LEN=*), INTENT(IN) :: HFMT
 CHARACTER(LEN=*), INTENT(IN) :: HREC
 CHARACTER(LEN=*), INTENT(IN) :: HSURF
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PTAB
 CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: HREC2
!
 CHARACTER(LEN=1) :: YREC2
 CHARACTER (LEN=100) :: YFMT     ! format for writing
 CHARACTER(LEN=12)   :: YRECFM   ! Name of the article to be read
 CHARACTER(LEN=4) :: YNLAYER     !Format depending on the number of layers
INTEGER :: JLAYER, IRESP
!
IF (PRESENT(HREC2)) THEN
  YREC2=TRIM(HREC2)
ELSE
  YREC2=""
ENDIF
!
YNLAYER='I1.1'
DO JLAYER = 1,KNL
  !
  IF (JLAYER==10) YNLAYER='I2.2'
  YFMT = TRIM(HFMT)//','//YNLAYER//')'
  !
  WRITE(YRECFM,YFMT) TRIM(HREC),TRIM(YREC2)//'_',TRIM(HSURF),JLAYER
  IF (OVERSION) YRECFM=ADJUSTL(TRIM(HPREF)//YRECFM)
  !
  CALL READ_SURF(HPROGRAM,YRECFM,PTAB(:,JLAYER,:),IRESP,HDIR=HDIRIN)
  !
ENDDO
!
END SUBROUTINE READ_LAYERS
!
END SUBROUTINE READ_GR_SNOW
