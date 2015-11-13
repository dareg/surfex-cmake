!     #########
SUBROUTINE PREP_SNOW_EXTERN (&
                             HPROGRAM,HSURF,HFILE,HFILETYPE,HFILEPGD,HFILEPGDTYPE,&
                            KLUOUT,PFIELD,OSNOW_IDEAL,KLAYER,KTEB_PATCH)
!     #################################################################################
!
!
!!****  *PREP_SNOW_EXTERN*  
!!
!!    PURPOSE
!!    -------
!       Read and prepare initial snow fields from external files
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    
!!    AUTHOR
!!    ------
!!         * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    ?
!!       02/2014 E. Martin : cor. for passing from from multilayer to a single layer
!!      B. Decharme  04/2014, external init with FA files
!!                            improve vertical interpolation
!-------------------------------------------------------------------------------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO
USE MODD_TYPE_SNOW
USE MODD_PREP,           ONLY : CINGRID_TYPE, CINTERP_TYPE
USE MODD_PREP_SNOW,      ONLY : XGRID_SNOW, NGRID_LEVEL
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_CSTS,           ONLY : XTT
!
USE MODE_SNOW3L
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_TOWN_PRESENCE
USE MODI_ABOR1_SFX
USE MODI_PREP_GRID_EXTERN
USE MODI_OPEN_AUX_IO_SURF
USE MODI_CLOSE_AUX_IO_SURF
USE MODI_ALLOCATE_GR_SNOW
USE MODI_DEALLOC_GR_SNOW
USE MODI_INTERP_GRID_NAT
USE MODI_READ_GR_SNOW
USE MODI_READ_SURF
USE MODI_SNOW_T_WLIQ_TO_HEAT
USE MODI_READ_TEB_PATCH
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=10),  INTENT(IN)  :: HSURF     ! type of field
CHARACTER(LEN=28),  INTENT(IN)  :: HFILE     ! name of file
CHARACTER(LEN=6),   INTENT(IN)  :: HFILETYPE ! type of file
CHARACTER(LEN=28),  INTENT(IN)  :: HFILEPGD     ! name of file
CHARACTER(LEN=6),   INTENT(IN)  :: HFILEPGDTYPE ! type of file
INTEGER,            INTENT(IN)  :: KLUOUT    ! logical unit of output listing
REAL,DIMENSION(:,:,:), POINTER  :: PFIELD    ! field to interpolate horizontally
LOGICAL,            INTENT(INOUT)  :: OSNOW_IDEAL
INTEGER,            INTENT(IN)  :: KLAYER    ! Number of layer of output snow scheme
INTEGER,            INTENT(IN) :: KTEB_PATCH
!
!*      0.2    declarations of local variables
!
TYPE(SURF_SNOW)                    :: TZSNOW ! snow characteristics

REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZFIELD       ! work field on input snow grid
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTEMP        ! snow temperature
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWLIQ        ! liquid water snow pack content
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZD           ! total snow depth
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDEPTH       ! thickness of each layer (m)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZGRID        ! normalized input grid
!
LOGICAL                           :: GTOWN          ! town variables written in the file
CHARACTER(LEN=12)                 :: YRECFM         ! record name
INTEGER                           :: IRESP          ! error return code
INTEGER                           :: IVERSION       ! SURFEX version
LOGICAL                           :: GOLD_NAME      ! old name flag 
INTEGER                           :: IBUGFIX        ! SURFEX bug version
INTEGER                           :: IVEGTYPE       ! actual number of vegtypes
INTEGER                           :: JLAYER         ! loop on snow vertical grids
INTEGER                           :: JI             ! loop on pts
INTEGER                           :: INI
CHARACTER(LEN=8)                  :: YAREA          ! area treated ('ROOF','ROAD','VEG ')
CHARACTER(LEN=3)                  :: YPREFIX        ! prefix to identify patch
INTEGER                           :: IPATCH         ! number of input patch
INTEGER                           :: ITEB_PATCH     ! number of input patch for TEB
INTEGER                           :: JPATCH         ! loop on patch
CHARACTER(LEN=6)                  :: YMASK          ! type of tile mask
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!*      3.     Area being treated
!              ------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_SNOW_EXTERN',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
YAREA='        '
YAREA(1:4) = HSURF(7:10)
!
IF (YAREA(1:4)=='VEG ') THEN
  IVEGTYPE = NVEGTYPE
  YMASK = 'NATURE'
  YPREFIX = '   '  
ELSE
  IVEGTYPE = 1
  YMASK    = 'TOWN  '
  IPATCH   = 1
  YPREFIX = '   '  
END IF
!
!*      1.     Preparation of IO for reading in the file
!              -----------------------------------------
!
!* Note that all points are read, even those without physical meaning.
!  These points will not be used during the horizontal interpolation step.
!  Their value must be defined as XUNDEF.
!
!* reading of version of the file being read
CALL OPEN_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE,'FULL  ')
CALL READ_SURF(HFILEPGDTYPE,'VERSION',IVERSION,IRESP,HDIR='-')
CALL READ_SURF(HFILEPGDTYPE,'BUG',IBUGFIX,IRESP,HDIR='-')
GOLD_NAME=(IVERSION<7 .OR. (IVERSION==7 .AND. IBUGFIX<3))
!
 CALL PREP_GRID_EXTERN(HFILEPGDTYPE,KLUOUT,CINGRID_TYPE,CINTERP_TYPE,INI)
!
 CALL TOWN_PRESENCE(HFILEPGDTYPE,GTOWN,HDIR='-')
!
CALL CLOSE_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE)
!
CALL OPEN_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE,YMASK)
!
IF (YAREA(1:4)=='VEG ') THEN
  YRECFM = 'PATCH_NUMBER'
  CALL READ_SURF(HFILEPGDTYPE,YRECFM,IPATCH,IRESP,HDIR='-')
ELSE
  IF (.NOT.GOLD_NAME) THEN
    IF (YAREA(1:4)=='ROOF') YAREA(1:4) = 'RF  '
    IF (YAREA(1:4)=='ROAD') YAREA(1:4) = 'RD  '
  ENDIF
  IF (GTOWN) THEN   
    CALL READ_TEB_PATCH(HFILEPGD,HFILEPGDTYPE,IVERSION,IBUGFIX,ITEB_PATCH,HDIR='-')
  ELSE
    ITEB_PATCH = 1
  ENDIF    
  IF (ITEB_PATCH>1) THEN
    WRITE(YPREFIX,FMT='(A,I1,A)') 'T',MIN(KTEB_PATCH,ITEB_PATCH),'_'
  END IF  
END IF
!
 CALL CLOSE_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE)
!-------------------------------------------------------------------------------------
!
!*      2.     Reading of grid
!              ---------------
!
IF (NRANK/=NPIO) INI = 0
!
!-------------------------------------------------------------------------------------
!
!*      4.     Reading of snow data
!              ---------------------
!
IF (YAREA(1:2)=='RO' .OR. YAREA(1:2)=='GA' .OR. YAREA(1:2)=='RF' .OR. YAREA(1:2)=='RD') THEN
  IF (.NOT. GTOWN) THEN
    TZSNOW%SCHEME='1-L'
    TZSNOW%NLAYER=1
    CALL ALLOCATE_GR_SNOW(TZSNOW,INI,IPATCH)
  ELSE
    CALL OPEN_AUX_IO_SURF(HFILE,HFILETYPE,YMASK)
    CALL READ_GR_SNOW(HFILETYPE,TRIM(YAREA),YPREFIX,INI,IPATCH,TZSNOW, &
                      HDIR='E',KVERSION=IVERSION,KBUGFIX=IBUGFIX)
    CALL CLOSE_AUX_IO_SURF(HFILE,HFILETYPE)
  ENDIF
ELSE
  CALL OPEN_AUX_IO_SURF(HFILE,HFILETYPE,YMASK)
  CALL READ_GR_SNOW(HFILETYPE,TRIM(YAREA),YPREFIX,INI,IPATCH,TZSNOW, &
                    HDIR='E',KVERSION=IVERSION,KBUGFIX=IBUGFIX)
  CALL CLOSE_AUX_IO_SURF(HFILE,HFILETYPE)
ENDIF
!
IF (TZSNOW%NLAYER.LT.KLAYER) THEN
   CALL DEALLOC_GR_SNOW(TZSNOW)
  CALL ABOR1_SFX("PREP_SNOW_EXTERN: SNOW NLAYER IN EXTERN FILE MUST BE GROWER THAN CURRENT NLAYER")
ELSEIF (TZSNOW%NLAYER.EQ.KLAYER) THEN
  OSNOW_IDEAL = .TRUE.
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      5.     Total snow content
!              ------------------
!
IF (NRANK==NPIO) THEN
  !
  SELECT CASE (HSURF(1:3))
    CASE ('WWW')
      IF (OSNOW_IDEAL) THEN
        ALLOCATE(PFIELD(INI,KLAYER,IPATCH))
        PFIELD(:,:,:) = TZSNOW%WSNOW(:,1:KLAYER,:)              
      ELSE
        ALLOCATE(PFIELD(INI,1,IPATCH))
        PFIELD = 0.
        DO JLAYER=1,TZSNOW%NLAYER
          PFIELD(:,1,:) = PFIELD(:,1,:) + TZSNOW%WSNOW(:,JLAYER,:)
        END DO 
        WHERE ( PFIELD(:,1,:)>XUNDEF ) PFIELD(:,1,:)=XUNDEF
      ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      6.     Albedo
!              ------
!
  CASE ('ALB')
    ALLOCATE(PFIELD(INI,1,IPATCH))
    PFIELD(:,1,:) = TZSNOW%ALB(:,:)
!
!-------------------------------------------------------------------------------------
!
!*      7.     Total depth to snow grid
!              ------------------------
!
  CASE ('DEP')
    ALLOCATE(PFIELD(INI,KLAYER,IPATCH))
    IF (OSNOW_IDEAL.OR.TZSNOW%NLAYER==KLAYER) THEN    
      PFIELD(:,:,:) = TZSNOW%WSNOW(:,1:KLAYER,:)/TZSNOW%RHO(:,1:KLAYER,:)
      WHERE(TZSNOW%WSNOW(:,1:KLAYER,:)==XUNDEF) PFIELD(:,:,:)=XUNDEF
    ELSE     
      ALLOCATE(ZD(INI,IPATCH))
      ZD(:,:) = 0.0
      DO JPATCH=1,IPATCH
        DO JLAYER=1,TZSNOW%NLAYER
          ZD(:,JPATCH) = ZD(:,JPATCH) + &
                        TZSNOW%WSNOW(:,JLAYER,JPATCH)/TZSNOW%RHO(:,JLAYER,JPATCH)
        END DO
        WHERE(TZSNOW%WSNOW(:,1,JPATCH)==XUNDEF) ZD(:,JPATCH)=XUNDEF
      ENDDO
      !* fine grid
      DO JPATCH=1,IPATCH
        CALL SNOW3LGRID(PFIELD(:,:,JPATCH),ZD(:,JPATCH))
      ENDDO  
      DEALLOCATE(ZD)
    ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      8.     Density or heat profile
!              -----------------------
!
  CASE ('RHO','HEA','SG1','SG2','HIS','AGE')
!
    SELECT CASE (TZSNOW%SCHEME)
      CASE ('D95','1-L','EBA')
        ALLOCATE(PFIELD(INI,1,IPATCH))   
        !* computes output physical variable
        IF (HSURF(1:3)=='RHO') PFIELD(:,1,:) = TZSNOW%RHO(:,1,:)
        IF (HSURF(1:3)=='HEA') THEN
          ALLOCATE(ZTEMP(INI,TZSNOW%NLAYER,IPATCH))
          ALLOCATE(ZWLIQ(INI,TZSNOW%NLAYER,IPATCH))
          IF (TZSNOW%SCHEME=='D95'.OR.TZSNOW%SCHEME=='EBA') ZTEMP(:,1,:) = XTT-2.
          IF (TZSNOW%SCHEME=='1-L') ZTEMP(:,1,:) = TZSNOW%T(:,1,:)
          ZWLIQ(:,:,:) = 0.0
          CALL SNOW_T_WLIQ_TO_HEAT(PFIELD,TZSNOW%RHO,ZTEMP,ZWLIQ)
          DEALLOCATE(ZTEMP)
          DEALLOCATE(ZWLIQ)
        END IF
        IF (HSURF(1:3)=='SG1') PFIELD(:,1,:) = -20.0
        IF (HSURF(1:3)=='SG2') PFIELD(:,1,:) = 80.0
        IF (HSURF(1:3)=='HIS') PFIELD(:,1,:) = 0.0
        IF (HSURF(1:3)=='AGE') PFIELD(:,1,:) = 3.0

      CASE ('3-L','CRO')
        ALLOCATE(ZFIELD(INI,TZSNOW%NLAYER,IPATCH))
        !* input physical variable
        IF (HSURF(1:3)=='RHO') ZFIELD(:,:,:) = TZSNOW%RHO (:,1:TZSNOW%NLAYER,:)
        IF (HSURF(1:3)=='HEA') ZFIELD(:,:,:) = TZSNOW%HEAT(:,1:TZSNOW%NLAYER,:)
        IF (HSURF(1:3)=='AGE') ZFIELD(:,:,:) = TZSNOW%AGE(:,1:KLAYER,:)
        IF (TZSNOW%SCHEME=='CRO')THEN
          IF (HSURF(1:3)=='SG1') ZFIELD(:,:,:) = TZSNOW%GRAN1(:,1:KLAYER,:)
          IF (HSURF(1:3)=='SG2') ZFIELD(:,:,:) = TZSNOW%GRAN2(:,1:KLAYER,:)
          IF (HSURF(1:3)=='HIS') ZFIELD(:,:,:) = TZSNOW%HIST(:,1:KLAYER,:)
        ELSE
         IF (HSURF(1:3)=='SG1') ZFIELD(:,:,:) = -20.0
         IF (HSURF(1:3)=='SG2') ZFIELD(:,:,:) = 80.0
         IF (HSURF(1:3)=='HIS') ZFIELD(:,:,:) = 0.0                  
        ENDIF    
        !
        ALLOCATE(PFIELD(INI,KLAYER,IPATCH))
        IF (OSNOW_IDEAL.OR.TZSNOW%NLAYER==KLAYER) THEN 
          PFIELD(:,:,:) = ZFIELD(:,:,:)
        ELSE
          !
          !* input snow layer thickness
          ALLOCATE(ZDEPTH(INI,TZSNOW%NLAYER,IPATCH))
          DO JPATCH=1,IPATCH
              ZDEPTH(:,:,JPATCH) = TZSNOW%WSNOW(:,:,JPATCH)/TZSNOW%RHO(:,:,JPATCH)
          END DO
          !
          !* total depth
          ALLOCATE(ZD(INI,IPATCH))
          ZD(:,:) = 0.
          DO JLAYER=1,TZSNOW%NLAYER
            ZD(:,:) = ZD(:,:) + ZDEPTH(:,JLAYER,:)
          ENDDO
          !
          !* input normalized grid
          ALLOCATE(ZGRID(INI,TZSNOW%NLAYER,IPATCH))
          DO JPATCH=1,IPATCH
             DO JI=1,INI
                IF(ZD(JI,JPATCH)==0.0)THEN
                   DO JLAYER = 1,TZSNOW%NLAYER
                      ZGRID(JI,JLAYER,JPATCH)=REAL(JLAYER)/REAL(TZSNOW%NLAYER)
                   ENDDO
                ELSE
                   DO JLAYER = 1,TZSNOW%NLAYER
                      IF(JLAYER==1)THEN
                        ZGRID(JI,JLAYER,JPATCH)=ZDEPTH(JI,JLAYER,JPATCH)/ ZD(JI,JPATCH)
                      ELSE
                        ZGRID(JI,JLAYER,JPATCH) = ZGRID(JI,JLAYER-1,JPATCH) + ZDEPTH(JI,JLAYER,JPATCH)/ZD(JI,JPATCH)
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
          DEALLOCATE(ZDEPTH)
          DEALLOCATE(ZD)
          !    
          !* interpolation of profile onto fine normalized snow grid
          DO JPATCH=1,IPATCH
            CALL INTERP_GRID_NAT(ZGRID(:,:,JPATCH),ZFIELD(:,:,JPATCH),    &
                             XGRID_SNOW(:), PFIELD(:,:,JPATCH))
          END DO
          DEALLOCATE(ZGRID)
        ENDIF
        DEALLOCATE(ZFIELD)

      END SELECT
      !* put field form patch to all vegtypes    
  END SELECT
!
ENDIF
!
CALL DEALLOC_GR_SNOW(TZSNOW)
!
!-------------------------------------------------------------------------------------
!
!*      9.     End of IO
!              ---------
!
IF (LHOOK) CALL DR_HOOK('PREP_SNOW_EXTERN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
END SUBROUTINE PREP_SNOW_EXTERN
