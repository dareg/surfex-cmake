!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE PREP_HOR_TEB_GARDEN_FIELD (DTCO, UG, U, USS, GCP, IO, S, K, P, PEK, TG, TOP,  &
                                      HPROGRAM,HSURF,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!     #################################################################################
!
!!****  *PREP_HOR_TEB_GARDEN_FIELD* - reads, interpolates and prepares an ISBA field
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      P. Le Moigne 10/2005, Phasage Arome
!!      P. Le Moigne 03/2007, Ajout initialisation par ascllv
!!      B. Decharme  01/2009, Optional Arpege deep soil temperature initialization
!!      B. Decharme  03/2014, external init with FA files
!!                            new vertical interpol
!!      P. Marguinaud10/2014, Support for a 2-part PREP
!!------------------------------------------------------------------
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_GRID_CONF_PROJ_n, ONLY : GRID_CONF_PROJ_t
!
USE MODD_ISBA_n, ONLY : ISBA_S_t, ISBA_K_t, ISBA_P_t, ISBA_PE_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
!
USE MODD_SFX_GRID_n, ONLY : GRID_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
!
USE MODD_TYPE_SNOW, ONLY: SURF_SNOW, NSURF_SNOW, TYPE_SNOW_INIT
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
USE MODD_GRID_GRIB, ONLY : CINMODEL
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM, NPROC, LSFX_MPI
USE MODD_PREP,            ONLY : CINGRID_TYPE, CINTERP_TYPE, XZS_LS, LINTERP, CMASK

USE MODD_PREP_TEB_GARDEN, ONLY : XGRID_SOIL, NGRID_LEVEL,                  &
                                 XWSNOW_GD, XRSNOW_GD, XTSNOW_GD, XLWCSNOW_GD, &
                                 XAGESNOW_GD, XASNOW_GD, LSNOW_IDEAL_GD

USE MODD_ISBA_PAR,        ONLY : XWGMIN
USE MODD_DATA_COVER_PAR,  ONLY : NVEGTYPE
USE MODD_SURF_PAR,        ONLY : XUNDEF
USE MODD_PREP_SNOW, ONLY : NIMPUR
!
USE MODI_PREP_GRIB_GRID
USE MODI_READ_PREP_TEB_GARDEN_CONF
USE MODI_READ_PREP_GARDEN_SNOW
USE MODI_PREP_TEB_GARDEN_ASCLLV
USE MODI_PREP_TEB_GARDEN_GRIB
USE MODI_PREP_TEB_GARDEN_UNIF
USE MODI_PREP_TEB_GARDEN_BUFFER
USE MODI_HOR_INTERPOL
USE MODI_PUT_ON_ALL_VEGTYPES
USE MODI_VEGTYPE_GRID_TO_PATCH_GRID
USE MODI_PREP_HOR_SNOW_FIELDS
USE MODI_GET_LUOUT
USE MODI_PREP_TEB_GARDEN_EXTERN
USE MODI_ABOR1_SFX
USE MODI_ALLOCATE_GR_SNOW
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
!*      0.1    declarations of arguments
!
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SSO_t), INTENT(INOUT) :: USS
TYPE(GRID_CONF_PROJ_t),INTENT(INOUT) :: GCP
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_K_t), INTENT(INOUT) :: K
TYPE(ISBA_P_t), INTENT(INOUT) :: P
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
TYPE(GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=7),   INTENT(IN)  :: HSURF     ! type of field
 CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
!
INTEGER,            INTENT(IN)  :: KPATCH
!
!*      0.2    declarations of local variables
!
 CHARACTER(LEN=6)              :: YFILETYPE ! type of input file
 CHARACTER(LEN=28)             :: YFILE     ! name of file
 CHARACTER(LEN=6)              :: YFILEPGDTYPE ! type of input file
 CHARACTER(LEN=28)             :: YFILEPGD     ! name of file
 CHARACTER(LEN=6)              :: YFILETYPE_SNOW ! type of input file
 CHARACTER(LEN=28)             :: YFILE_SNOW     ! name of file
 CHARACTER(LEN=6)              :: YFILEPGDTYPE_SNOW ! type of input file
 CHARACTER(LEN=28)             :: YFILEPGD_SNOW     ! name of file 
REAL, POINTER,     DIMENSION(:,:,:) :: ZFIELDIN=>NULL()  ! field to interpolate horizontally
!
TYPE(NSURF_SNOW) :: TNPSNOW
!
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZFIELDOUTP ! field interpolated   horizontally
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZFIELDOUTV !
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZVEGTYPE_PATCH ! vegtype for each patch
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZW        ! work array (x, fine   soil grid)
REAL, ALLOCATABLE, DIMENSION(:)     :: ZSUM
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZF        ! work array (x, output soil grid)
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZDG       ! out T grid (x, output soil grid)
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZPATCH    ! work array for patches
REAL, ALLOCATABLE, DIMENSION(:)     :: ZSG1SNOW, ZSG2SNOW, ZHISTSNOW
REAL, ALLOCATABLE, DIMENSION(:,:)      :: ZIMPURSNOW
INTEGER                             :: ILUOUT    ! output listing logical unit
!
TYPE (DATE_TIME)                :: TZTIME_GRIB    ! current date and time
LOGICAL                             :: GUNIF     ! flag for prescribed uniform field
LOGICAL                             :: GUNIF_SNOW! flag for prescribed uniform field
INTEGER                             :: JVEGTYPE, JPATCH  ! loop on vegtypes
INTEGER                             :: JLAYER    ! loop on layers
INTEGER                             :: JI, INP, INL, INI
INTEGER                             :: IWORK     ! Work integer
INTEGER :: INFOMPI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
!
!*      1.     Reading of input file name and type
!
IF (LHOOK) CALL DR_HOOK('PREP_HOR_TEB_GARDEN_FIELD',0,ZHOOK_HANDLE)
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
 CALL READ_PREP_TEB_GARDEN_CONF(HPROGRAM,HSURF,YFILE,YFILETYPE,YFILEPGD,YFILEPGDTYPE,&
                               HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,ILUOUT,GUNIF)
!
CMASK = 'TOWN  '
!
INI=SIZE(TG%XLAT)
!-------------------------------------------------------------------------------------
!
!*      2.     Snow variables case?
!
IF (HSURF=='SN_VEG ') THEN
  CALL READ_PREP_GARDEN_SNOW(HPROGRAM,PEK%TSNOW%SCHEME,PEK%TSNOW%NLAYER,YFILE_SNOW,&
        YFILETYPE_SNOW,YFILEPGD_SNOW,YFILEPGDTYPE_SNOW,GUNIF_SNOW)
  !
  IF(.NOT.GUNIF_SNOW.AND.LEN_TRIM(YFILE_SNOW)==0.AND.LEN_TRIM(YFILETYPE_SNOW)==0)THEN
    !IF(LEN_TRIM(YFILE)/=0.AND.LEN_TRIM(YFILETYPE)/=0)THEN
    IF (YFILETYPE=='GRIB') THEN
      YFILE_SNOW    =YFILE
      YFILETYPE_SNOW=YFILETYPE
      YFILEPGD_SNOW    =YFILEPGD
      YFILEPGDTYPE_SNOW=YFILEPGDTYPE       
    ELSE          
      GUNIF_SNOW=.TRUE.
      IF(ALL(XWSNOW_GD==XUNDEF))XWSNOW_GD=0.0 
    ENDIF 
  ENDIF    
  !
  ALLOCATE(ZSG1SNOW(SIZE(XWSNOW_GD)))
  ALLOCATE(ZSG2SNOW(SIZE(XWSNOW_GD)))
  ALLOCATE(ZHISTSNOW(SIZE(XWSNOW_GD)))
  ALLOCATE(ZIMPURSNOW(SIZE(XWSNOW_GD),NIMPUR))
  !
  ALLOCATE(TNPSNOW%AL(1))
  TNPSNOW%AL(1)%SCHEME = PEK%TSNOW%SCHEME
  TNPSNOW%AL(1)%NLAYER = PEK%TSNOW%NLAYER
  !
  CALL PREP_HOR_SNOW_FIELDS(DTCO, TG, U, GCP, HPROGRAM,HSURF, &
                            YFILE,YFILETYPE,                &
                            YFILEPGD, YFILEPGDTYPE,         &
                            ILUOUT,GUNIF_SNOW,1,KPATCH,     &
                            INI, TNPSNOW, TOP%TTIME,        &
                            XWSNOW_GD, XRSNOW_GD, XTSNOW_GD,&
                            XLWCSNOW_GD, XASNOW_GD,         &
                            LSNOW_IDEAL_GD, ZSG1SNOW,       &
                            ZSG2SNOW, ZHISTSNOW, XAGESNOW_GD, &
                            ZIMPURSNOW, &
                            PVEGTYPE_PATCH=S%XVEGTYPE_PATCH, PPATCH=S%XPATCH )
  !
  CALL ALLOCATE_GR_SNOW(PEK%TSNOW,INI)
  PEK%TSNOW%WSNOW = TNPSNOW%AL(1)%WSNOW
  PEK%TSNOW%RHO   = TNPSNOW%AL(1)%RHO
  PEK%TSNOW%ALB   = TNPSNOW%AL(1)%ALB
  IF (PEK%TSNOW%SCHEME/='D95') PEK%TSNOW%HEAT = TNPSNOW%AL(1)%HEAT
  IF (PEK%TSNOW%SCHEME=='CRO'.OR.PEK%TSNOW%SCHEME=='3-L') &
    PEK%TSNOW%AGE = TNPSNOW%AL(1)%AGE
  IF (PEK%TSNOW%SCHEME=='CRO') THEN
    PEK%TSNOW%GRAN1 = TNPSNOW%AL(1)%GRAN1
    PEK%TSNOW%GRAN2 = TNPSNOW%AL(1)%GRAN2
    PEK%TSNOW%HIST = TNPSNOW%AL(1)%HIST
  ENDIF
  !
  CALL TYPE_SNOW_INIT(TNPSNOW%AL(1))
  DEALLOCATE(TNPSNOW%AL)
  DEALLOCATE(ZSG1SNOW)
  DEALLOCATE(ZSG2SNOW)
  DEALLOCATE(ZHISTSNOW)  
  DEALLOCATE(ZIMPURSNOW)
  IF (LHOOK) CALL DR_HOOK('PREP_HOR_TEB_GARDEN_FIELD',1,ZHOOK_HANDLE)
  RETURN
END IF
!
!-------------------------------------------------------------------------------------
!
!*      3.     Reading of input  configuration (Grid and interpolation type)
!
IF (GUNIF) THEN
  CALL PREP_TEB_GARDEN_UNIF(ILUOUT,HSURF,ZFIELDIN)
ELSE IF (YFILETYPE=='ASCLLV') THEN
  CALL PREP_TEB_GARDEN_ASCLLV(DTCO, UG, U, USS, HPROGRAM,HSURF,ILUOUT,ZFIELDIN)
ELSE IF (YFILETYPE=='GRIB  ') THEN
  CALL PREP_GRIB_GRID(YFILE,ILUOUT,CINMODEL,CINGRID_TYPE,CINTERP_TYPE,TZTIME_GRIB)            
   IF (NRANK==NPIO) CALL PREP_TEB_GARDEN_GRIB(HPROGRAM,HSURF,YFILE,ILUOUT,ZFIELDIN)        
ELSE IF (YFILETYPE=='MESONH' .OR. YFILETYPE=='ASCII ' .OR. YFILETYPE=='LFI   '&
        .OR.YFILETYPE=='FA    '.OR. YFILETYPE=='AROME '.OR.YFILETYPE=='NC    ') THEN
   CALL PREP_TEB_GARDEN_EXTERN(DTCO, IO, U, GCP, &
                               HPROGRAM,HSURF,YFILE,YFILETYPE,YFILEPGD,YFILEPGDTYPE,ILUOUT,KPATCH,ZFIELDIN)
ELSE IF (YFILETYPE=='BUFFER') THEN
   CALL PREP_TEB_GARDEN_BUFFER(HPROGRAM,HSURF,ILUOUT,ZFIELDIN)
ELSE
   CALL ABOR1_SFX('PREP_HOR_TEB_GARDEN_FIELD: data file type not supported : '//YFILETYPE)
END IF
!
!-------------------------------------------------------------------------------------
!
!*      5.     Horizontal interpolation
!
IF (NRANK==NPIO) THEN
INL = SIZE(ZFIELDIN,2)
INP = SIZE(ZFIELDIN,3)
ELSE
 IF (.NOT.ASSOCIATED(ZFIELDIN)) ALLOCATE(ZFIELDIN(0,0,0))
ENDIF
!
IF (NPROC>1) THEN
#ifdef SFX_MPI
  IF (LSFX_MPI) THEN
    CALL MPI_BCAST(INL,KIND(INL)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(INP,KIND(INP)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
  ENDIF
#endif
ENDIF
!
ALLOCATE(ZFIELDOUTP(INI,INL,INP))
!
DO JPATCH = 1, INP
  IF (INP==NVEGTYPE) LINTERP = (S%XVEGTYPE(:,JPATCH) > 0.)
  CALL HOR_INTERPOL(DTCO, U, GCP, ILUOUT,ZFIELDIN(:,:,JPATCH),ZFIELDOUTP(:,:,JPATCH))
  LINTERP = .TRUE.
END DO
!
DEALLOCATE(ZFIELDIN )
!
ALLOCATE(ZW (INI,INL))
ZW = 0.
!
IF (1/=INP) THEN
!
ALLOCATE(ZFIELDOUTV(INI,INL,NVEGTYPE))
CALL PUT_ON_ALL_VEGTYPES(INI,INL,INP,NVEGTYPE,ZFIELDOUTP,ZFIELDOUTV)
!
  ALLOCATE(ZSUM (INI))
DO JLAYER=1,SIZE(ZW,2)
    ZSUM(:) = SUM(S%XVEGTYPE(:,:),2,ZFIELDOUTV(:,JLAYER,:)/=XUNDEF)
  DO JVEGTYPE=1,NVEGTYPE
    WHERE (ZFIELDOUTV(:,JLAYER,JVEGTYPE)/=XUNDEF) 
        ZW(:,JLAYER) = ZW(:,JLAYER) + S%XVEGTYPE(:,JVEGTYPE) * ZFIELDOUTV(:,JLAYER,JVEGTYPE) / ZSUM(:)
    END WHERE
  END DO
  DO JI=1,SIZE(ZW,1)
    IF (ALL(ZFIELDOUTV(JI,JLAYER,:)==XUNDEF)) ZW(JI,JLAYER) = XUNDEF
  ENDDO
END DO
  DEALLOCATE(ZFIELDOUTV)
  DEALLOCATE(ZSUM)
  !
ELSE
  !
  ZW(:,:) = ZFIELDOUTP(:,:,1)
  !
ENDIF
!
DEALLOCATE(ZFIELDOUTP)
!
!-------------------------------------------------------------------------------------
!
!*      7.     Return to historical variable
!
!
SELECT CASE (HSURF)
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('WG     ') 
  ALLOCATE(ZF (INI,IO%NGROUND_LAYER))
  !
  !* interpolates on output levels
  CALL INIT_FROM_REF_GRID(XGRID_SOIL,ZW,P%XDG(:,:),ZF)
  !
  !* retrieves soil water content from soil relative humidity
  ALLOCATE(PEK%XWG(INI,IO%NGROUND_LAYER))
  PEK%XWG(:,:) = K%XWWILT + ZF(:,:) * (K%XWFC-K%XWWILT)
  PEK%XWG(:,:) = MAX(MIN(PEK%XWG(:,:),K%XWSAT),XWGMIN)
  !
  WHERE(ZF(:,:)==XUNDEF)PEK%XWG(:,:)=XUNDEF
  !
  DEALLOCATE(ZF)
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('WGI    ')
  ALLOCATE(ZF (INI,IO%NGROUND_LAYER))
  !
  !* interpolates on output levels
  CALL INIT_FROM_REF_GRID(XGRID_SOIL,ZW,P%XDG(:,:),ZF) 
  !
  !* retrieves soil ice content from soil relative humidity
  ALLOCATE(PEK%XWGI(INI,IO%NGROUND_LAYER))
  PEK%XWGI(:,:) = ZF(:,:) * K%XWSAT
  PEK%XWGI(:,:) = MAX(MIN(PEK%XWGI(:,:),K%XWSAT),0.)
  !
  WHERE(ZF(:,:)==XUNDEF)PEK%XWGI(:,:)=XUNDEF
  !
  DEALLOCATE(ZF)
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('TG     ') 
  IWORK=IO%NGROUND_LAYER
  ALLOCATE(PEK%XTG(INI,IWORK))
  ALLOCATE(ZDG(SIZE(P%XDG,1),IWORK))
  IF (IO%CISBA=='2-L'.OR.IO%CISBA=='3-L') THEN
    ZDG(:,1) = 0.01
    ZDG(:,2) = 0.40   ! deep temperature for force-restore taken at 20cm
    IF(IO%CISBA=='3-L') ZDG(:,3) = 5.00   ! climatological temperature, usually not used
  ELSE
    !* diffusion method, the soil grid is the same as for humidity
    ZDG(:,:) = P%XDG(:,:)
  END IF
  CALL INIT_FROM_REF_GRID(XGRID_SOIL,ZW,ZDG,PEK%XTG(:,:))
  DEALLOCATE(ZDG)
  !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('WR     ') 
  ALLOCATE(PEK%XWR(INI))
  PEK%XWR(:) = ZW(:,1)
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('LAI    ') 
  !* LAI is updated only if present and pertinent (evolutive LAI) in input file

   WHERE (ZW(:,1)/=XUNDEF) PEK%XLAI(:) = ZW(:,1)
  !
END SELECT
!
DEALLOCATE(ZW)
!-------------------------------------------------------------------------------------
!
!*      8.     Deallocations
!
!
IF (LHOOK) CALL DR_HOOK('PREP_HOR_TEB_GARDEN_FIELD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!
CONTAINS
!
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!
SUBROUTINE INIT_FROM_REF_GRID(PGRID1,PT1,PD2,PT2)
!
USE MODI_INTERP_GRID_NAT
!
REAL, DIMENSION(:,:), INTENT(IN)  :: PT1    ! variable profile
REAL, DIMENSION(:),   INTENT(IN)  :: PGRID1 ! normalized grid
REAL, DIMENSION(:,:), INTENT(IN)  :: PD2    ! output layer thickness
REAL, DIMENSION(:,:), INTENT(OUT) :: PT2    ! variable profile
!
INTEGER                                  :: JI, JL  ! loop counter
REAL, DIMENSION(SIZE(PT1),SIZE(PT1,2)) :: ZD1 ! input grid
!
INTEGER :: ILAYER1, ILAYER2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (LHOOK) CALL DR_HOOK('INIT_FROM_REF_GRID',0,ZHOOK_HANDLE)
!
IF (SIZE(PT1,2)==3) THEN
!
!* 1. case with only 3 input levels (typically coming from 'UNIF')
!     -----------------------------
!
  IF (IO%CISBA=='2-L' .OR. IO%CISBA=='3-L') THEN
    !* Possible LTEMP_ARP case
    IF(SIZE(PT2,2)>3)THEN
       ILAYER1=3
       ILAYER2=SIZE(PT2,2)
    ELSE
       ILAYER1=SIZE(PT2,2)
       ILAYER2=0
    ENDIF
    !* historical 2L or 3L ISBA version
    PT2(:,1:ILAYER1) = PT1(:,1:ILAYER1) 
    !* Possible LTEMP_ARP case
    IF(ILAYER2>0)THEN
       DO JL=ILAYER1+1,ILAYER2
         PT2(:,JL) = PT2(:,ILAYER1)
       ENDDO
    ENDIF
!    
  ELSEIF(IO%CISBA=='DIF')THEN
       !surface layer (generally 0.01m imposed)
       PT2(:,1) = PT1(:,1) 
       !deep layers
       DO JL=2,IO%NGROUND_LAYER
          PT2(:,JL) = PT1(:,3)
       END DO
       !if root layers
       DO JI=1,SIZE(PT1,1)
          DO JL=2,IO%NGROUND_LAYER
             IF(P%XROOTFRAC(JI,JL)<=1.0)THEN 
                PT2(JI,JL) = PT1(JI,2)
                EXIT
             ENDIF
          END DO
       END DO 
  END IF    
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ELSE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!* 2. case with fine grid as input (general case)
!     ----------------------------
!
  DO JL=1,SIZE(PT1,2)
    ZD1(:,JL) = PGRID1(JL)
  END DO
!
  CALL INTERP_GRID_NAT(ZD1,PT1(:,:),PD2,PT2(:,:))
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (LHOOK) CALL DR_HOOK('INIT_FROM_REF_GRID',1,ZHOOK_HANDLE)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
END SUBROUTINE INIT_FROM_REF_GRID
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_HOR_TEB_GARDEN_FIELD
