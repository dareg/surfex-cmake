!     #########
SUBROUTINE PREP_HOR_ISBA_FIELD(HPROGRAM,HSURF,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,OKEY)
!     #################################################################################
!
!!****  *PREP_HOR_ISBA_FIELD* - reads, interpolates and prepares an ISBA field
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
!!      M. Lafaysse  07/2012, allow netcdf input files
!!      B. Decharme  07/2012, Bug init uniform snow
!!      M. Lafaysse 11/2012,  snow liquid water content
!!      B. Decharme  03/2014, external init with FA files
!!                            new vertical interpolation
!!      P Samuelsson 10/2014, MEB
!!------------------------------------------------------------------
!
USE MODD_PREP,     ONLY : XZS_LS, LINTERP, CMASK

USE MODD_PREP_ISBA, ONLY : XGRID_SOIL, NGRID_LEVEL, LSNOW_IDEAL,    &
                           XWSNOW, XRSNOW, XTSNOW, XLWCSNOW, XASNOW,          &
                           XSG1SNOW, XSG2SNOW, XHISTSNOW, XAGESNOW

USE MODD_ISBA_n,    ONLY : TTIME, XWG, XWGI, XTG, XWR, XLAI,         &
                           NGROUND_LAYER, NPATCH, NWG_LAYER,         &
                           XVEGTYPE_PATCH, XDG, XWWILT, XWFC, XPATCH,&
                           CISBA, XDG2, XWSAT, TSNOW, XVEGTYPE,      &
                           LTEMP_ARP, NTEMPLAYER_ARP, XICE_STO,      &
                           XVEGTYPE,                                 &
                           XWRV, XWRVN, XTV, XTC, XQC

USE MODD_ISBA_GRID_n,    ONLY : XLAT, XLON
USE MODD_ISBA_PAR,       ONLY : XWGMIN
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF,NUNDEF
!
USE MODI_READ_PREP_ISBA_CONF
USE MODI_READ_PREP_ISBA_SNOW
USE MODI_PREP_ISBA_ASCLLV
USE MODI_PREP_ISBA_GRIB
USE MODI_PREP_ISBA_UNIF
USE MODI_PREP_ISBA_BUFFER
USE MODI_ABOR1_SFX
USE MODI_HOR_INTERPOL
USE MODI_PUT_ON_ALL_VEGTYPES
USE MODI_VEGTYPE_GRID_TO_PATCH_GRID
USE MODI_PREP_HOR_SNOW_FIELDS
USE MODI_GET_LUOUT
USE MODI_PREP_ISBA_EXTERN
USE MODI_PREP_ISBA_NETCDF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=7),   INTENT(IN)  :: HSURF     ! type of field
 CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
!
LOGICAL, OPTIONAL,  INTENT(INOUT):: OKEY
!
!*      0.2    declarations of local variables
!
 CHARACTER(LEN=6)              :: YFILETYPE ! type of input file
 CHARACTER(LEN=28)             :: YFILE     ! name of file
 CHARACTER(LEN=6)              :: YFILETYPE_SNOW ! type of input file
 CHARACTER(LEN=28)             :: YFILE_SNOW     ! name of file
 CHARACTER(LEN=6)              :: YFILEPGDTYPE_SNOW ! type of input file
 CHARACTER(LEN=28)             :: YFILEPGD_SNOW     ! name of file 
 CHARACTER(LEN=6)              :: YFILEPGDTYPE ! type of input file
 CHARACTER(LEN=28)             :: YFILEPGD     ! name of file
REAL, POINTER, DIMENSION(:,:,:)     :: ZFIELDIN  ! field to interpolate horizontally
REAL, POINTER, DIMENSION(:,:)       :: ZFIELD ! field to interpolate horizontally
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZFIELDOUTP ! field interpolated   horizontally
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZFIELDOUTV !
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZW        ! work array (x, fine   soil grid, npatch)
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZF        ! work array (x, output soil grid, npatch)
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZDG       ! out T grid (x, output soil grid, npatch)
INTEGER                       :: ILUOUT    ! output listing logical unit
!
LOGICAL                       :: GUNIF     ! flag for prescribed uniform field
LOGICAL                       :: GUNIF_SNOW! flag for prescribed uniform field
INTEGER                       :: JPATCH    ! loop on patches
INTEGER                       :: JVEGTYPE  ! loop on vegtypes
INTEGER                       :: INI, INL, INP, JJ, JL! Work integer
INTEGER, DIMENSION(SIZE(XDG,1),SIZE(XDG,3)) :: IWORK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
!
!*      1.     Reading of input file name and type
!
IF (LHOOK) CALL DR_HOOK('PREP_HOR_ISBA_FIELD',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
 CALL READ_PREP_ISBA_CONF(HPROGRAM,HSURF,YFILE,YFILETYPE,YFILEPGD,YFILEPGDTYPE,&
                         HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,ILUOUT,GUNIF)
!
CMASK = 'NATURE'
!
INI=SIZE(XLAT)
!
!-------------------------------------------------------------------------------------
!
!*      2.     Snow variables case?
!
IF (HSURF=='SN_VEG ') THEN
  CALL READ_PREP_ISBA_SNOW(HPROGRAM,TSNOW%SCHEME,TSNOW%NLAYER,YFILE_SNOW,YFILETYPE_SNOW,&
                                YFILEPGD_SNOW,YFILEPGDTYPE_SNOW,GUNIF_SNOW)
  IF(.NOT.GUNIF_SNOW.AND.LEN_TRIM(YFILE_SNOW)==0.AND.LEN_TRIM(YFILETYPE_SNOW)==0)THEN
    IF(LEN_TRIM(YFILE)/=0.AND.LEN_TRIM(YFILETYPE)/=0)THEN
       YFILE_SNOW    =YFILE
       YFILETYPE_SNOW=YFILETYPE
       YFILEPGD_SNOW    =YFILEPGD
       YFILEPGDTYPE_SNOW=YFILEPGDTYPE       
    ELSE
       GUNIF_SNOW=.TRUE.
       IF(ALL(XWSNOW==XUNDEF))XWSNOW=0.0
    ENDIF
  ENDIF
  CALL PREP_HOR_SNOW_FIELDS(HPROGRAM, HSURF,                     &
                            YFILE_SNOW, YFILETYPE_SNOW,          &
                            YFILEPGD_SNOW, YFILEPGDTYPE_SNOW,    &
                            ILUOUT, GUNIF_SNOW, NPATCH,          &
                            INI,TSNOW, TTIME,                    &
                            XWSNOW, XRSNOW, XTSNOW, XLWCSNOW,    &
                            XASNOW, LSNOW_IDEAL, XSG1SNOW,       &
                            XSG2SNOW, XHISTSNOW, XAGESNOW,       &
                            XVEGTYPE, XVEGTYPE_PATCH, XPATCH,    &
                            OKEY                                 )
  DEALLOCATE(XWSNOW)
  DEALLOCATE(XRSNOW)
  DEALLOCATE(XTSNOW)
  DEALLOCATE(XLWCSNOW)
  DEALLOCATE(XSG1SNOW)
  DEALLOCATE(XSG2SNOW)
  DEALLOCATE(XHISTSNOW)
  DEALLOCATE(XAGESNOW)
  IF (LHOOK) CALL DR_HOOK('PREP_HOR_ISBA_FIELD',1,ZHOOK_HANDLE)
  RETURN
END IF
!
!-------------------------------------------------------------------------------------
!
!*      3.     Reading of input  configuration (Grid and interpolation type)
!
IF (GUNIF) THEN
  CALL PREP_ISBA_UNIF(ILUOUT,HSURF,ZFIELDIN)
ELSE IF (YFILETYPE=='ASCLLV') THEN
  CALL PREP_ISBA_ASCLLV(HPROGRAM,HSURF,ILUOUT,ZFIELDIN)
ELSE IF (YFILETYPE=='GRIB  ') THEN
  CALL PREP_ISBA_GRIB(HPROGRAM,HSURF,YFILE,ILUOUT,ZFIELDIN,OKEY)
ELSE IF (YFILETYPE=='MESONH' .OR. YFILETYPE=='ASCII ' .OR. YFILETYPE=='LFI   '.OR.YFILETYPE=='FA    ') THEN
   CALL PREP_ISBA_EXTERN(HPROGRAM,HSURF,YFILE,YFILETYPE,YFILEPGD,YFILEPGDTYPE,ILUOUT,ZFIELDIN,OKEY)
ELSE IF (YFILETYPE=='BUFFER') THEN
   CALL PREP_ISBA_BUFFER(HPROGRAM,HSURF,ILUOUT,ZFIELDIN)
ELSE IF (YFILETYPE=='NETCDF') THEN
   CALL PREP_ISBA_NETCDF(HPROGRAM,HSURF,YFILE,ILUOUT,ZFIELDIN)
ELSE
   CALL ABOR1_SFX('PREP_HOR_ISBA_FIELD: data file type not supported : '//YFILETYPE)
END IF
!
!-------------------------------------------------------------------------------------
!
!*      5.     Horizontal interpolation
!
INL = SIZE(ZFIELDIN,2)
INP = SIZE(ZFIELDIN,3)
!
ALLOCATE(ZFIELDOUTP(INI,INL,INP))
ALLOCATE(ZFIELD(SIZE(ZFIELDIN,1),INL))
!
DO JPATCH = 1, INP
  ZFIELD=ZFIELDIN(:,:,JPATCH)
  IF (INP==NVEGTYPE) LINTERP = (XVEGTYPE(:,JPATCH) > 0.)
  CALL HOR_INTERPOL(ILUOUT,ZFIELD,ZFIELDOUTP(:,:,JPATCH))
  LINTERP = .TRUE.
END DO
!
DEALLOCATE(ZFIELD)
!
ALLOCATE(ZFIELDOUTV(INI,INL,NVEGTYPE))
!
CALL PUT_ON_ALL_VEGTYPES(INI,INL,INP,NVEGTYPE,ZFIELDOUTP,ZFIELDOUTV)
!
DEALLOCATE(ZFIELDOUTP)
!
!-------------------------------------------------------------------------------------
!
!*      6.     Transformation from vegtype grid to patch grid
!
ALLOCATE(ZW (INI,SIZE(ZFIELDOUTV,2),NPATCH))
!
ZW = 0.
CALL VEGTYPE_GRID_TO_PATCH_GRID(NPATCH,XVEGTYPE_PATCH,XPATCH,ZFIELDOUTV,ZW)
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
 CASE('ZS     ') 
  ALLOCATE(XZS_LS(INI))
  XZS_LS(:) = ZFIELDOUTV(:,1,1)
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('WG     ') 
  ALLOCATE(ZF (INI,NGROUND_LAYER,NPATCH))
  !
  !* interpolates on output levels
  CALL INIT_FROM_REF_GRID(XGRID_SOIL,ZW,XDG,ZF)
  !
  !* retrieves soil water content from soil relative humidity
  ALLOCATE(XWG(INI,NGROUND_LAYER,NPATCH))
  XWG(:,:,:)=XUNDEF
  IF(CISBA=='DIF')THEN
     IWORK(:,:)=NWG_LAYER(:,:)
  ELSE
     IWORK(:,:)=SIZE(XWG,2)
  ENDIF
  DO JPATCH=1,NPATCH
    DO JJ=1,INI
       IF(IWORK(JJ,JPATCH)==NUNDEF)CYCLE
       INL=IWORK(JJ,JPATCH)
       DO JL=1,INL
          XWG(JJ,JL,JPATCH) = XWWILT(JJ,JL) + ZF(JJ,JL,JPATCH) * (XWFC(JJ,JL)-XWWILT(JJ,JL))
          XWG(JJ,JL,JPATCH) = MAX(MIN(XWG(JJ,JL,JPATCH),XWSAT(JJ,JL)),XWGMIN)
       ENDDO
    ENDDO
  ENDDO
  !
  WHERE(ZF(:,:,:)==XUNDEF)XWG(:,:,:)=XUNDEF
  !
  DEALLOCATE(ZF)
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('WGI    ')
  ALLOCATE(ZF (INI,NGROUND_LAYER,NPATCH))
  !
  !* interpolates on output levels
  CALL INIT_FROM_REF_GRID(XGRID_SOIL,ZW,XDG,ZF)
  !
  !* retrieves soil ice content from soil relative humidity
  ALLOCATE(XWGI(INI,NGROUND_LAYER,NPATCH))
  XWGI(:,:,:)=XUNDEF
  IF(CISBA=='DIF')THEN
     IWORK(:,:)=NWG_LAYER(:,:)
  ELSE
     IWORK(:,:)=2
  ENDIF  
  DO JPATCH=1,NPATCH
    DO JJ=1,INI
       IF(IWORK(JJ,JPATCH)==NUNDEF)CYCLE
       INL=IWORK(JJ,JPATCH)
       DO JL=1,INL
          XWGI(JJ,JL,JPATCH) = ZF(JJ,JL,JPATCH) * XWSAT(JJ,JL)
          XWGI(JJ,JL,JPATCH) = MAX(MIN(XWGI(JJ,JL,JPATCH),XWSAT(JJ,JL)),0.)
       ENDDO
    ENDDO
  END DO
  !
  WHERE(ZF  (:,:,:)==XUNDEF )XWGI(:,:,:)=XUNDEF
  WHERE(XWGI(:,:,:)<=1.0E-10)XWGI(:,:,:)=0.0
  !
  DEALLOCATE(ZF)
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('TG     ') 
  IF(LTEMP_ARP)THEN
    INL=NTEMPLAYER_ARP
  ELSE
    INL=NGROUND_LAYER
  ENDIF
  ALLOCATE(XTG(INI,INL,NPATCH))
  ALLOCATE(ZDG(SIZE(XDG,1),INL,SIZE(XDG,3)))
  IF (CISBA=='2-L'.OR.CISBA=='3-L') THEN
     DO JPATCH=1,NPATCH
        ZDG(:,1,JPATCH) = 0.01
        ZDG(:,2,JPATCH) = 0.40                    ! deep temperature for force-restore taken at 20cm
        IF(CISBA=='3-L') ZDG(:,3,JPATCH) = 5.00   ! climatological temperature, usually not used
     ENDDO         
     IF(LTEMP_ARP)THEN
       DO JPATCH=1,NPATCH
          ZDG(:,3,JPATCH) = 1.0
          DO JL=4,INL
             ZDG(:,JL,JPATCH) = ZDG(:,JL-1,JPATCH)+1.0
          ENDDO
       ENDDO
     ENDIF
  ELSE
    !* diffusion method, the soil grid is the same as for humidity
    ZDG(:,:,:) = XDG(:,:,:)
  END IF
  CALL INIT_FROM_REF_GRID(XGRID_SOIL,ZW,ZDG,XTG)
  DEALLOCATE(ZDG)
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('WR     ') 
  ALLOCATE(XWR(INI,NPATCH))
  DO JPATCH=1,NPATCH
    XWR(:,JPATCH) = ZW(:,1,JPATCH)
  END DO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('WRV    ') 
  ALLOCATE(XWRV(INI,NPATCH))
  DO JPATCH=1,NPATCH
    XWRV(:,JPATCH) = ZW(:,1,JPATCH)
  END DO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('WRVN   ') 
  ALLOCATE(XWRVN(INI,NPATCH))
  DO JPATCH=1,NPATCH
    XWRVN(:,JPATCH) = ZW(:,1,JPATCH)
  END DO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('TV     ') 
  ALLOCATE(XTV(INI,NPATCH))
  DO JPATCH=1,NPATCH
    XTV(:,JPATCH) = ZW(:,1,JPATCH)
  END DO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('TC     ') 
  ALLOCATE(XTC(INI,NPATCH))
  DO JPATCH=1,NPATCH
    XTC(:,JPATCH) = ZW(:,1,JPATCH)
  END DO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('QC     ') 
  ALLOCATE(XQC(INI,NPATCH))
  DO JPATCH=1,NPATCH
    XQC(:,JPATCH) = ZW(:,1,JPATCH)
  END DO
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('LAI    ') 
  !* LAI is updated only if present and pertinent (evolutive LAI) in input file
  IF (ANY(ZW(:,:,:)/=XUNDEF)) THEN
    DO JPATCH=1,NPATCH
      XLAI(:,JPATCH) = ZW(:,1,JPATCH)
    END DO
  END IF
  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
 CASE('ICE_STO') 
  ALLOCATE(XICE_STO(INI,NPATCH))
  DO JPATCH=1,NPATCH
    XICE_STO(:,JPATCH) = ZW(:,1,JPATCH)
  END DO
  !
END SELECT
!
DEALLOCATE(ZW)
!-------------------------------------------------------------------------------------
!
!*      8.     Deallocations
!
DEALLOCATE(ZFIELDIN )
DEALLOCATE(ZFIELDOUTV)
!
IF (LHOOK) CALL DR_HOOK('PREP_HOR_ISBA_FIELD',1,ZHOOK_HANDLE)
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
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PT1    ! variable profile
REAL, DIMENSION(:),     INTENT(IN)  :: PGRID1 ! normalized grid
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PD2    ! output layer thickness
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PT2    ! variable profile
!
INTEGER                                  :: JI,JL  ! loop counter
REAL, DIMENSION(SIZE(PT1,1),SIZE(PT1,2)) :: ZD1 ! input grid
REAL, DIMENSION(SIZE(PD2,1),SIZE(PD2,2)) :: ZD2 ! output grid
!
INTEGER :: ILAYER1, ILAYER2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (LHOOK) CALL DR_HOOK('INIT_FROM_REF_GRID',0,ZHOOK_HANDLE)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
IF (SIZE(PT1,2)==3) THEN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!* 1. case with only 3 input levels (typically coming from 'UNIF')
!     -----------------------------
!
  IF (CISBA=='2-L' .OR. CISBA=='3-L') THEN
    !* Possible LTEMP_ARP case
    IF(SIZE(PT2,2)>3)THEN
      ILAYER1=3
      ILAYER2=SIZE(PT2,2)
    ELSE
      ILAYER1=SIZE(PT2,2)
      ILAYER2=0
    ENDIF
    !* historical 2L or 3L ISBA version
    DO JPATCH=1,NPATCH
      PT2(:,1:ILAYER1,JPATCH) = PT1(:,1:ILAYER1,JPATCH) 
      !* Possible LTEMP_ARP case
      IF(ILAYER2>0)THEN
        DO JL=ILAYER1+1,ILAYER2
          PT2(:,JL,JPATCH) = PT2(:,ILAYER1,JPATCH)
        ENDDO
      ENDIF
    END DO
  ELSEIF(CISBA=='DIF')THEN
    DO JPATCH=1,NPATCH
       !surface layer (generally 0.01m imposed)
       PT2(:,1,JPATCH) = PT1(:,1,JPATCH) 
       !second layer
       PT2(:,2,JPATCH) = 0.25*PT1(:,1,JPATCH)+0.75*PT1(:,2,JPATCH)
       !others layers
       DO JI=1,SIZE(PT1,1)
         DO JL=3,NGROUND_LAYER
           IF(PD2(JI,JL,JPATCH)<=XDG2(JI,JPATCH))THEN 
            !root layers
             PT2(JI,JL,JPATCH) = PT1(JI,2,JPATCH)
           ELSE
            !deep layers
             PT2(JI,JL,JPATCH) = PT1(JI,3,JPATCH)
           ENDIF
         END DO
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
  ENDDO
!
  DO JPATCH=1,NPATCH
     ZD2(:,:) = PD2(:,:,JPATCH)
     CALL INTERP_GRID_NAT(ZD1,PT1(:,:,JPATCH),ZD2,PT2(:,:,JPATCH))
  ENDDO
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
END SUBROUTINE PREP_HOR_ISBA_FIELD
