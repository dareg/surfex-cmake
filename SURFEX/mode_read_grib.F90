!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #####################
MODULE MODE_READ_GRIB
!     #####################
!-------------------------------------------------------------------
!
USE MODI_ABOR1_SFX
USE GRIB_API
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
CONTAINS
!-------------------------------------------------------------------
!     ####################
      SUBROUTINE MAKE_GRIB_INDEX(HGRIB)
!     ####################
!
USE MODD_GRID_GRIB, ONLY : CGRIB_FILE, NIDX,NIDX2
!
IMPLICIT NONE
!
 CHARACTER(LEN=*), INTENT(IN) :: HGRIB
!
INTEGER(KIND=kindOfInt) :: IRET
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:MAKE_GRIB_INDEX',0,ZHOOK_HANDLE)
!
IF (CGRIB_FILE==HGRIB .AND. LHOOK) CALL DR_HOOK('MODE_READ_GRIB:MAKE_GRIB_INDEX',1,ZHOOK_HANDLE)
IF (CGRIB_FILE==HGRIB) RETURN
!
CGRIB_FILE=HGRIB
!
 CALL GRIB_INDEX_CREATE(NIDX,HGRIB,'paramId',IRET)
 CALL GRIB_INDEX_CREATE(NIDX2,HGRIB,'shortName',IRET)

IF (IRET/=0) CALL ABOR1_SFX("MODE_READ_GRIB:MAKE_GRIB_INDEX: error while creating the grib index")
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:MAKE_GRIB_INDEX',1,ZHOOK_HANDLE)
!
END SUBROUTINE MAKE_GRIB_INDEX
!-------------------------------------------------------------------
!     ####################
      SUBROUTINE CLEAR_GRIB_INDEX
!     ####################
!
USE MODD_GRID_GRIB, ONLY : CGRIB_FILE, NIDX,NIDX2
!
IMPLICIT NONE
!
INTEGER(KIND=kindOfInt) :: IRET
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:CLEAR_GRIB_INDEX',0,ZHOOK_HANDLE)
!
IF (CGRIB_FILE.NE."") THEN
  CGRIB_FILE=""
  CALL GRIB_INDEX_RELEASE(NIDX,IRET)
  IF (IRET/=0) CALL ABOR1_SFX("MODE_READ_GRIB:MAKE_GRIB_INDEX: error while deleting the grib index NIDX")
  CALL GRIB_INDEX_RELEASE(NIDX2,IRET)
  IF (IRET/=0) CALL ABOR1_SFX("MODE_READ_GRIB:MAKE_GRIB_INDEX: error while deleting the grib index NIDX2")
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:CLEAR_GRIB_INDEX',1,ZHOOK_HANDLE)
!
END SUBROUTINE CLEAR_GRIB_INDEX
!-------------------------------------------------------------------
!     ####################
      SUBROUTINE GET_GRIB_MESSAGE(KLUOUT,KLTYPE,KLEV1,KLEV2,KGRIB,KFOUND)
!     ####################
!
USE MODD_GRID_GRIB, ONLY : NIDX
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KLUOUT
INTEGER, INTENT(INOUT)  :: KLTYPE
INTEGER, INTENT(INOUT)  :: KLEV1
INTEGER, INTENT(INOUT)  :: KLEV2
INTEGER(KIND=kindOfInt), INTENT(INOUT) :: KGRIB
INTEGER, INTENT(OUT) :: KFOUND
!
INTEGER :: ILTYPE
INTEGER :: ILEV1
INTEGER :: ILEV2
INTEGER(KIND=kindOfInt) :: IRET
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:GET_GRIB_MESSAGE',0,ZHOOK_HANDLE)
!
IRET = 0
KFOUND=0
!
DO WHILE (IRET /= GRIB_END_OF_INDEX .AND. KFOUND/=3)
  !
  IRET = 0
  KFOUND=0
  !
  IF (KLTYPE/=-2) THEN
    CALL GRIB_GET(KGRIB,'indicatorOfTypeOfLevel',ILTYPE,IRET)
    CALL TEST_IRET(KLUOUT,ILTYPE,KLTYPE,IRET)
  ENDIF
  !
  IF (IRET.EQ.0) THEN
    !
    KFOUND = KFOUND + 1
    !
    IF (KLEV1/=-2) THEN
      CALL GRIB_GET(KGRIB,'topLevel',ILEV1,IRET)
      CALL TEST_IRET(KLUOUT,ILEV1,KLEV1,IRET)
    ENDIF
    !
    IF (IRET.EQ.0) THEN
      !
      KFOUND = KFOUND + 1
      !
      IF (KLEV2/=-2) THEN
        CALL GRIB_GET(KGRIB,'bottomLevel',ILEV2,IRET)
        CALL TEST_IRET(KLUOUT,ILEV2,KLEV2,IRET)
      ENDIF
      !
      IF (IRET.EQ.0) KFOUND = KFOUND + 1
      !
    ENDIF
    !
  ENDIF
  !
  IF (KFOUND.NE.3) THEN
    CALL GRIB_RELEASE(KGRIB)
    CALL GRIB_NEW_FROM_INDEX(NIDX,KGRIB,IRET)
  ENDIF
  !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:GET_GRIB_MESSAGE',1,ZHOOK_HANDLE)
!
CONTAINS
!
!       ##############
        SUBROUTINE TEST_IRET(KLUOUT,VAL1,VAL0,KRET)
!       ##############
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KLUOUT ! logical unit of output listing
INTEGER, INTENT(IN) :: VAL1
INTEGER, INTENT(INOUT) :: VAL0
INTEGER(KIND=kindOfInt), INTENT(INOUT) :: KRET   ! number of the message researched
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:TEST_IRET',0,ZHOOK_HANDLE)
!
IF (KRET > 0) THEN
  WRITE (KLUOUT,'(A)')' | Error encountered in the Grib file, skipping field'
ELSE IF (KRET == -6) THEN
  WRITE (KLUOUT,'(A)')' | ECMWF pseudo-Grib data encountered, skipping field'
ELSEIF (VAL1 /= VAL0) THEN
  IF (VAL0 == -1) THEN
    VAL0 = VAL1
  ELSE
    KRET=1
  ENDIF
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:TEST_IRET',1,ZHOOK_HANDLE)
END SUBROUTINE TEST_IRET
!
END SUBROUTINE GET_GRIB_MESSAGE
!-------------------------------------------------------------------
!     ####################
      SUBROUTINE READ_GRIB(HGRIB,KLUOUT,KPARAM,KRET,PFIELD, &
                         & KLTYPE,KLEV1,KLEV2,KPARAM2,      &
                         & SHORTNAME,PLATS,PLONS)
!     ####################
!
USE MODD_GRID_GRIB, ONLY : NIDX,NIDX2
!
IMPLICIT NONE
!
 CHARACTER(LEN=*), INTENT(IN)           :: HGRIB      ! name of the GRIB file
INTEGER, INTENT(IN)                     :: KLUOUT
INTEGER,INTENT(IN)                      :: KPARAM ! Parameter to read
INTEGER(KIND=kindOfInt), INTENT(OUT)    :: KRET
REAL, DIMENSION(:), POINTER             :: PFIELD
INTEGER,INTENT(INOUT), OPTIONAL         :: KLTYPE ! Level type
INTEGER,INTENT(INOUT), OPTIONAL         :: KLEV1  ! Level parameter 1
INTEGER,INTENT(INOUT), OPTIONAL         :: KLEV2  ! Level parameter 2
INTEGER, INTENT(INOUT), OPTIONAL        :: KPARAM2
 CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: SHORTNAME
REAL, DIMENSION(:), POINTER, OPTIONAL   :: PLATS, PLONS
!
INTEGER :: ILTYPE, ILEV1, ILEV2
INTEGER(KIND=kindOfInt) :: IGRIB
INTEGER :: ISIZE, IFOUND
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB',0,ZHOOK_HANDLE)
!
ILTYPE=-2
IF (PRESENT(KLTYPE)) ILTYPE=KLTYPE
ILEV1=-2
IF (PRESENT(KLEV1)) ILEV1=KLEV1
ILEV2=-2
IF (PRESENT(KLEV2)) ILEV2=KLEV2
!
 CALL MAKE_GRIB_INDEX(HGRIB)
!
IFOUND=0
KRET=0
!
IF (PRESENT(SHORTNAME)) THEN
  CALL GRIB_INDEX_SELECT(NIDX2,'shortName',SHORTNAME,KRET)
  CALL GRIB_NEW_FROM_INDEX(NIDX2,IGRIB,KRET)

  IF (.NOT.ASSOCIATED(PFIELD)) THEN
     CALL GRIB_GET_SIZE(IGRIB,'values',ISIZE,KRET)
     IF (KRET.NE.0) CALL ABOR1_SFX("MODE_READ_GRIB:READ_GRIB: Problem getting size of values")
     ALLOCATE(PFIELD(ISIZE))
  ENDIF
  !
  CALL GRIB_GET(IGRIB,'values',PFIELD,KRET)
  IF (KRET.NE.0) CALL ABOR1_SFX("MODE_READ_GRIB:READ_GRIB: Problem getting values")
  CALL GRIB_RELEASE(IGRIB,KRET)
  IF (KRET.NE.0) CALL ABOR1_SFX("MODE_READ_GRIB:READ_GRIB: Problem releasing memory")
   !
ELSE
  CALL GRIB_INDEX_SELECT(NIDX,'paramId',KPARAM,KRET)
  CALL GRIB_NEW_FROM_INDEX(NIDX,IGRIB,KRET)
  IF (KRET.EQ.0) CALL GET_GRIB_MESSAGE(KLUOUT,ILTYPE,ILEV1,ILEV2,IGRIB,IFOUND)
  !
  IF (PRESENT(KPARAM2)) THEN
    IF (IFOUND/=3) THEN
      CALL GRIB_INDEX_SELECT(NIDX,'paramId',KPARAM2,KRET)
      CALL GRIB_NEW_FROM_INDEX(NIDX,IGRIB,KRET)
      IF (KRET.EQ.0) THEN
        ILTYPE=-2
        IF (PRESENT(KLTYPE)) ILTYPE=KLTYPE
        CALL GET_GRIB_MESSAGE(KLUOUT,ILTYPE,ILEV1,ILEV2,IGRIB,IFOUND)
      ENDIF
    ELSE
      KPARAM2 = KPARAM
    ENDIF
  ENDIF
  !
  IF (IFOUND==3) THEN
    !
    IF (PRESENT(KLTYPE)) KLTYPE = ILTYPE
    IF (PRESENT(KLEV1))  KLEV1  = ILEV1
    IF (PRESENT(KLEV2))  KLEV2  = ILEV2
    !
    IF (.NOT.ASSOCIATED(PFIELD)) THEN
      CALL GRIB_GET_SIZE(IGRIB,'values',ISIZE,KRET)
      IF (KRET.NE.0) CALL ABOR1_SFX("MODE_READ_GRIB:READ_GRIB: Problem getting size of values")
      ALLOCATE(PFIELD(ISIZE))
    ENDIF
    !
  
    IF (PRESENT(PLATS).AND.PRESENT(PLONS)) THEN
        ALLOCATE( PLATS(ISIZE),PLONS(ISIZE))
        CALL GRIB_GET_DATA(IGRIB,PLATS,PLONS,PFIELD,KRET)
        IF (KRET /= 0) CALL ABOR1_SFX("MODE_READ_GRIB:READ_GRIB: Problem getting grid or field")
    ELSE
        CALL GRIB_GET(IGRIB,'values',PFIELD,KRET)
        IF (KRET.NE.0) CALL ABOR1_SFX("MODE_READ_GRIB:READ_GRIB: Problem getting values")
    ENDIF
  
    CALL GRIB_RELEASE(IGRIB,KRET)
    IF (KRET.NE.0) CALL ABOR1_SFX("MODE_READ_GRIB:READ_GRIB: Problem releasing memory")
    !
  ELSE
    !
    KRET = 1
    !
  ENDIF
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB
!
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!     ####################
      SUBROUTINE READ_GRIB_LAND_MASK(HGRIB,KLUOUT,HINMODEL,PMASK,PLFRAC)
!     ####################
!
USE MODD_PREP_ISBA, ONLY : XRM_WM_ECMWF        
!
IMPLICIT NONE
!
CHARACTER(LEN=*),   INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,            INTENT(IN)    :: KLUOUT    ! logical unit of output listing
CHARACTER(LEN=6),   INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), POINTER       :: PMASK     ! Land mask
REAL, DIMENSION(:), POINTER, OPTIONAL       :: PLFRAC    ! Land fraction
!
INTEGER(KIND=kindOfInt)                 :: IRET      ! return code
INTEGER                           :: ILTYPE    ! leveltype
INTEGER                           :: ILEV      ! level
REAL                              :: ZFR_SMALL, ZFR_LARGE
INTEGER                           :: IPOINT    ! horizontal grid-box number
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_LAND_MASK',0,ZHOOK_HANDLE)
WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_LAND_MASK: | Reading land mask from ',HINMODEL
!
!

SELECT CASE (HINMODEL)
  CASE ('ECMWF ','RACMO ')
    CALL READ_GRIB(HGRIB,KLUOUT,172,IRET,PMASK)
  CASE ('ARPEGE','ALADIN','MOCAGE')
    CALL READ_GRIB(HGRIB,KLUOUT,81,IRET,PMASK)          
  CASE ('HIRLAM')        
    ILTYPE=105
    ILEV  =0    
    CALL READ_GRIB(HGRIB,KLUOUT,81,IRET,PMASK,KLTYPE=ILTYPE,KLEV1=ILEV)            
  CASE DEFAULT
    CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_LAND_MASK: OPTION NOT SUPPORTED '//HINMODEL)
END SELECT
!
IF (IRET /= 0) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB: LAND SEA MASK MISSING (READ_GRIB_LAND_MASK)')
END IF

IF (PRESENT(PLFRAC)) THEN
  IF (HINMODEL == 'ECMWF ') THEN
    ALLOCATE(PLFRAC(SIZE(PMASK)))
    PLFRAC=PMASK  
    IF (XRM_WM_ECMWF /= 0.0 ) THEN
      ZFR_SMALL=XRM_WM_ECMWF
      ZFR_LARGE=1.-ZFR_SMALL
      DO IPOINT=1,SIZE(PLFRAC)
        IF (PLFRAC(IPOINT).GT.ZFR_LARGE) THEN
          PLFRAC(IPOINT) = 1.
        END IF
        IF (PLFRAC(IPOINT).LT.ZFR_SMALL) THEN
          PLFRAC(IPOINT) = 0.
        END IF       
      END DO
    END IF
  ELSE
    CALL ABOR1_SFX('MODE_READ_GRIB: LAND FRACTION IS IMPLEMENTED ONLY FOR ECMWF INIT. DATA (READ_GRIB_LAND_MASK)') 
  END IF
END IF
!
WHERE (PMASK>0.5)
  PMASK = 1.
ELSEWHERE
  PMASK = 0.
END WHERE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_LAND_MASK',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_LAND_MASK
!-------------------------------------------------------------------
!     ############################
      SUBROUTINE READ_GRIB_ZS(HGRIB,KLUOUT,HINMODEL,PZS)
!     ############################
!
USE MODD_CSTS,       ONLY : XG
!
IMPLICIT NONE
!
 CHARACTER(LEN=*),   INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,            INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),   INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), POINTER       :: PZS       ! 
!
INTEGER(KIND=kindOfInt)                           :: IRET      ! return code
INTEGER :: ILTYPE,ILEV1
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------
!* Read orography
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_ZS',0,ZHOOK_HANDLE)
WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_ZS: | Reading orography from ',HINMODEL
!
SELECT CASE (HINMODEL)
  CASE ('ECMWF ','RACMO ')
    CALL READ_GRIB(HGRIB,KLUOUT,129,IRET,PZS) 
    IF ( IRET /= 0 ) THEN
      WRITE(KLUOUT,'(A)')' | 129 0 105 not found, try z 1 109 instead'
      ILTYPE=109
      ILEV1=1
      CALL READ_GRIB(HGRIB,KLUOUT,129,IRET,PZS, &
           & KLTYPE=ILTYPE,KLEV1=ILEV1,SHORTNAME='z') 
    ENDIF
  CASE ('ARPEGE','MOCAGE')
    CALL READ_GRIB(HGRIB,KLUOUT,8,IRET,PZS)               
  CASE ('HIRLAM','ALADIN')
    CALL READ_GRIB(HGRIB,KLUOUT,6,IRET,PZS)                  
  CASE DEFAULT
    CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_ZS:OPTION NOT SUPPORTED '//HINMODEL)
END SELECT
!
IF (IRET /= 0) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB: OROGRAPHY MISSING (READ_GRIB_ZS_LAND)')
END IF
!
! Datas given in archives are multiplied by the gravity acceleration
PZS = PZS / XG
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_ZS',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_ZS
!-------------------------------------------------------------------
!     ############################
      SUBROUTINE READ_GRIB_ZS_LAND(HGRIB,KLUOUT,HINMODEL,PMASK,PZSL)
!     ############################
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
 CHARACTER(LEN=*),   INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,            INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),   INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:), POINTER       :: PZSL      ! 
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_ZS_LAND',0,ZHOOK_HANDLE)
!
 CALL READ_GRIB_ZS(HGRIB,KLUOUT,HINMODEL,PZSL)
!
IF (SIZE(PMASK)==SIZE(PZSL)) &
  WHERE (PMASK(:)/=1.) PZSL = XUNDEF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_ZS_LAND',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_ZS_LAND
!-------------------------------------------------------------------
!     ############################
      SUBROUTINE READ_GRIB_ZS_SEA(HGRIB,KLUOUT,HINMODEL,PMASK,PZSS)
!     ############################
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
 CHARACTER(LEN=*),   INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,            INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),   INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:), POINTER       :: PZSS      ! 
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_ZS_SEA',0,ZHOOK_HANDLE)
!
IF (HINMODEL=='HIRLAM') THEN
  CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_ZS_SEA:OPTION NOT SUPPORTED '//HINMODEL)
ELSE
  CALL READ_GRIB_ZS(HGRIB,KLUOUT,HINMODEL,PZSS)
ENDIF
!
IF (SIZE(PMASK)==SIZE(PZSS)) &
  WHERE (PMASK(:)/=0.) PZSS = XUNDEF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_ZS_SEA',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_ZS_SEA
!-------------------------------------------------------------------
!     ###########################
      SUBROUTINE READ_GRIB_T(HGRIB,KLUOUT,HINMODEL,PT)
!     ###########################
!
IMPLICIT NONE
!
 CHARACTER(LEN=*),   INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,            INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),   INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), POINTER       :: PT        ! 
!
INTEGER(KIND=kindOfInt)                           :: IRET      ! return code
INTEGER                           :: ILTYPE    ! type of level (Grib code table 3)
INTEGER                           :: ILEV      ! level definition
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------
!* Read surface temperature
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_T',0,ZHOOK_HANDLE)
WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_T: | Reading surface temperature'
!
SELECT CASE (HINMODEL)
  CASE ('ECMWF ','RACMO ')
    WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_T: | Reading upper soil layer temperature'
    CALL READ_GRIB(HGRIB,KLUOUT,139,IRET,PT)

  CASE ('ARPEGE','ALADIN','MOCAGE')
    ILEV=0
    ILTYPE=111
    CALL READ_GRIB(HGRIB,KLUOUT,11,IRET,PT,KLTYPE=ILTYPE,KLEV1=ILEV) 
    IF (IRET /= 0) THEN
      ILTYPE=1
      CALL READ_GRIB(HGRIB,KLUOUT,11,IRET,PT,KLTYPE=ILTYPE) 
      IF (IRET /= 0) THEN
         ILTYPE=105
        CALL READ_GRIB(HGRIB,KLUOUT,11,IRET,PT,KLTYPE=ILTYPE,KLEV1=ILEV)   
      ENDIF        
    END IF

  CASE ('HIRLAM ')
    WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_T: | Reading surface temperature tile 4'
     ILTYPE=105
     ILEV=904
    CALL READ_GRIB(HGRIB,KLUOUT,11,IRET,PT,KLTYPE=ILTYPE,KLEV1=ILEV)

  CASE DEFAULT
    CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_T:OPTION NOT SUPPORTED '//HINMODEL)
END SELECT
!
IF (IRET /= 0) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB: SURFACE TEMPERATURE MISSING (READ_GRIB_T)')
END IF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_T',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_T
!-------------------------------------------------------------------
!     ###########################
      SUBROUTINE READ_GRIB_TS(HGRIB,KLUOUT,HINMODEL,PMASK,PTS,PLFRAC)
!     ###########################
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_PREP_ISBA,  ONLY : XRM_WM_ECMWF
!
IMPLICIT NONE
!
 CHARACTER(LEN=*),   INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,            INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),   INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:), POINTER       :: PTS       ! 
REAL, DIMENSION(:),   INTENT(IN), OPTIONAL :: PLFRAC ! grib land fraction
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TS',0,ZHOOK_HANDLE)
!
 CALL READ_GRIB_T(HGRIB,KLUOUT,HINMODEL,PTS)
IF (PRESENT(PLFRAC)) THEN
  IF (SIZE(PLFRAC)==SIZE(PTS)) THEN
    WHERE (PLFRAC(:).LT.XRM_WM_ECMWF) PTS(:) = XUNDEF
  ELSE
    CALL ABOR1_SFX('MODE_READ_GRIB: TS AND LAND FRACTION DIFFER (READ_GRIB_TS)')
  END IF       
ELSE
  IF (SIZE(PMASK)==SIZE(PTS)) &
    WHERE (PMASK(:)/=1.) PTS = XUNDEF
END IF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TS',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_TS
!-------------------------------------------------------------------
!     ###########################
      SUBROUTINE READ_GRIB_SST(HGRIB,KLUOUT,HINMODEL,PMASK,PSST)
!     ###########################
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
 CHARACTER(LEN=*),   INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,            INTENT(IN)     :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),   INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), INTENT(IN)     :: PMASK     ! grib land mask
REAL, DIMENSION(:), POINTER        :: PSST      ! 
!
!* local variables
!  ---------------
INTEGER                           :: IL, JL 
INTEGER(KIND=kindOfInt)           :: IRET
REAL, DIMENSION(:), ALLOCATABLE   :: ZTS       ! surface temperature
!
REAL  :: ZMIN, ZMAX  ! Min and max values
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_SST',0,ZHOOK_HANDLE)
WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_SST: | Reading sea surface temperature from ',HINMODEL
!
SELECT CASE (HINMODEL)
  CASE ('HIRLAM')
    CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_SST:OPTION NOT SUPPORTED '//HINMODEL)
  CASE ('ECMWF ')
!
  CALL READ_GRIB(HGRIB,KLUOUT,235,IRET,PSST)
  IF (IRET /= 0) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB_SST: SURFACE TEMPERATURE (235) MISSING ')
  END IF
  ALLOCATE(ZTS(SIZE(PSST)))
  ZTS = PSST
!
  CALL READ_GRIB(HGRIB,KLUOUT,34,IRET,PSST)
  IF (IRET /= 0) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB_SST: SEA SURFACE TEMPERATURE MISSING ')
  END IF
!
! Calculates min and max values of ZTS and PSST 
  ZMIN = 500.
  ZMAX = 0. 
  DO JL=1,SIZE(PSST)
    ZMIN = MIN(ZMIN,  PSST(JL))
    ZMAX = MAX(ZMAX,  PSST(JL))
  ENDDO
  WRITE(KLUOUT,*) ' mode_read_grib  PSSTmin =', ZMIN,' PSSTmax =', ZMAX
!
  ZMIN = 500.
  ZMAX = 0. 
  DO JL=1,SIZE(ZTS)
     ZMIN = MIN(ZMIN,  ZTS(JL))
     ZMAX = MAX(ZMAX,  ZTS(JL))
  ENDDO
  WRITE(KLUOUT,*) ' mode_read_grib   ZTSmin =', ZMIN,' ZTSmax =', ZMAX
!
! Replace SST freezing temperature over sea ice and missing values over land with surface skin temperatures
  WHERE ( PSST(:).LT.272. .OR. PSST(:).GT.9990. ) 
    PSST(:) = ZTS(:)
  ENDWHERE
  DEALLOCATE(ZTS)
!
  CASE ('ARPEGE','ALADIN','MOCAGE')
    CALL READ_GRIB_T(HGRIB,KLUOUT,HINMODEL,PSST)
  CASE ('RACMO ')
    CALL READ_GRIB(HGRIB,KLUOUT,34,IRET,PSST)
    IF (IRET /= 0) THEN
      CALL ABOR1_SFX('MODE_READ_GRIB_SST: SEA SURFACE TEMPERATURE (11-102-0) MISSING ')
    END IF
  CASE DEFAULT
     CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_SST:OPTION NOT SUPPORTED '//HINMODEL)    
!
END SELECT
!
IF (SIZE(PMASK)==SIZE(PSST)) WHERE (PMASK(:)/=0.) PSST = XUNDEF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_SST',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_SST
!-------------------------------------------------------------------
!     ###########################
      SUBROUTINE READ_GRIB_SIC(HGRIB,KLUOUT,HINMODEL,PMASK,PSIC)
!     ###########################
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODI_ABOR1_SFX

USE NN_EXTRAPOLATE
!
IMPLICIT NONE
!
CHARACTER(LEN=*),     INTENT(IN)  :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)  :: KLUOUT    ! logical unit of output listing
CHARACTER(LEN=6),     INTENT(IN)  :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)  :: PMASK     ! grib land mask
REAL, DIMENSION(:),   POINTER     :: PSIC      !
!
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)           :: IRET
REAL, POINTER                     :: ZLATS_RAW(:),  &
                                     ZLONS_RAW(:)
!
REAL  :: ZMIN, ZMAX  ! Min and max values
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_SIC',0,ZHOOK_HANDLE)
WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_SIC: | Reading sea ice fraction from ',HINMODEL

IF (HINMODEL=='ECMWF' .OR. HINMODEL=='ALADIN' .OR. HINMODEL=='RACMO ') THEN
  CALL READ_GRIB(HGRIB,KLUOUT,31,IRET,PSIC, PLATS = ZLATS_RAW, PLONS = ZLONS_RAW)
  IF (IRET /= 0) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB_SIC: SEA ICE CONCENTRATION MISSING ')
  END IF
ELSE
  CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_SIC:OPTION NOT SUPPORTED '//HINMODEL)
ENDIF

IF( ANY( ABS(PSIC)>1.1 ) ) THEN
  CALL NEAREST_NEIGHBOUR_EXTRAPOLATE( PSIC, ZLONS_RAW, ZLATS_RAW, PSIC >= 0. .AND. PSIC <= 1. )
END IF
IF (SIZE(PMASK)==SIZE(PSIC)) WHERE (PMASK(:)/=0.) PSIC = XUNDEF

IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_SIC',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_SIC
!-------------------------------------------------------------------
!     ###########################
      SUBROUTINE READ_GRIB_TSWATER(HGRIB,KLUOUT,HINMODEL,PMASK,PTS)
!     ###########################
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
 CHARACTER(LEN=*),   INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,            INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),   INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:), POINTER       :: PTS     ! 
!
INTEGER :: IRET
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TSWATER',0,ZHOOK_HANDLE)
WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_TSWATER: | Reading water temperature from ',HINMODEL
!
SELECT CASE (HINMODEL)
  CASE ('ECMWF ')
     CALL READ_GRIB(HGRIB,KLUOUT,3080,IRET,PTS)
     IF (IRET /= 0) CALL READ_GRIB_T2(HGRIB,KLUOUT,HINMODEL,PMASK,PTS)
  CASE ('ARPEGE','ALADIN','MOCAGE','HIRLAM','RACMO ')
    CALL READ_GRIB_T2(HGRIB,KLUOUT,HINMODEL,PMASK,PTS)
   CASE DEFAULT
     CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_TSWATER:OPTION NOT SUPPORTED '//HINMODEL)    
END SELECT
!
IF (SIZE(PMASK)==SIZE(PTS)) WHERE (PMASK(:)/=0.) PTS = XUNDEF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TSWATER',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_TSWATER
!-------------------------------------------------------------------
!     ###########################
      SUBROUTINE READ_GRIB_T2(HGRIB,KLUOUT,HINMODEL,PMASK,PT2)
!     ###########################
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
 CHARACTER(LEN=*),   INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,            INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),   INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:), POINTER       :: PT2       ! 
!
INTEGER(KIND=kindOfInt)                           :: IRET
INTEGER                           :: ILTYPE, ILEV    ! type of level (Grib code table 3)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------
!* Read deep soil temperature
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_T2',0,ZHOOK_HANDLE)
WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_T2: | Reading deep soil temperature'
!
SELECT CASE (HINMODEL)
  CASE ('ECMWF ')
     CALL READ_GRIB(HGRIB,KLUOUT,170,IRET,PT2)   
  CASE ('ARPEGE','ALADIN','MOCAGE')
     ILTYPE=111
     CALL READ_GRIB(HGRIB,KLUOUT,11,IRET,PT2,KLTYPE=ILTYPE) 
    IF (IRET /= 0) THEN
       ILTYPE=105
      CALL READ_GRIB(HGRIB,KLUOUT,11,IRET,PT2,KLTYPE=ILTYPE)   
    ENDIF  
  CASE ('HIRLAM ')
     ILTYPE=105
     ILEV=954
     CALL READ_GRIB(HGRIB,KLUOUT,11,IRET,PT2,KLTYPE=ILTYPE,KLEV1=ILEV)
  CASE ('RACMO ')
     ! Read third soil layer
     CALL READ_GRIB(HGRIB,KLUOUT,183,IRET,PT2)
  CASE DEFAULT
     CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_T2:OPTION NOT SUPPORTED '//HINMODEL)
END SELECT
!
IF (IRET /= 0) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB: DEEP SOIL TEMPERATURE MISSING (READ_GRIB_T2)')
END IF
!
!IF (SIZE(PMASK)==SIZE(PT2)) &
!  WHERE (PMASK(:)/=1.) PT2 = XUNDEF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_T2',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_T2
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!     ###########################
      SUBROUTINE READ_GRIB_T2_LAND(HGRIB,KLUOUT,HINMODEL,PMASK,ZFIELD)
!     ###########################
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
 CHARACTER(LEN=*),   INTENT(IN)   :: HGRIB     ! Grib file name
INTEGER,            INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),   INTENT(IN)   :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:), POINTER       :: ZFIELD    ! 
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------
!* Read deep soil temperature
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_T2_LAND',0,ZHOOK_HANDLE)
!
 CALL READ_GRIB_T2(HGRIB,KLUOUT,HINMODEL,PMASK,ZFIELD)
!
IF (SIZE(PMASK)==SIZE(ZFIELD)) WHERE (PMASK(:)/=1.) ZFIELD = XUNDEF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_T2_LAND',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_T2_LAND
!-------------------------------------------------------------------
!-------------------------------------------------------------------
SUBROUTINE PUT_LAYER_DEPTH(KLUOUT,KLEV,HROUT,KLTYPE,KLEV1,KLEV2,KNLAYERDEEP,PV4,PV,PD)
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KLUOUT
INTEGER, INTENT(IN) :: KLEV
 CHARACTER(LEN=*), INTENT(IN) :: HROUT
INTEGER, INTENT(INOUT) :: KLTYPE
INTEGER, INTENT(IN) :: KLEV1
INTEGER, INTENT(IN) :: KLEV2
INTEGER, INTENT(IN) :: KNLAYERDEEP
REAL, INTENT(IN) :: PV4
REAL, INTENT(IN) :: PV
REAL, INTENT(OUT) :: PD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:PUT_LAYER_DEPTH',0,ZHOOK_HANDLE)
!
IF (KLEV2 == -1) KLTYPE = 0
IF (KLTYPE==112) THEN
  PD = (KLEV2 - KLEV1) / 100.
ELSE
  IF (KNLAYERDEEP == 4) THEN
    PD = PV4
  ELSE
    PD = PV
  END IF
  WRITE (KLUOUT,'(A,i1,A,f5.2,A)') 'MODE_READ_GRIB:'//TRIM(HROUT)//': | Level ',&
                                    KLEV,' height not available, assuming ',PD,' m'
END IF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:PUT_LAYER_DEPTH',1,ZHOOK_HANDLE)
END SUBROUTINE PUT_LAYER_DEPTH
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE FILL_PFIELD(KLUOUT,HROUT,KNLAYERDEEP,PDIN,PFIELDIN,PMASK,PFIELDOUT,PDOUT,PLFRAC,HVAR)
!     #######################
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_PREP_ISBA,  ONLY : XRM_WM_ECMWF
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KLUOUT
 CHARACTER(LEN=*), INTENT(IN) :: HROUT
INTEGER, INTENT(IN) :: KNLAYERDEEP
REAL, DIMENSION(:), INTENT(IN) :: PDIN
REAL, DIMENSION(:,:), INTENT(IN) :: PFIELDIN
REAL, DIMENSION(:), INTENT(IN) :: PMASK
REAL, DIMENSION(:,:), POINTER :: PFIELDOUT
REAL, DIMENSION(:,:), POINTER :: PDOUT
REAL, DIMENSION(:),   INTENT(IN), OPTIONAL :: PLFRAC ! grib land fraction
CHARACTER(3), INTENT(IN), OPTIONAL :: HVAR ! which field to fill
!
 CHARACTER(LEN=20) :: FMT0
INTEGER :: JL, JP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:FILL_PFIELD',0,ZHOOK_HANDLE)
!--------------------------------------------------------------------------------
! 1.  Display the number of layer found
!     -----------------------
WRITE(FMT0,FMT='(A8,I1,A11)') '(A,I1,A,',KNLAYERDEEP,'(F5.2,","))'
WRITE (KLUOUT,FMT=FMT0) 'MODE_READ_GRIB:'//TRIM(HROUT)//': | ',KNLAYERDEEP,&
                        ' deep layers, heights are : ',PDIN(1:KNLAYERDEEP)
!--------------------------------------------------------------------------------
! 2.  Set temperature profile and layer thicknesses
!     -----------------------------------------------

ALLOCATE(PFIELDOUT(SIZE(PFIELDIN,1),KNLAYERDEEP))
ALLOCATE(PDOUT(SIZE(PFIELDIN,1),KNLAYERDEEP))
!
DO JL=1,KNLAYERDEEP
  PDOUT(:,JL)=SUM(PDIN(1:JL))
  PFIELDOUT(:,JL)=PFIELDIN(:,JL)
  IF(PRESENT(PLFRAC)) THEN
     IF (SIZE(PLFRAC)==SIZE(PFIELDOUT,1)) THEN
        IF (PRESENT(HVAR)) THEN
          SELECT CASE(HVAR)
          CASE('TG ')
            WHERE (PLFRAC(:).LT.XRM_WM_ECMWF)
              PFIELDOUT(:,JL) = XUNDEF
              PDOUT(:,JL) = XUNDEF
            END WHERE
          CASE('WG ')
             DO JP=1,SIZE(PFIELDOUT,1)
                IF(PLFRAC(JP).LT.XRM_WM_ECMWF) THEN
                   PFIELDOUT(JP,JL) = XUNDEF
                   PDOUT(JP,JL) = XUNDEF
                ELSE
                   PFIELDOUT(JP,JL) = MAX(PFIELDOUT(JP,JL),0.196) ! <- EK: Removing unrealistically low WG values which come from interpolation.
                                                                  ! Options: 0.059 is the minimum between all soil types wilting point (see READ_GRIB_WG_ECMWF) 
                                                                  !          0.196 is the average between all soil types wilting point
                END IF              
             END DO
          END SELECT
        ELSE
          CALL ABOR1_SFX('MODE_READ_GRIB: WHICH FIELD IS TO TREAT WITH LAND FRACTION? (FILL_PFIELD)')
        END IF
     ELSE
        CALL ABOR1_SFX('MODE_READ_GRIB: FIELD AND LAND FRACTION DIFFER (FILL_PFIELD)')
     END IF
  ELSE
     IF (SIZE(PMASK)==SIZE(PFIELDOUT,1)) &
       WHERE (PMASK(:)/=1.) PFIELDOUT(:,JL) = XUNDEF
  END IF  
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:FILL_PFIELD',1,ZHOOK_HANDLE)
END SUBROUTINE FILL_PFIELD
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_TG_ECMWF(HGRIB,KLUOUT,HINMODEL,PMASK,PTG,PD,PLFRAC)
!     #######################
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PTG       ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PD        ! thickness of each layer
REAL, DIMENSION(:),   INTENT(IN), OPTIONAL :: PLFRAC ! grib land fraction
!
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)                           :: IRET      ! return code
INTEGER                           :: ILTYPE    ! type of level (Grib code table 3)
INTEGER                           :: ILEV1     ! level definition
INTEGER                           :: ILEV2     ! level definition
INTEGER                           :: JL         ! layer loop counter
INTEGER                           :: INLAYERDEEP! number of deep moisture layers
REAL,    DIMENSION(:), POINTER    :: ZFIELD => NULL()  ! first layer temperature
REAL,  DIMENSION(:,:), ALLOCATABLE:: ZTG      ! first layer temperature
REAL, DIMENSION(:)   , ALLOCATABLE:: ZD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TG_ECMWF',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_TG_ECMWF: | Reading soil temperature'
!
ALLOCATE(ZD(5))
ZD=0.
!
! 1.  Search and read level 1 (and its depth)
!     --------------------------------------
ILTYPE= -1
ILEV1 = -1
ILEV2 = -1
 CALL READ_GRIB(HGRIB,KLUOUT,139,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2)
!
IF (IRET== 0) THEN
  CALL PUT_LAYER_DEPTH(KLUOUT,1,'READ_GRIB_TG_ECMWF',ILTYPE,ILEV1,ILEV2,4,0.07,0.07,ZD(1))
  ALLOCATE(ZTG(SIZE(ZFIELD),5))
  ZTG=0.
  ZTG(:,1)=ZFIELD
ELSE
  CALL ABOR1_SFX('MODE_READ_GRIB: SOIL TEMPERATURE LEVEL 1 MISSING (READ_GRIB_TG_ECMWF)')
ENDIF
!
! 2.  Search and read level 4 (and its depth) This level is optionnal
!     ---------------------------------------------------------------
ILTYPE= -1
ILEV1 = -1
ILEV2 = -1
 CALL READ_GRIB(HGRIB,KLUOUT,236,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2)
!
IF (IRET == 0) THEN
  INLAYERDEEP = 4
  CALL PUT_LAYER_DEPTH(KLUOUT,4,'READ_GRIB_TG_ECMWF',ILTYPE,ILEV1,ILEV2,INLAYERDEEP,1.89,1.89,ZD(4))
  ZTG(:,4)=ZFIELD
ELSE
  INLAYERDEEP = 3
  ZD(4) = 0.
ENDIF
!
! 3.  Search and read level 3 (and its depth) This level is optionnal
!     ---------------------------------------------------------------
ILTYPE= -1
ILEV1 = -1
ILEV2 = -1
 CALL READ_GRIB(HGRIB,KLUOUT,183,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2)
!
IF (IRET == 0) THEN      
  CALL PUT_LAYER_DEPTH(KLUOUT,3,'READ_GRIB_TG_ECMWF',ILTYPE,ILEV1,ILEV2,INLAYERDEEP,0.72,0.42,ZD(3))
  ZTG(:,3)=ZFIELD 
ELSE
  INLAYERDEEP = 2
  ZD(3) = 0.        
ENDIF
!
! 4.  Search and read level 2 (and its depth)
!     ---------------------------------------
ILTYPE= -1
ILEV1 = -1
ILEV2 = -1
 CALL READ_GRIB(HGRIB,KLUOUT,170,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2)
!
IF (IRET== 0) THEN         
  CALL PUT_LAYER_DEPTH(KLUOUT,2,'READ_GRIB_TG_ECMWF',ILTYPE,ILEV1,ILEV2,INLAYERDEEP,0.21,0.42,ZD(2))
  ZTG(:,2)=ZFIELD
  DEALLOCATE(ZFIELD)    
ELSE
  CALL ABOR1_SFX('MODE_READ_GRIB: SOIL TEMPERATURE LEVEL 2 MISSING (READ_GRIB_TG_ECMWF)')
ENDIF
!--------------------------------------------------------------------------------
! 5.  Assumes uniform temperature profile up to 3m depth
!     -------------------------------------------------
!
IF(SUM(ZD(1:INLAYERDEEP)) < 3.) THEN
  !We add a temperature layer
  INLAYERDEEP=INLAYERDEEP+1
  ZD(INLAYERDEEP)=3.-SUM(ZD(1:INLAYERDEEP-1))
  ZTG(:,INLAYERDEEP)=ZTG(:,INLAYERDEEP-1)
ENDIF
!
!--------------------------------------------------------------------------------
! 6.  Set temperature profile and layer thicknesses
!     ----------------------------------------------

IF (PRESENT(PLFRAC)) THEN
   CALL FILL_PFIELD(KLUOUT,'READ_GRIB_TG_ECMWF',INLAYERDEEP,ZD,ZTG,PMASK,PTG,PD,PLFRAC,'TG ')
ELSE
   CALL FILL_PFIELD(KLUOUT,'READ_GRIB_TG_ECMWF',INLAYERDEEP,ZD,ZTG,PMASK,PTG,PD)
END IF

DEALLOCATE(ZD)
DEALLOCATE(ZTG)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TG_ECMWF',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_TG_ECMWF
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_WG_ECMWF_1(HGRIB,KLUOUT,HINMODEL,PMASK,PWG,PD,PLFRAC)
!     #######################
!
! This tasks is divided in the following steps :
!  - computing the MesoNH constants
!  - reading the grib datas according to the type of file (ECMWF/Arpege/Aladin)
!  - converting from specific humidity to relative humidity
!  - interpolation with land mask
!  - converting back from relative humidity to specific humidity with MesoNH constants
! Five different models are supported :
!  - ECMWF with 2 layers (untested)
!  - ECMWF with 3 layers (archive before 1991 - Blondin model)
!  - ECMWF with 4 layers (archive after 1991 - Viterbo model)
!  - Arpege/Aladin before ISBA (I don't know the name of this model)
!  - Arpege/Aladin with ISBA model
! The available model is detect according to the fields which are presents :
!  - ECMWF archive : loads as many layers as possible
!  - Arpege/Aladin archive : ISBA model needs Clay and Sans fraction fields, if they
!    are present, they are used and the model is declared to be ISBA.
! To detect the height of the layers, two methods are used :
!  - if level type is not 112, a default value is assumed and a warning message is
!    displayed
!  - if level type is ID 112, then the position of the top and bottom surface may be
!    given. If they are present, they are used, if not the default value is taken and
!    a warning message is issued.
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PWG       ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PD        ! thickness of each layer
REAL, DIMENSION(:),   INTENT(IN), OPTIONAL :: PLFRAC ! grib land fraction
!
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)                           :: IRET      ! return code
INTEGER                           :: IPAR      ! parameter number for field reading
INTEGER                           :: ILTYPE    ! type of level (Grib code table 3)
INTEGER                           :: ILEV1     ! level definition
INTEGER                           :: ILEV2     ! level definition
INTEGER                           :: INLAYERDEEP! number of deep moisture layers
REAL,    DIMENSION(:), POINTER    :: ZFIELD => NULL()  ! first layer temperature
REAL,  DIMENSION(:,:), ALLOCATABLE:: ZWG      ! first layer temperature
REAL, DIMENSION(:)   , ALLOCATABLE:: ZD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WG_ECMWF_1',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_WG_ECMWF_1: | Reading soil moisture'
!
ALLOCATE(ZD(4))
ZD=0.
!
! 1.  Search and read level 1 (and its depth)
!     --------------------------------------
ILTYPE= -1
ILEV1 = -1
ILEV2 = -1
IPAR=39
 CALL READ_GRIB(HGRIB,KLUOUT,140,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2,KPARAM2=IPAR)
!
IF (IRET == 0) THEN 
  CALL PUT_LAYER_DEPTH(KLUOUT,1,'READ_GRIB_WG_ECMWF_1',ILTYPE,ILEV1,ILEV2,4,0.07,0.07,ZD(1))
  ALLOCATE(ZWG(SIZE(ZFIELD,1),4))
  ZWG=0.
  ZWG(:,1)=ZFIELD
  !
  IF (IPAR==140) ZWG(:,1)=ZWG(:,1) / ZD(1)
ELSE
  CALL ABOR1_SFX('MODE_READ_GRIB: SOIL MOISTURE LEVEL 1 MISSING (READ_GRIB_WG_ECMWF_1)')
ENDIF
!
! 2.  Search and read level 4 (and its depth) This level is optionnal
!     ---------------------------------------------------------------
ILTYPE= -1
ILEV1 = -1
ILEV2 = -1
IPAR=42
 CALL READ_GRIB(HGRIB,KLUOUT,237,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2,KPARAM2=IPAR)
!
IF (IRET == 0) THEN   
  INLAYERDEEP = 4
  CALL PUT_LAYER_DEPTH(KLUOUT,4,'READ_GRIB_WG_ECMWF_1',ILTYPE,ILEV1,ILEV2,INLAYERDEEP,1.89,1.89,ZD(4))
  ZWG(:,4)=ZFIELD   
  !
  IF (IPAR==237) ZWG(:,4)=ZWG(:,4) / ZD(1) 
ELSE
  INLAYERDEEP = 3  
  ZD(4) = 0.
ENDIF
!
! 3.  Search and read level 3 (and its depth) This level is optionnal
!     ---------------------------------------------------------------
ILTYPE= -1
ILEV1 = -1
ILEV2 = -1
IPAR=41
 CALL READ_GRIB(HGRIB,KLUOUT,184,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2,KPARAM2=IPAR)
!
IF (IRET == 0) THEN
  CALL PUT_LAYER_DEPTH(KLUOUT,3,'READ_GRIB_WG_ECMWF_1',ILTYPE,ILEV1,ILEV2,INLAYERDEEP,0.72,0.42,ZD(3))
  ZWG(:,3)=ZFIELD 
  !
  IF (IPAR==184) ZWG(:,3)=ZWG(:,3) / ZD(1)
ELSE
  INLAYERDEEP = 2  
  ZD(3) = 0.  
ENDIF
!
! 4.  Search and read level 2 (and its depth)
!     ---------------------------------------
ILTYPE= -1
ILEV1 = -1
ILEV2 = -1
IPAR=40
 CALL READ_GRIB(HGRIB,KLUOUT,171,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2,KPARAM2=IPAR)
!
IF (IRET == 0) THEN
  CALL PUT_LAYER_DEPTH(KLUOUT,2,'READ_GRIB_WG_ECMWF_1',ILTYPE,ILEV1,ILEV2,INLAYERDEEP,0.21,0.42,ZD(2))
  ZWG(:,2)=ZFIELD
  DEALLOCATE(ZFIELD)  
  !
  IF (IPAR==171) ZWG(:,2)=ZWG(:,2) / ZD(1)
ELSE
  CALL ABOR1_SFX('MODE_READ_GRIB: SOIL MOISTURE LEVEL 2 MISSING (READ_GRIB_WG_ECMWF_1)')
ENDIF
!
!--------------------------------------------------------------------------------
!
IF(PRESENT(PLFRAC)) THEN
   CALL FILL_PFIELD(KLUOUT,'READ_GRIB_WG_ECMWF_1',INLAYERDEEP,ZD,ZWG,PMASK,PWG,PD,PLFRAC,'WG ')
ELSE
   CALL FILL_PFIELD(KLUOUT,'READ_GRIB_WG_ECMWF_1',INLAYERDEEP,ZD,ZWG,PMASK,PWG,PD)
END IF
DEALLOCATE(ZD)
DEALLOCATE(ZWG)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WG_ECMWF_1',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_WG_ECMWF_1
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE ECMWF_WGI(PTG,PWG,PWGI,PLFRAC)
!     #######################
!
! ECMWF grib only contain (ice+water) content.
! This routine computes iced part and water part according to the formula
! given in ECMWF documentation. But we use real water content instead of
! (CL+CH) times saturation capacity.
!
USE MODD_CSTS,        ONLY : XTT, XPI
USE MODD_SURF_PAR,    ONLY : XUNDEF
USE MODD_PREP_ISBA,   ONLY : XRM_WM_ECMWF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
REAL, DIMENSION(:,:), INTENT(IN)    :: PTG       ! Temperature profil
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWG       ! INPUT contains (water+ice) profil
                                                 ! OUTPUT contains only water profil
REAL, DIMENSION(:,:), INTENT(OUT)   :: PWGI      ! ice profil
REAL, DIMENSION(:),   INTENT(IN), OPTIONAL :: PLFRAC ! grib land fraction
!
!* local variables
!  ---------------
REAL  :: ZT1, ZT2  ! Temperature threshold
INTEGER :: JP,JL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:ECMWF_WGI',0,ZHOOK_HANDLE)
!
ZT1=XTT + 1.
ZT2=XTT - 3.

!
IF (PRESENT(PLFRAC)) THEN
   DO JL=1,SIZE(PWG,2)
      DO JP=1,SIZE(PWG,1)
         IF(PLFRAC(JP).LT.XRM_WM_ECMWF) THEN
            PWGI(JP,JL) = XUNDEF
         ELSE
            IF(PTG(JP,JL).GT.ZT1) THEN
               PWGI(JP,JL) = 0.
            ELSE
               IF(PTG(JP,JL).LT.ZT2) THEN
                  PWGI(JP,JL) = PWG(JP,JL)
                  PWG(JP,JL) = 0.
               ELSE
                  PWGI(JP,JL)=PWG(JP,JL) * 0.5* (1 - sin(XPI * (PTG(JP,JL) - 0.5*ZT1 - 0.5*ZT2) / &
                                                               (ZT1 - ZT2                   )   ))
                  PWG(JP,JL) = PWG(JP,JL) - PWGI(JP,JL)
               END IF
            END IF            
         END IF         
      END DO
   END DO
ELSE
   WHERE(PTG(:,:) > ZT1)
     PWGI(:,:) = 0.
   ELSEWHERE(PTG(:,:) < ZT2)
     PWGI(:,:) = PWG(:,:)
     PWG(:,:) = 0.
   ELSEWHERE
     PWGI(:,:)=PWG(:,:) * 0.5* (1 - sin(XPI * (PTG(:,:) - 0.5*ZT1 - 0.5*ZT2) / &
                                              (ZT1 - ZT2                   )   ))
     PWG(:,:) = PWG(:,:) - PWGI(:,:)
   ENDWHERE
END IF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:ECMWF_WGI',1,ZHOOK_HANDLE)
END SUBROUTINE ECMWF_WGI
!--------------------------------------------------------------------------------
!     #######################
      SUBROUTINE HARMONIZE_GRIB_WG_WGI_ECMWF(HGRIB,KLUOUT,HINMODEL,PMASK,PLFRAC,PWG,PD,PWGI)
!     #######################
!
! ECMWF/RACMO grib only contain (ice+water) content.
! This routine computes iced part and water part according to the formula
! given in ECMWF documentation. But we use real water content instead of
! (CL+CH) times saturation capacity.
!
USE MODD_CSTS,        ONLY : XTT, XPI
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:),   INTENT(IN), OPTIONAL :: PLFRAC ! grib land fraction 
REAL, DIMENSION(:,:), OPTIONAL, POINTER :: PWG   ! INPUT contains (water+ice) profil
                                                 ! OUTPUT contains only water profil
REAL, DIMENSION(:,:), OPTIONAL, POINTER :: PD    ! thickness of each layer                          
REAL, DIMENSION(:,:), OPTIONAL, POINTER :: PWGI  ! ice profil
!* local variables
!  ---------------
REAL,  DIMENSION(:,:), POINTER :: ZWG => NULL()          ! profile of soil water contents
REAL,  DIMENSION(:,:), POINTER :: ZD => NULL()           ! thickness of each layer
REAL,  DIMENSION(:,:), POINTER :: ZTG => NULL()          ! profile of temperature
REAL,  DIMENSION(:,:), POINTER :: ZDT => NULL()          ! thickness of each temperature layer
REAL,  DIMENSION(:,:), ALLOCATABLE:: ZWGI      ! profile of soil ice contents
REAL  :: ZT1, ZT2  ! Temperature threshold
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:HARMONIZE_GRIB_WG_WGI_ECMWF',0,ZHOOK_HANDLE)
!
SELECT CASE (TRIM(HINMODEL))
  CASE ('ECMWF')
     IF(PRESENT(PLFRAC)) THEN
       CALL READ_GRIB_TG_ECMWF(HGRIB,KLUOUT,HINMODEL,PMASK,ZTG,ZDT,PLFRAC=PLFRAC)
       CALL READ_GRIB_WG_ECMWF_1(HGRIB,KLUOUT,HINMODEL,PMASK,ZWG,ZD,PLFRAC=PLFRAC)
     ELSE
       CALL READ_GRIB_TG_ECMWF(HGRIB,KLUOUT,HINMODEL,PMASK,ZTG,ZDT)
       CALL READ_GRIB_WG_ECMWF_1(HGRIB,KLUOUT,HINMODEL,PMASK,ZWG,ZD)
     END IF
  CASE ('RACMO')
    CALL READ_GRIB_TG_RACMO(HGRIB,KLUOUT,HINMODEL,PMASK,ZTG,ZDT)
    CALL READ_GRIB_WG_RACMO_1(HGRIB,KLUOUT,HINMODEL,PMASK,ZWG,ZD)
  CASE DEFAULT
    CALL ABOR1_SFX("MODE_READ_GRIB:HARMONIZE_GRIB_WG_WGI_ECMWF: UNSUPPORTED MODEL " // TRIM(HINMODEL))
END SELECT
!
IF (SIZE(ZTG,2) .LT. SIZE(ZWG,2)) THEN
  WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:HARMONIZE_GRIB_WG_WGI_ECMWF: '
  WRITE  (KLUOUT,'(A)') 'ERROR, YOU HAVE NOT THE SAME NUMBER OF LEVELS '
  WRITE  (KLUOUT,'(A)') 'IN SOIL FOR TEMPERATURE AND HUMIDITY '
  WRITE  (KLUOUT,'(A)') 'VERIFY GRIB FILE '
  CALL ABOR1_SFX("MODE_READ_GRIB:HARMONIZE_GRIB_WG_WGI_ECMWF: VERIFY NUMBER OF LEVELS IN GRIB FILE")
ENDIF
!
IF (PRESENT(PD)) THEN
  ALLOCATE(PD(SIZE(ZD,1),SIZE(ZD,2)))
  PD(:,:)=ZD(:,:)
ENDIF
IF (PRESENT(PWGI)) THEN
  ALLOCATE(PWGI(SIZE(ZWG,1),SIZE(ZWG,2)))
  PWGI(:,:)=0.
ENDIF
!
!If same vertical grids are taken into account for WG and TG we can
!compute ice content and new water content
IF(ALL(ZDT(:,1:SIZE(ZWG,2))==ZD(:,1:SIZE(ZWG,2)))) THEN     
  ALLOCATE(ZWGI(SIZE(ZWG,1),SIZE(ZWG,2)))
  IF(PRESENT(PLFRAC)) THEN
    CALL ECMWF_WGI(ZTG(:,1:SIZE(ZWG,2)),ZWG,ZWGI,PLFRAC=PLFRAC)
  ELSE
    CALL ECMWF_WGI(ZTG(:,1:SIZE(ZWG,2)),ZWG,ZWGI)
  END IF
  IF (PRESENT(PWGI)) PWGI(:,:)=ZWGI(:,:)
  DEALLOCATE(ZWGI)
ENDIF
!
IF (PRESENT(PWG)) THEN
  ALLOCATE(PWG(SIZE(ZWG,1),SIZE(ZWG,2)))
  PWG(:,:)=ZWG(:,:)
ENDIF
!
DEALLOCATE(ZWG)
DEALLOCATE(ZD)
DEALLOCATE(ZTG)
DEALLOCATE(ZDT)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:HARMONIZE_GRIB_WG_WGI_ECMWF',1,ZHOOK_HANDLE)
END SUBROUTINE HARMONIZE_GRIB_WG_WGI_ECMWF
!--------------------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_WG_ECMWF(HGRIB,KLUOUT,HINMODEL,PMASK,PFIELD,PD,PLFRAC)
!     #######################
!
USE MODD_GRID_GRIB,  ONLY : NNI
USE MODD_SURF_PAR,   ONLY : XUNDEF, NUNDEF
USE MODD_PREP_ISBA,  ONLY : XRM_WM_ECMWF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PFIELD    ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PD        ! thickness of each layer     
REAL, DIMENSION(:),   INTENT(IN), OPTIONAL :: PLFRAC ! grib land fraction 
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)                           :: IRET       ! return code
REAL,    DIMENSION(:), POINTER    :: ZSLT => NULL()       ! soil type
REAL,    DIMENSION(:), ALLOCATABLE:: ZWWILT     ! ECMWF wilting point
REAL,    DIMENSION(:), ALLOCATABLE:: ZWFC       ! ECMWF field capacity
INTEGER                           :: JL         ! loop counter on layers
INTEGER, DIMENSION(:), POINTER    :: ISLT => NULL()       ! soil type integer
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WG_ECMWF',0,ZHOOK_HANDLE)
!
IF (PRESENT(PLFRAC)) THEN
   CALL HARMONIZE_GRIB_WG_WGI_ECMWF(HGRIB,KLUOUT,HINMODEL,PMASK,PLFRAC=PLFRAC,PWG=PFIELD,PD=PD)
ELSE
   CALL HARMONIZE_GRIB_WG_WGI_ECMWF(HGRIB,KLUOUT,HINMODEL,PMASK,PWG=PFIELD,PD=PD)
END IF 
!
! 1.  Get soil type to compute SWI
!     ----------------------------
SELECT CASE (TRIM(HINMODEL))
  CASE ('ECMWF')
    CALL READ_GRIB(HGRIB,KLUOUT,43,IRET,ZSLT)
  CASE ('RACMO')
    CALL READ_GRIB(HGRIB,KLUOUT,43,IRET,ZSLT)
    IF (IRET /= 0) THEN
      CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_WG_ECMWF: COULD NOT READ SOIL TYPE')
    END IF
    IF (NINT(MINVAL(ZSLT)) < 1 .OR. NINT(MAXVAL(ZSLT)) > 7) THEN
      CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_WG_ECMWF: SOIL TYPE OUTSIDE EXPECTED RANGE (1-7)')
    END IF
  CASE DEFAULT
    CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_WG_ECMWF: UNSUPPORTED INPUT MODEL ' // TRIM(HINMODEL))
END SELECT
!--------------------------------------------------------------------------------
ALLOCATE (ZWFC(SIZE(PFIELD,1)))
ALLOCATE (ZWWILT(SIZE(PFIELD,1)))
ZWFC  (:) = 0.
ZWWILT(:) = 0.
!
IF (IRET == 0) THEN
!        
! 2.1 Convert from specific humidity to relative humidity using soil types
!     --------------------------------------------------------------------
  ALLOCATE(ISLT(SIZE(ZSLT)))
  ISLT(:) = NINT(ZSLT(:))   
  IF (PRESENT(PLFRAC)) THEN
    WHERE(PLFRAC(:).LT.XRM_WM_ECMWF)
       ISLT(:) = NUNDEF
    ENDWHERE
  END IF 
  WHERE(ISLT(:).EQ.0) ! <- EK: 0 is a dummy value; set some realistic (most common) value here
    ISLT(:) = 1
  ENDWHERE
  WHERE (ISLT(:)== NUNDEF)
    ZWFC(:) = XUNDEF
    ZWWILT(:) = XUNDEF
  ELSEWHERE (ISLT(:)==1)
    ZWFC(:)         = 0.242
    ZWWILT(:)       = 0.059
  ELSEWHERE (ISLT(:)==2)
    ZWFC(:)         = 0.346
    ZWWILT(:)       = 0.151
  ELSEWHERE (ISLT(:)==3)
    ZWFC(:)         = 0.382
    ZWWILT(:)       = 0.133
  ELSEWHERE (ISLT(:)==4)
    ZWFC(:)         = 0.448
    ZWWILT(:)       = 0.279
  ELSEWHERE (ISLT(:)==5)
    ZWFC(:)         = 0.541
    ZWWILT(:)       = 0.335
  ELSEWHERE (ISLT(:)==6)
    ZWFC(:)         = 0.662
    ZWWILT(:)       = 0.267
  ELSEWHERE (ISLT(:)==7)
    ! Tropical organic soils
    ZWFC(:)         = 0.346
    ZWWILT(:)       = 0.151     
 ENDWHERE
 DEALLOCATE(ISLT)
  !
ELSE
!
! 2.2 Convert from specific humidity to relative humidity single soil type
!     --------------------------------------------------------------------
  ! Compute model's constants
  IF (SIZE(PFIELD,2)==4) THEN
    ZWFC(:)   = 0.323
    ZWWILT(:) = 0.171
  ELSE
    ZWFC(:)   = 0.171
    ZWWILT(:) = 0.086
  END IF
  IF (PRESENT(PLFRAC)) THEN
    WHERE(PLFRAC(:).LT.XRM_WM_ECMWF)
      ZWFC(:)   =  XUNDEF
      ZWWILT(:) =  XUNDEF
    ENDWHERE 
  END IF
  !
ENDIF
!
DO JL=1,SIZE(PFIELD,2)
  IF (PRESENT(PLFRAC)) THEN
    WHERE(PLFRAC(:).LT.XRM_WM_ECMWF)
      PFIELD(:,JL) = XUNDEF
    ELSEWHERE
      PFIELD(:,JL) = (PFIELD(:,JL) - ZWWILT(:)) / (ZWFC(:) - ZWWILT(:))
    ENDWHERE 
  ELSE
    WHERE ( PFIELD(:,JL).NE.XUNDEF .AND. ZWFC(:).NE.0. ) 
      PFIELD(:,JL) = (PFIELD(:,JL) - ZWWILT(:)) / (ZWFC(:) - ZWWILT(:))
    ELSEWHERE
      PFIELD(:,JL) = 0.
    ENDWHERE
  END IF 
ENDDO
!
IF (ASSOCIATED(ZSLT)) DEALLOCATE(ZSLT)
DEALLOCATE(ZWFC)
DEALLOCATE(ZWWILT)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WG_ECMWF',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_WG_ECMWF
!----------------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_WGI_ECMWF(HGRIB,KLUOUT,HINMODEL,PMASK,PFIELD,PD,PLFRAC)
!     #######################
!
USE MODD_GRID_GRIB,  ONLY : NNI
USE MODD_SURF_PAR,   ONLY : XUNDEF, NUNDEF
USE MODD_PREP_ISBA,  ONLY : XRM_WM_ECMWF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PFIELD    ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PD        ! thickness of each layer
REAL, DIMENSION(:),   INTENT(IN), OPTIONAL :: PLFRAC ! grib land fraction
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)                           :: IRET      ! return code
REAL,    DIMENSION(:), POINTER    :: ZSLT => NULL()      ! soil type
REAL,  DIMENSION(:)  , ALLOCATABLE:: ZWSAT     ! ECMWF saturation
INTEGER                           :: JL        ! loop counter on layers
INTEGER, DIMENSION(:), POINTER    :: ISLT => NULL()       ! soil type integer
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WGI_ECMWF',0,ZHOOK_HANDLE)
!
IF (PRESENT(PLFRAC)) THEN
   CALL HARMONIZE_GRIB_WG_WGI_ECMWF(HGRIB,KLUOUT,HINMODEL,PMASK,PLFRAC=PLFRAC,PD=PD,PWGI=PFIELD)
ELSE
   CALL HARMONIZE_GRIB_WG_WGI_ECMWF(HGRIB,KLUOUT,HINMODEL,PMASK,PD=PD,PWGI=PFIELD)
END IF
!
! 1.  Get soil type to compute WSAT
!----------------------------
SELECT CASE (TRIM(HINMODEL))
  CASE ('ECMWF')
    CALL READ_GRIB(HGRIB,KLUOUT,43,IRET,ZSLT)
  CASE ('RACMO')
    CALL READ_GRIB(HGRIB,KLUOUT,43,IRET,ZSLT)
    IF (IRET /= 0) THEN
      CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_WGI_ECMWF: COULD NOT READ SOIL TYPE')
    END IF
    IF (NINT(MINVAL(ZSLT)) < 1 .OR. NINT(MAXVAL(ZSLT)) > 7) THEN
      CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_WGI_ECMWF: SOIL TYPE OUTSIDE EXPECTED RANGE (1-7)')
    END IF
  CASE DEFAULT
    CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_WGI_ECMWF: UNSUPPORTED INPUT MODEL ' // TRIM(HINMODEL))
END SELECT
!--------------------------------------------------------------------------------
ALLOCATE (ZWSAT(SIZE(PFIELD,1)))
ZWSAT(:)=0.
!
IF (IRET == 0) THEN
!        
! 2.1 Convert from specific humidity to relative humidity using soil types
   !     --------------------------------------------------------------------
  ALLOCATE(ISLT(SIZE(ZSLT)))
  ISLT(:) = NINT(ZSLT(:))
  IF (PRESENT(PLFRAC)) THEN
    WHERE(PLFRAC(:).LT.XRM_WM_ECMWF)
      ISLT(:) = NUNDEF
    ENDWHERE
  END IF
  WHERE(ISLT(:).EQ.0) ! <- EK: 0 is a dummy value; set some realistic (most common) value here
    ISLT(:) = 1
  ENDWHERE
  WHERE (ISLT(:)== NUNDEF)
    ZWSAT(:) = XUNDEF
  ELSEWHERE (ISLT(:)==1)
    ZWSAT(:) = 0.403
  ELSEWHERE (ISLT(:)==2)
    ZWSAT(:) = 0.439
  ELSEWHERE (ISLT(:)==3)
    ZWSAT(:) = 0.430
  ELSEWHERE (ISLT(:)==4)
    ZWSAT(:) = 0.520
  ELSEWHERE (ISLT(:)==5)
    ZWSAT(:) = 0.614
  ELSEWHERE (ISLT(:)==6)
    ZWSAT(:) = 0.766
  ELSEWHERE (ISLT(:)==7)
    ! Tropical organic soils
    ZWSAT(:) = 0.439
 ENDWHERE
 DEALLOCATE(ISLT)
!
ELSE
!
! 2.2 Convert from specific humidity to relative humidity single soil type
!     --------------------------------------------------------------------
  ! Compute model's constants
  IF (SIZE(PFIELD,2)==4) THEN
    ZWSAT(:)  = 0.472
  ELSE
    ZWSAT(:)  = 0.286
  END IF
  IF (PRESENT(PLFRAC)) THEN
    WHERE(PLFRAC(:).LT.XRM_WM_ECMWF)
      ZWSAT(:) = XUNDEF
    ENDWHERE
  END IF
  !
ENDIF
!
! Then perform conversion
DO JL=1,SIZE(PFIELD,2)
  IF (PRESENT(PLFRAC)) THEN
    WHERE(PLFRAC(:).LT.XRM_WM_ECMWF)
      PFIELD(:,JL) = XUNDEF
    ELSEWHERE
      PFIELD(:,JL) = PFIELD(:,JL) / ZWSAT(:)
    ENDWHERE
  ELSE
    WHERE ( PFIELD(:,JL).NE.XUNDEF .AND. ZWSAT(:).NE.0. ) 
      PFIELD(:,JL) = PFIELD(:,JL) / ZWSAT(:)
    ELSEWHERE
      PFIELD(:,JL) = 0.
    ENDWHERE     
  END IF
ENDDO
!
IF (ASSOCIATED(ZSLT)) DEALLOCATE(ZSLT)
DEALLOCATE(ZWSAT)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WGI_ECMWF',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_WGI_ECMWF
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_TG_METEO_FRANCE(HGRIB,KLUOUT,HINMODEL,PMASK,PTG,PDT)
!     #######################
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PTG       ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PDT       ! thickness of each layer
!* local variables
!  ---------------
REAL,    DIMENSION(:), POINTER    :: ZFIELD => NULL()    ! field to read
INTEGER                           :: JL         ! layer loop counter
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TG_METEO_FRANCE',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_TG_METEO_FRANCE: | Reading soil temperature'
!--------------------------------------------------------------------------------
! 1.  Allocate soil temperature profile
!     ---------------------------------
!--------------------------------------------------------------------------------
! 2.  Search and read level 1 (and its depth)
!     ---------------------------------------
 CALL READ_GRIB_TS(HGRIB,KLUOUT,HINMODEL,PMASK,ZFIELD)
!
ALLOCATE(PTG(SIZE(ZFIELD),3))
ALLOCATE(PDT(SIZE(ZFIELD),3))
!
PTG(:,1) = ZFIELD(:)
PDT(:,1) = 0.01
!--------------------------------------------------------------------------------
! 3.  Deep soil temperature
!     ---------------------
 CALL READ_GRIB_T2_LAND(HGRIB,KLUOUT,HINMODEL,PMASK,ZFIELD)
!
PTG(:,2) = ZFIELD(:)
PDT(:,2) = 0.4         ! deep temperature layer depth assumed equal to 0.4m
DEALLOCATE(ZFIELD)
!--------------------------------------------------------------------------------
! 4.  Assumes uniform temperature profile below
!     -----------------------------------------
PTG(:,3) = PTG(:,2)
PDT(:,3) = 5.          ! temperature profile down to 5m
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TG_METEO_FRANCE',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_TG_METEO_FRANCE
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_SAND_CLAY_METEO_FRANCE(HGRIB,KLUOUT,HINMODEL,PSAND,PCLAY,GISBA)
!     ######################
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), POINTER         :: PSAND     ! field to initialize
REAL, DIMENSION(:), POINTER         :: PCLAY     ! thickness of each layer
LOGICAL, INTENT(OUT)                :: GISBA     ! T: surface scheme in file is ISBA
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)                           :: IRET      ! return code
INTEGER                           :: IPAR      ! parameter number for field reading
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
! 1.  Search and read clay fraction if available
!     ------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_SAND_CLAY_METEO_FRANCE',0,ZHOOK_HANDLE)
!
IF (HINMODEL == 'ARPEGE' .OR. HINMODEL == 'MOCAGE') IPAR=171
IF (HINMODEL == 'ALADIN') IPAR=128
 CALL READ_GRIB(HGRIB,KLUOUT,IPAR,IRET,PCLAY)
!
! if not available, the model is not ISBA (IWMODE=1)
IF (IRET /= 0) THEN
  GISBA = .FALSE.
ELSE
  GISBA = .TRUE.
  PCLAY(:) = PCLAY(:) / 100. ! this field is given in percent
  WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_SAND_CLAY_METEO_FRANCE: | The soil model is ISBA'
END IF
!-------------------------------------------------------------------------------
! 2.  Search and read sand fraction if available
!     ------------------------------------------
IF (HINMODEL == 'ARPEGE' .OR. HINMODEL == 'MOCAGE') IPAR=172
IF (HINMODEL == 'ALADIN') IPAR=129
 CALL READ_GRIB(HGRIB,KLUOUT,IPAR,IRET,PSAND)
!
! if not available, the model is not ISBA (IWMODE=1)
IF (GISBA) THEN
  IF (IRET /= 0) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB: SAND FRACTION MISSING (READ_GRIB_SAND_CLAY_METEO_FRANCE)')
  ELSE
    PSAND(:) = PSAND(:) / 100. ! this field is given in percent
  END IF
END IF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_SAND_CLAY_METEO_FRANCE',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_SAND_CLAY_METEO_FRANCE
!-----------------------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_WG_METEO_FRANCE(HGRIB,KLUOUT,HINMODEL,PMASK,PFIELD,PD)
!     #######################
!
USE MODD_GRID_GRIB,  ONLY : NNI
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PFIELD    ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PD        ! thickness of each layer
!
!* local variables
!  ---------------
LOGICAL                           :: GISBA     ! T: surface scheme in file is ISBA
INTEGER(KIND=kindOfInt)                           :: IRET      ! return code
INTEGER                           :: ILTYPE    ! type of level (Grib code table 3)
INTEGER                           :: ILEV1     ! level definition
INTEGER                           :: ILEV2     ! level definition
REAL,  DIMENSION(:),   POINTER    :: ZCLAY => NULL()     ! clay fraction
REAL,  DIMENSION(:),   POINTER    :: ZSAND => NULL()     ! sand fraction
REAL,  DIMENSION(:),   POINTER    :: ZFIELD => NULL()
REAL,  DIMENSION(:),   ALLOCATABLE:: ZWWILT     ! wilting point
REAL,  DIMENSION(:),   ALLOCATABLE:: ZWFC       ! field capacity
REAL,  DIMENSION(:),   ALLOCATABLE:: ZWSAT      ! saturation
INTEGER                           :: JL         ! layer loop counter
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WG_METEO_FRANCE',0,ZHOOK_HANDLE)
!
! 1.  Search and read clay and sand fractions if available
!     ----------------------------------------------------
 CALL READ_GRIB_SAND_CLAY_METEO_FRANCE(HGRIB,KLUOUT,HINMODEL,ZSAND,ZCLAY,GISBA)
!-------------------------------------------------------------------------------
IF (GISBA) THEN
  ALLOCATE(PFIELD(SIZE(ZSAND),3))
  ALLOCATE(PD(SIZE(ZSAND),3))
ELSE
  ALLOCATE(PFIELD(NNI,3))
  ALLOCATE(PD(NNI,3))
ENDIF
!
PD(:,1) = 0.01
PD(:,2) = 0.20
!-------------------------------------------------------------------------------
! 2.  Read layer 1 moisture
!     ---------------------
ILEV1   = 0
IF (HINMODEL == 'ARPEGE' .OR. HINMODEL=='MOCAGE') THEN
  ILTYPE  = 112     
  ILEV2   = 1        
  CALL READ_GRIB(HGRIB,KLUOUT,153,IRET,ZFIELD,KLEV1=ILEV1,KLEV2=ILEV2)
ELSE
  ILTYPE  = 105
  ILEV2   = 0
  CALL READ_GRIB(HGRIB,KLUOUT,86,IRET,ZFIELD,KLEV1=ILEV1,KLEV2=ILEV2)
ENDIF
IF (IRET /= 0) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB: SOIL MOISTURE LEVEL 1 MISSING (READ_GRIB_WG_METEO_FRANCE)')
END IF
!
PFIELD(:,1) = ZFIELD(:)
!-------------------------------------------------------------------------------
! 3.  Read layer 2 moisture
!     ---------------------
IF (HINMODEL == 'ARPEGE' .OR. HINMODEL=='MOCAGE') THEN
  ILTYPE  = 112  
  ILEV1   = 0  
  ILEV2   = 250             
  CALL READ_GRIB(HGRIB,KLUOUT,153,IRET,ZFIELD,KLEV1=ILEV1,KLEV2=ILEV2)
ELSE
  ILTYPE  = 111
  ILEV1   = -1
  ILEV2   = -1        
  CALL READ_GRIB(HGRIB,KLUOUT,86,IRET,ZFIELD,KLEV1=ILEV1,KLEV2=ILEV2)
ENDIF
IF (IRET /= 0) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB: SOIL MOISTURE LEVEL 2 MISSING (READ_GRIB_WG_METEO_FRANCE)')
END IF
!        
PFIELD(:,2) = ZFIELD(:)
!-------------------------------------------------------------------------------
! 4.  Read layer 2 depth (ISBA only)
!     -----------------------------
!* note that soil water reservoir is considered uniform between 0.2m and GRIB soil depth
IF (GISBA) THEN
  IF (HINMODEL == 'ARPEGE' .OR. HINMODEL == 'MOCAGE') THEN
    CALL READ_GRIB(HGRIB,KLUOUT,173,IRET,ZFIELD)
  ELSE
    CALL READ_GRIB(HGRIB,KLUOUT,130,IRET,ZFIELD)
  ENDIF
  IF (IRET /= 0) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB: SOIL MOISTURE LEVEL 2 DEPTH MISSING (READ_GRIB_WG_METEO_FRANCE)')          
  END IF
  PD(:,3) = ZFIELD(:)
  DEALLOCATE(ZFIELD)
ELSE
  PD(:,3) = 2.
END IF
!-------------------------------------------------------------------------------
! 5.  Compute relative humidity from units kg/m^2
!     -------------------------------------------
! Compute ISBA model constants (if needed)
IF (GISBA) THEN
  !
  !* updates Wg in m3/m3
  PFIELD(:,1) = PFIELD(:,1) / 10.
  PFIELD(:,2) = PFIELD(:,2) /(1000. * PD(:,3))
  !
  ALLOCATE (ZWSAT (SIZE(ZSAND)))
  ZWSAT (:) = (-1.08*100. * ZSAND(:) + 494.305) * 0.001
  PFIELD(:,1) = MAX(MIN(PFIELD(:,1),ZWSAT(:)),0.)
  PFIELD(:,2) = MAX(MIN(PFIELD(:,2),ZWSAT(:)),0.)
  DEALLOCATE(ZWSAT)
  DEALLOCATE (ZSAND)
  !
  ALLOCATE (ZWWILT(SIZE(ZCLAY)))
  ALLOCATE (ZWFC  (SIZE(ZCLAY)))
  ZWWILT(:) = 37.1342E-3 * SQRT( 100. * ZCLAY(:) )
  ZWFC  (:) = 89.0467E-3 * (100. * ZCLAY(:) )**0.3496
  PFIELD(:,1) = (PFIELD(:,1) - ZWWILT(:)) / (ZWFC(:) - ZWWILT(:))
  PFIELD(:,2) = (PFIELD(:,2) - ZWWILT(:)) / (ZWFC(:) - ZWWILT(:))
  DEALLOCATE (ZWWILT)
  DEALLOCATE (ZWFC)
  DEALLOCATE (ZCLAY)
  !
ELSE ! Non ISBA
 !
  PFIELD(:,2) = (PFIELD(:,1)+PFIELD(:,2)) / (20. + 100.) 
  PFIELD(:,1) =  PFIELD(:,1)           /  20.
  !
END IF
!
PFIELD(:,3) = PFIELD(:,2)
!--------------------------------------------------------------------------------
! 6.  Apply land mask
!     ---------------
IF (SIZE(PMASK)==SIZE(PFIELD,1)) THEN
  DO JL=1,SIZE(PFIELD,2)
    WHERE (PMASK(:)/=1.) PFIELD(:,JL) = XUNDEF
  END DO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WG_METEO_FRANCE',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_WG_METEO_FRANCE
!--------------------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_WGI_METEO_FRANCE(HGRIB,KLUOUT,HINMODEL,PMASK,PFIELD,PD)
!     #######################
!
USE MODD_GRID_GRIB,  ONLY : NNI
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PFIELD    ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PD        ! thickness of each layer
!
!* local variables
!  ---------------
LOGICAL                           :: GISBA     ! T: surface scheme in file is ISBA
INTEGER(KIND=kindOfInt)                           :: IRET      ! return code
INTEGER                           :: ILTYPE    ! type of level (Grib code table 3)
INTEGER                           :: ILEV1     ! level definition
INTEGER                           :: ILEV2     ! level definition
REAL,  DIMENSION(:),   POINTER    :: ZCLAY => NULL()     ! clay fraction
REAL,  DIMENSION(:),   POINTER    :: ZSAND => NULL()     ! sand fraction
REAL,  DIMENSION(:),   POINTER    :: ZFIELD => NULL()
REAL,  DIMENSION(:),   ALLOCATABLE:: ZWSAT      ! saturation
INTEGER                           :: JL         ! layer loop counter
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WGI_METEO_FRANCE',0,ZHOOK_HANDLE)
!
! 1.  Search and read clay fraction if available
!     ------------------------------------------
 CALL READ_GRIB_SAND_CLAY_METEO_FRANCE(HGRIB,KLUOUT,HINMODEL,ZSAND,ZCLAY,GISBA)
!-------------------------------------------------------------------------------
IF (GISBA) THEN
  ALLOCATE(PFIELD(SIZE(ZSAND),2))
  ALLOCATE(PD(SIZE(ZSAND),2))
ELSE
  ALLOCATE(PFIELD(NNI,2))
  ALLOCATE(PD(NNI,2))
ENDIF
!
PD(:,1) = 0.01
!-------------------------------------------------------------------------------
! 2.  Read layer 1 soil ice
!     ---------------------
ILEV1   = 0  
IF (HINMODEL == 'ARPEGE' .OR. HINMODEL=='MOCAGE') THEN
  ILTYPE  = 112     
  ILEV2   = 1
  CALL READ_GRIB(HGRIB,KLUOUT,152,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2)
ELSE
  ILTYPE  = 105      
  ILEV2   = 0        
  CALL READ_GRIB(HGRIB,KLUOUT,139,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2)
END IF
!
IF (IRET == 0) THEN
  PFIELD(:,1) = ZFIELD(:)
  WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_WGI_METEO_FRANCE: -> Soil ice level 1 is present'
ELSE
  PFIELD(:,1) = 0.
END IF
!-------------------------------------------------------------------------------
! 3.  Read layer 2 soil ice
!     ---------------------
IF (HINMODEL == 'ARPEGE' .OR. HINMODEL=='MOCAGE') THEN
  ILTYPE  = 112
  ILEV1   = 0        
  ILEV2   = 250
  CALL READ_GRIB(HGRIB,KLUOUT,152,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2)
ELSE
  ILTYPE  = 111
  ILEV1   = -1        
  ILEV2   = -1        
  CALL READ_GRIB(HGRIB,KLUOUT,139,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2)
END IF
!
IF (IRET == 0) THEN
  PFIELD(:,2) = ZFIELD(:)     
  WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_WGI_METEO_FRANCE: -> Soil ice level 2 is present'
ELSE
  PFIELD(:,2) = 0.
END IF
!-------------------------------------------------------------------------------
! 4.  Read layer 2 depth (ISBA only)
!     ------------------------------
IF (GISBA) THEN
  IF (HINMODEL == 'ARPEGE' .OR. HINMODEL=='MOCAGE') THEN 
    CALL READ_GRIB(HGRIB,KLUOUT,173,IRET,ZFIELD)
  ELSE
    CALL READ_GRIB(HGRIB,KLUOUT,130,IRET,ZFIELD)
  ENDIF
  IF (IRET /= 0) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB: SOIL ICE LEVEL 2 MISSING (READ_GRIB_WGI_METEO_FRANCE)')
  END IF
  PD(:,2) = ZFIELD(:)
  DEALLOCATE(ZFIELD)
ELSE
  PD(:,2) = 2.
END IF
!-------------------------------------------------------------------------------
! 5.  Compute relative humidity from units kg/m^2
!     -------------------------------------------
IF (GISBA) THEN
  !
  !* updates Wgi in m3/m3
  PFIELD(:,1) = PFIELD(:,1) / 10.
  PFIELD(:,2) = PFIELD(:,2) /(1000. * PD(:,2))
  !        
  ALLOCATE (ZWSAT (NNI))
  ZWSAT (:) = (-1.08*100. * ZSAND(:) + 494.305) * 0.001
  PFIELD(:,1) = PFIELD(:,1) / ZWSAT(:)
  PFIELD(:,2) = PFIELD(:,2) / ZWSAT(:)
  DEALLOCATE (ZWSAT) 
  DEALLOCATE (ZSAND)
  DEALLOCATE (ZCLAY)
  !
ELSE ! Non ISBA
  !
  PFIELD(:,1) = 0.
  PFIELD(:,2) = 0.
  !
END IF
!--------------------------------------------------------------------------------
! 6.  Apply land mask
!     ---------------
IF (SIZE(PMASK)==SIZE(PFIELD,1)) THEN
  DO JL=1,SIZE(PFIELD,2)
    WHERE (PMASK(:)/=1.) PFIELD(:,JL) = XUNDEF
  END DO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WGI_METEO_FRANCE',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_WGI_METEO_FRANCE
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_TG_HIRLAM(HGRIB,KLUOUT,HINMODEL,PMASK,PTG,PDT)
!     #######################
!
USE MODD_GRID_GRIB,  ONLY : NNI
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PTG       ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PDT       ! thickness of each layer
!
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)                           :: IRET      ! return code
INTEGER                           :: ILEV1     ! level definition
INTEGER                           :: ILEV2     ! level definition
REAL,  DIMENSION(:),   POINTER    :: ZFIELD => NULL()
INTEGER                           :: JL         ! layer loop counter
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TG_HIRLAM',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_TG_HIRLAM: | Reading soil temperature'
!--------------------------------------------------------------------------------
! 1.  Search and read level 1 (and its depth)
!     -----------------------
ILEV1 = 904
ILEV2 = -1
 CALL READ_GRIB(HGRIB,KLUOUT,11,IRET,ZFIELD,KLEV1=ILEV1,KLEV2=ILEV2)
IF (IRET /= 0 ) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB: SOIL TEMPERATURE LEVEL 1 MISSING (READ_GRIB_TG_HIRLAM)')
END IF
!
ALLOCATE(PTG(SIZE(ZFIELD),3))
ALLOCATE(PDT(SIZE(ZFIELD),3))
PTG(:,1)= ZFIELD(:)
PDT(:,1) = 0.01
!--------------------------------------------------------------------------------
! 2.  Deep soil temperature
!     ---------------------
WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_TG_HIRLAM: | Reading deep soil temperature'
!
ILEV1 = 954
ILEV2 = -1
 CALL READ_GRIB(HGRIB,KLUOUT,11,IRET,ZFIELD,KLEV1=ILEV1,KLEV2=ILEV2)
IF (IRET /= 0) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB: DEEP SOIL TEMPERATURE MISSING (READ_GRIB_TG_HIRLAM)')
END IF
!
PTG(:,2)= ZFIELD(:)
DEALLOCATE(ZFIELD)
PDT(:,2) = 0.4         ! deep temperature layer depth assumed equal to 0.40m
!--------------------------------------------------------------------------------
! 4.  Assumes uniform temperature profile below
!     -----------------------------------------
PTG(:,3) = PTG(:,2)
PDT(:,3) = 5.          ! temperature profile down to 5m
!--------------------------------------------------------------------------------
! 5.  Apply land mask
!     ---------------
IF (SIZE(PMASK)==SIZE(PTG,1)) THEN
  DO JL=1,SIZE(PTG,2)
    WHERE (PMASK(:)/=1.) PTG(:,JL) = XUNDEF
  END DO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TG_HIRLAM',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_TG_HIRLAM
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_WG_HIRLAM(HGRIB,KLUOUT,HINMODEL,PMASK,PFIELD,PD)
!     #######################
!
USE MODD_GRID_GRIB,  ONLY : NNI
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), INTENT(IN)      :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PFIELD    ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PD        ! thickness of each layer
!
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)                           :: IRET      ! return code
INTEGER                           :: ILTYPE    ! type of level (Grib code table 3)
INTEGER                           :: ILEV1     ! level definition
INTEGER                           :: ILEV2     ! level definition
REAL, DIMENSION(:), POINTER       :: ZFIELD => NULL()
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWG        ! first water reservoir
REAL, DIMENSION(:), ALLOCATABLE   :: ZD         ! Height of each layer
INTEGER                           :: INLAYERDEEP! number of deep moisture layers
REAL                              :: ZWWILT     ! ECMWF wilting point
REAL                              :: ZWFC       ! ECMWF field capacity
REAL                              :: ZWSAT      ! ECMWF saturation
INTEGER                           :: JL         ! loop counter on layers
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WG_HIRLAM',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_WG_HIRLAM: | Reading soil moisture'
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_WG_HIRLAM: | WARNING READING LOW VEGETATION TILE (NR 4) ONLY'
!
ALLOCATE(ZD(2))
!
! 1.  Search and read level 1 (and its depth)
!     -----------------------
ILTYPE=105
ILEV1=904
ILEV2=-1
 CALL READ_GRIB(HGRIB,KLUOUT,86,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2)
IF (IRET /= 0 ) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB: SOIL MOISTURE LEVEL 1 MISSING (READ_GRIB_WG_HIRLAM)')
END IF
!
ALLOCATE(ZWG(SIZE(ZFIELD),2))
ZWG(:,1)=ZFIELD
!
ZD(1) = 0.01
WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_WG_HIRLAM: | Level 1 height set to 0.01 m '
!
ZWG(:,1) = ZWG(:,1) / ZD(1)      ! convert units to m3/m3  
!--------------------------------------------------------------------------------
! 2.  Search and read level 2 (and its depth) This level is optionnal
!     -----------------------
ILTYPE=105
ILEV1=954
ILEV2=-1
 CALL READ_GRIB(HGRIB,KLUOUT,86,IRET,ZFIELD,KLTYPE=ILTYPE,KLEV1=ILEV1,KLEV2=ILEV2)
IF (IRET /= 0 ) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB: SOIL MOISTURE LEVEL 2 MISSING (READ_GRIB_WG_HIRLAM)')
END IF
!
ZWG(:,2)=ZFIELD ! already units m3/m3 
DEALLOCATE(ZFIELD)
!
INLAYERDEEP = 2
ZD(2) = 0.42
WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_WG_HIRLAM: | Level 2 height set to 0.42 m '
!
WRITE  (KLUOUT,'(A)') 'WARNING MODE_READ_GRIB: ZWG3 AND ZWG4 SET TO 0. (READ_GRIB_WG_HIRLAM)'
!--------------------------------------------------------------------------------
! 3.  Set water content profile and layer thicknesses
!     -----------------------------------------------
 CALL FILL_PFIELD(KLUOUT,'READ_GRIB_WG_HIRLAM',INLAYERDEEP,ZD,ZWG,PMASK,PFIELD,PD)
DEALLOCATE(ZD)
DEALLOCATE(ZWG)
!--------------------------------------------------------------------------------
! 4.  Convert from specific humidity to relative humidity
!     ---------------------------------------------------
! Compute model's constants
ZWFC   = 0.171
ZWWILT = 0.086
!
! Then perform conversion
DO JL=1,INLAYERDEEP
  WHERE (PFIELD(:,JL).NE.XUNDEF) PFIELD(:,JL) = (PFIELD(:,JL) - ZWWILT) / (ZWFC - ZWWILT)
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WG_HIRLAM',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_WG_HIRLAM
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_WGI_HIRLAM(HGRIB,KLUOUT,PFIELD,PD)
!     #######################
!
USE MODD_GRID_GRIB,  ONLY : NNI
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
REAL, DIMENSION(:,:), POINTER       :: PFIELD    ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PD        ! thickness of each layer
!
!* local variables
!  ---------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WGI_HIRLAM',0,ZHOOK_HANDLE)
!
ALLOCATE (PFIELD(NNI,2))
ALLOCATE (PD    (NNI,2))
PFIELD(:,:) = 0.
!
PD    (:,1) = 0.01
PD    (:,2) = 1.
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WGI_HIRLAM',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_WGI_HIRLAM
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE SURF_INQ(KLUOUT,HINMODEL,KNCSS,PTHICK)
!     #######################
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
INTEGER,              INTENT(IN)              :: KLUOUT    ! logical unit of output listing
CHARACTER(LEN=6),     INTENT(IN)              :: HINMODEL  ! Grib originating model
INTEGER,              INTENT(OUT)  , OPTIONAL :: KNCSS     ! number of soil levels
REAL, DIMENSION(:),   POINTER      , OPTIONAL :: PTHICK    ! thickness of each layer
!
!* local variables
!  ---------------
INTEGER                                       :: IERR      ! error handling
INTEGER                                       :: INCSS     ! number of soil levels
REAL(KIND=JPHOOK)                               :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:SURF_INQ',0,ZHOOK_HANDLE)
!
SELECT CASE (TRIM(HINMODEL))
  CASE ('RACMO')
    ! Number of soil layers
    INCSS = 4
    IF (PRESENT(KNCSS)) THEN
      KNCSS = INCSS
    END IF
    
    ! Thickness of each soil layer
    IF (PRESENT(PTHICK)) THEN
      IF (ASSOCIATED(PTHICK)) THEN
        CALL ABOR1_SFX('MODE_READ_GRIB:SURF_INQ: PTHICK ALREADY ASSOCIATED, CHECK CODE')
      END IF
      
      ALLOCATE(PTHICK(INCSS),STAT=IERR)
      IF (IERR /= 0) THEN
        CALL ABOR1_SFX('MODE_READ_GRIB:SURF_INQ: FAILED TO ALLOCATE MEMORY FOR PTHICK')
      END IF
      PTHICK(1) = 0.07
      PTHICK(2) = 0.21
      PTHICK(3) = 0.72
      PTHICK(4) = 1.89
    END IF
  CASE DEFAULT
    CALL ABOR1_SFX('MODE_READ_GRIB:SURF_INQ: UNSUPPORTED INPUT MODEL ' // TRIM(HINMODEL))
END SELECT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:SURF_INQ',1,ZHOOK_HANDLE)
END SUBROUTINE SURF_INQ
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_TG_RACMO(HGRIB,KLUOUT,HINMODEL,PMASK,PTG,PD)
!     #######################
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PTG       ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PD        ! depth of each layer
!
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)             :: IRET      ! return code
INTEGER                             :: JL        ! layer loop counter
INTEGER                             :: INCSS     ! number of soil layers
REAL, DIMENSION(:),   POINTER       :: ZFIELD => NULL() ! temperature for 1 layer
REAL, DIMENSION(:),   POINTER       :: ZTHICK => NULL() ! soil layer thickness per layer
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZTG       ! temperature per layer
CHARACTER(LEN=1)                    :: YLEV      ! printing string for level
REAL(KIND=JPHOOK)                     :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TG_RACMO',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_TG_RACMO: | Reading soil temperature'
!
! 1.  Get number of levels and their thickness
!     ----------------------------------------
CALL SURF_INQ(KLUOUT,HINMODEL,KNCSS=INCSS,PTHICK=ZTHICK)
!
! 2.  Search and read temperature per level, use shortName as paramId is not consecutive
!     ----------------------------------------------------------------------------------
DO JL = 1, INCSS, 1
  WRITE(YLEV,FMT='(I1)') JL
  CALL READ_GRIB(HGRIB,KLUOUT,-999,IRET,ZFIELD,SHORTNAME='stl'//YLEV)
  IF (IRET /= 0) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB: SOIL TEMPERATURE FOR LEVEL ' // YLEV // ' MISSING (READ_GRIB_TG_RACMO)')
  ENDIF
  IF (.NOT. ALLOCATED(ZTG)) THEN
    ALLOCATE(ZTG(SIZE(ZFIELD),INCSS))
  ENDIF
  ZTG(:,JL) = ZFIELD
ENDDO
!
! 3.  Set temperature profile and layer thicknesses
!     ---------------------------------------------
CALL FILL_PFIELD(KLUOUT,'READ_GRIB_TG_RACMO',INCSS,ZTHICK,ZTG,PMASK,PTG,PD)
!
! 4.  Cleanup
!     -------
DEALLOCATE(ZTHICK)
DEALLOCATE(ZTG)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TG_RACMO',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_TG_RACMO
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_WG_RACMO_1(HGRIB,KLUOUT,HINMODEL,PMASK,PWG,PD)
!     #######################
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PWG       ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PD        ! depth of each layer
!
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)             :: IRET      ! return code
INTEGER                             :: JL        ! layer loop counter
INTEGER                             :: INCSS     ! number of soil layers
REAL, DIMENSION(:),   POINTER       :: ZFIELD => NULL() ! moisture for 1 layer
REAL, DIMENSION(:),   POINTER       :: ZTHICK => NULL() ! soil layer thickness per layer
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZWG       ! moisture per layer
CHARACTER(LEN=1)                    :: YLEV      ! printing string for level
REAL(KIND=JPHOOK)                     :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WG_RACMO_1',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_WG_RACMO_1: | Reading soil moisture'
!
! 1.  Get number of levels and their thickness
!     ----------------------------------------
CALL SURF_INQ(KLUOUT,HINMODEL,KNCSS=INCSS,PTHICK=ZTHICK)
!
! 2.  Search and read moisture per level, use shortName
!     -------------------------------------------------
DO JL = 1, INCSS, 1
  WRITE(YLEV,FMT='(I1)') JL
  CALL READ_GRIB(HGRIB,KLUOUT,-999,IRET,ZFIELD,SHORTNAME='swvl'//YLEV)
  IF (IRET /= 0) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB: SOIL MOISTURE FOR LEVEL ' // YLEV // ' MISSING (READ_GRIB_WG_RACMO_1)')
  ENDIF
  IF (.NOT. ALLOCATED(ZWG)) THEN
    ALLOCATE(ZWG(SIZE(ZFIELD),INCSS))
  ENDIF
  ZWG(:,JL) = ZFIELD
ENDDO
!
! 3.  Set moisture profile and layer thicknesses
!     ------------------------------------------
CALL FILL_PFIELD(KLUOUT,'READ_GRIB_WG_RACMO_1',INCSS,ZTHICK,ZWG,PMASK,PWG,PD)
!
! 4.  Cleanup
!     -------
DEALLOCATE(ZTHICK)
DEALLOCATE(ZWG)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WG_RACMO_1',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_WG_RACMO_1
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!     #######################

      SUBROUTINE READ_GRIB_SNOW_VEG_AND_DEPTH(HGRIB,KLUOUT,HINMODEL,PMASK,PLFRAC,PSNV,PSNVD)

!     #######################
!!
!!    MODIFICATIONS
!!    -------------
!!    C Ardilouze 07/2013 : possibility to read snow density (ERAI-land)
!!
!-------------------------------------------------------------------------------
!
USE MODD_GRID_GRIB,  ONLY : NNI
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_SNOW_PAR,   ONLY : XRHOSMAX
USE MODD_PREP_ISBA,  ONLY : XRM_WM_ECMWF
USE MODD_PREP_SNOW,  ONLY : XRM_LITTLE_SNOW
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:),   INTENT(IN), OPTIONAL :: PLFRAC ! grib land fraction
REAL, DIMENSION(:), OPTIONAL, POINTER :: PSNV    ! field to initialize
REAL, DIMENSION(:), OPTIONAL, POINTER :: PSNVD   ! field to initialize
!
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)                           :: IRET      ! return code
REAL, DIMENSION(:), POINTER       :: ZFIELD => NULL()    ! field to initialize
REAL, DIMENSION(:), POINTER       :: ZFIELD2 => NULL()    ! field to initialize
REAL, DIMENSION(:), POINTER       :: ZRHO => NULL() ! local field for snow density
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_SNOW_VEG_AND_DEPTH',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_SNOW_VEG_AND_DEPTH: | Reading snow depth and density (if present)'
!
IF (PRESENT(PLFRAC).AND.HINMODEL.NE.'ECMWF ') THEN
   CALL ABOR1_SFX('MODE_READ_GRIB: OPTION IS IMPLEMENTED ONLY FOR ECMWF (READ_GRIB_SNOW_VEG_AND_DEPTH)')
END IF
!
SELECT CASE(HINMODEL)
  CASE('ECMWF ','RACMO ')
    CALL READ_GRIB(HGRIB,KLUOUT,141,IRET,ZFIELD)
  CASE('ARPEGE','ALADIN','MOCAGE','HIRLAM')
    CALL READ_GRIB(HGRIB,KLUOUT,66,IRET,ZFIELD)          
  CASE DEFAULT
    CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_SNOW_VEG_AND_DEPTH: OPTION NOT SUPPORTED '//HINMODEL)
END SELECT
!
IF (IRET /= 0 ) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB: SNOW AND VEG DEPTH MISSING (READ_GRIB_SNOW_VEG_AND_DEPTH)')
END IF
!
IF (PRESENT(PLFRAC)) THEN
  CALL READ_GRIB_SNOW_DEN(HGRIB,KLUOUT,HINMODEL,PMASK,ZRHO,PLFRAC)
ELSE
  CALL READ_GRIB_SNOW_DEN(HGRIB,KLUOUT,HINMODEL,PMASK,ZRHO)
END IF
!
IF (SIZE(ZFIELD).NE.SIZE(ZRHO)) THEN
   CALL ABOR1_SFX('MODE_READ_GRIB: SNOW AND SNOW DENS. DIFFER (READ_GRIB_SNOW_VEG_AND_DEPTH)')
END IF
IF (PRESENT(PLFRAC)) THEN
   IF(SIZE(PLFRAC).NE.SIZE(ZFIELD)) THEN
      CALL ABOR1_SFX('MODE_READ_GRIB: SNOW AND LAND FRACTION DIFFER (READ_GRIB_SNOW_VEG_AND_DEPTH)')
   END IF
END IF
!
IF (PRESENT(PSNV)) THEN
  ALLOCATE(PSNV(SIZE(ZFIELD)))
  PSNV(:)=ZFIELD(:)
!  IF (HINMODEL=='ECMWF ') PSNV(:) = PSNV(:) * ZRHO(:)                                                                                                            
! Transform SWE from tonn/m2 to kg/m2                                                                                                                             
  IF (HINMODEL=='ECMWF '.OR. HINMODEL=='RACMO ') THEN
     PSNV(:) = PSNV(:) * 1000.
  END IF
  IF (XRM_LITTLE_SNOW.GT.0.) THEN
     WHERE(PSNV(:).LT.XRM_LITTLE_SNOW)
        PSNV(:)=0.
     END WHERE
  END IF
  IF (PRESENT(PLFRAC)) THEN
     WHERE (PLFRAC(:).LT.XRM_WM_ECMWF)
        PSNV(:) = XUNDEF
     END WHERE
  ELSE
    IF (SIZE(PMASK)==SIZE(PSNV)) &
      WHERE (PMASK(:)/=1.) PSNV(:) = XUNDEF
  END IF
ENDIF
!
IF (PRESENT(PSNVD)) THEN
  ALLOCATE(PSNVD(SIZE(ZFIELD)))
  PSNVD(:)=ZFIELD(:)  
  IF (HINMODEL=='ECMWF '.OR. HINMODEL=='RACMO ') THEN
     ! convert snow amount from m. water equivalent to kg/m2, then divide by density to get m. snow     
     PSNVD(:)=PSNVD(:)*1000.
  END IF
  IF (XRM_LITTLE_SNOW.GT.0.) THEN
     WHERE(PSNVD(:).LT.XRM_LITTLE_SNOW)
        PSNVD(:)=0.
     END WHERE
  END IF
  PSNVD(:) = PSNVD(:) / ZRHO(:) 
  IF (PRESENT(PLFRAC)) THEN
     WHERE (PLFRAC(:).LT.XRM_WM_ECMWF)
        PSNVD(:) = XUNDEF
     END WHERE
  ELSE
    IF (SIZE(PMASK)==SIZE(PSNVD)) THEN
      WHERE (PMASK(:)/=1.) &
        PSNVD(:) = XUNDEF
    END IF
  END IF 
ENDIF
!
DEALLOCATE(ZFIELD)
DEALLOCATE(ZRHO)
!
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_SNOW_VEG_AND_DEPTH',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_SNOW_VEG_AND_DEPTH
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_SNOW_ALB(HGRIB,KLUOUT,HINMODEL,PMASK,PSNVA,PLFRAC)
!     #######################
!!
!!    AUTHOR
!!    -------------
!!    C Ardilouze 07/2013 : possibility to read snow albedo (ERAI-land)
!!
!-------------------------------------------------------------------------------
!
USE MODD_GRID_GRIB,  ONLY : NNI
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_SNOW_PAR,   ONLY : XANSMIN, XANSMAX
USE MODD_PREP_ISBA,  ONLY : XRM_WM_ECMWF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB    ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT   ! logical unit of output listing
CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK    ! grib land mask
REAL, DIMENSION(:), POINTER         :: PSNVA    ! field to initialize
REAL, DIMENSION(:),   INTENT(IN), OPTIONAL :: PLFRAC ! grib land fraction
!
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)           :: IRET      ! return code
REAL, DIMENSION(:), POINTER       :: ZFIELD => NULL()    ! field to initialize 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_SNOW_ALB',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_SNOW_ALB: | Reading snow albedo'
!
ALLOCATE(PSNVA(NNI))
IF (HINMODEL == 'ECMWF') THEN
  CALL READ_GRIB(HGRIB,KLUOUT,32,IRET,ZFIELD)
  IF (IRET == 0 ) THEN
    IF (SIZE(ZFIELD).NE.NNI) THEN
      CALL ABOR1_SFX('MODE_READ_GRIB: SNOW ALB DIFFERS FROM OTHER FIELDS (READ_GRIB_SNOW_DEN)')
    END IF    
    IF (PRESENT(PLFRAC)) THEN
      WHERE (PLFRAC(:).LT.XRM_WM_ECMWF) ZFIELD(:) = XUNDEF
      PSNVA(:) = ZFIELD(:)
    ELSE
      WHERE (PMASK(:)/=1.) ZFIELD(:) = XUNDEF
      PSNVA(:) = ZFIELD(:) 
    END IF
    DEALLOCATE(ZFIELD)
  ELSE
    PSNVA(:) =  0.5 * ( XANSMIN + XANSMAX )
    WHERE (PMASK(:)/=1.) PSNVA(:) = XUNDEF 
  END IF

ELSEIF (HINMODEL == 'RACMO ') THEN
  ! read albedo
  DEALLOCATE(PSNVA)
  CALL READ_GRIB(HGRIB,KLUOUT,32,IRET,PSNVA)
  IF (IRET /= 0 ) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB: SNOW ALBEDO MISSING (READ_GRIB_SNOW_ALB)')
  END IF
!
  ! read snow amount
  CALL READ_GRIB(HGRIB,KLUOUT,141,IRET,ZFIELD)
  IF (IRET /= 0 ) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB: SNOW AMOUNT MISSING (READ_GRIB_SNOW_ALB)')
  END IF
!
  IF (SIZE(PMASK) /= SIZE(PSNVA)) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_SNOW_ALB: ARRAY SIZES OF MASK AND SNOW ALBEDO DIFFER')
  ENDIF
  IF (SIZE(ZFIELD) /= SIZE(PSNVA)) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_SNOW_ALB: ARRAY SIZES OF SNOW AMOUNT AND SNOW ALBEDO DIFFER')
  ENDIF
!
  ! Filter out sea points, and points without snow
  WHERE (PMASK(:)/=1. .OR. ZFIELD(:)<=0.)
    PSNVA(:) = XUNDEF
  END WHERE
  DEALLOCATE(ZFIELD)  
 
ELSE 
 PSNVA(:) = 0.5 * ( XANSMIN + XANSMAX )
 WHERE (PMASK(:)/=1.) PSNVA(:) = XUNDEF  
END IF 
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_SNOW_ALB',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_SNOW_ALB
!!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_SNOW_DEN(HGRIB,KLUOUT,HINMODEL,PMASK,PSNV,PLFRAC)
!     #######################
!!
!!    AUTHOR
!!    -------------
!!    C Ardilouze 08/2013 : possibility to read snow density (ERAI-land)
!!
!-------------------------------------------------------------------------------
!
USE MODD_GRID_GRIB,  ONLY : NNI
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_SNOW_PAR,   ONLY : XRHOSMAX
USE MODD_PREP_ISBA,  ONLY : XRM_WM_ECMWF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB    ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT   ! logical unit of output listing
CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL ! Grib originating model
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK    ! grib land mask
REAL, DIMENSION(:), POINTER         :: PSNV    ! field to initialize
REAL, DIMENSION(:),   INTENT(IN), OPTIONAL :: PLFRAC ! grib land fraction 
!
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)           :: IRET      ! return code
REAL, DIMENSION(:), POINTER       :: ZFIELD => NULL()    ! field to initialize 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_SNOW_DEN',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_SNOW_DEN: | Reading snow density'
!
ALLOCATE(PSNV(NNI))
IF (HINMODEL == 'ECMWF') THEN
  CALL READ_GRIB(HGRIB,KLUOUT,33,IRET,ZFIELD)
  IF (IRET == 0 ) THEN
    IF (SIZE(ZFIELD).NE.NNI) THEN
      CALL ABOR1_SFX('MODE_READ_GRIB: SNOW DENS DIFFERS FROM OTHER FIELDS (READ_GRIB_SNOW_DEN)')
    END IF    
    IF (PRESENT(PLFRAC)) THEN
      WHERE (PLFRAC(:).LT.XRM_WM_ECMWF) ZFIELD(:) = XUNDEF
      PSNV(:) = ZFIELD(:)
    ELSE
      WHERE (PMASK(:)/=1.) ZFIELD(:) = XUNDEF
      PSNV(:) = ZFIELD(:) 
    END IF
    DEALLOCATE(ZFIELD)
  ELSE
    PSNV(:) = XRHOSMAX
    WHERE (PMASK(:)/=1.) PSNV(:) = XUNDEF 
  END IF

ELSEIF (HINMODEL == 'RACMO ') THEN
  ! read snow density
  DEALLOCATE(PSNV)
  CALL READ_GRIB(HGRIB,KLUOUT,33,IRET,PSNV)
  IF (IRET /= 0 ) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB: SNOW DENSITY MISSING (READ_GRIB_SNOW_DEN)')
  END IF
!
  ! read snow amount
  CALL READ_GRIB(HGRIB,KLUOUT,141,IRET,ZFIELD)
  IF (IRET /= 0 ) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB: SNOW AMOUNT MISSING (READ_GRIB_SNOW_DEN)')
  END IF
!
  IF (SIZE(PMASK) /= SIZE(PSNV)) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_SNOW_DEN: ARRAY SIZES OF MASK AND SNOW DENSITY DIFFER')
  ENDIF
  IF (SIZE(ZFIELD) /= SIZE(PSNV)) THEN
    CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_SNOW_DEN: ARRAY SIZES OF SNOW AMOUNT AND SNOW DENSITY DIFFER')
  ENDIF
!
  ! Filter out sea points, and points without snow
  WHERE (PMASK(:)/=1. .OR. ZFIELD(:)<=0.)
    PSNV(:) = XUNDEF
  END WHERE
  DEALLOCATE(ZFIELD)
  
ELSE
 PSNV(:) = XRHOSMAX
 WHERE (PMASK(:)/=1.) PSNV(:) = XUNDEF  
END IF 
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_SNOW_DEN',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_SNOW_DEN
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_WR(HGRIB,KLUOUT,HINMODEL,PMASK,PWR)
!     #######################
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
CHARACTER(LEN=*),   INTENT(IN)        :: HGRIB     ! Grib file name
INTEGER,            INTENT(IN)        :: KLUOUT    ! logical unit of output listing
CHARACTER(LEN=6),   INTENT(IN)        :: HINMODEL  ! Grib originating model
REAL, DIMENSION(:), INTENT(IN)        :: PMASK     ! Grib land mask
REAL, DIMENSION(:), POINTER           :: PWR       ! field to initialize, skin reservoir content
!
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)               :: IRET      ! return code
REAL(KIND=JPHOOK)                       :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_WR',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_WR: | Reading skin reservoir content'
!
SELECT CASE(TRIM(HINMODEL))
  CASE('RACMO')
    CALL READ_GRIB(HGRIB,KLUOUT,198,IRET,PWR)
  CASE DEFAULT
    CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_WR: OPTION NOT SUPPORTED '//HINMODEL)
END SELECT
!
IF (IRET /= 0 ) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB: SKIN RESERVOIR CONTENT MISSING (READ_GRIB_WR)')
END IF
!
IF (SIZE(PMASK) /= SIZE(PWR)) THEN
  CALL ABOR1_SFX('MODE_READ_GRIB:READ_GRIB_WR: ARRAY SIZES OF MASK AND SKIN RESERVOIR CONTENT DIFFER')
ENDIF
!
SELECT CASE(TRIM(HINMODEL))
  CASE('RACMO')
    ! go from m. water equivalent to kg/m2
    PWR(:) = PWR(:) * 1000.
END SELECT
!
IF (SIZE(PMASK)==SIZE(PWR)) &
  WHERE (PMASK(:)/=1.) PWR = XUNDEF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_WR',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_WR
!-------------------------------------------------------------------
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_T_TEB(HGRIB,KLUOUT,HINMODEL,PTI,PMASK,PT,PD)
!     #######################
!
USE MODD_GRID_GRIB,  ONLY : NNI
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL,                 INTENT(IN)    :: PTI       ! internal temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PT        ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PD        ! normalized grid
!
!* local variables
!  ---------------
REAL, DIMENSION(:), POINTER       :: ZFIELD => NULL()    ! field to initialize
INTEGER                           :: JL         ! layer loop counter
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_T_TEB',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_T_TEB: | Reading temperature for buildings'
!
 CALL READ_GRIB_TS(HGRIB,KLUOUT,HINMODEL,PMASK,ZFIELD)
!
ALLOCATE(PT(SIZE(ZFIELD),3))
ALLOCATE(PD(SIZE(ZFIELD),3))
!
PT(:,1) = ZFIELD
DEALLOCATE(ZFIELD)
PD(:,1) = 0.
!
PT(:,2) = PTI
PD(:,2) = 0.5         ! deep temperature depth assumed at half of wall or roof
!
PT(:,3) = PTI
PD(:,3) = 1.          ! temperature at building interior
!
IF (SIZE(PMASK)==SIZE(PT,1)) THEN
  DO JL=1,SIZE(PT,2)
    WHERE (PMASK(:)/=1.) PT(:,JL) = XUNDEF
  END DO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_T_TEB',1,ZHOOK_HANDLE)
END SUBROUTINE READ_GRIB_T_TEB
!-------------------------------------------------------------------
!     #######################
      SUBROUTINE READ_GRIB_TF_TEB(HGRIB,KLUOUT,HINMODEL,PTI,PMASK,PTF,PD)
!     #######################
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
 CHARACTER(LEN=*),     INTENT(IN)    :: HGRIB     ! Grib file name
INTEGER,              INTENT(IN)    :: KLUOUT    ! logical unit of output listing
 CHARACTER(LEN=6),     INTENT(IN)    :: HINMODEL  ! Grib originating model
REAL,                 INTENT(IN)    :: PTI       ! internal temperature
REAL, DIMENSION(:),   INTENT(IN)    :: PMASK     ! grib land mask
REAL, DIMENSION(:,:), POINTER       :: PTF       ! field to initialize
REAL, DIMENSION(:,:), POINTER       :: PD        ! thickness of each layer
!
!* local variables
!  ---------------
INTEGER(KIND=kindOfInt)           :: IRET      ! return code
REAL,    DIMENSION(:), POINTER    :: ZFIELD => NULL()    ! field to read
INTEGER                           :: JL         ! layer loop counter
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!--------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TF_TEB',0,ZHOOK_HANDLE)
WRITE  (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_TF_TEB: | Reading temperature for building floor'
!
! 1.  Deep soil temperature
!     ---------------------
!
WRITE (KLUOUT,'(A)') 'MODE_READ_GRIB:READ_GRIB_TF_TEB: | Reading deep soil temperature'
!
 CALL READ_GRIB_T2_LAND(HGRIB,KLUOUT,HINMODEL,PMASK,ZFIELD)
!
ALLOCATE(PTF(SIZE(ZFIELD),3))
ALLOCATE(PD (SIZE(ZFIELD),3))
!
PTF(:,2) = ZFIELD(:)
PD (:,2) = 0.5           ! deep temperature depth assumed at half of the floor
!
DEALLOCATE(ZFIELD)
!
! 2.  level 1 is the internal building temperature
!     -----------------------
!
PTF(:,1) = PTI
PD (:,1) = 0.
!
! 3.  Assumes uniform temperature profile below
!     -----------------------------------------
!
PTF(:,3) = PTF(:,2)
PD (:,3) = 1.          ! deep temperature value
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_GRIB:READ_GRIB_TF_TEB',1,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
END SUBROUTINE READ_GRIB_TF_TEB
!-------------------------------------------------------------------
END MODULE MODE_READ_GRIB
