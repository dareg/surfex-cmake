!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!
!#############################################
SUBROUTINE WRITE_HEADER_FA (GCP, HGRID, PGRID_PAR, CFILETYPE, HWRITE)
!#############################################
!
!!    PURPOSE
!!    -------
!!    Create and write a header for an ARPEGE FA file
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!                           
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      A. Voldoire          Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2007
!!         F. Taillefer 06/2008 : add Gauss and Conf Proj cases
!!         B. Decharme  01/2009 : FA can be used only if NDIM_FULL >=289 in LATLON
!!         A. Alias     10/2010 : FA header modified
!!         R. El Khatib 30-Mar-2012 fanmsg with 2 arguments
!!         A. Mary      07/2018 : width of I and E zones
!!         A. Mary      05/2019 : Gauss KNOZPA != 0
!!         A. Mary      06/2019 : absolute value of RPK
!!         S. Riette    01/2021 : add cartesian support
!!         O. Vignes    12/2022 : spectral truncation factor
!!         A.NApoly     08/2024 : impose ITRONC to 0 as it is meaningless for Surfex grid point files 
!!                                 and it causes issues for cubic grids                            
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_GRID_CONF_PROJ_n,  ONLY : GRID_CONF_PROJ_t
!
USE MODD_IO_SURF_FA
!
USE MODD_CSTS,  ONLY : XPI 
!
USE MODE_GRIDTYPE_CONF_PROJ
USE MODE_GRIDTYPE_LONLAT_REG
USE MODE_GRIDTYPE_GAUSS
USE MODE_GRIDTYPE_CARTESIAN
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
USE MODI_ABOR1_SFX
!
USE MODI_IO_BUFF_CLEAN
!
USE MODI_GET_LUOUT
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
TYPE(GRID_CONF_PROJ_t),INTENT(INOUT) :: GCP
!
 CHARACTER(LEN=*), INTENT(IN) :: HGRID
REAL, DIMENSION(:), INTENT(IN) :: PGRID_PAR
!
 CHARACTER(LEN=3),    INTENT(IN)  :: HWRITE       ! 'PGD' : only physiographic fields are written
 CHARACTER(LEN=6),    INTENT(IN)  :: CFILETYPE    ! 'FA' could also be 'LFI' in future developments
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: IL
INTEGER :: ILON
INTEGER :: ILAT
!
REAL :: ZLONMIN
REAL :: ZLONMAX
REAL :: ZLATMIN
REAL :: ZLATMAX
!
REAL :: ZSLAPO
REAL :: ZCLOPO
REAL :: ZSLOPO
REAL :: ZCODIL
REAL :: ZPRPK
REAL :: ZBETA
!
REAL :: ZLAPO
REAL :: ZLOPO
!
REAL :: ZRAD
!
REAL, DIMENSION(:), ALLOCATABLE :: ZSINLA, ZAHYBR, ZBHYBR
!
REAL, DIMENSION(:), ALLOCATABLE :: ZLAT_XY, ZDX, ZDY
!
REAL, DIMENSION(0:1), PARAMETER :: ZNIVA = (/0.,0./)
!
REAL, DIMENSION(0:1), PARAMETER :: ZNIVB = (/0.,1./)
!
REAL, PARAMETER :: ZREFER = 101325.
!      
INTEGER, DIMENSION(11) :: IDATE
INTEGER, DIMENSION(:), ALLOCATABLE :: INLOPA, INOZPA
INTEGER :: ITYPTR
INTEGER :: INB ! number of articles in the file
INTEGER :: IRET
INTEGER :: ITRONC
INTEGER :: INLATI
INTEGER :: INXLON
INTEGER :: IWORK
INTEGER :: ICOUNT
INTEGER :: JLAT
INTEGER :: ILATE
INTEGER :: ILONE
INTEGER :: IWIDTH_I_X
INTEGER :: IWIDTH_I_Y
REAL    :: ZTRUNC
!
INTEGER :: ILUOUT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_HEADER_FA',0,ZHOOK_HANDLE)
!
#ifdef SFX_FA
!
 CALL IO_BUFF_CLEAN
!
ZRAD=XPI/180.0
!
ZSLAPO=0.0
ZCLOPO=0.0
ZSLOPO=0.0
ZCODIL=0.0
!
IF (HGRID=="CONF PROJ ") THEN
!
  CALL GET_GRIDTYPE_CONF_PROJ(PGRID_PAR,ZLAPO,ZLOPO,ZPRPK,ZBETA, &
                              ZLATMIN,ZLONMIN,ILON,ILAT          )
!
  ICOUNT=ILON*ILAT
  ALLOCATE(ZDX(ICOUNT))
  ALLOCATE(ZDY(ICOUNT))
!
  CALL GET_GRIDTYPE_CONF_PROJ(PGRID_PAR,PDX=ZDX,PDY=ZDY,KLONE=ILONE,KLATE=ILATE, &
                              KWIDTH_I_X=IWIDTH_I_X,KWIDTH_I_Y=IWIDTH_I_Y, &
                              PTRUNC=ZTRUNC)
!
  ALLOCATE(ZSINLA(18))
  ALLOCATE(INLOPA(8))
  ALLOCATE(INOZPA((1+ILAT)/2))
!
  ZSINLA(:)=0.0
  INLOPA(:)=0
  INOZPA(:)=0
!
  ZSINLA(1) = -1.0
  ZSINLA(2) = ABS(ZPRPK)
  ZSINLA(3) = ZLOPO*ZRAD
  ZSINLA(4) = ZLAPO*ZRAD
  ZSINLA(5) = GCP%XLONC*ZRAD
  ZSINLA(6) = GCP%XLATC*ZRAD  
  ZSINLA(7) = ZDX(1)
  ZSINLA(8) = ZDY(1)
  ZSINLA(9) = ZDX(1)*(ILON+ILONE)
  ZSINLA(10) = ZDY(1)*(ILAT+ILATE)
  ZSINLA(13) = 0.0
  ZSINLA(14) = 0.0
!
  INLOPA(1) = 10
  INLOPA(2) = -1
  INLOPA(3) = 1
  INLOPA(4) = ILON
  INLOPA(5) = 1
  INLOPA(6) = ILAT
  INLOPA(7) = IWIDTH_I_X
  INLOPA(8) = IWIDTH_I_Y
!
!
  INLATI = ILAT+ILATE
  INXLON = ILON+ILONE
  !ITRONC = INT(REAL(INLATI-2)/ZTRUNC)
  ITRONC=0 !ANTMPTEST : impose ITRONC to 0 as it is meaningless for Surfex grid point files
  ITYPTR = -INT(REAL(INXLON-2)/ZTRUNC)

  IF (ILONE == 0 .AND. ILATE == 0) THEN
    INLOPA(2) = 0
  ENDIF
!
ELSEIF (HGRID=="CARTESIAN ") THEN
!
!
  CALL GET_GRIDTYPE_CARTESIAN(PGRID_PAR,ZLAPO,ZLOPO,ILON,ILAT)  
!
  ICOUNT=ILON*ILAT
  ALLOCATE(ZDX(ICOUNT))
  ALLOCATE(ZDY(ICOUNT))
!
  CALL GET_GRIDTYPE_CARTESIAN(PGRID_PAR,PDX=ZDX,PDY=ZDY)
!
  ALLOCATE(ZSINLA(18))
  ALLOCATE(INLOPA(8))
  ALLOCATE(INOZPA((1+ILAT)/2))
!
  ZSINLA(:)=0.0
  INLOPA(:)=0
  INOZPA(:)=1
!
  !In FA header, 4 lon/lat are provided
  !Only one lon/lat is avaliable in GRIDTYPE_CARTESIAN
  ZSINLA (1)  = 0.
  ZSINLA (2)  = ZLOPO*ZRAD
  ZSINLA (3)  = ZLAPO*ZRAD
  ZSINLA (4)  = ZLOPO*ZRAD
  ZSINLA (5)  = ZLAPO*ZRAD
  ZSINLA (6)  = ZLOPO*ZRAD
  ZSINLA (7)  = ZLAPO*ZRAD
  ZSINLA (8)  = ZLOPO*ZRAD
  ZSINLA (9)  = ZLAPO*ZRAD
  ZSINLA (10) = 1.
  ZSINLA (11) = 1.
  ZSINLA (12) = 0.
  ZSINLA (13) = ZDX(1)*ILON
  ZSINLA (14) = ZDY(1)*ILAT
  ZSINLA (15) = ZDX(1)
  ZSINLA (16) = ZDY(1)
  ZSINLA (17) = 0.
  ZSINLA (18) = 0.
!
  INLOPA  (1)  = 0
  INLOPA  (2)  = 1
  INLOPA  (3)  = 1
  INLOPA  (4)  = 1
  INLOPA  (5)  = ILON
  INLOPA  (6)  = ILAT
  INLOPA  (7)  = 0
  INLOPA  (8)  = 0
!
  ITYPTR = -INT(REAL(ILON-1)/2.)
  !ITRONC = INT(REAL(ILAT-1)/2.)
  ITRONC=0 !ANTMPTEST : impose ITRONC to 0 as it is meaningless for Surfex grid point files
!
  INLATI = ILAT
  INXLON = ILON
!
  ZCODIL = -1.
!
  IF(ILON .NE. 1 .OR. ILAT .NE. 4)THEN
    !Code has been written by taking example on acadfa_sueframe subroutine
    !Some parts may be missing *or wrong* for other geometries than the 1D one
    CALL ABOR1_SFX('WRITE_HEADER_FA: CARTESIAN NOT YET FULLY IMPLEMENTED')
  ENDIF
!
ELSEIF (HGRID=="LONLAT REG") THEN
!
  CALL GET_GRIDTYPE_LONLAT_REG(PGRID_PAR,ZLONMIN,ZLONMAX, &
                                 ZLATMIN,ZLATMAX,ILON,ILAT  )  
!
  CALL GET_LUOUT(CFILETYPE,ILUOUT)
  IL=ILON*ILAT
  IF(IL<289)THEN
    WRITE(ILUOUT,*)' When Fa is used, NDIM_FULL must be >= 289, here NDIM_FULL = ',IL
    CALL ABOR1_SFX(' WRITE_HEADER_FA: LONLAT REG, With Fa, NDIM_FULL must be >= 289')
  ENDIF
!
  ALLOCATE(ZSINLA(18))
  ALLOCATE(INLOPA(8))
  ALLOCATE(INOZPA((1+ILAT)/2))
!
  !ITRONC= MIN(INT((REAL(ILAT-2)/2.0)),21)
  ITRONC=0 !ANTMPTEST : impose ITRONC to 0 as it is meaningless for Surfex grid point files
  ITYPTR=-MIN(INT((REAL(ILON-2)/2.0)),21)
  INLATI=ILAT
  INXLON=ILON
!
  ZSINLA(:)=0.
  INLOPA(:)=0
  INOZPA(:)=0
!
  ZSINLA(1) =-1.
  ZSINLA(2) =-9.
  ZSINLA(5) =(ZLONMIN+(ZLONMAX-ZLONMIN)/2.)*ZRAD
  ZSINLA(6) =(ZLATMIN+(ZLATMAX-ZLATMIN)/2.)*ZRAD
  ZSINLA(7) =((ZLONMAX-ZLONMIN)/REAL(ILON))*ZRAD
  ZSINLA(8) =((ZLATMAX-ZLATMIN)/REAL(ILAT))*ZRAD
  ZSINLA(9) =(ZLONMAX-ZLONMIN)*ZRAD
  ZSINLA(10)=(ZLATMAX-ZLATMIN)*ZRAD
  ZSINLA(13)=ZLONMIN*ZRAD
  ZSINLA(14)=ZLATMIN*ZRAD
  ZSINLA(15)=ZLONMAX*ZRAD
  ZSINLA(16)=ZLATMAX*ZRAD
!
  INLOPA(1) = 10
  INLOPA(2) = -1
  INLOPA(3) = 1
  INLOPA(4) = ILON
  INLOPA(5) = 1
  INLOPA(6) = ILAT
  INLOPA(7) = 8
  INLOPA(8) = 8
!
ELSEIF (HGRID=="GAUSS     ") THEN
!
  CALL GET_GRIDTYPE_GAUSS(PGRID_PAR,KNLATI=INLATI,KL=IL)
!
  ALLOCATE(INLOPA(INLATI))
  ALLOCATE(ZSINLA(INLATI))
  ALLOCATE(INOZPA(INLATI))
!
  ALLOCATE(ZLAT_XY(IL))
!
  CALL GET_GRIDTYPE_GAUSS(PGRID_PAR,PLAPO=ZLAPO,PLOPO=ZLOPO,          &
                            PCODIL=ZCODIL,KNLOPA=INLOPA,PLAT_XY=ZLAT_XY )  
!
! voir plus tard si ce parametre n'est pas deja dans un module !
  IF (ZLAPO>89.99 .AND. ABS(ZLOPO)<0.00001) THEN
    ITYPTR=1
  ELSE
    ITYPTR=2
  ENDIF
!
  ZSLAPO=SIN(ZLAPO*ZRAD)
  ZCLOPO=COS(ZLOPO*ZRAD)
  ZSLOPO=SIN(ZLOPO*ZRAD)
!
  IWORK = INT(REAL(INLATI)/2.0)
  INXLON=INLOPA(IWORK)
!
  !IF (ITYPTR==1) THEN
  !  ITRONC=INT(REAL(INXLON-1)/2.)
  !ELSE
  !  ITRONC=INT(REAL(INXLON-3)/2.)
  !ENDIF
!
  ITRONC=0 !ANTMPTEST : impose ITRONC to 0 as it is meaningless for Surfex grid point files
  INOZPA(:)=0
  DO JLAT = 1,INLATI
    INOZPA(JLAT) = INT(FLOOR(REAL(INLOPA(JLAT) - 1) / 2.)) ! 2:linear truncation
  ENDDO
  WHERE(INOZPA(:) == INOZPA(INLATI/2)) INOZPA(:) = ITRONC ! might be ITRONC + 1  
!
  ICOUNT=1
  DO JLAT = 1,INLATI
    ZSINLA(JLAT)=SIN(ZLAT_XY(ICOUNT)*ZRAD)
    ICOUNT=ICOUNT+INLOPA(JLAT)
  ENDDO
!
  DEALLOCATE(ZLAT_XY)
!
ELSEIF (HGRID=="IGN       ") THEN
!
  CALL ABOR1_SFX('WRITE_HEADER_FA: IGN NOT YET IMPLEMENTED')
!
ELSEIF (HGRID=="LONLATVAL ") THEN
!
  CALL ABOR1_SFX('WRITE_HEADER_FA: LONLATVAL NOT YET IMPLEMENTED')
!
END IF
!
ALLOCATE(ZAHYBR(0:1))
ALLOCATE(ZBHYBR(0:1))
ZAHYBR(0:1)=ZNIVA(0:1)
ZBHYBR(0:1)=ZNIVB(0:1)
!
! Reduce verbosity (in case it is not already done)
 CALL FANMSG(0,NLUOUT)
 CALL FACADE(CDNOMC,ITYPTR,ZSLAPO,ZCLOPO,ZSLOPO,ZCODIL,ITRONC,INLATI,INXLON, &
              INLOPA,INOZPA,ZSINLA,1,ZREFER,ZAHYBR,ZBHYBR,.TRUE.)  
!
 CALL FAITOU(IRET,NUNIT_FA,.TRUE.,CFILEOUT_FA,'UNKNOWN', &
              .TRUE.,.FALSE.,IVERBFA,0,INB,CDNOMC)  
!
IDATE(:)=0
IDATE(1)=1992
IDATE(2)=1
IDATE(3)=1
IDATE(6)=1
 CALL FANDAR(IRET,NUNIT_FA,IDATE)
!
DEALLOCATE(ZSINLA)
DEALLOCATE(INLOPA)
DEALLOCATE(INOZPA)
!
DEALLOCATE(ZAHYBR)
DEALLOCATE(ZBHYBR)
!
#endif
!
IF (LHOOK) CALL DR_HOOK('WRITE_HEADER_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_HEADER_FA
