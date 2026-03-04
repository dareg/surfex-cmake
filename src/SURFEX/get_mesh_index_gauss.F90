!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###############################################################
      SUBROUTINE GET_MESH_INDEX_GAUSS(KNBLINES,KSSO,PGRID_PAR,PLAT,PLON,&
                                      KINDEX,KISSOX,KISSOY,KFSSO,KFISSOX,KFISSOY, &
                                      PVALUE,PNODATA)
!     ###############################################################
!
!!**** *GET_MESH_INDEX_GAUSS* get the grid mesh where point (lat,lon) is located
!!
!!    PURPOSE
!!    -------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson         Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    12/09/95
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE EGGANGLES,ONLY : LOLA, ANGLE_DOMAIN
!
USE MODD_GET_MESH_INDEX_GAUSS, ONLY : XYINF, XYSUP, XXINF, XXSUP,       &
                                      NNLATI, NNLOPA, XLAPO, XLOPO, XCODIL, XSINTS,   &
                                      XLON, XLAT, XCOST, XSINTC, XCOSN, XSINN, XLONP, &
                                      XLATP, XCOSP, XSINP, XPI, X1, X2, XDR,          &
                                      NFRACDX, XYCEN, LROTSTRETCH
!
USE MODE_GRIDTYPE_GAUSS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
INTEGER,                       INTENT(IN)   :: KNBLINES
INTEGER,                       INTENT(IN)   :: KSSO        ! number of subgrid mesh in each direction
INTEGER,                       INTENT(IN)   :: KFSSO     ! number of fractional subgrid mesh in each direction
REAL,    DIMENSION(:),         INTENT(IN)   :: PGRID_PAR ! grid parameters
REAL,    DIMENSION(:),         INTENT(IN)   :: PLAT      ! latitude of the point  (degrees)
REAL,    DIMENSION(:),         INTENT(IN)   :: PLON      ! longitude of the point (degrees)
INTEGER, DIMENSION(:,:),       INTENT(OUT)  :: KINDEX    ! index of the grid mesh where the point is
INTEGER, DIMENSION(:,:),       INTENT(OUT)  :: KISSOX    ! X index of the subgrid mesh
INTEGER, DIMENSION(:,:),       INTENT(OUT)  :: KISSOY    ! Y index of the subgrid mesh
INTEGER, DIMENSION(:,:),       INTENT(OUT)  :: KFISSOX   ! X index of the fractional subgrid mesh
INTEGER, DIMENSION(:,:),       INTENT(OUT)  :: KFISSOY   ! Y index of the fractional subgrid mesh
!
REAL, DIMENSION(:), OPTIONAL, INTENT(IN)    :: PVALUE  ! value of the point to add
REAL, OPTIONAL, INTENT(IN) :: PNODATA

!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
REAL, DIMENSION(SIZE(PLAT))       :: ZX, ZY       ! pseudo longitude of input point
REAL, DIMENSION(SIZE(PLAT))       :: ZVALUE
REAL :: ZPC2
REAL :: ZNODATA
REAL :: ZSSOX, ZSSOY, ZSSO, ZFSSO
REAL :: ZRVRS90, RVRS360, ZNLATI
REAL(KIND=8) :: RVRSDP360
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: RDPNLOPA
!
INTEGER :: ILGRID, ISIZE_LON, ISIZE_DLAT, INBLINES  ! number of grid points
INTEGER :: JI, JL       ! loop counter in x
INTEGER :: IOFF, IINDEX, ILAT, ILON, IGUESS

LOGICAL :: LLREPRO_SINGLE_PGD_GAUSS ! .TRUE. to reproduce the results of the former algorithm in single precision.
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE2, ZHOOK_HANDLE4
!----------------------------------------------------------------------------
!
LLREPRO_SINGLE_PGD_GAUSS=(KIND(RVRSDP360)/=KIND(RVRS360))
!LLREPRO_SINGLE_PGD_GAUSS=.FALSE.

IF (LHOOK) CALL DR_HOOK('GET_MESH_INDEX_GAUSS_1',0,ZHOOK_HANDLE)
!
IF (PRESENT(PVALUE) .AND. PRESENT(PNODATA)) THEN
  ZVALUE(:) = PVALUE(:)
  ZNODATA = PNODATA
ELSE
  ZVALUE(:) = 1
  ZNODATA = 0
ENDIF
!
IF (KSSO /=0) ZSSO =FLOAT(KSSO)
IF (KFSSO/=0) ZFSSO=FLOAT(KFSSO)

IF (.NOT. ALLOCATED(NNLOPA)) THEN
  !
  IF (LHOOK) CALL DR_HOOK('GET_MESH_INDEX_GAUSS_2',0,ZHOOK_HANDLE2)

  !*    1.     Gets parameters of the projection
  !            ---------------------------------
  !
  XPI = 4.*ATAN(1.)
  XDR = XPI / 180.
  !
  CALL GET_GRIDTYPE_GAUSS(PGRID_PAR,NNLATI,KL=ILGRID)
  !
  ALLOCATE(NNLOPA(NNLATI))
  CALL GET_GRIDTYPE_GAUSS(PGRID_PAR,NNLATI,XLAPO,XLOPO,XCODIL,NNLOPA,ILGRID)
  !
  ALLOCATE(NFRACDX(NNLATI)) ! recycling : offset for each gaussian latitude
  NFRACDX(1)=0
  DO JI=2,NNLATI
    NFRACDX(JI)= NFRACDX(JI-1)+NNLOPA(JI-1)
  ENDDO

  ZPC2  = XCODIL*XCODIL
  X1 = 1.0-ZPC2
  X2 = 1.0+ZPC2
  !
  XLONP = ANGLE_DOMAIN(XLOPO,DOM='0+',UNIT='D') * XDR
  XLATP = XLAPO * XDR
  !
  XCOSP = COS(XLATP)
  XSINP = SIN(XLATP)
  !
  !*    2.     Limits of grid meshes in x and y
  !            --------------------------------
  !
  ALLOCATE(XXINF(ILGRID))
  ALLOCATE(XYINF(ILGRID))
  ALLOCATE(XXSUP(ILGRID))
  ALLOCATE(XYSUP(ILGRID))
  !
  CALL GAUSS_GRID_LIMITS(NNLATI,NNLOPA,XXINF,XXSUP,XYINF,XYSUP)

  ALLOCATE(XYCEN(NNLATI)) ! recycling : lower latitude bound for each latitude
  DO JI=1,NNLATI
    XYCEN(JI) = XYINF(NFRACDX(JI)+1)
  ENDDO

  !
  !*    3.     Find if rotated pole and/or stretching to improve CPU time
  !            ----------------------------------------------------------
  !
  LROTSTRETCH = .TRUE.
  IF (XCODIL==1.0.AND.XLAPO==90.0.AND.XLOPO==0.0) LROTSTRETCH = .FALSE.
  !
  IF (LHOOK) CALL DR_HOOK('GET_MESH_INDEX_GAUSS_2',1,ZHOOK_HANDLE2)

ENDIF
!
! case where the grid is not regular: all points are considered independently
IF (KNBLINES==0) THEN
  INBLINES = SIZE(PLAT)
ELSE
  INBLINES = KNBLINES
ENDIF
!
ISIZE_DLAT = SIZE(PLAT)/INBLINES
!
! case where the grid is not regular: all points are considered independently
IF (KNBLINES==0) THEN
  ISIZE_LON = INBLINES
ELSE
  ISIZE_LON = ISIZE_DLAT
ENDIF
!
IF (ALLOCATED(XLON)) THEN
  IF ( ISIZE_LON/=SIZE(XLON) .OR. INBLINES/=SIZE(XLAT) ) THEN
    DEALLOCATE(XLON)
    DEALLOCATE(XLAT)
    DEALLOCATE(XCOST)
    DEALLOCATE(XSINTC)
    DEALLOCATE(XSINTS)
    DEALLOCATE(XCOSN)
    DEALLOCATE(XSINN)
  ENDIF
ENDIF
!
IF (.NOT.ALLOCATED(XLON)) THEN
  !
  ALLOCATE(XLON(ISIZE_LON))
  ALLOCATE(XLAT(INBLINES))
  !
  ALLOCATE(XCOST (SIZE(XLAT)))
  ALLOCATE(XSINTC(SIZE(XLAT)))
  ALLOCATE(XSINTS(SIZE(XLAT)))
  ALLOCATE(XCOSN (SIZE(XLON)))
  ALLOCATE(XSINN (SIZE(XLON)))
  !
  XLON(:) = ANGLE_DOMAIN(PLON(1:ISIZE_LON),DOM='0+',UNIT='D') * XDR
  XCOSN(:) = COS(XLON(:)-XLONP)
  XSINN(:) = SIN(XLON(:)-XLONP)
  !
ENDIF
!
!
!*    3.     Projection of input points into pseudo coordinates
!            --------------------------------------------------
!
DO JI=1,INBLINES
  XLAT  (JI) = PLAT(JI*ISIZE_DLAT) * XDR
  XSINTC(JI) = SIN(XLAT(JI)) * XCOSP
  XSINTS(JI) = SIN(XLAT(JI)) * XSINP
  XCOST (JI) = COS(XLAT(JI))
ENDDO
!
IF (LROTSTRETCH) THEN
  CALL XY_GAUSS(XCODIL,ISIZE_DLAT,ISIZE_LON,ZNODATA,ZVALUE,ZY,ZX)
ELSE
  ZX(:) = PLON(:)
  ZY(:) = PLAT(:) 
ENDIF
!
IF (LHOOK) CALL DR_HOOK('GET_MESH_INDEX_GAUSS_1',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*    5.     Localisation of the data points on (x,y) grid
!            ---------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('GET_MESH_INDEX_GAUSS_4',0,ZHOOK_HANDLE4)

ZRVRS90=1./90.
ZNLATI=REAL(NNLATI)

IF (LLREPRO_SINGLE_PGD_GAUSS) THEN
  RVRSDP360=1._8/360._8
  ALLOCATE(RDPNLOPA(NNLATI))
  DO JI=1,NNLATI
    RDPNLOPA(JI)=REAL(NNLOPA(JI),KIND=8)
  ENDDO
ELSE
  RVRS360=1./360.
ENDIF

!$OMP PARALLEL DO PRIVATE(JL,IGUESS,JI,ILAT,IOFF,ILON,IINDEX,ZSSOX,ZSSOY)
DO JL=1,SIZE(PLAT)
  !
  IF (ZVALUE(JL)/=ZNODATA) THEN
  !
    ZY(JL)=MIN(MAX(ZY(JL),-90.),90.)

    IGUESS=MAX(1,NINT((ZNLATI*(1.-(ZY(JL)*ZRVRS90))*0.5))) ! rather smaller
    DO JI=IGUESS,NNLATI
      IF (ZY(JL)>=XYCEN(JI)) THEN
        ILAT = JI
        IOFF = NFRACDX(JI)
        EXIT
      ENDIF
    ENDDO
    !
    IF (ZX(JL) < 0.) ZX(JL)=ZX(JL)+360.

    IF (LLREPRO_SINGLE_PGD_GAUSS) THEN
      ILON=MOD(NINT((REAL(ZX(JL),KIND=8)*RDPNLOPA(ILAT))*RVRSDP360),NNLOPA(ILAT))+1
      IF (ZX(JL)==XXSUP(IOFF+ILON)) THEN
        ILON=ILON+1
        IF (ILON==NNLOPA(ILAT)) ILON=1
      ENDIF
    ELSE
      ILON=MOD(NINT((ZX(JL)*REAL(NNLOPA(ILAT)))*RVRS360),NNLOPA(ILAT))+1
    ENDIF
     
    IINDEX =  IOFF + ILON

  !*    6.     Localisation of the data points on in the subgrid of this mesh
  !            --------------------------------------------------------------
    IF (IINDEX/=0 .AND. (KSSO/=0 .OR. KFSSO/=0)) THEN
      IF (ZX(JL) > XXSUP(IINDEX)) THEN
        ZSSOX = (ZX(JL)-360.-XXINF(IINDEX))/(XXSUP(IINDEX)-XXINF(IINDEX))
      ELSE
        ZSSOX = (ZX(JL)-XXINF(IINDEX))/(XXSUP(IINDEX)-XXINF(IINDEX))
      ENDIF
      IF (ZY(JL) >= XYSUP(IINDEX)) THEN
        ! shift slightly to the south =>
        ZSSOY = (XYSUP(IINDEX)-1000.*EPSILON(XYSUP(IINDEX))-XYINF(IINDEX))/(XYSUP(IINDEX)-XYINF(IINDEX))
      ELSE
        ZSSOY = (ZY(JL)-XYINF(IINDEX))/(XYSUP(IINDEX)-XYINF(IINDEX))
      ENDIF
      IF (KSSO/=0) THEN
        KISSOX(1,JL) = 1 + INT( ZSSO * ZSSOX )
        KISSOY(1,JL) = 1 + INT( ZSSO * ZSSOY ) 
      ENDIF
      IF (KFSSO/=0) THEN
        KFISSOX(1,JL) = 1 + INT( ZFSSO * ZSSOX )
        KFISSOY(1,JL) = 1 + INT( ZFSSO * ZSSOY ) 
      ENDIF
    ENDIF

  ELSE

    IINDEX=0

  ENDIF
  !
  KINDEX(1,JL)=IINDEX
ENDDO
!$OMP END PARALLEL DO

IF (LHOOK) CALL DR_HOOK('GET_MESH_INDEX_GAUSS_4',1,ZHOOK_HANDLE4)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_MESH_INDEX_GAUSS
