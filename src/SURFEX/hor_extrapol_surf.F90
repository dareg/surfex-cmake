!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE HOR_EXTRAPOL_SURF(KLUOUT,HCOORTYPE, KILEN,PILA1,PILA2,PILO1,PILO2,&
                                   KINLA,KINLO,KP,PFIELD_IN,PLAT,PLON,PFIELD,OINTERP,&
                                   PILATARRAY)
!     ###################################################################
!
!!**** *HOR_EXTRAPOL_SURF* extrapolate a surface field
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!       For each point to interpolate, the nearest valid point value is set.
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
!!
!!    V. Masson          Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     01/12/98
!!     V. Masson    01/2004 extrapolation in latitude and longitude
!!     M. Jidane    11/2013 add OpenMP directives
!!     A. Napoly    10/2022 add OpenMP directives and optimisations in loops
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NPIO, NCOMM, LSFX_MPI
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_CSTS,       ONLY : XPI
USE MODN_PREP_SURF_ATM, ONLY : NHALO_PREP
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
!*    0.1    Declaration of arguments
!            ------------------------
!

INTEGER,               INTENT(IN)     :: KLUOUT   ! output listing logical unit
CHARACTER(LEN=4),      INTENT(IN)     :: HCOORTYPE! type of coordinate
INTEGER,               INTENT(IN)     :: KILEN    ! size of input grid ("NGPTOTG")
REAL,                      INTENT(IN)  :: PILA1   ! Lat. (y) of first input point ("NDGSA=1")
REAL,                      INTENT(IN)  :: PILA2   ! Lat. (y) of last  input point ("NDGEN=NDGLG")
REAL,                      INTENT(IN)  :: PILO1   ! Lon. (x) of first input point
REAL,                      INTENT(IN)  :: PILO2   ! Lon. (x) of last  input point
INTEGER,                   INTENT(IN)  :: KINLA   ! Number of parallels of input grid ("NDGLG")
INTEGER, DIMENSION(:), INTENT(IN)      :: KINLO   ! Number of point along a parallel of input grid ("NLOENG")
INTEGER, DIMENSION(:,:), INTENT(IN)    :: KP
REAL,   DIMENSION(:,:),  INTENT(IN), TARGET    :: PFIELD_IN! input field on grid mesh
REAL,   DIMENSION(:),  INTENT(IN)      :: PLAT     ! latitude of output points
REAL,   DIMENSION(:),  INTENT(IN)      :: PLON     ! longitude of output points
REAL,   DIMENSION(:,:),  INTENT(INOUT) :: PFIELD   ! output field
LOGICAL,DIMENSION(:),  INTENT(IN)      :: OINTERP  ! .true. where a physical value is present
REAL, DIMENSION(KINLA), INTENT(IN), OPTIONAL :: PILATARRAY ! input latitudes (if not regular)
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL :: ZLAT  ! latitude of point to define
REAL :: ZLON  ! longitude of point to define
REAL :: ZDIST ! current distance to valid point (in lat/lon grid)
REAL, DIMENSION(SIZE(PFIELD,2)) :: ZNDIST ! smallest distance to valid point for each layer
REAL :: ZCOSLA! cosine of latitude
REAL :: ZLONSC! longitude of valid point
REAL :: ZRAD ! conversion degrees to radians
!
INTEGER :: ID0, IST, II
INTEGER :: INFOMPI
INTEGER  :: JI, JL, JLAT, JLON, JRANK   ! loop index on points
#ifdef SFX_MPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
LOGICAL  :: GLALO ! flag true is second coordinate is a longitude or pseudo-lon.
                  !      false if metric coordinates

REAL, DIMENSION(KINLA) :: ZIDLA ! distance between two consecutive latitudes
REAL, DIMENSION(KINLA) :: ZSTLAT ! start latitude value of each latitude
INTEGER, DIMENSION(KINLA) :: IPOS ! start address of each latitude in input grid
REAL, DIMENSION(KINLA) :: ZIDLOA ! distance between two consecutive longitudes
REAL :: ZIDLO, ZIDLOMAX, ZIDLOMIN, ZIDLAMAX, ZIDLAMIN
REAL, DIMENSION(KILEN) :: ZLA       ! input "latitude"  coordinate
REAL, DIMENSION(KILEN) :: ZLO       ! input "longitude" coordinate

INTEGER :: ICPT ! actual number of points to extrapolate
INTEGER :: INO ! Number of gripoint in output array, or maximum number of points to extrapolate
INTEGER :: INL ! Number of field layers
INTEGER, DIMENSION(SIZE(PFIELD,1)) :: IADDR ! addresses of the actual points to extrapolate (overdimensionned)
! ZTLONMIN, ZTLONMAX, ZTLATMIN, ZTLATMAX contain for each point to extrapolate
! the limits of the domain where to search for the valid points, according to NHALO_PREP:
REAL, DIMENSION(:), ALLOCATABLE :: ZTLONMIN, ZTLONMAX, ZTLATMIN, ZTLATMAX
INTEGER, DIMENSION(:), ALLOCATABLE :: IMASKP
REAL :: ZLOMIN, ZLOMAX, ZLAMIN, ZLAMAX, ZLONWIDTH, ZLATWIDTH

REAL, DIMENSION(:,:), ALLOCATABLE, TARGET :: ZBUF ! input field on grid mesh recieved from NPIO task
REAL, DIMENSION(:,:), POINTER :: ZFIELD_IN ! pointer on input field on grid mesh
INTEGER :: ISIZE_IN ! size of input field on grid mesh
INTEGER :: ITAG ! communication tag
INTEGER :: IRCV ! recv request handler
INTEGER, DIMENSION(0:NPROC-1) :: ISND ! send requests handlers
REAL :: ZZLO,ZZLA ! longitude and latitude in cache

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE1
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HOR_EXTRAPOL_SURF',0,ZHOOK_HANDLE)
!
INL=SIZE(PFIELD,2)
ISIZE_IN=KILEN*INL
ITAG=1

! Distribute input grid data:
IF (NPROC > 1) THEN
#ifdef SFX_MPI    
  IF (LSFX_MPI) THEN
    IF (NRANK==NPIO) THEN
      DO JRANK=0,NPROC-1
        IF (JRANK/=NPIO) CALL MPI_ISEND(PFIELD_IN,ISIZE_IN*KIND(PFIELD_IN)/4,MPI_REAL,JRANK,ITAG,NCOMM,ISND(JRANK),INFOMPI)
      ENDDO
    ELSE
      ALLOCATE(ZBUF(KILEN,INL))
      CALL MPI_IRECV(ZBUF,ISIZE_IN*KIND(ZBUF)/4,MPI_REAL,NPIO,ITAG,NCOMM,IRCV,INFOMPI)
    ENDIF
  ENDIF
#endif
ENDIF

!-------------------------------------------------------------------------------
!
! Compute the distance between two consecutive latitudes:
IF (PRESENT(PILATARRAY)) THEN
  ZIDLA(1) = 0.
  DO JLAT=2,KINLA
    ZIDLA(JLAT) = PILATARRAY(JLAT) - PILATARRAY(JLAT-1)
  ENDDO
ELSE
  DO JLAT=1,KINLA
    ZIDLA(JLAT) = (PILA2 - PILA1) / (KINLA - 1)
  ENDDO
ENDIF
ZIDLAMAX = MAXVAL(ABS(ZIDLA(:)))
ZIDLAMIN = MINVAL(ABS(ZIDLA(2:KINLA)))
!
! Compute start latitude value of each latitude:
! Compute start address of each latitude in input grid
ZSTLAT(1) = PILA1
IPOS(1)=0
DO JLAT=2,KINLA
  ZSTLAT(JLAT) = ZSTLAT(JLAT-1) + ZIDLA(JLAT)
  IPOS(JLAT)=IPOS(JLAT-1)+KINLO(JLAT-1)
END DO

! Compute the distance between two consecutive longitude:
GLALO = HCOORTYPE=='LALO'
IF (GLALO) THEN
  DO JLAT=1,KINLA
    ZIDLOA(JLAT) = (PILO2-PILO1) / KINLO(JLAT)
  END DO
ELSE
  DO JLAT=1,KINLA
    ZIDLOA(JLAT) = (PILO2-PILO1) / (KINLO(JLAT)-1)
  END DO
ENDIF
ZIDLOMAX = 0.
ZIDLOMIN = XUNDEF
DO JLAT=1,KINLA
  ZIDLO = ABS(ZIDLOA(JLAT))
  IF (ZIDLO>ZIDLOMAX) ZIDLOMAX = ZIDLO
  IF (ZIDLO<ZIDLOMIN) ZIDLOMIN = ZIDLO
END DO

!$OMP PARALLEL DO PRIVATE(IST,JLON)
DO JLAT=1,KINLA
  IST=IPOS(JLAT)
  DO JLON=1,KINLO(JLAT)
    ZLA(IST+JLON) = ZSTLAT(JLAT)
    ZLO(IST+JLON) = PILO1 + (JLON-1) * ZIDLOA(JLAT) 
  END DO
END DO
!$OMP END PARALLEL DO

!-------------------------------------------------------------------------------
!
!*      4.   Loop on points to define
!            ------------------------
!
INO = SIZE(PFIELD,1)

ALLOCATE(ZTLONMIN(INO),ZTLONMAX(INO),ZTLATMIN(INO),ZTLATMAX(INO))
!
ZLOMIN=MINVAL(ZLO(:))
ZLOMAX=MAXVAL(ZLO(:))
ZLAMIN=MINVAL(ZLA(:))
ZLAMAX=MAXVAL(ZLA(:))
ZLONWIDTH=ZIDLOMAX*NHALO_PREP
ZLATWIDTH=ZIDLAMAX*NHALO_PREP
!
ALLOCATE(IMASKP(SIZE(KP,2)))

ICPT = 0
DO JI=1,INO
  IF (ALL(PFIELD(JI,:)/=XUNDEF)) CYCLE
  IF (.NOT. OINTERP(JI))  CYCLE
      ICPT = ICPT + 1
      IADDR(ICPT) = JI
      IF (NHALO_PREP>0) THEN
        IMASKP(:)=KP(JI,:)
        ZTLONMIN(JI) = MAX(ZLOMIN,MINVAL(ZLO(IMASKP))-ZLONWIDTH)
        ZTLONMAX(JI) = MIN(ZLOMAX,MAXVAL(ZLO(IMASKP))+ZLONWIDTH)
        ZTLATMIN(JI) = MAX(ZLAMIN,MINVAL(ZLA(IMASKP))-ZLATWIDTH)
        ZTLATMAX(JI) = MIN(ZLAMAX,MAXVAL(ZLA(IMASKP))+ZLATWIDTH)
      ELSE
        ZTLONMIN(JI) = ZLOMIN
        ZTLONMAX(JI) = ZLOMAX
        ZTLATMIN(JI) = ZLAMIN
        ZTLATMAX(JI) = ZLAMAX
      ENDIF
ENDDO
DEALLOCATE(IMASKP)

!2: loop on the points 

! Receive input grid data:

IF (NPROC > 1) THEN
#ifdef SFX_MPI    
  IF (LSFX_MPI) THEN
    IF (NRANK/=NPIO) THEN
      CALL MPI_WAIT(IRCV,ISTATUS,INFOMPI)
      ZFIELD_IN => ZBUF
    ELSE
      ZFIELD_IN => PFIELD_IN
    ENDIF
  ENDIF
#endif
ELSE
  ZFIELD_IN => PFIELD_IN
ENDIF

! Compute:

IF (ICPT /= 0) THEN
ZRAD=XPI/180.0
!$OMP PARALLEL PRIVATE(JI,II,ZNDIST,ZZLO,ZZLA,ZCOSLA,JLAT,ZLAT,JLON,ZLON,ID0,ZLONSC,ZDIST,JL)
!$OMP DO SCHEDULE(DYNAMIC,1)
  DO JI=1,ICPT
    II=IADDR(JI)
    ZNDIST(:) = XUNDEF
    ZZLO=PLON(II)
    ZZLA=PLAT(II)
    IF (GLALO) ZCOSLA=COS(ZZLA*ZRAD)
    DO JLAT = 1,KINLA
      ZLAT = ZSTLAT(JLAT)
      IF (ZLAT>=ZTLATMIN(II) .AND. ZLAT<=ZTLATMAX(II)) THEN
        DO JLON = 1,KINLO(JLAT)
          ZLON = ZLO(IPOS(JLAT)+JLON)
          IF (ZLON>=ZTLONMIN(II) .AND. ZLON<=ZTLONMAX(II)) THEN
            !index in the whole grid of the point used to interpolate
            ID0 = IPOS(JLAT) + JLON
            IF (ANY(ZFIELD_IN(ID0,:)/=XUNDEF)) THEN
              ZLONSC = ZLO(ID0)
              IF (GLALO) THEN
                IF (ZLONSC-ZZLO> 180.) ZLONSC = ZLONSC - 360.
                IF (ZLONSC-ZZLO<-180.) ZLONSC = ZLONSC + 360.
                ZDIST= (ZLA(ID0)-ZZLA) ** 2 + ((ZLONSC-ZZLO)*ZCOSLA) ** 2
              ELSE
                ZDIST= (ZLA(ID0)-ZZLA) ** 2 + (ZLONSC-ZZLO) ** 2
              END IF
              DO JL=1,INL
                IF (ZDIST<=ZNDIST(JL)) THEN
                  IF (ZFIELD_IN(ID0,JL)/=XUNDEF) THEN
                    PFIELD(II,JL) = ZFIELD_IN(ID0,JL)
                    ZNDIST(JL) = ZDIST
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
ENDIF

! Check validity of computation:

DO JL=1,INL
  IF (ANY(PFIELD(:,JL)==XUNDEF .AND. OINTERP(:))) THEN
    WRITE(*,*) 'LAYER ',JL,': NO EXTRAPOLATION : INCREASE YOUR HALO_PREP IN NAM_PREP_SURF_ATM'
    CALL ABOR1_SFX('NO EXTRAPOLATION : INCREASE YOUR HALO_PREP IN NAM_PREP_SURF_ATM')
  ENDIF
  WHERE (.NOT. OINTERP(:)) PFIELD(:,JL) = XUNDEF
ENDDO

! Finalize communications:

IF (NPROC > 1) THEN
#ifdef SFX_MPI    
  IF (LSFX_MPI) THEN
    IF (NRANK==NPIO) THEN
      DO JRANK=0,NPROC-1
        IF (JRANK/=NPIO) CALL MPI_WAIT(ISND(JRANK),ISTATUS,INFOMPI)
      ENDDO
    ENDIF
  ENDIF
#endif
ENDIF

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HOR_EXTRAPOL_SURF',1,ZHOOK_HANDLE)
!
END SUBROUTINE HOR_EXTRAPOL_SURF
