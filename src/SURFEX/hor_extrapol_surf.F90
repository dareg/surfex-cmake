!     #########
      SUBROUTINE HOR_EXTRAPOL_SURF(KLUOUT,HCOORTYPE,KILEN,PILA1,PILA2,PILO1,PILO2,&
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
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NPIO, NCOMM, IDX_I
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_CSTS,       ONLY : XPI
USE MODN_PREP_SURF_ATM, ONLY : NHALO_PREP
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*    0.1    Declaration of arguments
!            ------------------------
!
INTEGER,               INTENT(IN)     :: KLUOUT   ! output listing logical unit
 CHARACTER(LEN=4),      INTENT(IN)     :: HCOORTYPE! type of coordinate
 INTEGER, INTENT(IN) :: KILEN
REAL, INTENT(IN) :: PILA1
REAL, INTENT(IN) :: PILA2
REAL, INTENT(IN) :: PILO1
REAL, INTENT(IN) :: PILO2
INTEGER, INTENT(IN) :: KINLA
INTEGER, DIMENSION(:), INTENT(IN) :: KINLO
INTEGER, DIMENSION(:,:), INTENT(IN) :: KP
REAL,   DIMENSION(:,:),  INTENT(IN)     :: PFIELD_IN! input field on grid mesh
REAL,   DIMENSION(:),  INTENT(IN)     :: PLAT     ! latitude of each grid mesh.
REAL,   DIMENSION(:),  INTENT(IN)     :: PLON     ! longitude of each grid mesh.
REAL,   DIMENSION(:,:),  INTENT(INOUT)  :: PFIELD   ! field on grid mesh
LOGICAL,DIMENSION(:),  INTENT(IN)     :: OINTERP  ! .true. where physical value is needed
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PILATARRAY
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL, DIMENSION(:), ALLOCATABLE :: ZTLONMIN, ZTLONMAX, ZTLATMIN, ZTLATMAX
REAL, DIMENSION(:,:), ALLOCATABLE :: ZFIELD
REAL :: ZLAT  ! latitude of point to define
REAL :: ZLON  ! longitude of point to define
REAL :: ZDIST ! current distance to valid point (in lat/lon grid)
REAL :: ZNDIST! smallest distance to valid point
REAL :: ZCOSLA! cosine of latitude
REAL :: ZLONSC! longitude of valid point
REAL :: ZIDLA, ZIDLO, ZIDLOMAX, ZIDLOMIN
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCOOR, ZCO2
REAL,    DIMENSION(:), ALLOCATABLE :: ZLA       ! input "latitude"  coordinate
REAL,    DIMENSION(:), ALLOCATABLE :: ZLO       ! input "longitude" coordinate
REAL(KIND=JPRB) :: ZRAD ! conversion degrees to radians
!
INTEGER, DIMENSION(:), ALLOCATABLE :: IMASK, IMASKR
INTEGER, DIMENSION(:,:), ALLOCATABLE :: IVAL_EXT, ITAB
INTEGER, DIMENSION(NPROC) :: INO_TAB
INTEGER  :: INO     ! output array size
INTEGER, DIMENSION(2) :: ITSIZE, ITDIM
INTEGER, DIMENSION(2,0:NPROC-1) :: IBOR
INTEGER :: ISIZE, ISIZE_MAX, J, ID0, ICOMPT, ICPT
INTEGER :: INFOMPI, IDX, INL
INTEGER  :: JI, JL, JLAT, JLON, JIPOS   ! loop index on points
INTEGER  :: JISC  ! loop index on valid points
#ifndef NOMPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
LOGICAL  :: GLALO ! flag true is second coordinate is a longitude or pseudo-lon.
                  !      false if metric coordinates
REAL(KIND=JPRB) :: ZHOOK_HANDLE, ZHOOK_HANDLE_OMP
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HOR_EXTRAPOL_SURF',0,ZHOOK_HANDLE)
!
INO = SIZE(PFIELD,1)
INL = SIZE(PFIELD,2)
!
!-------------------------------------------------------------------------------
!
GLALO = HCOORTYPE=='LALO'
!
!-------------------------------------------------------------------------------
!
ALLOCATE(ZLA(KILEN))
ALLOCATE(ZLO(KILEN))
!
ZIDLOMAX = 0.
ZIDLOMIN = XUNDEF
JIPOS = 0
IF (PRESENT(PILATARRAY)) THEN
  ZIDLA = PILATARRAY(2) - PILATARRAY(1)
ELSE
  ZIDLA = (PILA2 - PILA1) / (KINLA - 1)
ENDIF
!
DO JLAT=1,KINLA
  IF (GLALO) THEN
    ZIDLO = (PILO2-PILO1) / KINLO(JLAT)
  ELSE
    ZIDLO = (PILO2-PILO1) / (KINLO(JLAT)-1)
  ENDIF
  DO JLON=1,KINLO(JLAT)
    JIPOS = JIPOS + 1
    ZLA(JIPOS) = PILA1 + (JLAT-1) * ZIDLA
    ZLO(JIPOS) = PILO1 + (JLON-1) * ZIDLO 
  END DO
  ZIDLO = ABS(ZIDLO)
  IF (ZIDLO>ZIDLOMAX) ZIDLOMAX = ZIDLO
  IF (ZIDLO<ZIDLOMIN) ZIDLOMIN = ZIDLO
END DO
!
!-------------------------------------------------------------------------------
!
!*      4.   Loop on points to define
!            ------------------------
!
ALLOCATE(ZTLONMIN(INO),ZTLONMAX(INO),ZTLATMIN(INO),ZTLATMAX(INO))
ZTLONMIN(:) = 0.
ZTLONMAX(:) = 0.
ZTLATMIN(:) = 0.
ZTLATMAX(:) = 0.
!
ZRAD=XPI/180.0_JPRB
!
!1: ZTLONMIN, ZTLONMAX, ZTLATMIN, ZTLATMAX contain for each point to extrapol 
! the limits of the domain where to search for the valid points, according
! to NHALO_PREP
ICPT = 0
ISIZE_MAX = 0
ISIZE = 0
DO JI=1,INO
  !
  IF (ALL(PFIELD(JI,:)/=XUNDEF)) CYCLE
  IF (.NOT. OINTERP(JI))  CYCLE
  ICPT = ICPT + 1
  !
  ZTLONMIN(JI) = MINVAL(ZLO(KP(JI,:)))-ZIDLOMAX*NHALO_PREP
  ZTLONMAX(JI) = MAXVAL(ZLO(KP(JI,:)))+ZIDLOMAX*NHALO_PREP
  ZTLATMIN(JI) = MINVAL(ZLA(KP(JI,:)))-ABS(ZIDLA)*NHALO_PREP
  ZTLATMAX(JI) = MAXVAL(ZLA(KP(JI,:)))+ABS(ZIDLA)*NHALO_PREP
  ISIZE = CEILING((ZTLONMAX(JI)-ZTLONMIN(JI)+1)/ZIDLOMIN)*&
          CEILING((ZTLATMAX(JI)-ZTLATMIN(JI)+1)/ABS(ZIDLA))
  IF (ISIZE>ISIZE_MAX) ISIZE_MAX = ISIZE
  !
ENDDO
!
ALLOCATE(IMASK(ICPT))
IMASK(:) = 0
ALLOCATE(IVAL_EXT(ICPT,ISIZE_MAX))
IVAL_EXT(:,:) = 0
ALLOCATE(ZCOOR(ICPT,2))
!
!number of points to extrapol and maximum number of points to sound
ITSIZE(1) = ICPT
ITSIZE(2) = ISIZE_MAX
!
ICPT = 0
!2: loop on the points 
DO JI=1,INO
  !
  IF (ALL(PFIELD(JI,:)/=XUNDEF)) CYCLE
  IF (.NOT. OINTERP(JI))  CYCLE
  !
  ICPT = ICPT + 1
  !imask associated the number of the point to extrapolate to its real index in
  !the field
  IMASK(ICPT) = JI
  !
  ICOMPT = 0
  JISC = 0
  !
  !coordinates of the point in the grid 
  ZCOOR(ICPT,1) = PLAT(JI)
  ZCOOR(ICPT,2) = PLON(JI)
  !
  !loop on the grid
  DO JLAT = 1,KINLA
    ZLAT = PILA1 + (JLAT-1) * ZIDLA
    IF (ZLAT>=ZTLATMIN(JI) .AND. ZLAT<=ZTLATMAX(JI)) THEN
      IF (GLALO) THEN
        ZIDLO = (PILO2-PILO1) / KINLO(JLAT)
      ELSE
        ZIDLO = (PILO2-PILO1) / (KINLO(JLAT)-1)
      ENDIF
      IF (JLAT>1) JISC = SUM(KINLO(1:JLAT-1))
      DO JLON = 1,KINLO(JLAT)
        ZLON = PILO1 + (JLON-1) * ZIDLO
        IF (ZLON>=ZTLONMIN(JI) .AND. ZLON<=ZTLONMAX(JI)) THEN
          ICOMPT = ICOMPT + 1
          !ival_ext contains the indexes of the points needed to interpolate
          !in the complete grid
          IVAL_EXT(ICPT,ICOMPT) = JISC + JLON
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  !
ENDDO
!
!NPIO knows the numbers of points to extrapolate for all tasks
IF (NPROC>1) THEN
#ifndef NOMPI
  CALL MPI_GATHER(ITSIZE,2*KIND(ITSIZE)/4,MPI_INTEGER,&
                  IBOR,2*KIND(IBOR)/4,MPI_INTEGER,& 
                  NPIO,NCOMM,INFOMPI)
#endif
ELSE
  IBOR(:,0) = ITSIZE(:)
ENDIF
!
!
IF (NRANK/=NPIO) THEN

  !if some points of this task need to be extrapolated
  IF (SUM(ITSIZE)/=0) THEN

    !itab contains indexes and mask
    ALLOCATE(ITAB(ITSIZE(1),ITSIZE(2)+1))
    ITAB(:,1:ITSIZE(2)) = IVAL_EXT(:,:)
    ITAB(:,ITSIZE(2)+1) = IMASK(:)
    !zfield will contain the values of the field
    ALLOCATE(ZFIELD(ITSIZE(1),INL))

    IDX_I = IDX_I + 1
    !sends indexes to npio
    CALL MPI_SEND(ITAB,SIZE(ITAB)*KIND(ITAB)/4,MPI_INTEGER,NPIO,IDX_I,NCOMM,INFOMPI)
    DEALLOCATE(ITAB)

    IDX_I = IDX_I + 1
    !send coords of the points to extrapolate
    CALL MPI_SEND(ZCOOR,SIZE(ZCOOR)*KIND(ZCOOR)/4,MPI_REAL,NPIO,IDX_I,NCOMM,INFOMPI)

    IDX_I = IDX_I + 1
    !receives values of the field from NPIO
    CALL MPI_RECV(ZFIELD,SIZE(ZFIELD)*KIND(ZFIELD)/4,MPI_REAL,NPIO,IDX_I,NCOMM,ISTATUS,INFOMPI)

    DO JI=1,ITSIZE(1)
      PFIELD(IMASK(JI),:) = ZFIELD(JI,:)
    ENDDO
    DEALLOCATE(ZFIELD)

  ELSE
    IDX_I = IDX_I + 3
  ENDIF

ELSE
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(J,ITDIM,ITAB,ZFIELD, &
!$OMP             JI,ZNDIST,ZCOSLA,JISC,ID0,IDX,ZLONSC,ZDIST,ZCO2,ZHOOK_HANDLE_OMP)
  DO J=0,NPROC-1

    IF (SUM(IBOR(:,J))/=0) THEN

      ALLOCATE(ITAB(IBOR(1,J),IBOR(2,J)+1))
      ALLOCATE(ZCO2(IBOR(1,J),2))
      ALLOCATE(ZFIELD(IBOR(1,J),INL))
      ZFIELD=XUNDEF

      IF (J/=NPIO) THEN

        !receives indexes and coordinaites
        CALL MPI_RECV(ITAB,SIZE(ITAB)*KIND(ITAB)/4,MPI_INTEGER,J,IDX_I+1,NCOMM,ISTATUS,INFOMPI)
        CALL MPI_RECV(ZCO2,SIZE(ZCO2)*KIND(ZCO2)/4,MPI_REAL,J,IDX_I+2,NCOMM,ISTATUS,INFOMPI)

      ELSE

        ITAB(:,1:IBOR(2,J)) = IVAL_EXT(:,:)
        ITAB(:,IBOR(2,J)+1) = IMASK(:)
        ZCO2(:,:) = ZCOOR(:,:)

      ENDIF
    
      DO JL=1,INL
        DO JI=1,IBOR(1,J)
          ZNDIST=XUNDEF
          IDX = IBOR(2,J)+1
          ZCOSLA=COS(ZCO2(JI,1)*ZRAD)
          DO JISC = 1,IBOR(2,J)
            !index in the whole grid of the point used to interpolate
            ID0 = ITAB(JI,JISC)
            IF (ID0==0) EXIT
            IF (PFIELD_IN(ID0,JL)/=XUNDEF) THEN
              ZLONSC = ZLO(ID0)
              IF (GLALO) THEN
                IF (ZLONSC-ZCO2(JI,2)> 180.) ZLONSC = ZLONSC - 360.
                IF (ZLONSC-ZCO2(JI,2)<-180.) ZLONSC = ZLONSC + 360.
                ZDIST= (ZLA(ID0)-ZCO2(JI,1)) ** 2 + ((ZLONSC-ZCO2(JI,2))*ZCOSLA) ** 2
              ELSE
                ZDIST= (ZLA(ID0)-ZCO2(JI,1)) ** 2 + (ZLONSC-ZCO2(JI,2)) ** 2
              END IF
              IF (ZDIST<=ZNDIST) THEN
                ZFIELD(JI,JL) = PFIELD_IN(ID0,JL)
                ZNDIST = ZDIST
              ENDIF
            ENDIF
          END DO   
        ENDDO
      ENDDO
      !
      IF (J/=NPIO) THEN
        !send values found to extrapolate
        CALL MPI_SEND(ZFIELD,SIZE(ZFIELD)*KIND(ZFIELD)/4,MPI_REAL,J,IDX_I+3,NCOMM,INFOMPI)
      ELSE
        DO JI = 1,IBOR(1,J)
          PFIELD(IMASK(JI),:) = ZFIELD(JI,:)
        ENDDO
      ENDIF            
      !
      DEALLOCATE(ZCO2)
      DEALLOCATE(ITAB)
      DEALLOCATE(ZFIELD)
      !
    ENDIF
    !
  ENDDO
!$OMP END PARALLEL DO
!$OMP BARRIER
  !
  !
  IDX_I = IDX_I + 3
  !
ENDIF
!
DEALLOCATE(ZCOOR)
IF (ALLOCATED(ZLA)) DEALLOCATE(ZLA)
IF (ALLOCATED(ZLO)) DEALLOCATE(ZLO)
DEALLOCATE(IVAL_EXT)
DEALLOCATE(ZTLONMIN,ZTLONMAX,ZTLATMIN,ZTLATMAX)
DEALLOCATE(IMASK)
!

DO JL=1,INL
  IF (ANY(PFIELD(:,JL)==XUNDEF .AND. OINTERP(:))) THEN
    WRITE(*,*) 'LAYER ',JL,': NO EXTRAPOLATION : INCREASE YOUR HALO_PREP IN NAM_PREP_SURF_ATM'
    CALL ABOR1_SFX('NO EXTRAPOLATION : INCREASE YOUR HALO_PREP IN NAM_PREP_SURF_ATM')
  ENDIF
  WHERE (.NOT. OINTERP(:)) PFIELD(:,JL) = XUNDEF
ENDDO
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HOR_EXTRAPOL_SURF',1,ZHOOK_HANDLE)
!
END SUBROUTINE HOR_EXTRAPOL_SURF
