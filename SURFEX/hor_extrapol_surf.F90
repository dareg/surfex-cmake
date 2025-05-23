!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
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
!!     A. Napoly    10/2022 add OpenMP directives and optimisations in loops
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NPIO, NCOMM, IDX_I, LSFX_MPI
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
REAL :: ZWORK
REAL, DIMENSION(:,:), ALLOCATABLE :: ZNDIST! smallest distance to valid point
REAL :: ZCOSLA! cosine of latitude
REAL :: ZLONSC! longitude of valid point
REAL :: ZIDLO, ZIDLOMAX, ZIDLOMIN, ZIDLAMAX, ZIDLAMIN
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCOOR
REAL, DIMENSION(:), ALLOCATABLE :: ZIDLA
REAL, DIMENSION(:), ALLOCATABLE :: ZLA       ! input "latitude"  coordinate
REAL, DIMENSION(:), ALLOCATABLE :: ZLO       ! input "longitude" coordinate
REAL :: ZRAD ! conversion degrees to radians
!
INTEGER, DIMENSION(:), ALLOCATABLE :: IMASK
INTEGER, DIMENSION(:,:), ALLOCATABLE :: IVAL_EXT
INTEGER, DIMENSION(NPROC) :: INO_TAB
INTEGER  :: INO     ! output array size
INTEGER, DIMENSION(2) :: ITSIZE, ITDIM
INTEGER, DIMENSION(2,0:NPROC-1) :: IBOR
INTEGER :: ISIZE, ISIZE_MAX, J, ID0, ICOMPT, ICPT
INTEGER :: INFOMPI, IDX, INL
INTEGER  :: JI, JL, JLAT, JLON, JIPOS, JP   ! loop index on points
INTEGER  :: JISC  ! loop index on valid points
INTEGER, DIMENSION(SIZE(KP,2)) :: IKP
#ifdef SFX_MPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
LOGICAL  :: GLALO ! flag true is second coordinate is a longitude or pseudo-lon.
                  !      false if metric coordinates
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE_OMP
REAL :: ZWORK1,ZWORK2,ZWORK3,ZWORK4,ZWORK5,ZWORK6
REAL, DIMENSION(KINLA) :: ZWORK7
INTEGER, DIMENSION(KINLA) :: IWORK1
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HOR_EXTRAPOL_SURF_1',0,ZHOOK_HANDLE)
!
INO = SIZE(PFIELD,1)
INL = SIZE(PFIELD,2)
GLALO = HCOORTYPE=='LALO'
!
!-------------------------------------------------------------------------------
!
ALLOCATE(ZIDLA(KINLA))
ALLOCATE(ZLA(KILEN))
ALLOCATE(ZLO(KILEN))
!
ZIDLOMAX = 0.
ZIDLOMIN = XUNDEF
JIPOS = 0
IF (PRESENT(PILATARRAY)) THEN
  ZIDLA(1) = 0.
  DO JLAT=2,KINLA
    ZIDLA(JLAT) = PILATARRAY(JLAT) - PILATARRAY(JLAT-1)
  ENDDO
ELSE
  ZIDLA(:) = (PILA2 - PILA1) / (KINLA - 1)
ENDIF
!
ZIDLAMAX = MAXVAL(ABS(ZIDLA))
ZIDLAMIN = MINVAL(ABS(ZIDLA(2:KINLA)))
!
ZWORK7(:)=0.
DO JLAT=2,KINLA
  ZWORK7(JLAT) = ZWORK7(JLAT-1) + ZIDLA(JLAT)
ENDDO
DO JLAT=1,KINLA
  ZWORK7(JLAT) = PILA1 + ZWORK7(JLAT)
ENDDO
!
DO JLAT=1,KINLA
  IF (GLALO) THEN
    ZIDLO = (PILO2-PILO1) / KINLO(JLAT)
  ELSE
    ZIDLO = (PILO2-PILO1) / (KINLO(JLAT)-1)
  ENDIF
  DO JLON=1,KINLO(JLAT)
    JIPOS = JIPOS + 1
    ZLA(JIPOS) = ZWORK7(JLAT)
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
ZRAD=XPI/180.0
!
!1: ZTLONMIN, ZTLONMAX, ZTLATMIN, ZTLATMAX contain for each point to extrapol 
! the limits of the domain where to search for the valid points, according
! to NHALO_PREP
ISIZE_MAX = 0
ISIZE = 0
!
ZWORK1=MINVAL(ZLO(:))
ZWORK2=MAXVAL(ZLO(:))
ZWORK3=MINVAL(ZLA(:))
ZWORK4=MAXVAL(ZLA(:))
ZWORK5=ZIDLOMAX*NHALO_PREP
ZWORK6=ZIDLAMAX*NHALO_PREP
!
ICPT = 0
ALLOCATE(IMASK(INO))
IMASK(:)=0
DO JI=1,INO
  IF (ALL(PFIELD(JI,:)/=XUNDEF)) CYCLE
  IF (.NOT. OINTERP(JI))  CYCLE
  ICPT = ICPT + 1
  IMASK(ICPT) = JI
ENDDO
ITSIZE(1) = ICPT

!$OMP PARALLEL DO PRIVATE(JP,ISIZE,IKP) REDUCTION(max:ISIZE_MAX)
DO JP=1,ITSIZE(1)
  !
  JI=IMASK(JP)
  !
  IF (NHALO_PREP>0) THEN
    IKP=KP(JI,:)
    ZTLONMIN(JI) = MAX(ZWORK1,MINVAL(ZLO(IKP))-ZWORK5)
    ZTLONMAX(JI) = MIN(ZWORK2,MAXVAL(ZLO(IKP))+ZWORK5)
    ZTLATMIN(JI) = MAX(ZWORK3,MINVAL(ZLA(IKP))-ZWORK6)
    ZTLATMAX(JI) = MIN(ZWORK4,MAXVAL(ZLA(IKP))+ZWORK6)
  ELSE
    ZTLONMIN(JI) = ZWORK1
    ZTLONMAX(JI) = ZWORK2
    ZTLATMIN(JI) = ZWORK3
    ZTLATMAX(JI) = ZWORK4
  ENDIF
  ISIZE = CEILING((ZTLONMAX(JI)-ZTLONMIN(JI)+1)/ZIDLOMIN)*&
          CEILING((ZTLATMAX(JI)-ZTLATMIN(JI)+1)/ZIDLAMIN)
  ISIZE_MAX = MAX(ISIZE_MAX,ISIZE)
  !
ENDDO
ITSIZE(2) = ISIZE_MAX
!
!NPIO knows the numbers of points to extrapolate for all tasks
IF (NPROC>1) THEN
#ifdef SFX_MPI
  IF (LSFX_MPI) THEN
    CALL MPI_GATHER(ITSIZE,2*KIND(ITSIZE)/4,MPI_INTEGER,&
                    IBOR,2*KIND(IBOR)/4,MPI_INTEGER,& 
                    NPIO,NCOMM,INFOMPI)
  ENDIF
#endif
ELSE
  IBOR(:,0) = ITSIZE(:)
ENDIF
!
IF (NRANK==NPIO) THEN
  ALLOCATE(IVAL_EXT(MAXVAL(IBOR(1,:)),MAXVAL(IBOR(2,:))))
  ALLOCATE(ZCOOR   (MAXVAL(IBOR(1,:)),2))
ELSE
  ALLOCATE(IVAL_EXT(ITSIZE(1),ITSIZE(2)))
  ALLOCATE(ZCOOR(ITSIZE(1),2))
ENDIF
IVAL_EXT(:,:) = 0
ZCOOR   (:,:) = 0.
!
!2: loop on the points 

IF (KINLA >= 1) THEN
  IWORK1(1) = 0
  DO JLAT = 1, KINLA - 1
    IWORK1(JLAT + 1) = IWORK1(JLAT) + KINLO(JLAT)
  ENDDO
ENDIF

!$OMP PARALLEL DO PRIVATE(JP,ICOMPT,JI,ICPT,JLAT,ZLAT,ZIDLO,JLON) 
DO JP=1,ITSIZE(1)
  !
  ICOMPT = 0
  JI=IMASK(JP)
  ICPT=JP
  !coordinates of the point in the grid 
  ZCOOR(ICPT,1) = PLAT(JI)
  ZCOOR(ICPT,2) = PLON(JI)
  !
  !loop on the grid
  DO JLAT = 1,KINLA
    ZLAT = ZWORK7(JLAT)
    IF (ZLAT>=ZTLATMIN(JI) .AND. ZLAT<=ZTLATMAX(JI)) THEN
      IF (GLALO) THEN
        ZIDLO = (PILO2-PILO1) / KINLO(JLAT)
      ELSE
        ZIDLO = (PILO2-PILO1) / (KINLO(JLAT)-1)
      ENDIF
      DO JLON = 1,KINLO(JLAT)
        ZLON = PILO1 + (JLON-1) * ZIDLO
        IF (ZLON>=ZTLONMIN(JI) .AND. ZLON<=ZTLONMAX(JI)) THEN
          ICOMPT = ICOMPT + 1
          !ival_ext contains the indexes of the points needed to interpolate
          !in the complete grid
          IVAL_EXT(ICPT,ICOMPT) = IWORK1(JLAT) + JLON
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  !
ENDDO
!$OMP END PARALLEL DO
!
DEALLOCATE(ZIDLA)
!
!
!
IF (NRANK/=NPIO) THEN
  !if some points of this task need to be extrapolated
  IF (SUM(ITSIZE)/=0) THEN

    !zfield will contain the values of the field
    ALLOCATE(ZFIELD(ITSIZE(1),INL))

    IDX_I = IDX_I + 1
    !sends indexes to npio
#ifdef SFX_MPI    
    IF (LSFX_MPI) CALL MPI_SEND(IVAL_EXT,SIZE(IVAL_EXT)*KIND(IVAL_EXT)/4,MPI_INTEGER,NPIO,IDX_I,NCOMM,INFOMPI)
#endif

    IDX_I = IDX_I + 1
    !send coords of the points to extrapolate
#ifdef SFX_MPI    
    IF (LSFX_MPI) CALL MPI_SEND(ZCOOR,SIZE(ZCOOR)*KIND(ZCOOR)/4,MPI_REAL,NPIO,IDX_I,NCOMM,INFOMPI)
#endif

    IDX_I = IDX_I + 1
    !receives values of the field from NPIO
#ifdef SFX_MPI    
    IF (LSFX_MPI) CALL MPI_RECV(ZFIELD,SIZE(ZFIELD)*KIND(ZFIELD)/4,MPI_REAL,NPIO,IDX_I,NCOMM,ISTATUS,INFOMPI)
#endif
    DO JI=1,ITSIZE(1)
      PFIELD(IMASK(JI),:) = ZFIELD(JI,:)
    ENDDO
    DEALLOCATE(ZFIELD)
  ELSE
    IDX_I = IDX_I + 3
  ENDIF
ELSE
  DO JP = NPIO,NPROC-1+NPIO
    J = JP
    IF (JP>NPROC-1) J = JP-NPROC
    IF (SUM(IBOR(:,J))/=0) THEN
      ALLOCATE(ZFIELD(IBOR(1,J),INL))
      ZFIELD=XUNDEF
      IF (J/=NPIO) THEN
        !receives indexes and coordinaites
#ifdef SFX_MPI        
        IF (LSFX_MPI) THEN
          CALL MPI_RECV(IVAL_EXT(1:IBOR(1,J),1:IBOR(2,J)), IBOR(1,J)*IBOR(2,J)*KIND(IVAL_EXT)/4, &
                          MPI_INTEGER, J, IDX_I+1, NCOMM, ISTATUS, INFOMPI)
          CALL MPI_RECV(ZCOOR(1:IBOR(1,J),:), IBOR(1,J)*SIZE(ZCOOR,2)*KIND(ZCOOR)/4,&
                          MPI_REAL, J, IDX_I+2, NCOMM, ISTATUS, INFOMPI)
        ENDIF
#endif

      ENDIF
      ALLOCATE(ZNDIST(IBOR(1,J),INL))
!$OMP PARALLEL DO PRIVATE(JI,IDX,ZCOSLA,JISC,JL,ID0,ZLONSC,ZDIST) 
        DO JI=1,IBOR(1,J)
          ZNDIST(JI,:) = XUNDEF
          IDX = IBOR(2,J)+1
          ZCOSLA=COS(ZCOOR(JI,1)*ZRAD)
          DO JISC = 1,IBOR(2,J)
            !index in the whole grid of the point used to interpolate
            ID0 = IVAL_EXT(JI,JISC)
            IF (ID0==0) EXIT
            IF (ANY(PFIELD_IN(ID0,:)/=XUNDEF)) THEN
              ZLONSC = ZLO(ID0)
              IF (GLALO) THEN
                IF (ZLONSC-ZCOOR(JI,2)> 180.) ZLONSC = ZLONSC - 360.
                IF (ZLONSC-ZCOOR(JI,2)<-180.) ZLONSC = ZLONSC + 360.
                ZDIST= (ZLA(ID0)-ZCOOR(JI,1)) ** 2 + ((ZLONSC-ZCOOR(JI,2))*ZCOSLA) ** 2
              ELSE
                ZDIST= (ZLA(ID0)-ZCOOR(JI,1)) ** 2 + (ZLONSC-ZCOOR(JI,2)) ** 2
              END IF
              DO JL=1,INL
                IF (ZDIST<=ZNDIST(JI,JL)) THEN
                  IF (PFIELD_IN(ID0,JL)/=XUNDEF) THEN
                    ZFIELD(JI,JL) = PFIELD_IN(ID0,JL)
                    ZNDIST(JI,JL) = ZDIST
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          END DO   
        ENDDO
!$OMP END PARALLEL DO
      !
      DEALLOCATE(ZNDIST)
      !
      IF (J/=NPIO) THEN
        !send values found to extrapolate
#ifdef SFX_MPI        
        IF (LSFX_MPI) CALL MPI_SEND(ZFIELD,SIZE(ZFIELD)*KIND(ZFIELD)/4,MPI_REAL,J,IDX_I+3,NCOMM,INFOMPI)
#endif
      ELSE
        DO JI = 1,IBOR(1,J)
          PFIELD(IMASK(JI),:) = ZFIELD(JI,:)
        ENDDO
      ENDIF
      !
      DEALLOCATE(ZFIELD)
      !
    ENDIF
    !
  ENDDO
  !
  IDX_I = IDX_I + 3
  !
ENDIF
!
DEALLOCATE(ZCOOR)
DEALLOCATE(IVAL_EXT)
DEALLOCATE(IMASK)
!
IF (ALLOCATED(ZLA)) DEALLOCATE(ZLA)
IF (ALLOCATED(ZLO)) DEALLOCATE(ZLO)
DEALLOCATE(ZTLONMIN,ZTLONMAX,ZTLATMIN,ZTLATMAX)
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
IF (LHOOK) CALL DR_HOOK('HOR_EXTRAPOL_SURF_1',1,ZHOOK_HANDLE)
!
END SUBROUTINE HOR_EXTRAPOL_SURF
