!     #########
      SUBROUTINE TREAT_FIELD (UG, U, USS, &
                              HPROGRAM,HSCHEME,HFILETYPE,    &
                              HSUBROUTINE,HFILENAME,HFIELD,   &
                              PPGDARRAY,HSFTYPE               )  
!     ##############################################################
!
!!**** *TREAT_FIELD* chooses which treatment subroutine to use
!!
!!    PURPOSE
!!    -------
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
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    11/09/95
!!
!!    Modification
!!    25/05/96    (V. Masson) remove useless case for HSUBROUTINE   
!!    29/11/2002  (D. Gazen)  add HSFTYPE argument + call to read_binllvfast routine
!!    03/2004     (V. MAsson) externalization
!!    04/2009     (B. Decharme) Special treatement for gaussian grid
!!    06/2009     (B. Decharme)  call Topographic index statistics calculation
!!    09/2010     (E. Kourzeneva) call reading of the lake database
!!    03/2012     (M. Lafaysse) NETCDF
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
!
!
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_PGDWORK, ONLY : NSIZE, NSIZE_ALL, X1D_ALL, X2D_ALL, X3D_ALL, N3D_ALL, &
                         XSUMVAL, XSUMVAL2, XEXTVAL2, XSUMVAL2, &
                         XSSQO, LSSQO, NSSO, XMIN_WORK, XMAX_WORK
!
USE MODD_SURFEX_OMP, ONLY : NBLOCKTOT
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NPROC, NCOMM, NREQ, NINDEX, IDX_R, &
                                NSIZE_TASK,NREQ, NSIZE_max=>NSIZE
!
USE MODD_DATA_LAKE,      ONLY : NGRADDEPTH_LDB, NGRADSTATUS_LDB 
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
!
USE MODI_INI_SSOWORK
USE MODI_GET_LUOUT
USE MODI_READ_DIRECT
USE MODI_READ_DIRECT_GAUSS
USE MODI_READ_LATLON
USE MODI_READ_BINLLV
USE MODI_READ_BINLLVFAST
USE MODI_READ_ASCLLV
USE MODI_READ_AND_SEND_MPI
USE MODI_READ_PGD_NETCDF
USE MODI_AVERAGE2_MESH
USE MODI_MAKE_LCOVER
USE MODI_ABOR1_SFX
USE MODI_AVERAGE2_COVER
USE MODI_AVERAGE2_CTI
USE MODI_AVERAGE2_LDB
USE MODI_AVERAGE2_OROGRAPHY
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
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
!
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
!
 CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM      ! Type of program
 CHARACTER(LEN=6),  INTENT(IN) :: HSCHEME       ! Scheme treated
 CHARACTER(LEN=6),  INTENT(IN) :: HFILETYPE     ! Type of the data file
 CHARACTER(LEN=6),  INTENT(IN) :: HSUBROUTINE   ! Name of the subroutine to call
 CHARACTER(LEN=28), INTENT(IN) :: HFILENAME     ! Name of the field file.
 CHARACTER(LEN=20), INTENT(IN) :: HFIELD        ! Name of the field.
REAL, DIMENSION(:), INTENT(INOUT), OPTIONAL :: PPGDARRAY ! field on MESONH grid
 CHARACTER(LEN=3),   INTENT(IN),    OPTIONAL :: HSFTYPE
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: I3D_ALL
INTEGER, DIMENSION(:), ALLOCATABLE :: IMASK
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZSUMVAL3
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSUMVAL2
REAL, DIMENSION(:), ALLOCATABLE :: ZEXTVAL
INTEGER, DIMENSION(:), ALLOCATABLE :: ISIZE
!
#ifdef SFX_MPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
INTEGER, DIMENSION(MPI_STATUS_SIZE,NPROC-1) :: ISTATUS2
#endif
INTEGER, DIMENSION(0:NPROC-1) :: ITCOV
INTEGER :: ILUOUT, IS2, INFOMPI, JP, ICPT, JCOV, JI, JL, IREQ, IDX,&
                        IDX_SAVE, ISIZE_OMP
REAL(KIND=JPRB) :: ZHOOK_HANDLE, ZHOOK_HANDLE_OMP
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TREAT_FIELD',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!*    1.     Selection of type of reading (and point by point treatment)
!            -----------------------------------------------------------
!
SELECT CASE (HFILETYPE)

  CASE ('DIRECT')
    IF(UG%CGRID=="GAUSS     " .OR. UG%CGRID=="IGN       " .OR. UG%CGRID=="LONLAT REG")THEN
      CALL READ_DIRECT_GAUSS(UG, U, USS, &
                             HPROGRAM,HSCHEME,HSUBROUTINE,HFILENAME,HFIELD)
    ELSE
      CALL READ_DIRECT(UG, U, USS, &
                       HPROGRAM,HSCHEME,HSUBROUTINE,HFILENAME,HFIELD)
    ENDIF

  CASE ('BINLLV')
    CALL INI_SSOWORK
    IF (NRANK==NPIO) CALL READ_BINLLV(UG, U, USS, &
                        HPROGRAM,HSUBROUTINE,HFILENAME)

  CASE ('BINLLF')
    CALL INI_SSOWORK
    IF (NRANK==NPIO) CALL READ_BINLLVFAST(UG, U, USS, &
                            HPROGRAM,HSUBROUTINE,HFILENAME)

  CASE ('ASCLLV')
    CALL INI_SSOWORK
    IF (NRANK==NPIO) CALL READ_ASCLLV(UG, U, USS, &
                        HPROGRAM,HSUBROUTINE,HFILENAME)

  CASE ('LATLON')
    CALL INI_SSOWORK
    IF (NRANK==NPIO) CALL READ_LATLON(UG, U, USS, &
                        HPROGRAM,HSCHEME,HSUBROUTINE,HFILENAME)

  CASE ('NETCDF')
    CALL INI_SSOWORK
    IF (NRANK==NPIO) CALL READ_PGD_NETCDF(UG, U, USS, &
                            HPROGRAM,HSCHEME,HSUBROUTINE,HFILENAME,HFIELD)

  CASE DEFAULT
    CALL ABOR1_SFX('TREAT_FIELD: FILE TYPE NOT SUPPORTED: '//HFILETYPE)

END SELECT
!
!-------------------------------------------------------------------------------
!
!nsize contains the number of points found for each of the domain, for each task
ALLOCATE(NSIZE(U%NSIZE_FULL))
!
IF (NPROC>1) THEN
  !
  ALLOCATE(ISIZE(NSIZE_max))
  !
  IDX_SAVE = IDX_R
  IDX = IDX_SAVE + NRANK
  !each task sends to each other task the part of NSIZE_ALL it got, stored in
  !isize
  CALL READ_AND_SEND_MPI(NSIZE_ALL,ISIZE(1:NSIZE_TASK(NRANK)),KPIO=NRANK,KDX=IDX)
  !
  NSIZE(:) = 0
  ISIZE_OMP = NPROC/NBLOCKTOT
  !for each task
  DO JP=0,NPROC-1
   !
    IF (JP/=NRANK) THEN
      !
      !each task receives each ISIZE from each task
      CALL MPI_RECV(ISIZE,NSIZE_max*KIND(ISIZE)/4,MPI_INTEGER,&
                JP,IDX_SAVE+1+JP,NCOMM,ISTATUS,INFOMPI)
      !
    ELSE
      !
      ICPT = 0
      DO JI = 1,SIZE(NINDEX)
        IF (NINDEX(JI)==JP) THEN
          ICPT = ICPT + 1
          ISIZE(ICPT) = NSIZE_ALL(JI)
        ENDIF
      ENDDO
      !
    ENDIF
    !
    !nsize is the sum of all parts isize
    NSIZE(:) = NSIZE(:) + ISIZE(1:NSIZE_TASK(NRANK))
    !
  ENDDO
  DEALLOCATE(ISIZE)
  CALL MPI_WAITALL(NPROC-1,NREQ(1:NPROC-1),ISTATUS2,INFOMPI)
ELSE
  NSIZE(:) = NSIZE_ALL(:)
ENDIF
!
!
DEALLOCATE(NSIZE_ALL)
!
!
IF (HSUBROUTINE=='A_COVR') THEN
  IS2 = SIZE(X3D_ALL,2)
ELSEIF (HSUBROUTINE=='A_LDBS') THEN
  IS2 = NGRADSTATUS_LDB
ELSEIF (HSUBROUTINE=='A_LDBD') THEN
  IS2 = NGRADDEPTH_LDB
ELSEIF (HSUBROUTINE=='A_OROG') THEN
  IS2 = 2
ELSEIF (HSUBROUTINE=='A_CTI ') THEN
  IS2 = 3 
ENDIF
!
!
SELECT CASE (HSUBROUTINE)

  CASE ('A_COVR')
    !
    !to gather the parts of LCOVER sparsed among the tasks
    CALL MAKE_LCOVER(U%LCOVER)
    !
    !contains the indexes of the covers in XCOVER, associated to their effective
    !num
    ALLOCATE(IMASK(SIZE(U%LCOVER)))
    IMASK(:) = 0
    ICPT = 0
    DO JCOV = 1,SIZE(U%LCOVER)
      IF (U%LCOVER(JCOV)) THEN
        ICPT = ICPT + 1
        IMASK(JCOV) = ICPT
      ENDIF
    ENDDO
    !
    !because each task did not necessarily meet the same number of coves, itcov 
    !contains the numbers of covers met for all tasks 
    IF (NPROC>1) THEN
#ifdef SFX_MPI
    CALL MPI_ALLGATHER(IS2,KIND(IS2)/4,MPI_INTEGER,&
                     ITCOV,KIND(ITCOV)/4,MPI_INTEGER,NCOMM,INFOMPI)
#endif
    ELSE
      ITCOV(:) = IS2
    ENDIF
    !
    !XSUMVAL2 needs to contain the numbers of times each cover is encountered
    !for the current task 
    ALLOCATE(XSUMVAL2(U%NSIZE_FULL,COUNT(U%LCOVER)))
    XSUMVAL2(:,:) = 0.
    !
    IF (NPROC>1) THEN
      ALLOCATE(ZSUMVAL3(U%NSIZE_FULL,MAXVAL(ITCOV),2))
      DO JP = 0,NPROC-1
        !the part of X3D_ALL concerning the current task is sent 
        !from each other task
        CALL READ_AND_SEND_MPI(X3D_ALL,ZSUMVAL3(:,1:ITCOV(JP),:),KPIO=JP)
        DO JL = 1,ITCOV(JP)
          DO JI=1,U%NSIZE_FULL
            JCOV = NINT(ZSUMVAL3(JI,JL,1))
            IF (JCOV/=0) THEN
              !xsumval2 is the sum of contributions zsumval3 coming from all tasks
              XSUMVAL2(JI,IMASK(JCOV)) = XSUMVAL2(JI,IMASK(JCOV)) + ZSUMVAL3(JI,JL,2)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(ZSUMVAL3)
    ELSE
      DO JL = 1,SIZE(X3D_ALL,2)
        DO JI=1,SIZE(X3D_ALL,1)
          JCOV = NINT(X3D_ALL(JI,JL,1))
          IF (JCOV/=0) THEN
            XSUMVAL2(JI,IMASK(JCOV)) = XSUMVAL2(JI,IMASK(JCOV)) + X3D_ALL(JI,JL,2)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    DEALLOCATE(X3D_ALL,IMASK)
    !
    !
  CASE ('A_LDBD','A_LDBS','A_OROG','A_CTI ')
    !
    !XSUMVAL2 needs to contain the numbers of times each quantity is encountered
    !for the current task    
    ALLOCATE(XSUMVAL2(U%NSIZE_FULL,IS2))
    XSUMVAL2(:,:) = 0.
    IF (NPROC>1) THEN
      ALLOCATE(ZSUMVAL2(U%NSIZE_FULL,IS2))
      DO JP = 0,NPROC-1
        CALL READ_AND_SEND_MPI(X2D_ALL,ZSUMVAL2,KPIO=JP)
        !sum of contributions coming from all tasks
        XSUMVAL2(:,:) = XSUMVAL2(:,:) + ZSUMVAL2(:,:)
      ENDDO
      DEALLOCATE(ZSUMVAL2)
    ELSE
      XSUMVAL2(:,:) = X2D_ALL(:,:)
    ENDIF
    DEALLOCATE(X2D_ALL)
    !
    !
    IF (HSUBROUTINE=='A_OROG' .OR. HSUBROUTINE=='A_CTI ') THEN
      !
      !special fields
      !
      IF (HSUBROUTINE=='A_CTI ') THEN
        !max and min
        IF (NPROC>1) THEN
          ALLOCATE(ZEXTVAL(U%NSIZE_FULL))
          DO JP = 0,NPROC-1
            CALL READ_AND_SEND_MPI(XEXTVAL2(:,1),ZEXTVAL,KPIO=JP)
            XMAX_WORK(:) = MAX(XMAX_WORK,ZEXTVAL)
          ENDDO
          DO JP = 0,NPROC-1
            CALL READ_AND_SEND_MPI(XEXTVAL2(:,2),ZEXTVAL,KPIO=JP)
            XMIN_WORK(:) = MIN(XMIN_WORK,ZEXTVAL)
          ENDDO   
          DEALLOCATE(ZEXTVAL)
        ELSE
          XMAX_WORK(:) = XEXTVAL2(:,1)
          XMIN_WORK(:) = XEXTVAL2(:,2)
        ENDIF
        !
      ELSEIF (HSUBROUTINE=='A_OROG') THEN
        !max and min
        IF (NPROC>1) THEN
          ALLOCATE(ZEXTVAL(U%NSIZE_FULL))
          DO JP = 0,NPROC-1
            CALL READ_AND_SEND_MPI(XEXTVAL2(:,1),ZEXTVAL,KPIO=JP)
            USS%XMAX_ZS(:) = MAX(USS%XMAX_ZS,ZEXTVAL)
          ENDDO 
          DO JP = 0,NPROC-1
            CALL READ_AND_SEND_MPI(XEXTVAL2(:,2),ZEXTVAL,KPIO=JP)
            USS%XMIN_ZS(:) = MIN(USS%XMIN_ZS,ZEXTVAL)
          ENDDO   
          DEALLOCATE(ZEXTVAL)
        ELSE
          USS%XMAX_ZS(:) = XEXTVAL2(:,1)
          USS%XMIN_ZS(:) = XEXTVAL2(:,2)
        ENDIF
        !
        !sso fields   
        ALLOCATE(XSSQO(U%NSIZE_FULL,NSSO,NSSO))
        XSSQO(:,:,:) = -XUNDEF
        IF (NPROC>1) THEN
          ALLOCATE(ZSUMVAL3(U%NSIZE_FULL,NSSO,NSSO))
          DO JP = 0,NPROC-1
            CALL READ_AND_SEND_MPI(X3D_ALL,ZSUMVAL3,KPIO=JP)
            XSSQO(:,:,:) = MAX(XSSQO(:,:,:),ZSUMVAL3)
          ENDDO
          DEALLOCATE(ZSUMVAL3)
        ELSE
          XSSQO(:,:,:) = X3D_ALL(:,:,:)
        ENDIF
        DEALLOCATE(X3D_ALL)
        !
        ALLOCATE(LSSQO(U%NSIZE_FULL,NSSO,NSSO))
        LSSQO(:,:,:) = .FALSE.
        IF (NPROC>1) THEN
          ALLOCATE(I3D_ALL(U%NSIZE_FULL,NSSO,NSSO)) 
          DO JP = 0,NPROC-1       
            CALL READ_AND_SEND_MPI(N3D_ALL,I3D_ALL,KPIO=JP)
            WHERE (I3D_ALL(:,:,:)==1) LSSQO(:,:,:) = .TRUE.
          ENDDO
          DEALLOCATE(I3D_ALL)
        ELSE
          WHERE(N3D_ALL(:,:,:)==1) LSSQO(:,:,:) = .TRUE.
        ENDIF
        DEALLOCATE(N3D_ALL)
      ENDIF
      DEALLOCATE(XEXTVAL2)
      !
    ENDIF
    !
    !
  CASE ('A_MESH')
    !most simple case
    ALLOCATE(XSUMVAL(U%NSIZE_FULL))
    IF (NPROC>1) THEN
      XSUMVAL(:) = 0.
      ALLOCATE(ZEXTVAL(U%NSIZE_FULL))
      DO JP = 0,NPROC-1
        CALL READ_AND_SEND_MPI(X1D_ALL,ZEXTVAL,KPIO=JP)
        XSUMVAL(:) = XSUMVAL(:) + ZEXTVAL(:)
      ENDDO
      DEALLOCATE(ZEXTVAL)
    ELSE
      XSUMVAL(:) = X1D_ALL(:)
    ENDIF
    DEALLOCATE(X1D_ALL)
    !
END SELECT
!
!-------------------------------------------------------------------------------
!
!*    2.     Call to the adequate subroutine (global treatment)
!            --------------------------------------------------
!
SELECT CASE (HSUBROUTINE)

  CASE ('A_COVR')
    CALL AVERAGE2_COVER(U, &
                        HPROGRAM)

  CASE ('A_OROG')
    CALL AVERAGE2_OROGRAPHY(USS)

  CASE ('A_CTI ')
    CALL AVERAGE2_CTI

  CASE ('A_LDBD')
    CALL AVERAGE2_LDB(PPGDARRAY,'D',1)

  CASE ('A_LDBS')
    CALL AVERAGE2_LDB(PPGDARRAY,'S',1)
    
  CASE ('A_MESH')
    IF (.NOT. PRESENT(PPGDARRAY)) THEN
      WRITE(ILUOUT,*) 'You asked to average a PGD field with A_MESH option,'
      WRITE(ILUOUT,*) 'but you did not give the array to store this field'
      CALL ABOR1_SFX('TREAT_FIELD: ARRAY IS MISSING')
    END IF
    CALL AVERAGE2_MESH(PPGDARRAY)

END SELECT
!
IF (LHOOK) CALL DR_HOOK('TREAT_FIELD',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE TREAT_FIELD
