MODULE MODE_READ_SURF_LFI
!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURF_LFI is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
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
!!
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
INTERFACE READ_SURF0_LFI
        MODULE PROCEDURE READ_SURFX0_LFI
        MODULE PROCEDURE READ_SURFN0_LFI
        MODULE PROCEDURE READ_SURFL0_LFI
        MODULE PROCEDURE READ_SURFC0_LFI
END INTERFACE
INTERFACE READ_SURFN_LFI
        MODULE PROCEDURE READ_SURFX1_LFI
        MODULE PROCEDURE READ_SURFN1_LFI
        MODULE PROCEDURE READ_SURFL1_LFI
        MODULE PROCEDURE READ_SURFX2_LFI
END INTERFACE
INTERFACE READ_SURFT_LFI
        MODULE PROCEDURE READ_SURFT0_LFI
        MODULE PROCEDURE READ_SURFT1_LFI
END INTERFACE
!
CONTAINS
!
!     #############################################################
      SUBROUTINE READ_SURFX0_LFI(HREC,PFIELD,KRESP,HCOMMENT)
!     #############################################################
!
USE MODD_IO_SURF_LFI,        ONLY : CFILE_LFI, CLUOUT_LFI
!
USE MODI_FMREAD
USE MODI_ERROR_READ_SURF_LFI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*), INTENT(IN)  :: HREC     ! name of the article to be read
REAL,              INTENT(OUT) :: PFIELD   ! the real scalar to be read
INTEGER,           INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
INTEGER          :: IGRID   ! position of data on grid
INTEGER          :: ILENCH  ! length of comment string
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFX0_LFI',0,ZHOOK_HANDLE)
!
KRESP=0
!
CALL FMREADX0(CFILE_LFI,HREC,CLUOUT_LFI,1,PFIELD,IGRID,ILENCH,HCOMMENT,KRESP)
!
CALL ERROR_READ_SURF_LFI(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFX0_LFI',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX0_LFI
!
!     #############################################################
      SUBROUTINE READ_SURFX1_LFI(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX1* - routine to fill a real 1D array for the externalised surface 
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ
!
USE MODD_IO_SURF_LFI, ONLY : CFILE_LFI, CLUOUT_LFI, NMASK, NFULL, &
                             LMNH_COMPATIBLE  
!
USE MODI_FMREAD
USE MODI_ERROR_READ_SURF_LFI
USE MODI_READ_AND_SEND_MPI
USE MODI_GET_SURF_UNDEF
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
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),   INTENT(IN)  :: HREC     ! name of the article to be read
REAL, DIMENSION(:),  INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,             INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),  INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),    INTENT(IN)  :: HDIR     ! type of field :
                                             ! 'H' : field with
                                             !       horizontal spatial dim.
                                             ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=18) :: YREC
REAL              :: ZUNDEF  ! default value
INTEGER           :: IGRID   ! position of data on grid
INTEGER           :: ILENCH  ! length of comment string
INTEGER           :: IVERSION, IBUGFIX
INTEGER           :: I, INFOMPI
!
#ifndef NOMPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
DOUBLE PRECISION   :: XTIME0
REAL, DIMENSION(:),   ALLOCATABLE :: ZWORK   ! work array read in the file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFX1_LFI',0,ZHOOK_HANDLE)
!
KRESP=0
!
#ifndef NOMPI
XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  !
  IF (HDIR=='H' .OR. HDIR=='A') THEN
    ALLOCATE(ZWORK(NFULL))
  ELSE
    ALLOCATE(ZWORK(SIZE(PFIELD)))
  ENDIF
  !
!$OMP SINGLE        
  !  
  YREC = HREC
  !
  !---------------------------------------------------------------------------
  !* patch to read some test files done before version 3.5
  !  this should be removed once all tests with reading lfi files done with 923
  !  configuration (with these early versions) are finished.
  !
  IF (HREC(1:2)=='D_') THEN
    CALL FMREADN0(CFILE_LFI,'VERSION',CLUOUT_LFI,1,IVERSION,IGRID,ILENCH,HCOMMENT,KRESP)
    CALL FMREADN0(CFILE_LFI,'BUG',CLUOUT_LFI,1,IBUGFIX,IGRID,ILENCH,HCOMMENT,KRESP)
    IF (IVERSION<=2 .OR. (IVERSION==3 .AND. IBUGFIX<=5)) YREC = 'DATA_'//HREC(3:12)
  END IF
  !---------------------------------------------------------------------------
  !
  IF (HDIR=='H' .OR. HDIR=='A') THEN
    IF (.NOT. LMNH_COMPATIBLE) THEN
      CALL FMREADX1(CFILE_LFI,YREC,CLUOUT_LFI,NFULL,ZWORK,IGRID,ILENCH,HCOMMENT,KRESP)
    ELSE
      CALL READ_IN_LFI_X1_FOR_MNH(YREC,ZWORK,KRESP,HCOMMENT,HDIR)
    END IF
  ELSE
    CALL FMREADX1(CFILE_LFI,YREC,CLUOUT_LFI,SIZE(ZWORK),ZWORK,IGRID,ILENCH,HCOMMENT,KRESP)
  END IF
  CALL ERROR_READ_SURF_LFI(YREC,KRESP)
  !
!$OMP END SINGLE COPYPRIVATE(ZWORK,HCOMMENT,KRESP)
  !  
ELSE
  ALLOCATE(ZWORK(0))
ENDIF
!
#ifndef NOMPI
XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
!
IF (HDIR=='A') THEN  ! no distribution on other tasks
  IF ( NRANK==NPIO ) THEN
#ifndef NOMPI          
    XTIME0 = MPI_WTIME()
#endif
    PFIELD(:) = ZWORK(1:SIZE(PFIELD))
#ifndef NOMPI
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif
  ENDIF
ELSEIF (HDIR=='-') THEN ! distribution of the total field on other tasks
!$OMP SINGLE    
  PFIELD(:) = ZWORK(1:SIZE(PFIELD))
  IF (NPROC>1) THEN
#ifndef NOMPI
    XTIME0 = MPI_WTIME()
    CALL MPI_BCAST(PFIELD,SIZE(PFIELD)*KIND(PFIELD)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif
  ENDIF
!$OMP END SINGLE COPYPRIVATE(PFIELD)    
ELSE
  CALL READ_AND_SEND_MPI(ZWORK,PFIELD,NMASK)
ENDIF
!
DEALLOCATE(ZWORK)
!
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFX1_LFI',1,ZHOOK_HANDLE)
!
CONTAINS
!
!     #############################################################
      SUBROUTINE READ_IN_LFI_X1_FOR_MNH(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to fill a read 2D array for the externalised surface 
!
USE MODD_IO_SURF_LFI,        ONLY : CFILE_LFI, CLUOUT_LFI, &
                                    NIU, NIB, NIE, NJU, NJB, NJE
!
USE MODI_FMREAD
USE MODI_ERROR_READ_SURF_LFI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),        INTENT(IN) :: HREC     ! name of the article to be read
REAL, DIMENSION(:),     INTENT(OUT):: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(OUT):: HCOMMENT ! comment string
CHARACTER(LEN=1),         INTENT(IN) :: HDIR     ! type of field :
                                                 ! 'H' : field with
                                                 !       horizontal spatial dim.
                                                 ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
! 
CHARACTER(LEN=4)   :: YREC1D
INTEGER :: JI, JJ
INTEGER :: ILEN
INTEGER :: IGRID, ILENCH
REAL, DIMENSION(:),   ALLOCATABLE :: ZWORK1D! 1D work array read in the file
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK2D ! work array read in a MNH file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFX1_LFI:READ_IN_LFI_X1_FOR_MNH',0,ZHOOK_HANDLE)
!
ALLOCATE(ZWORK2D(NIU,NJU))
ZWORK2D(:,:) = 999.
!
IF (HREC=='XX                 ' .OR. HREC=='DX                 ') THEN
  ALLOCATE(ZWORK1D(NIU))
  YREC1D = 'XHAT'
  ILEN   = NIU
ELSEIF (HREC=='YY                 ' .OR. HREC=='DY                 ') THEN
  ALLOCATE(ZWORK1D(NJU))
  YREC1D = 'YHAT'
  ILEN   = NJU
ELSEIF (NJB==NJE) THEN
  ALLOCATE(ZWORK1D(NIU))
  ZWORK1D(:) = 999.
ELSEIF (NIB==NIE) THEN
  ALLOCATE(ZWORK1D(NJU))
  ZWORK1D(:) = 999. 
ENDIF
!
IF (HREC=='XX' .OR. HREC=='YY'.OR. HREC=='DX' .OR. HREC=='DY') THEN
!
  CALL FMREADX1(CFILE_LFI,YREC1D,CLUOUT_LFI,ILEN,ZWORK1D,IGRID,ILENCH,HCOMMENT,KRESP)
  CALL ERROR_READ_SURF_LFI(YREC1D,KRESP)
!
  SELECT CASE(HREC)
    CASE('XX                  ')
      DO JJ = 1,SIZE(ZWORK2D,2)
        ZWORK2D(NIB:NIE,JJ) = 0.5 * ZWORK1D(NIB:NIE) + 0.5 * ZWORK1D(NIB+1:NIE+1)
      END DO
    CASE('DX                  ')
      DO JJ = 1,SIZE(ZWORK2D,2)
        ZWORK2D(NIB:NIE,JJ) = - ZWORK1D(NIB:NIE) + ZWORK1D(NIB+1:NIE+1)
      END DO
    CASE('YY                  ')
      DO JI = 1,SIZE(ZWORK2D,1)
        ZWORK2D(JI,NJB:NJE) = 0.5 * ZWORK1D(NJB:NJE) + 0.5 * ZWORK1D(NJB+1:NJE+1)
      END DO
    CASE('DY                  ')
      DO JI = 1,SIZE(ZWORK2D,1)
        ZWORK2D(JI,NJB:NJE) = - ZWORK1D(NJB:NJE) + ZWORK1D(NJB+1:NJE+1)
      END DO
  END SELECT
!
  DEALLOCATE(ZWORK1D)
!
ELSEIF (NJB==NJE) THEN
!
  CALL FMREADX1(CFILE_LFI,YREC,CLUOUT_LFI,SIZE(ZWORK1D),ZWORK1D,IGRID,ILENCH,HCOMMENT,KRESP)
  DO JJ = 1,SIZE(ZWORK2D,2)
    ZWORK2D(NIB:NIE,JJ) = ZWORK1D(NIB:NIE)
  END DO
!
  DEALLOCATE(ZWORK1D) 
!
ELSEIF (NIB==NIE) THEN
!
  CALL FMREADX1(CFILE_LFI,YREC,CLUOUT_LFI,SIZE(ZWORK1D),ZWORK1D,IGRID,ILENCH,HCOMMENT,KRESP) 
  DO JI = 1,SIZE(ZWORK2D,1)
    ZWORK2D(JI,NJB:NJE) = ZWORK1D(NJB:NJE)
  END DO
!
  DEALLOCATE(ZWORK1D)
!
ELSE
!
  CALL FMREADX2(CFILE_LFI,HREC,CLUOUT_LFI,SIZE(ZWORK2D),ZWORK2D,IGRID,ILENCH,HCOMMENT,KRESP)
!
ENDIF
!
DO JJ=1,NJE-NJB+1
  DO JI=1,NIE-NIB+1
    PFIELD(JI+(NIE-NIB+1)*(JJ-1)) = ZWORK2D(NIB+JI-1,NJB+JJ-1) 
  END DO
END DO
!
DEALLOCATE(ZWORK2D) 
!  
CALL ERROR_READ_SURF_LFI(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFX1_LFI:READ_IN_LFI_X1_FOR_MNH',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_IN_LFI_X1_FOR_MNH
!
END SUBROUTINE READ_SURFX1_LFI
!
!     #############################################################
      SUBROUTINE READ_SURFX2_LFI(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX2* - routine to fill a real 2D array for the externalised surface 
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ
!
USE MODD_IO_SURF_LFI, ONLY : CFILE_LFI, CLUOUT_LFI, NMASK, NFULL, &
                             LMNH_COMPATIBLE  
!
USE MODI_FMREAD
USE MODI_ERROR_READ_SURF_LFI
USE MODI_READ_AND_SEND_MPI
USE MODI_GET_SURF_UNDEF
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
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
REAL, DIMENSION(:,:),     INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),         INTENT(IN)  :: HDIR     ! type of field :
                                                  ! 'H' : field with
                                                  !       horizontal spatial dim.
                                                  ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
! 
CHARACTER(LEN=16) :: YREC
REAL              :: ZUNDEF  ! default value
INTEGER           :: IGRID   ! position of data on grid
INTEGER           :: ILENCH  ! length of comment string
INTEGER           :: IVERSION, IBUGFIX
INTEGER           :: I, INFOMPI
!
#ifndef NOMPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
DOUBLE PRECISION   :: XTIME0
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK   ! work array read in the file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFX2_LFI',0,ZHOOK_HANDLE)
!
KRESP=0
!
#ifndef NOMPI
XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  !
  IF (HDIR=='H' .OR. HDIR=='A') THEN
    ALLOCATE(ZWORK(NFULL,SIZE(PFIELD,2)))
  ELSE
    ALLOCATE(ZWORK(SIZE(PFIELD,1),SIZE(PFIELD,2)))
  ENDIF
  !
!$OMP SINGLE        
  !  
  YREC = HREC
  !
  !---------------------------------------------------------------------------
  !* patch to read some test files done before version 3.5
  !  this should be removed once all tests with reading lfi files done with 923
  !  configuration (with these early versions) are finished.
  !
  IF (HREC(1:2)=='D_') THEN
    CALL FMREADN0(CFILE_LFI,'VERSION',CLUOUT_LFI,1,IVERSION,IGRID,ILENCH,HCOMMENT,KRESP)
    CALL FMREADN0(CFILE_LFI,'BUG',CLUOUT_LFI,1,IBUGFIX,IGRID,ILENCH,HCOMMENT,KRESP)
    IF (IVERSION<=2 .OR. (IVERSION==3 .AND. IBUGFIX<=5)) YREC = 'DATA_'//HREC(3:12)
    IF (YREC(13:15)=='SOI') YREC=YREC(1:15)//'L'
    IF (YREC(12:14)=='SOI') YREC=YREC(1:14)//'L'
  END IF
  !---------------------------------------------------------------------------
  !
  IF (HDIR=='H' .OR. HDIR=='A') THEN
    IF (.NOT. LMNH_COMPATIBLE) THEN
      CALL FMREADX2(CFILE_LFI,YREC,CLUOUT_LFI,SIZE(ZWORK),ZWORK(:,:),IGRID,ILENCH,HCOMMENT,KRESP)
     ELSE
      CALL READ_IN_LFI_X2_FOR_MNH(YREC,ZWORK,KRESP,HCOMMENT,HDIR)
    END IF
  ELSE
    CALL FMREADX2(CFILE_LFI,YREC,CLUOUT_LFI,SIZE(ZWORK),ZWORK(:,:),IGRID,ILENCH,HCOMMENT,KRESP)
  END IF
  CALL ERROR_READ_SURF_LFI(YREC,KRESP)
  !
!$OMP END SINGLE COPYPRIVATE(ZWORK,HCOMMENT,KRESP)
  !  
ELSE
  ALLOCATE(ZWORK(0,0))
ENDIF
!
#ifndef NOMPI
XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
!
IF (HDIR=='A') THEN  ! no distribution on other tasks
  IF ( NRANK==NPIO ) THEN    
#ifndef NOMPI          
    XTIME0 = MPI_WTIME()
#endif    
    PFIELD(:,:) = ZWORK(1:SIZE(PFIELD,1),:)
#ifndef NOMPI    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif
  ENDIF
ELSEIF (HDIR=='-') THEN ! distribution of the total field on other tasks
!$OMP SINGLE          
  PFIELD(:,:) = ZWORK(1:SIZE(PFIELD,1),:)
  IF (NPROC>1) THEN
#ifndef NOMPI            
    XTIME0 = MPI_WTIME()
    CALL MPI_BCAST(PFIELD,SIZE(PFIELD)*KIND(PFIELD)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif    
  ENDIF
!$OMP END SINGLE COPYPRIVATE(PFIELD)    
ELSE
  CALL READ_AND_SEND_MPI(ZWORK,PFIELD,NMASK)
ENDIF
!
DEALLOCATE(ZWORK)
!
IF (HDIR=='H' .OR. HDIR=='A') THEN
  CALL GET_SURF_UNDEF(ZUNDEF)
  WHERE(PFIELD==999.) PFIELD=ZUNDEF
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFX2_LFI',1,ZHOOK_HANDLE)
!
CONTAINS
!
!     #############################################################
      SUBROUTINE READ_IN_LFI_X2_FOR_MNH(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to fill a read 2D array for the externalised surface 
!
USE MODD_IO_SURF_LFI,        ONLY : CFILE_LFI, CLUOUT_LFI, &
                                    NIU, NIB, NIE, NJU, NJB, NJE
!
USE MODI_FMREAD
USE MODI_ERROR_READ_SURF_LFI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),        INTENT(IN) :: HREC     ! name of the article to be read
REAL, DIMENSION(:,:), INTENT(OUT):: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(OUT):: HCOMMENT ! comment string
CHARACTER(LEN=1),         INTENT(IN) :: HDIR     ! type of field :
                                                 ! 'H' : field with
                                                 !       horizontal spatial dim.
                                                 ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
! 
INTEGER :: JI, JJ
INTEGER :: IGRID, ILENCH
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWORK3D ! work array read in a MNH file
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK2D
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFX2_LFI:READ_IN_LFI_X2_FOR_MNH',0,ZHOOK_HANDLE)
!
ALLOCATE(ZWORK3D(NIU,NJU,SIZE(PFIELD,2)))
ZWORK3D(:,:,:) = 999.
!
IF (NJB==NJE) THEN
  ALLOCATE(ZWORK2D(NIU,SIZE(PFIELD,2)))
  ZWORK2D(:,:) = 999.
ELSEIF (NIB==NIE) THEN
  ALLOCATE(ZWORK2D(NJU,SIZE(PFIELD,2)))
  ZWORK2D(:,:) = 999. 
ENDIF
!
IF (NJB==NJE) THEN
!
  CALL FMREADX2(CFILE_LFI,YREC,CLUOUT_LFI,SIZE(ZWORK2D),ZWORK2D,IGRID,ILENCH,HCOMMENT,KRESP)
  DO JJ = 1,SIZE(ZWORK3D,2)
    ZWORK3D(NIB:NIE,JJ,:) = ZWORK2D(NIB:NIE,:)
  END DO
!
  DEALLOCATE(ZWORK2D) 
!
ELSEIF (NIB==NIE) THEN
!
  CALL FMREADX2(CFILE_LFI,YREC,CLUOUT_LFI,SIZE(ZWORK2D),ZWORK2D,IGRID,ILENCH,HCOMMENT,KRESP)
  DO JI = 1,SIZE(ZWORK3D,1)
    ZWORK3D(JI,NIB:NIE,:) = ZWORK2D(NJB:NJE,:)
  END DO  
!
  DEALLOCATE(ZWORK2D)
!
ELSE
!
  CALL FMREADX3(CFILE_LFI,HREC,CLUOUT_LFI,SIZE(ZWORK3D),ZWORK3D,IGRID,ILENCH,HCOMMENT,KRESP)
!
ENDIF
!
DO JJ=1,NJE-NJB+1
  DO JI=1,NIE-NIB+1
    PFIELD(JI+(NIE-NIB+1)*(JJ-1),:) = ZWORK3D(NIB+JI-1,NJB+JJ-1,:) 
  END DO
END DO
DEALLOCATE(ZWORK3D)
!  
CALL ERROR_READ_SURF_LFI(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFX2_LFI:READ_IN_LFI_X2_FOR_MNH',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_IN_LFI_X2_FOR_MNH
!
END SUBROUTINE READ_SURFX2_LFI
!
!     #############################################################
      SUBROUTINE READ_SURFN0_LFI(HREC,KFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!!      B. Decharme 07/2011 : Grdid dimension only read in pgd file
!
USE MODD_IO_SURF_LFI, ONLY : CFILE_LFI, CLUOUT_LFI, CFILEPGD_LFI, &
                             LMNH_COMPATIBLE, NIU, NIB, NIE, NJU, NJB, NJE
!
USE MODD_SURFEX_OMP, ONLY : NBLOCK
!
USE MODI_FMREAD
USE MODI_ERROR_READ_SURF_LFI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,            INTENT(OUT) :: KFIELD   ! the integer to be read
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=40) :: YGRID
INTEGER           :: IGRID   ! position of data on grid
INTEGER           :: ILENCH  ! length of comment string
INTEGER           :: IIMAX, IJMAX
INTEGER           :: INB ! number of articles in the file
INTEGER           :: IRESP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFN0_LFI',0,ZHOOK_HANDLE)
!
KRESP=0
!
CALL FMREADN0(CFILE_LFI,HREC,CLUOUT_LFI,1,KFIELD,IGRID,ILENCH,HCOMMENT,KRESP)
!
CALL ERROR_READ_SURF_LFI(HREC,KRESP)
!
!* tests compatibility with MesoNH files
!
IF (HREC/='DIM_FULL' .AND. LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFN0_LFI',1,ZHOOK_HANDLE)
IF (HREC/='DIM_FULL') RETURN
!
!-----------------------------------------------------------------------------------------------------
! READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!
IF (CFILE_LFI/=CFILEPGD_LFI) THEN
  CALL FMOPEN(CFILEPGD_LFI,'OLD',CLUOUT_LFI,0,1,1,INB,IRESP)
ENDIF
!
CALL FMREADC0(CFILEPGD_LFI,'GRID_TYPE ',CLUOUT_LFI,1,YGRID,IGRID,ILENCH,HCOMMENT,KRESP)
CALL ERROR_READ_SURF_LFI('GRID_TYPE ',KRESP)
LMNH_COMPATIBLE = (YGRID=="CARTESIAN " .OR. YGRID=="CONF PROJ ")
!
IF (LMNH_COMPATIBLE) THEN
  CALL FMREADN0(CFILEPGD_LFI,'IMAX',CLUOUT_LFI,1,IIMAX,IGRID,ILENCH,HCOMMENT,KRESP)
  CALL ERROR_READ_SURF_LFI('IMAX',KRESP)
  NIU = IIMAX+2
  NIB = 2
  NIE = IIMAX+1
  CALL FMREADN0(CFILEPGD_LFI,'JMAX',CLUOUT_LFI,1,IJMAX,IGRID,ILENCH,HCOMMENT,KRESP)
  CALL ERROR_READ_SURF_LFI('JMAX',KRESP)
  NJU = IJMAX+2
  NJB = 2
  NJE = IJMAX+1
END IF
!
IF(CFILE_LFI/=CFILEPGD_LFI)THEN
  CALL FMCLOS(CFILEPGD_LFI,'KEEP',CLUOUT_LFI,IRESP)
ENDIF
!
!-----------------------------------------------------------------------------------------------------
! END READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFN0_LFI',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFN0_LFI
!
!     #############################################################
      SUBROUTINE READ_SURFN1_LFI(HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ
!
USE MODD_IO_SURF_LFI,        ONLY : CFILE_LFI, CLUOUT_LFI, NMASK, NFULL
!
USE MODI_FMREAD
USE MODI_ERROR_READ_SURF_LFI
USE MODI_READ_AND_SEND_MPI
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
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),      INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:),  INTENT(OUT) :: KFIELD   ! the integer to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
                                                ! 'H' : field with
                                                !       horizontal spatial dim.
                                                ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
INTEGER          :: IGRID   ! position of data on grid
INTEGER          :: ILENCH  ! length of comment string
INTEGER          :: I, INFOMPI
!
#ifndef NOMPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK  ! work array read in the file
DOUBLE PRECISION   :: XTIME0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFN1_LFI',0,ZHOOK_HANDLE)
!
KRESP=0
!
#ifndef NOMPI
XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  !
  IF (HDIR=='H' .OR. HDIR=='A') THEN
    ALLOCATE(IWORK(NFULL))
  ELSE
    ALLOCATE(IWORK(SIZE(KFIELD)))
  ENDIF
  !
!$OMP SINGLE        
  !  
  IF (HDIR=='H' .OR. HDIR=='A')  THEN
    CALL FMREADN1(CFILE_LFI,HREC,CLUOUT_LFI,NFULL,IWORK,IGRID,ILENCH,HCOMMENT,KRESP)
  ELSE
    CALL FMREADN1(CFILE_LFI,HREC,CLUOUT_LFI,SIZE(IWORK),IWORK(:),IGRID,ILENCH,HCOMMENT,KRESP)
  END IF
  !
!$OMP END SINGLE COPYPRIVATE(IWORK,HCOMMENT,KRESP)  
  !  
  CALL ERROR_READ_SURF_LFI(HREC,KRESP)
  !
ENDIF
!
#ifndef NOMPI
XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
!
IF (HDIR=='A') THEN  ! no distribution on other tasks
  IF ( NRANK==NPIO ) THEN
#ifndef NOMPI          
    XTIME0 = MPI_WTIME()
#endif    
    KFIELD(:) = IWORK(1:SIZE(KFIELD))
#ifndef NOMPI    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif 
  ENDIF
ELSEIF (HDIR=='-') THEN ! distribution of the total field on other tasks
!$OMP SINGLE     
  KFIELD(:) = IWORK(1:SIZE(KFIELD))
  IF (NPROC>1) THEN
#ifndef NOMPI           
    XTIME0 = MPI_WTIME()
    CALL MPI_BCAST(KFIELD,SIZE(KFIELD)*KIND(KFIELD)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif     
  ENDIF
!$OMP END SINGLE COPYPRIVATE(KFIELD)    
ELSE
  CALL READ_AND_SEND_MPI(IWORK,KFIELD,NMASK)
ENDIF
!
!$OMP SINGLE
DEALLOCATE(IWORK)
!$OMP END SINGLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFN1_LFI',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFN1_LFI
!
!     #############################################################
      SUBROUTINE READ_SURFC0_LFI(HREC,HFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READC0* - routine to read a character
!
USE MODD_IO_SURF_LFI,        ONLY : CFILE_LFI, CLUOUT_LFI
!
USE MODI_FMREAD
USE MODI_ERROR_READ_SURF_LFI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),  INTENT(IN)  :: HREC      ! name of the article to be read
CHARACTER(LEN=40),  INTENT(OUT) :: HFIELD    ! the integer to be read
INTEGER,            INTENT(OUT) :: KRESP     ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT  ! comment
!
!*      0.2   Declarations of local variables
!
INTEGER          :: IGRID   ! position of data on grid
INTEGER          :: ILENCH  ! length of comment string
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFC0_LFI',0,ZHOOK_HANDLE)
!
KRESP=0
!
CALL FMREADC0(CFILE_LFI,HREC,CLUOUT_LFI,1,HFIELD,IGRID,ILENCH,HCOMMENT,KRESP)
!
CALL ERROR_READ_SURF_LFI(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFC0_LFI',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFC0_LFI
!
!     #############################################################
      SUBROUTINE READ_SURFL0_LFI(HREC,OFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READL0* - routine to read a logical
!
USE MODD_IO_SURF_LFI,        ONLY : CFILE_LFI, CLUOUT_LFI
!
USE MODI_FMREAD
USE MODI_ERROR_READ_SURF_LFI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),  INTENT(IN)  :: HREC     ! name of the article to be read
LOGICAL,            INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
INTEGER          :: IGRID   ! position of data on grid
INTEGER          :: ILENCH  ! length of comment string
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFL0_LFI',0,ZHOOK_HANDLE)
!
KRESP=0
!
CALL FMREADL0(CFILE_LFI,HREC,CLUOUT_LFI,1,OFIELD,IGRID,ILENCH,HCOMMENT,KRESP)
!
CALL ERROR_READ_SURF_LFI(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFL0_LFI',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFL0_LFI
!
!     #############################################################
      SUBROUTINE READ_SURFL1_LFI(HREC,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READL1* - routine to read a logical array
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ
!
USE MODD_IO_SURF_LFI,        ONLY : CFILE_LFI, CLUOUT_LFI
!
USE MODI_FMREAD
USE MODI_ERROR_READ_SURF_LFI
USE MODI_ABOR1_SFX
USE MODI_GET_LUOUT
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
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),      INTENT(IN)  :: HREC     ! name of the article to be read
LOGICAL, DIMENSION(:), INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
                                                ! 'H' : field with
                                                !       horizontal spatial dim.
                                                ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
INTEGER :: ILUOUT
INTEGER :: IGRID   ! position of data on grid
INTEGER :: ILENCH  ! length of comment string
INTEGER :: INFOMPI
DOUBLE PRECISION   :: XTIME0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFL1_LFI',0,ZHOOK_HANDLE)
!
KRESP=0
!
#ifndef NOMPI
XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  !
!$OMP SINGLE
  !   
  IF (HDIR=='H') THEN
    CALL GET_LUOUT('LFI   ',ILUOUT)
    WRITE(ILUOUT,*) 'Error: 1D logical vector for reading on an horizontal grid:'
    WRITE(ILUOUT,*) 'this option is not coded in READ_SURFL1_LFI'
    CALL ABOR1_SFX('MODE_READ_SURF_LFI: 1D LOGICAL VECTOR FOR READING NOT CODED IN READ_SURFL1_LFI')
  END IF
  !
  CALL FMREADL1(CFILE_LFI,HREC,CLUOUT_LFI,SIZE(OFIELD),OFIELD,IGRID,ILENCH,HCOMMENT,KRESP)
  !
!$OMP END SINGLE COPYPRIVATE(OFIELD,HCOMMENT,KRESP)
  !
  CALL ERROR_READ_SURF_LFI(HREC,KRESP)
  !
ENDIF
!
#ifndef NOMPI
XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
!
IF (NPROC>1 .AND. HDIR/='A') THEN
#ifndef NOMPI        
  XTIME0 = MPI_WTIME()
!$OMP SINGLE
  CALL MPI_BCAST(OFIELD,SIZE(OFIELD),MPI_LOGICAL,NPIO,NCOMM,INFOMPI)
!$OMP END SINGLE COPYPRIVATE(OFIELD)   
  XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0) 
#endif
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFL1_LFI',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFL1_LFI
!
!     #############################################################
      SUBROUTINE READ_SURFT0_LFI(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT0* - routine to read a date
!
USE MODD_IO_SURF_LFI,  ONLY : CFILE_LFI, CLUOUT_LFI  
!
USE MODI_FMREAD
USE MODI_ERROR_READ_SURF_LFI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,            INTENT(OUT) :: KYEAR    ! year
INTEGER,            INTENT(OUT) :: KMONTH   ! month
INTEGER,            INTENT(OUT) :: KDAY     ! day
REAL,               INTENT(OUT) :: PTIME    ! year
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment

!*      0.2   Declarations of local variables
!
CHARACTER(LEN=18)     :: YREC     ! Name of the article to be read
INTEGER, DIMENSION(3) :: ITDATE
!
INTEGER          :: IGRID   ! position of data on grid
INTEGER          :: ILENCH  ! length of comment string
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFT0_LFI',0,ZHOOK_HANDLE)
!
KRESP=0
!
YREC=TRIM(HREC)//'%TDATE'
CALL FMREADN1(CFILE_LFI,YREC,CLUOUT_LFI,3,ITDATE,IGRID,ILENCH,HCOMMENT,KRESP)
CALL ERROR_READ_SURF_LFI(HREC,KRESP)
!
YREC=TRIM(HREC)//'%TIME'
CALL FMREADX0(CFILE_LFI,YREC,CLUOUT_LFI,1,PTIME,IGRID,ILENCH,HCOMMENT,KRESP)
CALL ERROR_READ_SURF_LFI(HREC,KRESP)
!
KYEAR  = ITDATE(1)
KMONTH = ITDATE(2)
KDAY   = ITDATE(3)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFT0_LFI',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT0_LFI
!
!     #############################################################
      SUBROUTINE READ_SURFT1_LFI(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT0* - routine to read a date
!
USE MODD_IO_SURF_LFI,        ONLY : CFILE_LFI, CLUOUT_LFI  
!
USE MODI_FMREAD
USE MODI_ERROR_READ_SURF_LFI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),     INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:), INTENT(OUT) :: KYEAR    ! year
INTEGER, DIMENSION(:), INTENT(OUT) :: KMONTH   ! month
INTEGER, DIMENSION(:), INTENT(OUT) :: KDAY     ! day
REAL,    DIMENSION(:), INTENT(OUT) :: PTIME    ! year
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment

!*      0.2   Declarations of local variables
!
CHARACTER(LEN=18)     :: YREC     ! Name of the article to be read
INTEGER          :: ILUOUT
INTEGER          :: IGRID   ! position of data on grid
INTEGER          :: ILENCH  ! length of comment string
INTEGER, DIMENSION(3,SIZE(KYEAR)) :: ITDATE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFT1_LFI',0,ZHOOK_HANDLE)
!
KRESP=0
!
YREC=TRIM(HREC)//'%TDATE'
CALL FMREADN2(CFILE_LFI,YREC,CLUOUT_LFI,SIZE(ITDATE),ITDATE,IGRID,ILENCH,HCOMMENT,KRESP)
CALL ERROR_READ_SURF_LFI(HREC,KRESP)
!
YREC=TRIM(HREC)//'%TIME'
CALL FMREADX1(CFILE_LFI,YREC,CLUOUT_LFI,SIZE(PTIME),PTIME,IGRID,ILENCH,HCOMMENT,KRESP)
CALL ERROR_READ_SURF_LFI(HREC,KRESP)
!
KYEAR (:) = ITDATE(1,:)
KMONTH(:) = ITDATE(2,:)
KDAY  (:) = ITDATE(3,:)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_LFI:READ_SURFT1_LFI',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT1_LFI
!
END MODULE MODE_READ_SURF_LFI
