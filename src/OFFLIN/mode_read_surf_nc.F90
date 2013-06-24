MODULE MODE_READ_SURF_NC
!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURF_NC is
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
!!      F. Habets      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
INTERFACE READ_SURF0_NC
      SUBROUTINE READ_SURFX0_NC(HREC,PFIELD,KRESP,HCOMMENT)
 CHARACTER(LEN=*),  INTENT(IN)  :: HREC     ! name of the article to be read
REAL,               INTENT(OUT) :: PFIELD   ! the real scalar to be read
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
END SUBROUTINE READ_SURFX0_NC
      SUBROUTINE READ_SURFN0_NC(HREC,KFIELD,KRESP,HCOMMENT)
 CHARACTER(LEN=*),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,            INTENT(OUT) :: KFIELD   ! the integer scalar to be read
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
END SUBROUTINE READ_SURFN0_NC
      SUBROUTINE READ_SURFC0_NC(HREC,HFIELD,KRESP,HCOMMENT)
 CHARACTER(LEN=*),   INTENT(IN)  :: HREC     ! name of the article to be read
 CHARACTER(LEN=40),   INTENT(OUT) :: HFIELD   ! the integer scalar to be read
INTEGER,             INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),  INTENT(OUT) :: HCOMMENT ! comment
END SUBROUTINE READ_SURFC0_NC
      SUBROUTINE READ_SURFL0_NC(HREC,OFIELD,KRESP,HCOMMENT)
 CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
LOGICAL,                  INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
END SUBROUTINE READ_SURFL0_NC
END INTERFACE
INTERFACE READ_SURFN_NC
      SUBROUTINE READ_SURFX1_NC(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
 CHARACTER(LEN=*),  INTENT(IN)  :: HREC     ! name of the article to be read
REAL, DIMENSION(:), INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),   INTENT(IN)  :: HDIR     ! type of field :
END SUBROUTINE READ_SURFX1_NC
      SUBROUTINE READ_SURFX2_NC(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
 CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
REAL, DIMENSION(:,:),     INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),         INTENT(IN)  :: HDIR     ! type of field :
END SUBROUTINE READ_SURFX2_NC
      SUBROUTINE READ_SURFN1_NC(HREC,KFIELD,KRESP,HCOMMENT,HDIR)
 CHARACTER(LEN=*),      INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:), INTENT(OUT) :: KFIELD   ! the integer scalar to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
END SUBROUTINE READ_SURFN1_NC
      SUBROUTINE READ_SURFN2_NC(HREC,KFIELD,KRESP,HCOMMENT,HDIR)
 CHARACTER(LEN=*),      INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KFIELD   ! the integer scalar to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
END SUBROUTINE READ_SURFN2_NC
      SUBROUTINE READ_SURFL1_NC(HREC,OFIELD,KRESP,HCOMMENT,HDIR)
 CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
LOGICAL, DIMENSION(:),   INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),         INTENT(IN)  :: HDIR     ! type of field :
END SUBROUTINE READ_SURFL1_NC
END INTERFACE
INTERFACE READ_SURFT_NC
      SUBROUTINE READ_SURFT0_NC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
 CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                  INTENT(OUT) :: KYEAR    ! year
INTEGER,                  INTENT(OUT) :: KMONTH   ! month
INTEGER,                  INTENT(OUT) :: KDAY     ! day
REAL,                     INTENT(OUT) :: PTIME    ! time
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
END SUBROUTINE READ_SURFT0_NC
      SUBROUTINE READ_SURFT1_NC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
 CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:),    INTENT(OUT) :: KYEAR    ! year
INTEGER, DIMENSION(:),    INTENT(OUT) :: KMONTH   ! month
INTEGER, DIMENSION(:),    INTENT(OUT) :: KDAY     ! day
REAL, DIMENSION(:),       INTENT(OUT) :: PTIME    ! time
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
END SUBROUTINE READ_SURFT1_NC
      SUBROUTINE READ_SURFT2_NC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
 CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:,:),  INTENT(OUT) :: KYEAR    ! year
INTEGER, DIMENSION(:,:),  INTENT(OUT) :: KMONTH   ! month
INTEGER, DIMENSION(:,:),  INTENT(OUT) :: KDAY     ! day
REAL, DIMENSION(:,:),     INTENT(OUT) :: PTIME    ! time
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
END SUBROUTINE READ_SURFT2_NC

END INTERFACE
!
END MODULE MODE_READ_SURF_NC
!
!     #############################################################
      SUBROUTINE READ_SURFX0_NC(HREC,PFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READX0* - routine to read a real scalar
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE MODD_IO_SURF_NC, ONLY : NID_NC
!
USE MODI_ERROR_READ_SURF_NC
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),  INTENT(IN)  :: HREC     ! name of the article to be read
REAL,               INTENT(OUT) :: PFIELD   ! the real scalar to be read
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
REAL*4 :: ZFIELD
 CHARACTER(LEN=100) :: YFILE          ! filename
INTEGER            :: IVAR_ID,JRET,IVAL,ITYPE,INDIMS
INTEGER,DIMENSION(4) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFX0_NC',0,ZHOOK_HANDLE)
!
KRESP=0
HCOMMENT = " "
!
IF (NID_NC.NE.0) THEN
  !       
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (NID_NC,HREC,IVAR_ID)
  IRET(1)=NF_INQ_VARTYPE (NID_NC,IVAR_ID,ITYPE)
  IRET(1)=NF_INQ_VARNDIMS(NID_NC,IVAR_ID,INDIMS)
  !  
  ! 2. Get variable
  !----------------------------
  IF (ITYPE==NF_DOUBLE) THEN
    IRET(2)=NF_GET_VAR_DOUBLE(NID_NC,IVAR_ID,PFIELD)
  ELSEIF (ITYPE==NF_FLOAT) THEN
    IRET(2)=NF_GET_VAR_REAL(NID_NC,IVAR_ID,ZFIELD)
    PFIELD = ZFIELD
  ENDIF
  !  
ENDIF
!
! 3. Check for errors
!--------------------
DO JRET=1,2
  IF ((PFIELD==XUNDEF).OR.(NID_NC==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    PFIELD=XUNDEF
    KRESP=1
  ENDIF
ENDDO
!     
IF (KRESP /=0) CALL ERROR_READ_SURF_NC(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFX0_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX0_NC
!
!     #############################################################
      SUBROUTINE READ_SURFX1_NC(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX1* - routine to fill a real 1D array for the externalised surface 
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, XTIME_NPIO_READ, NPROC, NCOMM, &
                                 XTIME_NPIO_READ, XTIME_COMM_READ
!
USE MODD_SURFEX_OMP, ONLY : XWORKD, NWORKB
!
USE MODD_IO_SURF_NC, ONLY: LMASK,NMASK,NID_NC
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE MODI_ERROR_READ_SURF_NC
USE MODI_READ_AND_SEND_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),  INTENT(IN)  :: HREC     ! name of the article to be read
REAL, DIMENSION(:), INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),   INTENT(IN)  :: HDIR     ! type of field :
                                            ! 'H' : field with
                                            !       horizontal spatial dim.
                                            ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=100) :: YFILE,YOUT          ! Filename
 CHARACTER(LEN=100)    :: YNAME
INTEGER :: IL1, IVAR_ID,JRET,JDIM,INDIMS, ITYPE, INFOMPI
INTEGER,DIMENSION(4) :: IDIMIDS,IDIMLEN
INTEGER,DIMENSION(4) :: IRET
!
REAL*4, DIMENSION(:), ALLOCATABLE :: ZTAB_1D4
DOUBLE PRECISION   :: XTIME0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFX1_NC',0,ZHOOK_HANDLE)
!
!$OMP SINGLE
NWORKB=0
!$OMP END SINGLE
!
HCOMMENT = " "
!
#ifndef NOMPI
XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  !
!$OMP SINGLE
  !  
  IF (NID_NC.NE.0) THEN
    !  
    ! 1. Find id of the variable
    !----------------------------
    IRET(1)=NF_INQ_VARID   (NID_NC,HREC,IVAR_ID)
    IRET(2)=NF_INQ_VARTYPE (NID_NC,IVAR_ID,ITYPE)
    IRET(3)=NF_INQ_VARNDIMS(NID_NC,IVAR_ID,INDIMS)
    IRET(4)=NF_INQ_VARDIMID(NID_NC,IVAR_ID,IDIMIDS(1:INDIMS))
    IDIMLEN(:) = 1.
    DO JDIM=1,INDIMS
      JRET=NF_INQ_DIMLEN(NID_NC,IDIMIDS(JDIM),IDIMLEN(JDIM))
    ENDDO
    !
    IRET(4)=NF_INQ_DIMNAME(NID_NC,IDIMIDS(1),YNAME)
    !
    IF (TRIM(YNAME).NE.'Number_of_points') THEN
      ALLOCATE(XWORKD(IDIMLEN(1)*IDIMLEN(2)))
    ELSE
      ALLOCATE(XWORKD(IDIMLEN(1)))
    ENDIF
    !
    ! 2. Get variable
    !----------------------------
    IF (ITYPE==NF_DOUBLE) THEN
      IRET(1)=NF_GET_VAR_DOUBLE(NID_NC,IVAR_ID,XWORKD)
    ELSEIF (ITYPE==NF_FLOAT) THEN
      ALLOCATE(ZTAB_1D4(SIZE(XWORKD)))
      IRET(2)=NF_GET_VAR_REAL(NID_NC,IVAR_ID,ZTAB_1D4)
      XWORKD(:) = ZTAB_1D4(:)
      DEALLOCATE(ZTAB_1D4)
    ENDIF            
  ENDIF
  !
  ! 3. Check for errors
  !--------------------
  DO JRET=1,1
    IF ((NID_NC==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
      XWORKD = XUNDEF
      NWORKB=1
    ENDIF
  ENDDO
  !
!$OMP END SINGLE
  !
  IF (NWORKB /=0) CALL ERROR_READ_SURF_NC(HREC,NWORKB)
  !
ELSE
!$OMP SINGLE
  ALLOCATE(XWORKD(0))
!$OMP END SINGLE
ENDIF
!
KRESP = NWORKB
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
    PFIELD(:) = XWORKD(:)
#ifndef NOMPI    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif    
  ENDIF
ELSEIF (HDIR=='-') THEN ! distribution of the total field on other tasks
!$OMP SINGLE    
#ifndef NOMPI
  IF (NPROC>1) THEN
    XTIME0 = MPI_WTIME()
    CALL MPI_BCAST(XWORKD,SIZE(XWORKD)*KIND(XWORKD)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)   
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
  ENDIF
#endif    
!$OMP END SINGLE
  PFIELD(:) = XWORKD(:)
ELSE
  IF (LMASK) THEN
    CALL READ_AND_SEND_MPI(XWORKD,PFIELD,NMASK)
  ELSE 
    CALL READ_AND_SEND_MPI(XWORKD,PFIELD)
  END IF
ENDIF
!
!$OMP BARRIER
!
!$OMP SINGLE
DEALLOCATE(XWORKD)
!$OMP END SINGLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFX1_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX1_NC
!
!     #############################################################
      SUBROUTINE READ_SURFX2_NC(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX2* - routine to fill a real 2D array for the externalised surface 
!
USE MODD_SURFEX_MPI, ONLY: NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ
!
USE MODD_SURFEX_OMP, ONLY : XWORKD2, NWORKB
!
USE MODD_IO_SURF_NC, ONLY: LMASK,NMASK,NID_NC
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE MODI_ERROR_READ_SURF_NC
USE MODI_READ_AND_SEND_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
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
 CHARACTER(LEN=100) :: YFILE,YOUT          ! filename
 CHARACTER(LEN=100)    :: YNAME 
INTEGER            :: IVAR_ID,JRET,JDIM,INDIMS,ITYPE, INFOMPI
INTEGER,DIMENSION(4) :: IDIMIDS,IDIMLEN
INTEGER,DIMENSION(4) :: IRET
REAL*4, DIMENSION(:,:), ALLOCATABLE :: ZTAB_2D4
DOUBLE PRECISION   :: XTIME0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFX2_NC',0,ZHOOK_HANDLE)
!
!$OMP SINGLE
NWORKB=0
!$OMP END SINGLE
!
HCOMMENT = " "
!
#ifndef NOMPI
XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  !
!$OMP SINGLE
  !  
  IF (NID_NC.NE.0) THEN
    !   
    ! 1. Find id of the variable
    !----------------------------
    IRET(1)=NF_INQ_VARID   (NID_NC,HREC,IVAR_ID)
    IRET(2)=NF_INQ_VARTYPE (NID_NC,IVAR_ID,ITYPE)
    IRET(3)=NF_INQ_VARNDIMS(NID_NC,IVAR_ID,INDIMS)
    IRET(4)=NF_INQ_VARDIMID(NID_NC,IVAR_ID,IDIMIDS(1:INDIMS))
    IDIMLEN(:) = 1.
    DO JDIM=1,INDIMS
      JRET=NF_INQ_DIMLEN(NID_NC,IDIMIDS(JDIM),IDIMLEN(JDIM))
    ENDDO
    !
    IRET(4)=NF_INQ_DIMNAME(NID_NC,IDIMIDS(1),YNAME)    
    ! 
    ! 2. Get variable
    !----------------------------
    IF (TRIM(YNAME).NE.'Number_of_points') THEN
      ALLOCATE(XWORKD2(IDIMLEN(1)*IDIMLEN(2),IDIMLEN(3)))
    ELSE
      ALLOCATE(XWORKD2(IDIMLEN(1),IDIMLEN(2)))
    ENDIF
    !
    IF (ITYPE==NF_DOUBLE) THEN
      IRET(2)=NF_GET_VAR_DOUBLE(NID_NC,IVAR_ID,XWORKD2)
    ELSEIF (ITYPE==NF_FLOAT) THEN
      ALLOCATE(ZTAB_2D4(PRODUCT(IDIMLEN(1:INDIMS-1)),IDIMLEN(INDIMS)))
      IRET(2)=NF_GET_VAR_REAL(NID_NC,IVAR_ID,ZTAB_2D4)
      XWORKD2(:,:) = ZTAB_2D4(:,:)
      DEALLOCATE(ZTAB_2D4)
    ENDIF      

  ENDIF

  ! 3. Check for errors
  !--------------------
  DO JRET=1,2
    IF ((NID_NC==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
      XWORKD2 = XUNDEF
      NWORKB=1
    ENDIF
  ENDDO
  !
!$OMP END SINGLE
  !
  IF (NWORKB /=0) CALL ERROR_READ_SURF_NC(HREC,NWORKB)
  !
ELSE
!$OMP SINGLE
  ALLOCATE(XWORKD2(0,0))
!$OMP END SINGLE
ENDIF
!
KRESP = NWORKB
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
    PFIELD(:,:) = XWORKD2(:,:)
#ifndef NOMPI    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif    
  ENDIF
ELSEIF (HDIR=='-') THEN ! distribution of the total field on other tasks
!$OMP SINGLE    
#ifndef NOMPI
  IF (NPROC>1) THEN
    XTIME0 = MPI_WTIME()
    CALL MPI_BCAST(XWORKD2,SIZE(XWORKD2)*KIND(XWORKD2)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)   
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
  ENDIF
#endif    
!$OMP END SINGLE
  PFIELD(:,:) = XWORKD2(:,:)
ELSE
  IF (LMASK) THEN
    CALL READ_AND_SEND_MPI(XWORKD2,PFIELD,NMASK)
  ELSE 
    CALL READ_AND_SEND_MPI(XWORKD2,PFIELD)
  END IF
ENDIF
!
!$OMP BARRIER
!
!$OMP SINGLE
DEALLOCATE(XWORKD2) 
!$OMP END SINGLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFX2_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX2_NC
!
!     #############################################################
      SUBROUTINE READ_SURFN0_NC(HREC,KFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!
USE MODD_SURF_PAR,   ONLY: NUNDEF
!
USE MODD_IO_SURF_NC, ONLY : NID_NC
!
USE MODI_ERROR_READ_SURF_NC
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,            INTENT(OUT) :: KFIELD   ! the integer scalar to be read
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
!
!
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=100):: YFILE          ! filename
INTEGER :: IVAR_ID,JRET,JDIM,INDIMS
INTEGER,DIMENSION(4) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFN0_NC',0,ZHOOK_HANDLE)
!
KRESP=0
HCOMMENT = " "
!
IF (NID_NC.NE.0) THEN
  !        
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (NID_NC,HREC,IVAR_ID)
  !  
  ! 2. Get variable
  !----------------------------
  IRET(2)=NF_GET_VAR_INT(NID_NC,IVAR_ID,KFIELD)
  !  
ENDIF
!
! 3. Check for errors
!--------------------
DO JRET=1,2
  IF ((KFIELD==NUNDEF).OR.(NID_NC==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    KFIELD=NUNDEF
    KRESP=1
  ENDIF
ENDDO
!
IF (KRESP /=0)  CALL ERROR_READ_SURF_NC(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFN0_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFN0_NC
!
!     #############################################################
      SUBROUTINE READ_SURFN1_NC(HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),      INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:), INTENT(OUT) :: KFIELD   ! the integer scalar to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
                                                ! 'H' : field with
                                                !       horizontal spatial dim.
                                                ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(KFIELD)) :: ZFIELD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFN1_NC',0,ZHOOK_HANDLE)
!
 CALL READ_SURFX1_NC(HREC,ZFIELD,KRESP,HCOMMENT,HDIR)
KFIELD = NINT(ZFIELD)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFN1_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFN1_NC
!
!     #############################################################
      SUBROUTINE READ_SURFN2_NC(HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),      INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KFIELD   ! the integer scalar to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
                                                ! 'H' : field with
                                                !       horizontal spatial dim.
                                                ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(KFIELD,1),SIZE(KFIELD,2)) :: ZFIELD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFN2_NC',0,ZHOOK_HANDLE)
!
 CALL READ_SURFX2_NC(HREC,ZFIELD,KRESP,HCOMMENT,HDIR)
KFIELD = NINT(ZFIELD)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFN2_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFN2_NC
!
!     #############################################################
      SUBROUTINE READ_SURFC0_NC(HREC,HFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READC0* - routine to read a STRING
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE MODD_IO_SURF_NC, ONLY : NID_NC
!
USE MODI_ERROR_READ_SURF_NC
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),   INTENT(IN)  :: HREC     ! name of the article to be read
 CHARACTER(LEN=40),   INTENT(OUT) :: HFIELD   ! the integer scalar to be read
INTEGER,             INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),  INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=100):: YFILE          ! filename
 CHARACTER(LEN=40):: YFIELD   
INTEGER :: IVAR_ID,JRET,JDIM,INDIMS
INTEGER,DIMENSION(4) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFC0_NC',0,ZHOOK_HANDLE)
!
KRESP=0
HCOMMENT = " "
!
IF (NID_NC.NE.0) THEN
  !       
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (NID_NC,HREC,IVAR_ID)
  !  
  ! 2. Get variable
  !----------------------------
  IRET(2)=NF_GET_VAR_TEXT(NID_NC,IVAR_ID,YFIELD)
  HFIELD=YFIELD(:LEN_TRIM(YFIELD))
  !  
ENDIF

! 3. Check for errors
!--------------------
DO JRET=1,2
  IF ((NID_NC==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    KRESP=1
  ENDIF
ENDDO  
!
IF (KRESP /=0) CALL ERROR_READ_SURF_NC(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFC0_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFC0_NC
!
!     #############################################################
      SUBROUTINE READ_SURFL0_NC(HREC,OFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READL0* - routine to read a logical
!    
USE MODD_IO_SURF_NC, ONLY : NID_NC
!
USE MODI_ERROR_READ_SURF_NC
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
LOGICAL,                  INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
INTEGER   :: IFIELD   ! work array read in the file
 CHARACTER(LEN=100) :: YFILE    ! Filename
INTEGER :: IVAR_ID,JRET
INTEGER,DIMENSION(2) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFL0_NC',0,ZHOOK_HANDLE)
!
KRESP=0
HCOMMENT = " "
!
IF (NID_NC.NE.0) THEN
  !       
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (NID_NC,HREC,IVAR_ID)
  !  
  ! 2. Get variable
  !----------------------------
  IRET(2)=NF_GET_VAR_INT(NID_NC,IVAR_ID,IFIELD)
  !  
  IF (IFIELD ==1) OFIELD=.TRUE.
  IF (IFIELD ==0) OFIELD=.FALSE.
  !
ENDIF
!
! 3. Check for errors
!--------------------
IF ((NID_NC==0).OR.IRET(1).NE.NF_NOERR) THEN 
  KRESP=1
ENDIF
!
IF (KRESP /=0)  CALL ERROR_READ_SURF_NC(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFL0_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFL0_NC
!
!     #############################################################
      SUBROUTINE READ_SURFL1_NC(HREC,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READL1* - routine to read a logical array
!    
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ
!
USE MODD_SURFEX_OMP, ONLY : LWORKD, NWORKB
!
USE MODD_IO_SURF_NC, ONLY : NID_NC
!
USE MODI_ERROR_READ_SURF_NC
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
LOGICAL, DIMENSION(:),   INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),         INTENT(IN)  :: HDIR     ! type of field :
                                                  ! 'H' : field with
                                                  !       horizontal spatial dim.
                                                  ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=100) :: YFILE          ! Filename
INTEGER, DIMENSION(:), ALLOCATABLE :: ITAB_1D  ! work array read in the file
!
INTEGER :: IVAR_ID,JRET,JDIM,INDIMS
INTEGER :: INFOMPI
INTEGER,DIMENSION(1) :: IDIMIDS,IDIMLEN
INTEGER,DIMENSION(2) :: IRET
INTEGER, DIMENSION(:),    POINTER     :: IMASK    ! 1D mask to read only interesting
DOUBLE PRECISION   :: XTIME0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFL1_NC',0,ZHOOK_HANDLE)
!
!$OMP BARRIER
!
!$OMP SINGLE
NWORKB=0
!$OMP END SINGLE
!
HCOMMENT = " "
!
#ifndef NOMPI
XTIME0 = MPI_WTIME()
#endif
!
!$OMP SINGLE
ALLOCATE(LWORKD(SIZE(OFIELD)))
!$OMP END SINGLE
!
IF (NRANK==NPIO) THEN
  !
!$OMP SINGLE
  ! 
  IF (NID_NC.NE.0) THEN
    !   
    ! 1. Find id of the variable
    !----------------------------
    IRET(1)=NF_INQ_VARID   (NID_NC,HREC,IVAR_ID)
    IRET(1)=NF_INQ_VARNDIMS(NID_NC,IVAR_ID,INDIMS)
    IRET(1)=NF_INQ_VARDIMID(NID_NC,IVAR_ID,IDIMIDS)
    DO JDIM=1,INDIMS
      JRET=NF_INQ_DIMLEN(NID_NC,IDIMIDS(JDIM),IDIMLEN(JDIM))
    ENDDO
    ALLOCATE(ITAB_1D(IDIMLEN(1)))
    !  
    ! 2. Get variable
    !----------------------------
    IRET(1)=NF_GET_VAR_INT(NID_NC,IVAR_ID,ITAB_1D)
    !
    DO JRET=1,SIZE(LWORKD)
      IF (ITAB_1D(JRET) ==1) LWORKD(JRET)=.TRUE.
      IF (ITAB_1D(JRET) ==0) LWORKD(JRET)=.FALSE.
    ENDDO
    !
  ENDIF
  !
  ! 3. Check for errors
  !--------------------
  DO JRET=1,1
    IF ((NID_NC==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
      NWORKB=1
    ENDIF
  ENDDO
  !
  DEALLOCATE(ITAB_1D)
  !
!$OMP END SINGLE
  !  
  IF (NWORKB /=0) CALL ERROR_READ_SURF_NC(HREC,NWORKB)
  !
ENDIF
!
KRESP = NWORKB
!
#ifndef NOMPI
XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
!
IF (NPROC>1) THEN
#ifndef NOMPI
  XTIME0 = MPI_WTIME()
!$OMP SINGLE  
  CALL MPI_BCAST(LWORKD,SIZE(LWORKD),MPI_LOGICAL,NPIO,NCOMM,INFOMPI)
!$OMP END SINGLE
  XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif
ENDIF
!
OFIELD = LWORKD
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFL1_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFL1_NC
!
!
!     #############################################################
      SUBROUTINE READ_SURFT0_NC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT0* - routine to read a NETCDF  date_time scalar
!
USE MODD_IO_SURF_NC, ONLY : NID_NC
!
USE MODI_ERROR_READ_SURF_NC
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                  INTENT(OUT) :: KYEAR    ! year
INTEGER,                  INTENT(OUT) :: KMONTH   ! month
INTEGER,                  INTENT(OUT) :: KDAY     ! day
REAL,                     INTENT(OUT) :: PTIME    ! time
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment

!
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=18)  :: YRECFM    ! Name of the article to be written
 CHARACTER(LEN=100) :: YFILE          ! Filename
INTEGER :: IVAR_ID,JRET,JDIM,INDIMS,JWRK
INTEGER, DIMENSION(1) :: IDIMIDS,IDIMLEN
INTEGER, DIMENSION(4) :: IRET
INTEGER, DIMENSION(:), POINTER :: IMASK    ! 1D mask to read only interesting
REAL:: ZTIME
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFT0_NC',0,ZHOOK_HANDLE)
!
KRESP=0
HCOMMENT = " "
!
DO JWRK=1,4
  !
  IF (JWRK == 1) THEN 
    YRECFM=TRIM(HREC)//'-YEAR'
  ELSEIF (JWRK == 2) THEN
    YRECFM = TRIM(HREC)//'-MONTH'
  ELSEIF (JWRK == 3) THEN
    YRECFM = TRIM(HREC)//'-DAY'
  ELSEIF (JWRK == 4) THEN    
    YRECFM=TRIM(HREC)//'-TIME'
  ENDIF
! 0. find filename
  !
  IF (NID_NC.NE.0) THEN
    !   
    ! 1. Find id of the variable
    !----------------------------
    JRET=NF_INQ_VARID   (NID_NC,YRECFM,IVAR_ID)
    !
    ! 2. Get variable
    !----------------------------
    IF (JWRK == 1) THEN 
      IRET(JWRK)=NF_GET_VAR_INT(NID_NC,IVAR_ID,KYEAR)
    ELSEIF (JWRK==2) THEN
      IRET(JWRK)=NF_GET_VAR_INT(NID_NC,IVAR_ID,KMONTH)
    ELSEIF (JWRK==3) THEN
      IRET(JWRK)=NF_GET_VAR_INT(NID_NC,IVAR_ID,KDAY)
    ELSEIF (JWRK==4) THEN      
      IRET(JWRK)=NF_GET_VAR_DOUBLE(NID_NC,IVAR_ID,PTIME)
    ENDIF
  ENDIF
ENDDO
!
! 3. Check for errors
!--------------------
DO JRET=1,4
  IF ((NID_NC==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    KRESP=1
  ENDIF
ENDDO
IF (KRESP /=0) CALL ERROR_READ_SURF_NC(YRECFM,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFT0_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT0_NC
!
!     #############################################################
      SUBROUTINE READ_SURFT1_NC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT0* - routine to read a NETCDF  date_time scalar
!
USE MODD_IO_SURF_NC, ONLY : NID_NC
!
USE MODI_ERROR_READ_SURF_NC
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:),    INTENT(OUT) :: KYEAR    ! year
INTEGER, DIMENSION(:),    INTENT(OUT) :: KMONTH   ! month
INTEGER, DIMENSION(:),    INTENT(OUT) :: KDAY     ! day
REAL, DIMENSION(:),       INTENT(OUT) :: PTIME    ! time
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment

!
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=18)  :: YRECFM    ! Name of the article to be written
 CHARACTER(LEN=100) :: YFILE          ! Filename
INTEGER :: IVAR_ID,JRET,JDIM,INDIMS,JWRK
INTEGER, DIMENSION(1) :: IDIMIDS,IDIMLEN
INTEGER, DIMENSION(4) :: IRET
INTEGER, DIMENSION(:), POINTER :: IMASK    ! 1D mask to read only interesting
REAL:: ZTIME
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFT1_NC',0,ZHOOK_HANDLE)
!
KRESP=0
HCOMMENT = " "
!
DO JWRK=1,4
  !
  IF (JWRK == 1) THEN 
    YRECFM=TRIM(HREC)//'-YEAR'
  ELSEIF (JWRK == 2) THEN
    YRECFM = TRIM(HREC)//'-MONTH'
  ELSEIF (JWRK == 3) THEN
    YRECFM = TRIM(HREC)//'-DAY'
  ELSEIF (JWRK == 4) THEN    
    YRECFM=TRIM(HREC)//'-TIME'
  ENDIF
! 0. find filename
  !
  IF (NID_NC.NE.0) THEN
    !   
    ! 1. Find id of the variable
    !----------------------------
    JRET=NF_INQ_VARID   (NID_NC,YRECFM,IVAR_ID)
    !
    ! 2. Get variable
    !----------------------------
    IF (JWRK == 1) THEN 
      IRET(JWRK)=NF_GET_VAR_INT(NID_NC,IVAR_ID,KYEAR)
    ELSEIF (JWRK==2) THEN
      IRET(JWRK)=NF_GET_VAR_INT(NID_NC,IVAR_ID,KMONTH)
    ELSEIF (JWRK==3) THEN
      IRET(JWRK)=NF_GET_VAR_INT(NID_NC,IVAR_ID,KDAY)
    ELSEIF (JWRK==4) THEN      
      IRET(JWRK)=NF_GET_VAR_DOUBLE(NID_NC,IVAR_ID,PTIME)
    ENDIF
  ENDIF
ENDDO
!
! 3. Check for errors
!--------------------
DO JRET=1,4
  IF ((NID_NC==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    KRESP=1
  ENDIF
ENDDO
IF (KRESP /=0) CALL ERROR_READ_SURF_NC(YRECFM,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFT1_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT1_NC
!
!     #############################################################
      SUBROUTINE READ_SURFT2_NC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT0* - routine to read a NETCDF  date_time scalar
!
USE MODD_IO_SURF_NC, ONLY : NID_NC
!
USE MODI_ERROR_READ_SURF_NC
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:,:),  INTENT(OUT) :: KYEAR    ! year
INTEGER, DIMENSION(:,:),  INTENT(OUT) :: KMONTH   ! month
INTEGER, DIMENSION(:,:),  INTENT(OUT) :: KDAY     ! day
REAL, DIMENSION(:,:),     INTENT(OUT) :: PTIME    ! time
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment

!
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=18)  :: YRECFM    ! Name of the article to be written
 CHARACTER(LEN=100) :: YFILE          ! Filename
INTEGER :: IVAR_ID,JRET,JDIM,INDIMS,JWRK
INTEGER, DIMENSION(1) :: IDIMIDS,IDIMLEN
INTEGER, DIMENSION(4) :: IRET
INTEGER, DIMENSION(:), POINTER :: IMASK    ! 1D mask to read only interesting
REAL:: ZTIME
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFT2_NC',0,ZHOOK_HANDLE)
!
KRESP=0
HCOMMENT = " "
!
DO JWRK=1,4
  !
  IF (JWRK == 1) THEN 
    YRECFM=TRIM(HREC)//'-YEAR'
  ELSEIF (JWRK == 2) THEN
    YRECFM = TRIM(HREC)//'-MONTH'
  ELSEIF (JWRK == 3) THEN
    YRECFM = TRIM(HREC)//'-DAY'
  ELSEIF (JWRK == 4) THEN    
    YRECFM=TRIM(HREC)//'-TIME'
  ENDIF
! 0. find filename
  !
  IF (NID_NC.NE.0) THEN
    !   
    ! 1. Find id of the variable
    !----------------------------
    JRET=NF_INQ_VARID   (NID_NC,YRECFM,IVAR_ID)
    !
    ! 2. Get variable
    !----------------------------
    IF (JWRK == 1) THEN 
      IRET(JWRK)=NF_GET_VAR_INT(NID_NC,IVAR_ID,KYEAR)
    ELSEIF (JWRK==2) THEN
      IRET(JWRK)=NF_GET_VAR_INT(NID_NC,IVAR_ID,KMONTH)
    ELSEIF (JWRK==3) THEN
      IRET(JWRK)=NF_GET_VAR_INT(NID_NC,IVAR_ID,KDAY)
    ELSEIF (JWRK==4) THEN      
      IRET(JWRK)=NF_GET_VAR_DOUBLE(NID_NC,IVAR_ID,PTIME)
    ENDIF
  ENDIF
ENDDO
!
! 3. Check for errors
!--------------------
DO JRET=1,4
  IF ((NID_NC==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    KRESP=1
  ENDIF
ENDDO
IF (KRESP /=0) CALL ERROR_READ_SURF_NC(YRECFM,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_NC:READ_SURFT2_NC',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT2_NC
!
