MODULE MODE_READ_SURF_FA
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURF_FA is
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
INTERFACE READ_SURF0_FA
        MODULE PROCEDURE READ_SURFX0_FA
        MODULE PROCEDURE READ_SURFN0_FA
        MODULE PROCEDURE READ_SURFL0_FA
        MODULE PROCEDURE READ_SURFC0_FA
END INTERFACE
INTERFACE READ_SURFN_FA
        MODULE PROCEDURE READ_SURFX1_FA
        MODULE PROCEDURE READ_SURFN1_FA
        MODULE PROCEDURE READ_SURFL1_FA
        MODULE PROCEDURE READ_SURFX2_FA
END INTERFACE
INTERFACE READ_SURFT_FA
        MODULE PROCEDURE READ_SURFT0_FA
        MODULE PROCEDURE READ_SURFT2_FA
END INTERFACE
!
CONTAINS
!
!     #############################################################
      SUBROUTINE READ_SURFX0_FA(HREC,PFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READX0* - routine to read a real scalar
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, CMASK
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF_n
USE MODI_ERROR_READ_SURF_FA
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
 CHARACTER(LEN=50) :: YCOMMENT
 CHARACTER(LEN=6)  :: YMASK
LOGICAL           :: GFOUND
LOGICAL           :: GKNOWN
 CHARACTER(LEN=18) :: YNAME ! Field Name
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX0_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
YMASK=CMASK
 CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
YNAME=TRIM(YMASK)//TRIM(HREC)
 CALL FALIT_R(KRESP,NUNIT_FA,YNAME,PFIELD)
IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
!
YCOMMENT = TRIM(YNAME)
HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX0_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX0_FA
!
!     #############################################################
      SUBROUTINE READ_SURFX1_FA(HREC,KL,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX1* - routine to fill a real 1D array for the externalised surface 
!
USE MODD_SURFEX_OMP, ONLY : NBLOCK
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ, &
                            WLOG_MPI
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, NMASK, NFULL, NFULL_EXT, &
                                   NDGL, NDLON, NDGUX, NDLUX  
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF_n
USE MODI_ERROR_READ_SURF_FA
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
 CHARACTER(LEN=*),   INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,             INTENT(IN)  :: KL       ! number of points
REAL, DIMENSION(:), INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,             INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),  INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),    INTENT(IN)  :: HDIR     ! type of field :
                                             ! 'H' : field with
                                             !       horizontal spatial dim.
                                             ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=50)          :: YCOMMENT
 CHARACTER(LEN=3)           :: YPREF
 CHARACTER(LEN=13)          :: YSUFF
LOGICAL                    :: GFOUND
!
INTEGER ::  I, J, INFOMPI
#ifndef NOMPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
!
REAL, DIMENSION(:), ALLOCATABLE :: ZWORK   ! work array read in the file
REAL, DIMENSION(:), ALLOCATABLE :: ZWORK2
DOUBLE PRECISION   :: XTIME0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX1_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
#ifndef NOMPI
XTIME0 = MPI_WTIME()
#endif
!
!
IF (NRANK==NPIO) THEN
  !
  ALLOCATE(ZWORK(NFULL))
  ALLOCATE(ZWORK2(NFULL_EXT))
  !
!$OMP SINGLE
  !
  YPREF=HREC(1:3)
  YSUFF=HREC(4:16) 
  !
  IF (YPREF=='CLS' .OR. YPREF=='SUR' .OR. YPREF=='PRO' .OR. YPREF=='ATM') THEN
    CALL FACILE(KRESP,NUNIT_FA,HREC(1:4),0,HREC(5:16),ZWORK2,.FALSE.)
    IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
    DO J=1,NDGUX
      DO I=1,NDLUX
        ZWORK((J-1)*NDLUX + I) = ZWORK2((J-1)*NDLON + I)
      ENDDO
    ENDDO
    YCOMMENT = TRIM(HREC)
  ELSE
    CALL FACILE(KRESP,NUNIT_FA,'S1D_',0,HREC,ZWORK,.FALSE.)
    IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)  
    YCOMMENT = 'S1D_'//TRIM(HREC)
  ENDIF
  !
  HCOMMENT = YCOMMENT
  !
!$OMP END SINGLE COPYPRIVATE(ZWORK,HCOMMENT,KRESP)
  !
ELSE
  ALLOCATE(ZWORK(0))
  ALLOCATE(ZWORK2(0))
ENDIF
!
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
    PFIELD(:) = ZWORK(1:KL)
#ifndef NOMPI    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif    
  ENDIF
ELSEIF (HDIR=='-') THEN ! distribution of the total field on other tasks
!$OMP SINGLE
    IF (NRANK==NPIO) PFIELD(:) = ZWORK(1:KL)
#ifndef NOMPI          
    IF (NPROC>1) THEN
      XTIME0 = MPI_WTIME()
      CALL MPI_BCAST(PFIELD,SIZE(PFIELD)*KIND(PFIELD)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)
      XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
    ENDIF
#endif   
!$OMP END SINGLE COPYPRIVATE(PFIELD)
ELSE
  CALL READ_AND_SEND_MPI(ZWORK,PFIELD,NMASK)
ENDIF
!
DEALLOCATE(ZWORK)
DEALLOCATE(ZWORK2)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX1_FA',1,ZHOOK_HANDLE)
! 
END SUBROUTINE READ_SURFX1_FA
!
!     #############################################################
      SUBROUTINE READ_SURFX2_FA(HREC,KL1,KL2,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX2* - routine to fill a real 2D array for the externalised surface 
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ, &
                            WLOG_MPI
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, NMASK, NFULL
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF_n
USE MODI_ERROR_READ_SURF_FA
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
 CHARACTER(LEN=*),        INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                  INTENT(IN)  :: KL1      ! number of points
INTEGER,                  INTENT(IN)  :: KL2      ! 2nd dimension
REAL, DIMENSION(:,:), INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),         INTENT(IN)  :: HDIR     ! type of field :
                                                  ! 'H' : field with
                                                  !       horizontal spatial dim.
                                                  ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
! 
 CHARACTER(LEN=50)                     :: YCOMMENT
 CHARACTER(LEN=4)                      :: YSUFFIX
 CHARACTER(LEN=2)                      :: YPATCH
LOGICAL                               :: GFOUND
!
INTEGER :: JL, I, INFOMPI ! loop counter
#ifndef NOMPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK   ! work array read in the file
DOUBLE PRECISION   :: XTIME0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX2_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
#ifndef NOMPI
XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  !
  ALLOCATE(ZWORK(NFULL,KL2))
  !
!$OMP SINGLE
  !  
  DO JL=1,KL2
    WRITE(YPATCH,'(I2.2)')JL
    YSUFFIX='S'//YPATCH//'_'
    CALL FACILE(KRESP,NUNIT_FA,YSUFFIX,JL,HREC,ZWORK(:,JL),.FALSE.)
     IF (KRESP/=0) THEN
       YCOMMENT = YSUFFIX//TRIM(HREC)
       CALL ERROR_READ_SURF_FA(YCOMMENT,KRESP)
     ENDIF  
  END DO
  !
  YCOMMENT = 'PATCH_'//TRIM(HREC)
  HCOMMENT = YCOMMENT  
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
    PFIELD(:,:) = ZWORK(1:KL1,1:KL2)
#ifndef NOMPI    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif    
  ENDIF
ELSEIF (HDIR=='-') THEN ! distribution of the total field on other tasks
!$OMP SINGLE
  IF (NRANK==NPIO) PFIELD(:,:) = ZWORK(1:KL1,1:KL2)
#ifndef NOMPI
  IF (NPROC>1) THEN
    XTIME0 = MPI_WTIME()
    CALL MPI_BCAST(PFIELD(:,:),SIZE(PFIELD)*KIND(PFIELD)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
  ENDIF
#endif    
!$OMP END SINGLE COPYPRIVATE(PFIELD)
ELSE
  CALL READ_AND_SEND_MPI(ZWORK,PFIELD,NMASK)
ENDIF
!
DEALLOCATE(ZWORK)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX2_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX2_FA
!
!     #############################################################
      SUBROUTINE READ_SURFN0_FA(HREC,KFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, NMASK, CMASK
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF_n
USE MODI_ERROR_READ_SURF_FA
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
 CHARACTER(LEN=50) :: YCOMMENT
 CHARACTER(LEN=6)  :: YMASK
LOGICAL           :: GFOUND
LOGICAL           :: GKNOWN
 CHARACTER(LEN=18) :: YNAME ! Field Name
REAL(KIND=JPRB)  :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN0_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
YMASK=CMASK
 CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
YNAME=TRIM(YMASK)//TRIM(HREC)
 CALL FALIT_I(KRESP,NUNIT_FA,YNAME,KFIELD)
IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
!
YCOMMENT = YNAME
HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN0_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFN0_FA
!
!     #############################################################
      SUBROUTINE READ_SURFN1_FA(HREC,KL,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ, & 
                            WLOG_MPI
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, NMASK, NFULL, CMASK
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF_n
USE MODI_ERROR_READ_SURF_FA
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
INTEGER,                INTENT(IN)  :: KL       ! number of points
INTEGER, DIMENSION(:), INTENT(OUT) :: KFIELD   ! the integer to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
                                                ! 'H' : field with
                                                !       horizontal spatial dim.
                                                ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=50) :: YCOMMENT
 CHARACTER(LEN=6)  :: YMASK
 CHARACTER(LEN=18) :: YNAME ! Field Name
LOGICAL           :: GFOUND
LOGICAL           :: GKNOWN
!
INTEGER ::  I, INFOMPI
#ifndef NOMPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
!
INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK  ! work array read in the file
DOUBLE PRECISION   :: XTIME0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN1_FA',0,ZHOOK_HANDLE)
!
KRESP = 0
!
#ifndef NOMPI
XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  !
  ALLOCATE(IWORK(NFULL))
  !
!$OMP SINGLE        
  !
  YMASK = CMASK
  CALL IO_BUFF_n(HREC,'R',GKNOWN)
  IF (GKNOWN) YMASK = 'FULL  '
  !
  YNAME = TRIM(YMASK)//TRIM(HREC)
  IF (HDIR=="-") THEN
    CALL FALIT_I_D(KRESP,NUNIT_FA,YNAME,KL,IWORK(1:KL))
    IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
  ELSE
    CALL FALIT_I_D(KRESP,NUNIT_FA,YNAME,NFULL,IWORK)
    IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
  ENDIF
  !
  YCOMMENT = YNAME
  HCOMMENT = YCOMMENT
  !
!$OMP END SINGLE COPYPRIVATE(IWORK,HCOMMENT,KRESP)  
  !
ELSE
  ALLOCATE(IWORK(0))
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
    KFIELD(:) = IWORK(1:KL)
#ifndef NOMPI    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif    
  ENDIF
ELSEIF (HDIR=='-') THEN ! distribution of the total field on other tasks
!$OMP SINGLE     
  IF (NRANK==NPIO) KFIELD(:) = IWORK(1:KL)
#ifndef NOMPI
  IF (NPROC>1) THEN
    XTIME0 = MPI_WTIME()
    CALL MPI_BCAST(KFIELD,SIZE(KFIELD)*KIND(KFIELD)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
  ENDIF
#endif    
!$OMP END SINGLE COPYPRIVATE(KFIELD)
ELSE
  CALL READ_AND_SEND_MPI(IWORK,KFIELD,NMASK)
ENDIF
!
DEALLOCATE(IWORK)
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN1_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFN1_FA
!
!     #############################################################
      SUBROUTINE READ_SURFC0_FA(HREC,HFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READC0* - routine to read a character
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, CMASK
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF_n
USE MODI_ERROR_READ_SURF_FA
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
 CHARACTER(LEN=50) :: YCOMMENT
 CHARACTER(LEN=6)  :: YMASK
 CHARACTER(LEN=18) :: YNAME ! Field Name
 CHARACTER,DIMENSION(40) :: YFIELD
LOGICAL           :: GFOUND
LOGICAL           :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFC0_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
YMASK=CMASK
 CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
YNAME=TRIM(YMASK)//TRIM(HREC)
 CALL FALIT_C(KRESP,NUNIT_FA,YNAME,40,YFIELD)
IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
WRITE(HFIELD,'(40A1)') YFIELD(:)
!
YCOMMENT = YNAME
HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFC0_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFC0_FA
!
!     #############################################################
      SUBROUTINE READ_SURFL0_FA(HREC,OFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READL0* - routine to read a logical
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, CMASK
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF_n
USE MODI_ERROR_READ_SURF_FA
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
 CHARACTER(LEN=50) :: YCOMMENT
 CHARACTER(LEN=6)  :: YMASK
 CHARACTER(LEN=18) :: YNAME ! Field Name
LOGICAL           :: GFOUND
LOGICAL           :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFL0_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
YMASK=CMASK
 CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
YNAME=TRIM(YMASK)//TRIM(HREC)
 CALL FALIT_L(KRESP,NUNIT_FA,YNAME,OFIELD)
IF (KRESP/=0)CALL ERROR_READ_SURF_FA(HREC,KRESP)
!
YCOMMENT = YNAME
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFL0_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFL0_FA
!
!     #############################################################
      SUBROUTINE READ_SURFL1_FA(HREC,KL,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READL1* - routine to read a logical array
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPROC, NCOMM, NPIO, XTIME_NPIO_READ, XTIME_COMM_READ, &
                            WLOG_MPI
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, CMASK
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF_n
USE MODI_ERROR_READ_SURF_FA
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
INTEGER,                INTENT(IN)  :: KL       ! number of points
LOGICAL, DIMENSION(:), INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
 CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
                                                ! 'H' : field with
                                                !       horizontal spatial dim.
                                                ! '-' : no horizontal dim.
!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=50) :: YCOMMENT
 CHARACTER(LEN=6)  :: YMASK
 CHARACTER(LEN=18) :: YNAME ! Field Name
!LOGICAL, DIMENSION(KL) :: GCOVER
LOGICAL           :: GFOUND
LOGICAL           :: GKNOWN
INTEGER           :: INFOMPI
DOUBLE PRECISION  :: XTIME0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFL1_FA',0,ZHOOK_HANDLE)
!
#ifndef NOMPI
XTIME0 = MPI_WTIME()
#endif
!
IF (NRANK==NPIO) THEN
  !
!$OMP SINGLE
  !   
  KRESP=0
  !
  YMASK=CMASK
  CALL IO_BUFF_n(HREC,'R',GKNOWN)
  IF (GKNOWN) YMASK='FULL  '
  !
  YNAME=TRIM(YMASK)//TRIM(HREC)
  CALL FALIT_L_D(KRESP,NUNIT_FA,YNAME,KL,OFIELD)
  IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
  !
  YCOMMENT = YNAME
  HCOMMENT = YCOMMENT
  !
!$OMP END SINGLE COPYPRIVATE(OFIELD,HCOMMENT,KRESP)  
  !
  !OFIELD(:) = GCOVER(:)
  !
ENDIF
!
#ifndef NOMPI
XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
!
#ifndef NOMPI
IF (NPROC>1 .AND. HDIR/='A') THEN
!$OMP SINGLE         
  XTIME0 = MPI_WTIME()
  CALL MPI_BCAST(OFIELD,SIZE(OFIELD),MPI_LOGICAL,NPIO,NCOMM,INFOMPI)
  XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
!$OMP END SINGLE COPYPRIVATE(OFIELD)  
ENDIF
#endif
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFL1_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFL1_FA
!
!     #############################################################
      SUBROUTINE READ_SURFT0_FA(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT0* - routine to read a date
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, CMASK
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF_n
USE MODI_ERROR_READ_SURF_FA
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
 CHARACTER(LEN=50) :: YCOMMENT
 CHARACTER(LEN=6)  :: YMASK
 CHARACTER(LEN=18) :: YNAME ! Field Name
LOGICAL           :: GFOUND
LOGICAL          :: GKNOWN
INTEGER, DIMENSION(3) :: ITDATE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT0_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
YMASK=CMASK
 CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
YNAME=TRIM(YMASK)//TRIM(HREC)//'%TDATE'
 CALL FALIT_I_D(KRESP,NUNIT_FA,YNAME,3,ITDATE)
IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
!
KYEAR  = ITDATE(1)
KMONTH = ITDATE(2)
KDAY   = ITDATE(3)
!
YNAME=TRIM(YMASK)//TRIM(HREC)//'%TIME'
 CALL FALIT_R(KRESP,NUNIT_FA,YNAME,PTIME)
IF (KRESP/=0) CALL ERROR_READ_SURF_FA(HREC,KRESP)
!
YCOMMENT = TRIM(HREC)
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT0_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT0_FA
!
!     #############################################################
      SUBROUTINE READ_SURFT2_FA(HREC,KL1,KL2,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT2* - routine to read a date
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, CMASK
!
USE MODE_FASURFEX
!
USE MODI_IO_BUFF_n
USE MODI_ABOR1_SFX
USE MODI_ERROR_READ_SURF_FA
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=*),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER                                  :: KL1, KL2
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KYEAR    ! year
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KMONTH   ! month
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KDAY     ! day
REAL,    DIMENSION(:,:), INTENT(OUT) :: PTIME    ! year
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment

!*      0.2   Declarations of local variables
!
 CHARACTER(LEN=50) :: YCOMMENT
 CHARACTER(LEN=6)  :: YMASK
 CHARACTER(LEN=18) :: YNAME ! Field Name
LOGICAL           :: GFOUND
LOGICAL           :: GKNOWN
INTEGER, DIMENSION(3,SIZE(KYEAR,1),SIZE(KYEAR,2)) :: ITDATE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT2_FA',0,ZHOOK_HANDLE)
!
KRESP=0
!
YMASK=CMASK
 CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
YNAME=TRIM(CMASK)//TRIM(HREC)
WRITE(NLUOUT,*) ' READ_SURFT2_FA : time in 2 dimensions not yet implemented : YNAME=',YNAME
 CALL ABOR1_SFX('MODE_READ_SURF_FA:READ_SURFT2_FA: time in 2 dimensions not yet implemented')
!
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT2_FA',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFT2_FA
!
END MODULE MODE_READ_SURF_FA
