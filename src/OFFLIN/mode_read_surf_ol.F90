MODULE MODE_READ_SURF_OL
!
INTERFACE READ_SURF0_OL
        MODULE PROCEDURE READ_SURFX0_OL
        MODULE PROCEDURE READ_SURFN0_OL
        MODULE PROCEDURE READ_SURFC0_OL
        MODULE PROCEDURE READ_SURFL0_OL
END INTERFACE
INTERFACE READ_SURFN_OL
        MODULE PROCEDURE READ_SURFX1_OL
        MODULE PROCEDURE READ_SURFN1_OL
        MODULE PROCEDURE READ_SURFL1_OL
        MODULE PROCEDURE READ_SURFX2_OL
        MODULE PROCEDURE READ_SURFX3_OL
END INTERFACE
INTERFACE READ_SURFT_OL
        MODULE PROCEDURE READ_SURFT0_OL
END INTERFACE
!
CONTAINS
!
!     #############################################################
      SUBROUTINE READ_SURFX0_OL(HREC,PFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READX0* - routine to read a real scalar
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READX0 is
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
!*      0.    DECLARATIONS
!             ------------
!
USE MODI_OL_FIND_FILE_READ
USE MODI_ERROR_READ_SURF_OL
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),  INTENT(IN)  :: HREC     ! name of the article to be read
REAL,               INTENT(OUT) :: PFIELD   ! the real scalar to be read
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100):: YFILE          ! filename
INTEGER :: IVAR_ID,IFILE_ID,JRET,IVAL,ITYPE,INDIMS
INTEGER,DIMENSION(4) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFX0_OL',0,ZHOOK_HANDLE)
KRESP=0
HCOMMENT = " "
!
! 0. find filename
! -----------------
CALL OL_FIND_FILE_READ(HREC,IFILE_ID)
!
IF (IFILE_ID.NE.0) THEN
  !       
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (IFILE_ID,HREC,IVAR_ID)
  IRET(1)=NF_INQ_VARTYPE (IFILE_ID,IVAR_ID,ITYPE)
  IRET(1)=NF_INQ_VARNDIMS(IFILE_ID,IVAR_ID,INDIMS)
  !  
  ! 2. Get variable
  !----------------------------
  IRET(2)=NF_GET_VAR_DOUBLE(IFILE_ID,IVAR_ID,PFIELD)
  !  
ENDIF
!
! 3. Check for errors
!--------------------
DO JRET=1,2
  IF ((PFIELD==XUNDEF).OR.(IFILE_ID==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    PFIELD=XUNDEF
    KRESP=1
  ENDIF
ENDDO
!     
IF (KRESP /=0) CALL ERROR_READ_SURF_OL(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFX0_OL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX0_OL
!
!     #############################################################
      SUBROUTINE READ_SURFX1_OL(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX1* - routine to fill a real 1D array for the externalised surface 
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURFX1 is
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
!*      0.    DECLARATIONS
!             ------------
!    
USE MODI_OL_FIND_FILE_READ
USE MODI_ERROR_READ_SURF_OL
USE MODI_PACK_SAME_RANK
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
USE MODD_IO_SURF_OL, ONLY: LMASK,NMASK,XSTART,XCOUNT,XSTRIDE,LPARTR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),  INTENT(IN)  :: HREC     ! name of the article to be read
REAL, DIMENSION(:),INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),   INTENT(IN)  :: HDIR     ! type of field :
                                            ! 'H' : field with
                                            !       horizontal spatial dim.
                                            ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
REAL, DIMENSION(:), ALLOCATABLE :: ZTAB_1D  ! work array read in the file
CHARACTER(LEN=100):: YFILE,YOUT          ! Filename
INTEGER :: IVAR_ID,IFILE_ID,JRET,JDIM,INDIMS
INTEGER,DIMENSION(2) :: IDIMIDS,IDIMLEN
INTEGER,DIMENSION(2) :: IRET
INTEGER,DIMENSION(:),ALLOCATABLE :: ISTART,ICOUNT,ISTRIDE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!!

IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFX1_OL',0,ZHOOK_HANDLE)
KRESP=0
HCOMMENT = " "

! 0. find filename
! -----------------
CALL OL_FIND_FILE_READ(HREC,IFILE_ID)
 
IF (IFILE_ID.NE.0) THEN
    
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (IFILE_ID,HREC,IVAR_ID)
  IRET(1)=NF_INQ_VARNDIMS(IFILE_ID,IVAR_ID,INDIMS)
  IRET(1)=NF_INQ_VARDIMID(IFILE_ID,IVAR_ID,IDIMIDS(1:INDIMS))
  IDIMLEN(:) = 1.
  DO JDIM=1,INDIMS
    JRET=NF_INQ_DIMLEN(IFILE_ID,IDIMIDS(JDIM),IDIMLEN(JDIM))
  ENDDO
  ALLOCATE(ZTAB_1D(IDIMLEN(1)*IDIMLEN(2)))
    
  ! 2. Get variable
  !----------------------------
  IF  (LPARTR) THEN
    ! write partially a time-matrix. 
    ! Have to find which of the dimension is the time dimension
    ALLOCATE(ISTART(INDIMS))
    ALLOCATE(ICOUNT(INDIMS))
    ALLOCATE(ISTRIDE(INDIMS))
    DO  JDIM=1,INDIMS
      IRET=NF_INQ_DIMNAME(IFILE_ID,IDIMIDS(JDIM),YOUT)
      IF ((INDEX(YOUT,'time') > 0).OR.(INDEX(YOUT,'TIME') >0) &
        .OR.(INDEX(YOUT,'Time')>0.)) THEN  
        ISTART(JDIM)=XSTART
        ICOUNT(JDIM)=XCOUNT
        ISTRIDE(JDIM)=XSTRIDE
      ELSE
        ISTART(JDIM)=1
        ICOUNT(JDIM)=IDIMLEN(JDIM)
        ISTRIDE(JDIM)=1
      ENDIF
    ENDDO

    IRET(1)=NF_GET_VARS_DOUBLE(IFILE_ID,IVAR_ID,ISTART,ICOUNT,ISTRIDE,ZTAB_1D)
    DEALLOCATE(ISTART)
    DEALLOCATE(ICOUNT)
    DEALLOCATE(ISTRIDE)

  ELSE
    IRET(1)=NF_GET_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZTAB_1D)
  ENDIF
  !
ENDIF

! 3. Check for errors
!--------------------
DO JRET=1,1
  IF ((IFILE_ID==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    PFIELD=XUNDEF
    KRESP=1
  ELSE
    IF (MINVAL(ZTAB_1D)==XUNDEF) THEN 
      KRESP=1
      PFIELD=XUNDEF
   ENDIF
  ENDIF
ENDDO
!
IF (KRESP /=0) THEN
  CALL ERROR_READ_SURF_OL(HREC,KRESP)
ELSE
  IF (LMASK) THEN
    CALL PACK_SAME_RANK(NMASK,ZTAB_1D,PFIELD)
  ELSE 
    PFIELD(:)=ZTAB_1D(:) 
  ENDIF
END IF
!
IF (ALLOCATED(ZTAB_1D))   DEALLOCATE(ZTAB_1D)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFX1_OL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX1_OL
!
!     #############################################################
      SUBROUTINE READ_SURFX2_OL(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX2* - routine to fill a real 2D array for the externalised surface 
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURFX2 is
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
!*      0.    DECLARATIONS
!             ------------
USE MODI_OL_FIND_FILE_READ
USE MODI_ERROR_READ_SURF_OL
USE MODI_PACK_SAME_RANK
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
USE MODD_IO_SURF_OL, ONLY: LMASK,NMASK,XSTART,XCOUNT,XSTRIDE,LPARTR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),        INTENT(IN)  :: HREC     ! name of the article to be read
REAL, DIMENSION(:,:), INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),         INTENT(IN)  :: HDIR     ! type of field :
                                                  ! 'H' : field with
                                                  !       horizontal spatial dim.
                                                  ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTAB_2D  ! work array read in the file
INTEGER, DIMENSION(:), ALLOCATABLE :: ISTART,ISTRIDE,ICOUNT
CHARACTER(LEN=100):: YFILE,YOUT          ! filename
INTEGER :: IVAR_ID,IFILE_ID,JRET,JDIM,INDIMS
INTEGER,DIMENSION(3) :: IDIMIDS,IDIMLEN
INTEGER,DIMENSION(2) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFX2_OL',0,ZHOOK_HANDLE)
!
KRESP=0
HCOMMENT = " "
!
! 0. find filename
! -----------------
CALL OL_FIND_FILE_READ(HREC,IFILE_ID)
 
IF (IFILE_ID.NE.0) THEN
  !   
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (IFILE_ID,HREC,IVAR_ID)
  IRET(1)=NF_INQ_VARNDIMS(IFILE_ID,IVAR_ID,INDIMS)
  IRET(1)=NF_INQ_VARDIMID(IFILE_ID,IVAR_ID,IDIMIDS(1:INDIMS))
  IDIMLEN(:) = 1.
  DO JDIM=1,INDIMS
    JRET=NF_INQ_DIMLEN(IFILE_ID,IDIMIDS(JDIM),IDIMLEN(JDIM))
  ENDDO
  ! 
  ! 2. Get variable
  !----------------------------
  IF (LPARTR) THEN
    ! write partially a time-matrix. 
    ! Have to find which of the dimension is the time dimension
    ALLOCATE(ISTART(INDIMS))
    ALLOCATE(ICOUNT(INDIMS))
    ICOUNT(:) = 1.
    ALLOCATE(ISTRIDE(INDIMS))
    DO JDIM=1,INDIMS
      IRET=NF_INQ_DIMNAME(IFILE_ID,IDIMIDS(JDIM),YOUT)
      IF ((INDEX(YOUT,'time') > 0).OR.(INDEX(YOUT,'TIME') >0) &
        .OR.(INDEX(YOUT,'Time')>0.)) THEN  
        ISTART(JDIM)=XSTART
        ICOUNT(JDIM)=XCOUNT
        ISTRIDE(JDIM)=XSTRIDE
      ELSE
        ISTART(JDIM)=1
        ICOUNT(JDIM)=IDIMLEN(JDIM)
        ISTRIDE(JDIM)=1
      ENDIF
    ENDDO

    ALLOCATE(ZTAB_2D(PRODUCT(ICOUNT(1:INDIMS-1)),ICOUNT(INDIMS)))

    IRET(2)=NF_GET_VARS_DOUBLE(IFILE_ID,IVAR_ID,ISTART,ICOUNT,ISTRIDE,ZTAB_2D)
    DEALLOCATE(ISTART)
    DEALLOCATE(ICOUNT)
    DEALLOCATE(ISTRIDE)

  ELSE
    ALLOCATE(ZTAB_2D(PRODUCT(IDIMLEN(1:INDIMS-1)),IDIMLEN(INDIMS)))
    IRET(2)=NF_GET_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZTAB_2D)
  ENDIF

ENDIF

! 3. Check for errors
!--------------------
DO JRET=1,2
  IF ((IFILE_ID==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    PFIELD=XUNDEF
    KRESP=1
  ELSE
    IF (MINVAL(ZTAB_2D)==XUNDEF) THEN 
      KRESP=1
      PFIELD=XUNDEF
    ENDIF
  ENDIF
ENDDO
!
!
IF (KRESP /=0) THEN
  CALL ERROR_READ_SURF_OL(HREC,KRESP)
ELSE        
  IF (LMASK) THEN
    CALL PACK_SAME_RANK(NMASK,TRANSPOSE(ZTAB_2D),PFIELD)
  ELSE 
    PFIELD(:,:)=TRANSPOSE(ZTAB_2D(:,:))
  ENDIF
END IF
!  
IF (ALLOCATED(ZTAB_2D))   DEALLOCATE(ZTAB_2D) 
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFX2_OL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX2_OL
!
!     #############################################################
      SUBROUTINE READ_SURFX3_OL(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX3* - routine to fill a real 2D array for the externalised surface 
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURFX3 is
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
!*      0.    DECLARATIONS
!             ------------
USE MODI_OL_FIND_FILE_READ
USE MODI_ERROR_READ_SURF_OL
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
USE MODD_IO_SURF_OL, ONLY: LMASK,NMASK,XSTART,XCOUNT,XSTRIDE,LPARTR
USE MODI_PACK_SAME_RANK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
INCLUDE "netcdf.inc"

!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),        INTENT(IN)  :: HREC     ! name of the article to be read
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),         INTENT(IN)  :: HDIR     ! type of field :
                                                  ! 'H' : field with
                                                  !       horizontal spatial dim.
                                                  ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTAB_3D  ! work array read in the file
INTEGER, DIMENSION(:), ALLOCATABLE :: ISTART,ISTRIDE,ICOUNT
CHARACTER(LEN=100):: YFILE,YOUT          ! filename
INTEGER :: IVAR_ID,IFILE_ID,JRET,JDIM,INDIMS
INTEGER,DIMENSION(3) :: IDIMIDS,IDIMLEN
INTEGER,DIMENSION(2) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFX3_OL',0,ZHOOK_HANDLE)
KRESP=0
HCOMMENT = " "
!
! 0. find filename
! -----------------
CALL OL_FIND_FILE_READ(HREC,IFILE_ID)
! 
IF (IFILE_ID.NE.0) THEN
  !      
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (IFILE_ID,HREC,IVAR_ID)
  IRET(1)=NF_INQ_VARNDIMS(IFILE_ID,IVAR_ID,INDIMS)
  IRET(1)=NF_INQ_VARDIMID(IFILE_ID,IVAR_ID,IDIMIDS(1:INDIMS))
  DO JDIM=1,INDIMS
    JRET=NF_INQ_DIMLEN(IFILE_ID,IDIMIDS(JDIM),IDIMLEN(JDIM))
  ENDDO
  ! 
  ! 2. Get variable
  !----------------------------
  IF (LPARTR) THEN
    ! write partially a time-matrix. 
    ! Have to find which of the dimension is the time dimension
    ALLOCATE(ISTART(INDIMS))
    ALLOCATE(ICOUNT(INDIMS))
    ALLOCATE(ISTRIDE(INDIMS))
    DO  JDIM=1,INDIMS
      IRET=NF_INQ_DIMNAME(IFILE_ID,IDIMIDS(JDIM),YOUT)
      IF ((INDEX(YOUT,'time') > 0).OR.(INDEX(YOUT,'TIME') >0) &
          .OR.(INDEX(YOUT,'Time')>0.)) THEN  
        ISTART(JDIM)=XSTART
        ICOUNT(JDIM)=XCOUNT
        ISTRIDE(JDIM)=XSTRIDE
      ELSE
        ISTART(JDIM)=1
        ICOUNT(JDIM)=IDIMLEN(JDIM)
        ISTRIDE(JDIM)=1
      ENDIF
    ENDDO

    ALLOCATE(ZTAB_3D(ICOUNT(1),ICOUNT(2),ICOUNT(3)))

    IRET(2)=NF_GET_VARS_DOUBLE(IFILE_ID,IVAR_ID,ISTART,ICOUNT,ISTRIDE,ZTAB_3D)
    DEALLOCATE(ISTART)
    DEALLOCATE(ICOUNT)
    DEALLOCATE(ISTRIDE)
    !
  ELSE
    ALLOCATE(ZTAB_3D(IDIMLEN(1),IDIMLEN(2),IDIMLEN(3)))
    IRET(2)=NF_GET_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZTAB_3D)
  ENDIF
ENDIF
!
! 3. Check for errors
!--------------------
DO JRET=1,2
  IF ((IFILE_ID==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    PFIELD=XUNDEF
    KRESP=1
  ELSE
    IF (MINVAL(ZTAB_3D)==XUNDEF) THEN 
      KRESP=1
      PFIELD=XUNDEF
    ENDIF
  ENDIF
ENDDO
!
!
IF (KRESP /=0) THEN
  CALL ERROR_READ_SURF_OL(HREC,KRESP)
ELSE
  !         
  IF (LMASK) THEN
    CALL PACK_SAME_RANK(NMASK,ZTAB_3D,PFIELD)
  ELSE 
    PFIELD(:,:,:)=ZTAB_3D(:,:,:) 
  ENDIF
END IF
!  
IF (ALLOCATED(ZTAB_3D))   DEALLOCATE(ZTAB_3D)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFX3_OL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX3_OL
!
!     #############################################################
      SUBROUTINE READ_SURFN0_OL(HREC,KFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READN0 is
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
!*      0.    DECLARATIONS
!             ------------
USE MODI_OL_FIND_FILE_READ
USE MODI_ERROR_READ_SURF_OL
!
USE MODD_SURF_PAR,   ONLY: NUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,            INTENT(OUT) :: KFIELD   ! the integer scalar to be read
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
!
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100):: YFILE          ! filename
INTEGER :: IVAR_ID,IFILE_ID,JRET,JDIM,INDIMS
INTEGER,DIMENSION(4) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFN0_OL',0,ZHOOK_HANDLE)
KRESP=0
HCOMMENT = " "

! 0. find filename
! -----------------
CALL OL_FIND_FILE_READ(HREC,IFILE_ID)
 
IF (IFILE_ID.NE.0) THEN
  !        
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (IFILE_ID,HREC,IVAR_ID)
  !  
  ! 2. Get variable
  !----------------------------
  IRET(2)=NF_GET_VAR_INT(IFILE_ID,IVAR_ID,KFIELD)
  !  
ENDIF
!
! 3. Check for errors
!--------------------
DO JRET=1,2
  IF ((KFIELD==NUNDEF).OR.(IFILE_ID==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    KFIELD=NUNDEF
    KRESP=1
  ENDIF
ENDDO
!
IF (KRESP /=0)  CALL ERROR_READ_SURF_OL(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFN0_OL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFN0_OL
!     #############################################################
      SUBROUTINE READ_SURFN1_OL(HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READN0 is
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
!*      0.    DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),      INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:), INTENT(OUT) :: KFIELD   ! the integer scalar to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
                                                ! 'H' : field with
                                                !       horizontal spatial dim.
                                                ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(KFIELD)) :: ZFIELD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFN1_OL',0,ZHOOK_HANDLE)
CALL READ_SURFX1_OL(HREC,ZFIELD,KRESP,HCOMMENT,HDIR)
KFIELD = NINT(ZFIELD)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFN1_OL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFN1_OL
!
!     #############################################################
      SUBROUTINE READ_SURFC0_OL(HREC,HFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READC0* - routine to read a STRING
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READC0 is
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
!*      0.    DECLARATIONS
!             ------------
USE MODI_OL_FIND_FILE_READ
USE MODI_ERROR_READ_SURF_OL
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),   INTENT(IN)  :: HREC     ! name of the article to be read
CHARACTER(LEN=40),   INTENT(OUT) :: HFIELD   ! the integer scalar to be read
INTEGER,             INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),  INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100):: YFILE          ! filename
CHARACTER(LEN=100):: YFIELD   
INTEGER :: IVAR_ID,IFILE_ID,JRET,JDIM,INDIMS
INTEGER,DIMENSION(4) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFC0_OL',0,ZHOOK_HANDLE)
KRESP=0
HCOMMENT = " "
!
! 0. find filename
! -----------------
CALL OL_FIND_FILE_READ(HREC,IFILE_ID)
!
IF (IFILE_ID.NE.0) THEN
  !       
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (IFILE_ID,HREC,IVAR_ID)
  !  
  ! 2. Get variable
  !----------------------------
  IRET(2)=NF_GET_VAR_TEXT(IFILE_ID,IVAR_ID,YFIELD)
  HFIELD=YFIELD(:LEN_TRIM(YFIELD))
  !  
ENDIF

! 3. Check for errors
!--------------------
DO JRET=1,2
  IF ((IFILE_ID==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    KRESP=1
  ENDIF
ENDDO  
!
IF (KRESP /=0) CALL ERROR_READ_SURF_OL(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFC0_OL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFC0_OL
!
!     #############################################################
      SUBROUTINE READ_SURFL1_OL(HREC,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READL1* - routine to read a logical array
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READL1 is
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
!*      0.    DECLARATIONS
!             ------------
!    
USE MODI_OL_FIND_FILE_READ
USE MODI_ERROR_READ_SURF_OL
USE MODI_PACK_SAME_RANK
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
INCLUDE "netcdf.inc"

!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),        INTENT(IN)  :: HREC     ! name of the article to be read
LOGICAL, DIMENSION(:),   INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),         INTENT(IN)  :: HDIR     ! type of field :
                                                  ! 'H' : field with
                                                  !       horizontal spatial dim.
                                                  ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: YTAB_1D  ! work array read in the file
CHARACTER(LEN=100):: YFILE          ! Filename
INTEGER, DIMENSION(:),    POINTER     :: IMASK    ! 1D mask to read only interesting
INTEGER :: IVAR_ID,IFILE_ID,JRET,JDIM,INDIMS
INTEGER,DIMENSION(1) :: IDIMIDS,IDIMLEN
INTEGER,DIMENSION(2) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFL1_OL',0,ZHOOK_HANDLE)
KRESP=0
HCOMMENT = " "
!
! 0. find filename
! -----------------
CALL OL_FIND_FILE_READ(HREC,IFILE_ID)
! 
IF (IFILE_ID.NE.0) THEN
  !   
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (IFILE_ID,HREC,IVAR_ID)
  IRET(1)=NF_INQ_VARNDIMS(IFILE_ID,IVAR_ID,INDIMS)
  IRET(1)=NF_INQ_VARDIMID(IFILE_ID,IVAR_ID,IDIMIDS)
  DO JDIM=1,INDIMS
    JRET=NF_INQ_DIMLEN(IFILE_ID,IDIMIDS(JDIM),IDIMLEN(JDIM))
  ENDDO
  ALLOCATE(YTAB_1D(IDIMLEN(1)))
  !  
  ! 2. Get variable
  !----------------------------
  IRET(1)=NF_GET_VAR_TEXT(IFILE_ID,IVAR_ID,YTAB_1D)
  !
  DO JRET=1,IDIMLEN(1)
    IF (YTAB_1D(JRET) =='T') OFIELD(JRET)=.TRUE.
    IF (YTAB_1D(JRET) =='F') OFIELD(JRET)=.FALSE.
  ENDDO
  !
ENDIF
!
! 3. Check for errors
!--------------------
DO JRET=1,1
  IF ((IFILE_ID==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    KRESP=1
  ENDIF
ENDDO
!
IF (KRESP /=0) THEN
  CALL ERROR_READ_SURF_OL(HREC,KRESP)
ELSE
  DEALLOCATE(YTAB_1D)
END IF
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFL1_OL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFL1_OL
!
!
!     #############################################################
      SUBROUTINE READ_SURFL0_OL(HREC,OFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READL0* - routine to read a logical
!!
!!    PURPOSE
!!    -------
!
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
!*      0.    DECLARATIONS
!             ------------
!    
USE MODI_OL_FIND_FILE_READ
USE MODI_ERROR_READ_SURF_OL
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
INCLUDE "netcdf.inc"

!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),        INTENT(IN)  :: HREC     ! name of the article to be read
LOGICAL,                  INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=1)  :: YFIELD   ! work array read in the file
CHARACTER(LEN=100):: YFILE    ! Filename
INTEGER :: IVAR_ID,IFILE_ID, JRET
INTEGER,DIMENSION(2) :: IRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFL0_OL',0,ZHOOK_HANDLE)
KRESP=0
HCOMMENT = " "
!
! 0. find filename
! -----------------
CALL OL_FIND_FILE_READ(HREC,IFILE_ID)
!
IF (IFILE_ID.NE.0) THEN
  !       
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (IFILE_ID,HREC,IVAR_ID)
  !  
  ! 2. Get variable
  !----------------------------
  IRET(2)=NF_GET_VAR_TEXT(IFILE_ID,IVAR_ID,YFIELD)
  !  
  IF (YFIELD =='T') OFIELD=.TRUE.
  IF (YFIELD =='F') OFIELD=.FALSE.
  !
ENDIF
!
! 3. Check for errors
!--------------------
IF ((IFILE_ID==0).OR.IRET(1).NE.NF_NOERR) THEN 
  KRESP=1
ENDIF
!
IF (KRESP /=0)  CALL ERROR_READ_SURF_OL(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFL0_OL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFL0_OL
!
!     #############################################################
      SUBROUTINE READ_SURFT0_OL(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT0* - routine to read a NETCDF  date_time scalar
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READT0 is
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
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     18/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
USE MODI_OL_FIND_FILE_READ
USE MODI_ERROR_READ_SURF_OL
!
USE MODD_SURF_PAR,   ONLY: XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),        INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                  INTENT(OUT) :: KYEAR    ! year
INTEGER,                  INTENT(OUT) :: KMONTH   ! month
INTEGER,                  INTENT(OUT) :: KDAY     ! day
REAL,                     INTENT(OUT) :: PTIME    ! time
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment

!
!*      0.2   Declarations of local variables
!
INTEGER, DIMENSION(3) :: ITDATE  ! work array read in the file
REAL:: ZTIME
CHARACTER(LEN=100):: YFILE          ! Filename
INTEGER, DIMENSION(:),    POINTER     :: IMASK    ! 1D mask to read only interesting
INTEGER :: IVAR_ID,IFILE_ID,JRET,JDIM,INDIMS,JWRK
INTEGER,DIMENSION(1) :: IDIMIDS,IDIMLEN
INTEGER,DIMENSION(4) :: IRET
CHARACTER(LEN=16)              :: YRECFM    ! Name of the article to be written
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!

IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFT0_OL',0,ZHOOK_HANDLE)
KRESP=0
HCOMMENT = " "

DO JWRK=1,2
  IF (JWRK == 1) THEN 
    YRECFM=TRIM(HREC)//'-TDATE'
  ELSE
    YRECFM=TRIM(HREC)//'-TIME'
  ENDIF
! 0. find filename
! -----------------
  CALL OL_FIND_FILE_READ(YRECFM,IFILE_ID)
  !
  IF (IFILE_ID.NE.0) THEN
    !   
    ! 1. Find id of the variable
    !----------------------------
    JRET=NF_INQ_VARID   (IFILE_ID,YRECFM,IVAR_ID)
    !
    ! 2. Get variable
    !----------------------------
    IF (JWRK == 1) THEN 
      IRET(JWRK)=NF_GET_VAR_INT(IFILE_ID,IVAR_ID,ITDATE)
      KYEAR  = ITDATE(1)
      KMONTH = ITDATE(2)
      KDAY   = ITDATE(3)
    ELSE
      IRET(JWRK)=NF_GET_VAR_DOUBLE(IFILE_ID,IVAR_ID,PTIME)
    ENDIF
  ENDIF
ENDDO
!
! 3. Check for errors
!--------------------
DO JRET=1,2
  IF ((IFILE_ID==0).OR.IRET(JRET).NE.NF_NOERR) THEN 
    KRESP=1
  ENDIF
ENDDO
IF (KRESP /=0) CALL ERROR_READ_SURF_OL(YRECFM,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_OL:READ_SURFT0_OL',1,ZHOOK_HANDLE)
END SUBROUTINE READ_SURFT0_OL
!
END MODULE MODE_READ_SURF_OL
