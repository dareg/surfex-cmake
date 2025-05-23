!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.

MODULE MODE_READ_SURF_LAYERS
!
INTERFACE WRITE_READ_LAYERS
  MODULE PROCEDURE READ_SURF_LAYERS
END INTERFACE
!
CONTAINS
!
!     #############################################################
      SUBROUTINE READ_SURF_LAYERS (HPROGRAM,HREC,ODIM,PFIELD,KRESP,HCOMMENT,HDIR,KPATCH)
!     #############################################################
!
!
!
USE MODD_SURFEX_MPI, ONLY : NPROC, NPIO, NRANK, IDX_R, NREQ, NCOMM, NSIZE, NINDEX, LSFX_MPI
USE MODD_SURF_PAR,  ONLY : XUNDEF
!
#ifdef SFX_LFI
USE MODD_IO_SURF_LFI, ONLY : NMASK_lfi=>NMASK, NFULL_lfi=>NFULL
#endif
#ifdef SFX_NC
USE MODD_IO_SURF_NC, ONLY : NMASK_nc=>NMASK, NFULL_nc=>NFULL
#endif
#ifdef SFX_ASC
USE MODD_IO_SURF_ASC, ONLY : NMASK_asc=>NMASK, NFULL_asc=>NFULL
#endif
#ifdef SFX_FA
USE MODD_IO_SURF_FA, ONLY : NMASK_fa=>NMASK, NFULL_fa=>NFULL
#endif
#ifdef SFX_MNH
USE MODI_READ_SURF
#endif
!
USE MODI_ABOR1_SFX
USE MODI_PACK_SAME_RANK
USE MODI_READ_AND_SEND_MPI
USE MODI_MAKE_CHOICE_ARRAY
USE MODD_SURFEX_HOST
!
USE YOMHOOK ,ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
!*      0.1   Declarations of arguments
!
!
!
 CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM    ! calling program
 CHARACTER(LEN=*), INTENT(IN) :: HREC        ! name of the article to be read
LOGICAL, INTENT(IN) :: ODIM
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PFIELD ! array containing the data field
INTEGER, INTENT(OUT) :: KRESP               ! KRESP  : return-code if a problem appears
 CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: HCOMMENT ! name of the article to be read
 CHARACTER(LEN=1), OPTIONAL, INTENT(IN)  :: HDIR     ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
INTEGER, OPTIONAL, INTENT(IN) :: KPATCH
!
!*      0.2   Declarations of local variables
!
#ifdef SFX_MPI
INTEGER, DIMENSION(MPI_STATUS_SIZE,NPROC-1) :: ISTATUS
#endif
!INTEGER, DIMENSION(NPROC) :: ITREQ
REAL, DIMENSION(:,:),ALLOCATABLE :: ZWORKR
REAL, DIMENSION(:,:),ALLOCATABLE :: ZFIELD
INTEGER, DIMENSION(:), POINTER :: IMASKF
 CHARACTER(LEN=100) :: YCOMMENT
 CHARACTER(LEN=16)  :: YREC
 CHARACTER(LEN=1)   :: YDIR
 CHARACTER(LEN=4)  :: YLVL
INTEGER :: IFLAG, IPATCH, INPATCH
INTEGER :: IPIO_SAVE, IPAS, JP, IDEB, IFIN, JJ, JL
INTEGER :: JLAYER, JPROC, IPROC, IRESP
INTEGER :: IL1, IL2, IL3, IDX_SAVE, IDX, IVAL
INTEGER :: INFOMPI, IREQ, JPROC2, IFULL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE_OMP
!
IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS',0,ZHOOK_HANDLE)
!
IDX_SAVE = IDX_R
YREC = HREC
YCOMMENT="empty"
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
INPATCH = SIZE(PFIELD,3)
IPATCH = -1
IF (PRESENT(KPATCH).AND.INPATCH==1) IPATCH = KPATCH
!
IL1 = SIZE(PFIELD,1)
IL2 = SIZE(PFIELD,2)
IL3 = SIZE(PFIELD,3)
!
!IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_1',1,ZHOOK_HANDLE)
!
IF (HPROGRAM=='MESONH') THEN
#ifdef SFX_MNH
  DO JL=1,IL2
    WRITE(YLVL,'(I4)') JL
    YREC=TRIM(HREC)//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL MAKE_CHOICE_ARRAY(HPROGRAM, INPATCH, ODIM, YREC, PFIELD(:,JL,:), HDIR=YDIR, KPATCH=IPATCH)
  END DO
#endif
ELSE
  !
!  IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_2',0,ZHOOK_HANDLE)
  !
  !the mask to call read_and_send_mpi depends on the I/O type
  IF (HPROGRAM=='LFI   ') THEN
#ifdef SFX_LFI
    IFULL = NFULL_lfi
    ALLOCATE(ZFIELD(NFULL_lfi,IL3))
    IMASKF=>NMASK_lfi
#endif
  ELSEIF (HPROGRAM=='ASCII ') THEN
#ifdef SFX_ASC
    IFULL = NFULL_asc
    ALLOCATE(ZFIELD(NFULL_asc,IL3))
    IMASKF=>NMASK_asc
#endif
  ELSEIF (HPROGRAM=='FA     ') THEN
#ifdef SFX_FA
    IFULL = NFULL_fa
    ALLOCATE(ZFIELD(NFULL_fa,IL3))
    IMASKF=>NMASK_fa
#endif
  ELSEIF (HPROGRAM=='NC     ') THEN
#ifdef SFX_NC
    IFULL = NFULL_nc
    ALLOCATE(ZFIELD(NFULL_nc,IL3))
    IMASKF=>NMASK_nc
#endif
  ELSE
    ALLOCATE(ZFIELD(0,0))
  ENDIF
  !
  !if we want to get covers for the current task or for the whole domain
  IF (YDIR=='H') THEN
     !second dimension because the reading of covers is parallelized with MPI
    ALLOCATE(ZWORKR(NSIZE,IL3))
  ELSEIF (NRANK==NPIO) THEN
    ALLOCATE(ZWORKR(IFULL,IL3))
  ELSE 
    ALLOCATE(ZWORKR(0,0))
  ENDIF
  ZWORKR(:,:) = 0.
  !
  IF (NPROC>1 .AND. YDIR=='H') THEN
    IFLAG = 0
    !for the parallelization of reading, NINDEX must be known by all tasks
    IF (NRANK/=NPIO) THEN
      IF (ALLOCATED(NINDEX)) THEN
        IF (SIZE(NINDEX)==IFULL) IFLAG=1
        DEALLOCATE(NINDEX)
      ENDIF
      ALLOCATE(NINDEX(IFULL))
    ENDIF
#ifdef SFX_MPI
    IF (LSFX_MPI) THEN
      CALL MPI_BCAST(NINDEX,SIZE(NINDEX)*KIND(NINDEX)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    ENDIF
#endif
  ENDIF
  !
  IPIO_SAVE = NPIO
  !number of covers read by each task
  IPAS = CEILING(IL2*1./NPROC)
  !
  PFIELD(:,:,:) = 0.
  !
  !first cover number read by the current task
  IDEB = IPAS*NRANK
  !
!  IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_2',1,ZHOOK_HANDLE)
  !
  DO JP = 1,IPAS
    !
!    IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_3',0,ZHOOK_HANDLE)
    !
    !index of the cover read by this task at this loop index
    JLAYER = IDEB + JP
    !
    IF (JLAYER<=IL2) THEN
      !
      WRITE(YLVL,'(I4)') JLAYER
      YREC = TRIM(HREC)//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
      YCOMMENT='X_Y_'//YREC
      !
      !
      IF (HPROGRAM=='AROME ') THEN
         CALL YRSURFEX_HOST%READ_SURF_LAYER (HPROGRAM, INPATCH, ODIM, YREC, PFIELD(:,JLAYER,:), HDIR=YDIR, KPATCH=IPATCH)
      ELSE
        !
        !reads one cover by task
        !
        !number of the I/O task for this read 
        NPIO = NRANK
        !
        !reading of the whol cover (HDIR='A')
        CALL MAKE_CHOICE_ARRAY(HPROGRAM, INPATCH, ODIM, YREC, ZFIELD, HDIR='A', KPATCH=IPATCH)
        !
        !NPIO rebecomes the I/O task 
        NPIO = IPIO_SAVE
        !
        IDX = IDX_SAVE + JP
        IF (YDIR=='H') THEN
          !
          !send covers to other tasks
          CALL READ_AND_SEND_MPI(ZFIELD,PFIELD(:,JLAYER,:),IMASKF,NRANK,IDX)
          !
        ELSEIF (YDIR=='A' .OR. YDIR=='E') THEN
          !
          !NPIO needs to know all covers read
          IF (NRANK/=NPIO) THEN
            IDX = IDX + 1 
#ifdef SFX_MPI
            IF (LSFX_MPI) THEN
              CALL MPI_SEND(ZFIELD,SIZE(ZFIELD)*KIND(ZFIELD)/4,MPI_REAL,NPIO,IDX,NCOMM,INFOMPI)
            ENDIF
#endif
          ELSE
            CALL PACK_SAME_RANK(IMASKF,ZFIELD,PFIELD(:,JLAYER,:))
          ENDIF
          !
        ELSE
          CALL ABOR1_SFX("READ_SURF_LAYERS:HDIR MUST BE H OR A OR E")
        ENDIF
        !   
      ENDIF
      !
    ENDIF
    !
!    IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_3',1,ZHOOK_HANDLE)
    !
    IF (NRANK==NPIO .OR. YDIR=='H') THEN
      !
      !receives pieces of cover fields
      !ITREQ(:) = 0
      !
!$OMP PARALLEL PRIVATE(ZHOOK_HANDLE_OMP)
!IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_4',0,ZHOOK_HANDLE_OMP)
#ifdef SFX_MPI
!$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(JPROC,IDX,IVAL,ZWORKR,ISTATUS,INFOMPI)
#endif
      DO JPROC=0,NPROC-1
        !
        !the cover exists and was read
        IF (IPAS*JPROC + JP<=IL2) THEN
          !
          !IF (JPROC<NRANK) THEN
          !  ITREQ(JPROC+1) = JPROC+1
          !ELSE
          !  ITREQ(JPROC+1) = JPROC
          !ENDIF    
          !     
          IF (JPROC/=NRANK) THEN
            IDX = IDX_SAVE + JP + 1
            !each task receives the part of the cover read that concerns it 
            !only NPIO in cas of HDIR/=H
#ifdef SFX_MPI            
            IF (LSFX_MPI) THEN
              CALL MPI_RECV(ZWORKR(:,:),SIZE(ZWORKR)*KIND(ZWORKR)/4,&
                            MPI_REAL,JPROC,IDX,NCOMM,ISTATUS,INFOMPI)
            ENDIF
#endif
            IVAL = IPAS*JPROC + JP
            CALL PACK_SAME_RANK(IMASKF,ZWORKR(:,:),PFIELD(:,IVAL,:))
            !
          ENDIF
          !
        ENDIF
        !
      ENDDO
#ifdef SFX_MPI
!$OMP END DO
#endif
!IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_4',1,ZHOOK_HANDLE_OMP)
!$OMP END PARALLEL 

!
!      IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_5',0,ZHOOK_HANDLE)
      !
      !waits that all cover pieces are sent
#ifdef SFX_MPI
      IF (LSFX_MPI) THEN
        IF (YDIR=='H' .AND. IPAS*NRANK+JP<=IL2 .AND. NPROC>1) THEN
          CALL MPI_WAITALL(NPROC-1,NREQ(1:NPROC-1),ISTATUS,INFOMPI)
        ENDIF
      ENDIF
#endif
      !
!      IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_5',1,ZHOOK_HANDLE)
      !
!      IF (YDIR=='H' .OR. NRANK==NPIO) THEN
!        !packs data
!        IREQ = MAXVAL(ITREQ)
!!$OMP PARALLEL PRIVATE(ZHOOK_HANDLE_OMP)
!IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_6',0,ZHOOK_HANDLE_OMP)
!!$OMP DO SCHEDULE(DYNAMIC,1) PRIVATE(JPROC,IVAL)
!        DO JPROC=0,IREQ-1
!          IVAL = IPAS*JPROC + JP
!          IF (JPROC>=NRANK ) IVAL = IVAL + IPAS
!          CALL PACK_SAME_RANK(IMASKF,ZWORKR(:,:,JPROC+1),PFIELD(:,IVAL,:))
!        ENDDO
!!$OMP END DO
!IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_6',1,ZHOOK_HANDLE_OMP)
!!$OMP END PARALLEL
!      ENDIF
      !
    ENDIF
    !
  ENDDO
  !
!  IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_7',0,ZHOOK_HANDLE)
  !
  IF (NRANK/=NPIO .AND. YDIR=='H' .AND. IFLAG==0) THEN
    DEALLOCATE(NINDEX)
    ALLOCATE(NINDEX(0))
  ENDIF
  !
  IDX_R = IDX_R + IPAS + 1
  DEALLOCATE(ZWORKR)
  IF (HPROGRAM/="AROME ") DEALLOCATE(ZFIELD)
  IMASKF=>NULL()
  !  
!  IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_7',1,ZHOOK_HANDLE)
  !
ENDIF
!
!IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS_8',0,ZHOOK_HANDLE)
!
!RJ: what is a point of comment here? last field comment? Should be 'COVER_PACKED' status?
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('READ_SURF_LAYERS',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURF_LAYERS

END MODULE MODE_READ_SURF_LAYERS

