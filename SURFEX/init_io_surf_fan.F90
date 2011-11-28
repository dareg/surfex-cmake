!     ######################
      SUBROUTINE INIT_IO_SURF_FA_n(HPROGRAM,HMASK,HACTION)
!     ######################
!
!!****  *INIT_IO_SURF_FA* Keep in memory the output files
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      P. Le Moigne 04/2004: distinguish in and out file name
!!      P. Le Moigne 04/2006: special HACTION='GTMSK' to initialize
!!                            a mask different of 'FULL ' in order 
!!                            to read dimensions only.
!!      B. Decharme   2008  : Change to switch between offline and online run
!!                            In online run, the mask must be always global
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_SURF_ATM_n, ONLY: NDIM_FULL
!
USE MODD_IO_SURF_FA,ONLY: NUNIT_FA, CFILEIN_FA,CFILEOUT_FA,CDNOMC,IVERBFA,   &
                           NLUOUT,NFULL,NFULL_EXT, CMASK, LOPEN,             &
                           NDGL, NDLON, NDLUX, NDGUX, PERPK, PEBETA,         &
                           PELON0, PELAT0, PEDELX, PEDELY, PELON1, PELAT1 
USE MODI_GET_LUOUT
USE MODI_READ_SURF
USE MODD_CSTS, ONLY : XPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_SURF_MASK_n
USE MODI_GET_TYPE_DIM_n
USE MODI_ABOR1_SFX
USE MODI_GET_1D_MASK
!
IMPLICIT NONE
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM
CHARACTER(LEN=6),  INTENT(IN)  :: HMASK    
CHARACTER(LEN=5),  INTENT(IN)  :: HACTION    
!
INTEGER                        :: ILU,IRET, IL
CHARACTER(LEN=6)               :: YMASK
!
INTEGER                :: ISIZE
INTEGER                :: INB ! number of articles in the file
INTEGER                :: ITYPTR, ITRONC, INLATI, INXLON, INIVER
INTEGER, DIMENSION (1000) :: INLOPA, INOZPA
!
REAL, DIMENSION (1000)  :: ZSINLA
REAL, DIMENSION (200)   :: ZAHYBR, ZBHYBR
REAL                    :: ZSLAPO, ZCLOPO, ZSLOPO, ZCODIL, ZREFER
LOGICAL                 :: LOUTFAC
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_FA_N',0,ZHOOK_HANDLE)

IF(HPROGRAM/='FA    '.AND.HPROGRAM/='AROME ') THEN
  CALL ABOR1_SFX('INIT_IO_SURF_FA_N -- HPROGRAM should be FA or AROME')
ENDIF
CALL GET_LUOUT(HPROGRAM,NLUOUT)

LOPEN=.FALSE.
!
ILU  =0
YMASK=HMASK
!
IF (HACTION=='GTMSK') THEN
   CALL FAITOU(IRET,NUNIT_FA,.TRUE.,CFILEIN_FA,'OLD',.TRUE.,.FALSE.,IVERBFA,0,INB,CDNOMC)
   WRITE(NLUOUT,*)'HPROGRAM ',HPROGRAM,' IO_INIT HACTION==GTMSK',NUNIT_FA,CFILEIN_FA
   CMASK = YMASK
   LOPEN=.TRUE.
   IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_FA_N',1,ZHOOK_HANDLE)
   RETURN
ENDIF
!
IF (HACTION == 'READ ') THEN
   CALL FAITOU(IRET,NUNIT_FA,.TRUE.,CFILEIN_FA,'OLD',.TRUE.,.FALSE.,IVERBFA,0,INB,CDNOMC)
   WRITE(NLUOUT,*)'HPROGRAM ',HPROGRAM,' IO_INIT HACTION==READ',NUNIT_FA,CFILEIN_FA
   CALL FACAGE(CDNOMC,.TRUE.)

   IF (YMASK == 'FULL  ') THEN
      CMASK = YMASK
      CALL READ_SURF(HPROGRAM,'DIM_FULL',ILU,IRET)
      NFULL = ILU
      NFULL_EXT = ILU
      IF (HPROGRAM=='AROME ') THEN
         NDIM_FULL=ILU
      ENDIF
   ENDIF
   IF (YMASK == 'EXTZON') THEN
      CALL FACIES(CDNOMC, ITYPTR, ZSLAPO, ZCLOPO, ZSLOPO,       &
                           ZCODIL, ITRONC, INLATI, INXLON, INLOPA, &
                           INOZPA, ZSINLA, INIVER, ZREFER, ZAHYBR, &
                           ZBHYBR, LOUTFAC) 
      NFULL_EXT = INLATI*INXLON
      NDGL   = INLATI
      NDLON  = INXLON
      NFULL  = INLOPA(4)*INLOPA(6)
      NDLUX  = INLOPA(4)
      NDGUX  = INLOPA(6)
      PEBETA = ZSLAPO
      PERPK  = ZSINLA(2)
      PELON0 = ZSINLA(3)*180./XPI
      PELAT0 = ZSINLA(4)*180./XPI
      PEDELX = ZSINLA(7)
      PEDELY = ZSINLA(8)
      PELON1 = ZSINLA(13)*180./XPI
      PELAT1 = ZSINLA(14)*180./XPI
   ENDIF 
   LOPEN=.TRUE.
ENDIF
!
!------------------------------------------------------------------------------
CMASK=YMASK
!------------------------------------------------------------------------------
!
IF (HPROGRAM=='AROME ') THEN
  NFULL = NDIM_FULL
  ISIZE = NFULL
ELSE
  ISIZE = NFULL
  CALL GET_TYPE_DIM_n(YMASK,ISIZE)
ENDIF
!
CALL GET_MASK(ISIZE)
!
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_FA_N',1,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
!
CONTAINS
!
SUBROUTINE GET_MASK(KSIZE)
!
USE MODD_MASK,       ONLY: NMASK_FULL
USE MODD_IO_SURF_FA, ONLY: NMASK, NFULL
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KSIZE
!
REAL, DIMENSION(KSIZE) :: ZFULL
INTEGER, DIMENSION(KSIZE) :: IMASK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_FA_N:GET_MASK',0,ZHOOK_HANDLE)
!
IF (HPROGRAM=='AROME ') THEN
  ZFULL = 1.
  CALL GET_1D_MASK(NFULL,NFULL,ZFULL,IMASK)
ELSE
  CALL GET_SURF_MASK_n(YMASK,KSIZE,IMASK,NFULL,NLUOUT)
ENDIF
!
IF (.NOT.ALLOCATED(NMASK_FULL)) ALLOCATE(NMASK_FULL(NFULL))
NMASK_FULL(:)=0
!
NMASK => NMASK_FULL(1:KSIZE)
NMASK(:) = IMASK(:)
!
IF (LHOOK) CALL DR_HOOK('INIT_IO_SURF_FA_N:GET_MASK',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_MASK
!
END SUBROUTINE INIT_IO_SURF_FA_n

