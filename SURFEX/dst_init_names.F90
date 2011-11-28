!     #########
SUBROUTINE DST_INIT_NAMES (KLUOUT,HSV,KDSTEQ,KSV_DSTBEG,KSV_DSTEND,OVARSIG,ORGFIX,HDSTYN)
!!    ###########################################
!!
!!*** *DST_INIT_NAMES*
!!
!!    PURPOSE
!!    -------
!!      Read and filter all chemical species into the CSV array
!!     initialize NSV_CHSBEG and  NSV_CHSEND index for the begin and the ending chemical index
!!     
!!
!!    REFERENCE
!!    ---------
!!    Modified ch_init_names (february 2005)    
!!
!!    AUTHOR
!!    ------
!!    Alf Grini <alf.grini@cnrm.meteo.fr>
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
USE MODD_DST_SURF, ONLY : JPMODE_DST
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!

INTEGER,                         INTENT(IN)  :: KLUOUT   ! output listing channel
CHARACTER(LEN=*), DIMENSION(:),  INTENT(IN)  :: HSV      ! name of chemical species
                                                         ! with character # for chemistry
INTEGER,                         INTENT(OUT) :: KDSTEQ     ! number of dust related variables
INTEGER,                         INTENT(OUT) :: KSV_DSTBEG  ! first dust related scalar
INTEGER,                         INTENT(OUT) :: KSV_DSTEND  ! last  dust related scalar
LOGICAL,                         INTENT(OUT) :: OVARSIG  !type of standard deviation
LOGICAL,                         INTENT(OUT) :: ORGFIX   !type of mean radius
CHARACTER(LEN=*), INTENT(OUT), OPTIONAL      :: HDSTYN      ! dust or not
!
!*      0.2    declarations of local variables
INTEGER :: JSV  !! loop on scalar variables
CHARACTER(LEN=4) :: YRC1
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------

!Initialize output variables
IF (LHOOK) CALL DR_HOOK('DST_INIT_NAMES',0,ZHOOK_HANDLE)
KDSTEQ = 0
KSV_DSTBEG=0
KSV_DSTEND=0
OVARSIG = .FALSE.
ORGFIX  = .TRUE.
IF(PRESENT(HDSTYN))HDSTYN='N'

DO JSV=1, SIZE(HSV)
  YRC1= HSV(JSV)(1:4)
  IF (YRC1 == 'DSTM') THEN
  IF (HSV(JSV)(5:5) == '6') OVARSIG = .TRUE.
  IF (HSV(JSV)(5:5) == '0') ORGFIX  = .FALSE.
     KDSTEQ = KDSTEQ + 1
     IF (KDSTEQ == 1) THEN
        KSV_DSTBEG=JSV
        IF(PRESENT(HDSTYN))HDSTYN='Y'
     ENDIF !Check on first time
  ELSE !Not dust variables
     !DO NOTHING
  ENDIF
ENDDO
!
! Set the output list of scalar to the input list of scalars

! Get the index of the last dust relevant tracer
KSV_DSTEND = KSV_DSTBEG + KDSTEQ -1
!
! Get number of dust modes. Each mode represents
! 3 moments, so 9 dust tracers represents 3 modes.
! 3 dust tracers represents 1 mode
IF (OVARSIG) THEN
 JPMODE_DST = INT((KSV_DSTEND - KSV_DSTBEG +1) / 3.)
ELSE IF (ORGFIX) THEN
 JPMODE_DST = INT((KSV_DSTEND - KSV_DSTBEG +1) )
ELSE
 JPMODE_DST = INT((KSV_DSTEND - KSV_DSTBEG +1) / 2.)
END IF
IF (LHOOK) CALL DR_HOOK('DST_INIT_NAMES',1,ZHOOK_HANDLE)


END SUBROUTINE DST_INIT_NAMES
