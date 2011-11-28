SUBROUTINE INIT_SLT_n(HPROGRAM  &! Program calling unit
       )  

USE MODD_SLT_n
USE MODD_SLT_SURF
USE MODD_ISBA_n, ONLY:          &
       NSIZE_NATURE_P              &! Number of nature points in a patch
       ,NR_NATURE_P                &! Mask from patch --> nature vectors
       ,NPATCH                     &! Maximum number of patches
       ,XVEGTYPE_PATCH              ! fraction (in a nature point) of a vegtype for a patch  
USE MODD_DATA_COVER_PAR, ONLY: NVT_NO, NVT_ROCK
USE MODD_SURF_ATM_n, ONLY: NSIZE_NATURE   ! Number of nature points

!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_LUOUT
IMPLICIT NONE

!PASSED VARIABLES
CHARACTER(LEN=*)      :: HPROGRAM              !Passing unit

!LOCAL VARIABLES
INTEGER             :: JVEG                  ! Counter for vegetation classes
INTEGER             :: JVEG_IN               ! Vegetation index
INTEGER             :: JPATCH                ! Counter for patches
INTEGER             :: JMODE                 ! Counter for sea salt modes
INTEGER             :: JMODE_IDX             ! Index for sea salt modes
INTEGER             :: ILUOUT
INTEGER             :: IMI                   ! Current model
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!get output listing unit
IF (LHOOK) CALL DR_HOOK('INIT_SLT_N',0,ZHOOK_HANDLE)
CALL GET_LUOUT(HPROGRAM,ILUOUT)

!Allocate memory for the real values which will be used by the model
ALLOCATE(XEMISRADIUS_SLT(NSLTMDE))
ALLOCATE(XEMISSIG_SLT(NSLTMDE))

!Get initial size distributions. This is cut and pasted
!from dead routine dstpsd.F90
!Check for different source parameterizations
IF(CEMISPARAM_SLT.eq."Vig01")THEN
   CRGUNITS   = 'NUMB'
   XEMISRADIUS_INI_SLT(:) =  (/ 0.2,  2.0, 12. /)         ! [um]  Number median radius She84 p. 75 Table 1
   XEMISSIG_INI_SLT(:) = (/ 1.9     ,  2.0     ,  3.00    /)  ! [frc] Geometric standard deviation She84 p. 75 Table 1
ELSE  ! use default of Schultz et al, 2004
  CRGUNITS   = 'MASS'
  XEMISRADIUS_INI_SLT(:)=0.5*(/0.28, 2.25, 15.32/)! [um] Mass median radius
  XEMISSIG_INI_SLT(:) =   (/1.59, 2.00, 2.00/) ! [frc] Geometric standard deviation
ENDIF


DO JMODE=1,NSLTMDE
   JMODE_IDX=JORDER_SLT(JMODE)
   XEMISSIG_SLT(JMODE)=XEMISSIG_INI_SLT(JMODE_IDX)
IF (CRGUNITS=="MASS") THEN
   XEMISRADIUS_SLT(JMODE) =  XEMISRADIUS_INI_SLT(JMODE_IDX)  *&
               EXP(-3.d0 * (LOG(XEMISSIG_INI_SLT(JMODE_IDX)))**2)  
ELSE
   XEMISRADIUS_SLT(JMODE) =  XEMISRADIUS_INI_SLT(JMODE_IDX)    
ENDIF
ENDDO
IF (LHOOK) CALL DR_HOOK('INIT_SLT_N',1,ZHOOK_HANDLE)

END SUBROUTINE INIT_SLT_n
