!     #########
SUBROUTINE CH_INIT_NAMES (KLUOUT,HSV, YSV, OVARSIGI, OVARSIGJ)  
!!    ###########################################
!!
!!*** *CH_INIT_NAMES*
!!
!!    PURPOSE
!!    -------
!!      Read and filter all chemical species into the CSV array
!!     initialize NSV_CHSBEG and  NSV_CHSEND index for the begin and the ending chemical index
!!     
!!
!!    REFERENCE
!!    ---------
!!    
!!    AUTHOR
!!    ------
!!    P. Tulet    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 16/10/01
!!    01/12/03    (D.Gazen) change emissions handling for surf. externalization
!!    01/06/05    (P.Tulet) add aerosols list
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_SV_n, ONLY : SV_t
!
USE MODD_CHS_AEROSOL
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!

INTEGER,                         INTENT(IN)  :: KLUOUT ! output listing channel
 CHARACTER(LEN=*), DIMENSION(:),  INTENT(IN)  :: HSV    ! name of chemical species
                                                       ! with character # (gas chemistry )
                                                       ! and  character @ (aerosols)
TYPE(SV_t), INTENT(INOUT) :: YSV
!
LOGICAL,                         INTENT(OUT) :: OVARSIGI, OVARSIGJ ! type of standard deviation
!
!*      0.2    declarations of local variables
INTEGER :: JSV  !! loop  NBEQ
 CHARACTER        :: YRC1
 CHARACTER(LEN=5) :: YRC2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CH_INIT_NAMES',0,ZHOOK_HANDLE)
YSV%NBEQ = 0
YSV%NAEREQ = 0
YSV%NSV_CHSBEG = 0
YSV%NSV_AERBEG = 0
YSV%NSV_CHSEND = 0
YSV%NSV_AEREND = 0
OVARSIGI = .FALSE.
OVARSIGJ = .FALSE.
NSOA = 0


DO JSV=1, SIZE(HSV)

  YSV%CSV(JSV) = HSV(JSV)
  YRC1= HSV(JSV)(1:1)
  YRC2 = HSV(JSV)(2:)

  IF (YRC1 == '#') THEN
    YSV%CSV(JSV) = TRIM(YRC2)
    YSV%NBEQ = YSV%NBEQ + 1
    IF (YSV%NBEQ == 1) YSV%NSV_CHSBEG=JSV
  ELSE IF (YRC1 == '@') THEN
    YSV%CSV(JSV) = TRIM(YRC2)
    YSV%NAEREQ = YSV%NAEREQ + 1
    IF (YSV%NAEREQ == 1) YSV%NSV_AERBEG=JSV
    IF (YSV%CSV(JSV) == "M6I") OVARSIGI = .TRUE.
    IF (YSV%CSV(JSV) == "M6J") OVARSIGJ = .TRUE.
    IF (YSV%CSV(JSV) == "SOA1I") NSOA = 10
  ENDIF

ENDDO

YSV%NSV_CHSEND = YSV%NSV_CHSBEG + YSV%NBEQ -1
YSV%NSV_AEREND = YSV%NSV_AERBEG + YSV%NAEREQ -1

IF (YSV%NAEREQ .GT. 0) THEN
DO JSV=1, size(YSV%CSV)
   IF (TRIM(YSV%CSV(JSV)) == "M0I") JP_CH_M0i=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "M0J") JP_CH_M0j=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "M6I") JP_CH_M6i=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "M6J") JP_CH_M6j=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "H2OI") JP_CH_H2Oi=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "H2OJ") JP_CH_H2Oj=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "SO4I") JP_CH_SO4i=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "SO4J") JP_CH_SO4j=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "NO3I") JP_CH_NO3i=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "NO3J") JP_CH_NO3j=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "NH3I") JP_CH_NH3i=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "NH3J") JP_CH_NH3j=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "OCI") JP_CH_OCi=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "OCJ") JP_CH_OCj=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "BCI") JP_CH_BCi=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "BCJ") JP_CH_BCj=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "DSTI") JP_CH_DSTi=JSV-YSV%NSV_CHSEND
   IF (TRIM(YSV%CSV(JSV)) == "DSTJ") JP_CH_DSTj=JSV-YSV%NSV_CHSEND
END DO

END IF
IF (LHOOK) CALL DR_HOOK('CH_INIT_NAMES',1,ZHOOK_HANDLE)

!
END SUBROUTINE CH_INIT_NAMES
