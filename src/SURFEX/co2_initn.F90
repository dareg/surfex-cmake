!     #########
      SUBROUTINE CO2_INIT_n (IO, IP, PCO2, PGMES, PGC, PDMAX        )
!     #####################
!
!!****  *CO2_INIT_n* - routine to initialize ISBA-AGS variables
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2003 
!!      J.C. Calvet 01/2004 Externalization
!!      P Le Moigne 11/2004 cotwoinit changed into cotwoinit_n
!!      P Le Moigne 09/2005 AGS modifs of L. Jarlan
!!      S Lafont    09/2008 Add initialisation of POI and ABC (needed for TORI)
!!      A.L. Gibelin 04/2009 : TAU_WOOD for NCB option 
!!      A.L. Gibelin 04/2009 : Add carbon spinup
!!      A.L. Gibelin 07/2009 : Suppress RDK and transform GPP as a diagnostic
!!      A.L. Gibelin 07/2009 : Suppress PPST and PPSTF as outputs
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_PGD_INIT
!
USE MODD_SURFEX_MPI, ONLY : NRANK,NPIO
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
!
USE MODI_COTWOINIT_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
!
REAL, DIMENSION(:), INTENT(IN) :: PCO2 ! air CO2 concentration (kg/kg)
!
REAL, DIMENSION(:,:), INTENT(IN) :: PGMES
REAL, DIMENSION(:,:), INTENT(IN) :: PGC
REAL, DIMENSION(:,:), INTENT(IN) :: PDMAX
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
TYPE(ISBA_INIT_PGD_t) :: YP_IP
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_GMES           ! 
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_CO2            ! air CO2 concentration (kg/kg)
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_GC             !
REAL, DIMENSION(:),   ALLOCATABLE :: ZP_DMAX           !
!
INTEGER :: ILU   ! size of arrays
INTEGER :: IPATCH
INTEGER :: INBIOMASS
INTEGER :: JP    ! loop on tiles
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CO2_INIT_N',0,ZHOOK_HANDLE)
!
ILU    = SIZE(IP%XVEGTYPE_PATCH,1)
IPATCH = SIZE(IP%XVEGTYPE_PATCH,3)
INBIOMASS = SIZE(IP%XINCREASE,2)
!
DO JP=1,IPATCH
!
  IF (IP%NSIZE_NATURE_P(JP) == 0 ) CYCLE
!
  IF (MAXVAL(PGMES(:,JP)).NE.XUNDEF .OR. MINVAL(PGMES(:,JP)).NE.XUNDEF) THEN

     CALL PACK_CO2_INIT(IP%NR_NATURE_P(:,JP),IP%NSIZE_NATURE_P(JP),JP)
!
     CALL COTWOINIT_n(IO, ZP_GMES,ZP_CO2,ZP_GC, ZP_DMAX,YP_IP )  

     YP_IP%XINCREASE = 0.
     YP_IP%XTURNOVER = 0.
!
     CALL UNPACK_CO2_INIT(IP%NR_NATURE_P(:,JP),IP%NSIZE_NATURE_P(JP),JP)

  ENDIF

ENDDO
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CO2_INIT_N',1,ZHOOK_HANDLE)
CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE PACK_CO2_INIT(KMASK,KSIZE,KPATCH)
IMPLICIT NONE
INTEGER, INTENT(IN)               :: KSIZE, KPATCH
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
!
INTEGER JJ, JI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('PACK_CO2_INIT',0,ZHOOK_HANDLE)
!
ALLOCATE(ZP_GMES         (KSIZE))
ALLOCATE(ZP_CO2          (KSIZE))
ALLOCATE(ZP_GC           (KSIZE))
ALLOCATE(ZP_DMAX         (KSIZE))
!
ALLOCATE(YP_IP%XVEGTYPE_PATCH(KSIZE,NVEGTYPE,1))
ALLOCATE(YP_IP%XANMAX        (KSIZE,1))
ALLOCATE(YP_IP%XFZERO        (KSIZE,1))
ALLOCATE(YP_IP%XEPSO         (KSIZE,1))
ALLOCATE(YP_IP%XGAMM         (KSIZE,1))
ALLOCATE(YP_IP%XQDGAMM       (KSIZE,1))
ALLOCATE(YP_IP%XQDGMES       (KSIZE,1))
ALLOCATE(YP_IP%XT1GMES       (KSIZE,1))
ALLOCATE(YP_IP%XT2GMES       (KSIZE,1))
ALLOCATE(YP_IP%XAMAX         (KSIZE,1))
ALLOCATE(YP_IP%XQDAMAX       (KSIZE,1))
ALLOCATE(YP_IP%XT1AMAX       (KSIZE,1))
ALLOCATE(YP_IP%XT2AMAX       (KSIZE,1))
ALLOCATE(YP_IP%XAH           (KSIZE,1))
ALLOCATE(YP_IP%XBH           (KSIZE,1))
ALLOCATE(YP_IP%XTAU_WOOD     (KSIZE,1))
ALLOCATE(YP_IP%XINCREASE     (KSIZE,INBIOMASS,1))
ALLOCATE(YP_IP%XTURNOVER     (KSIZE,INBIOMASS,1))
!
! initialisation needed for TORI
ALLOCATE(YP_IP%XABC(SIZE(IP%XABC)))
ALLOCATE(YP_IP%XPOI(SIZE(IP%XPOI)))
YP_IP%XABC(:)=0.
YP_IP%XPOI(:)=0.
!
DO JJ=1,KSIZE
  JI                     =    KMASK(JJ)
  YP_IP%XVEGTYPE_PATCH(JJ,:,1) =    IP%XVEGTYPE_PATCH(JI,:,KPATCH)
  ZP_GMES         (JJ)   =    PGMES         (JI,KPATCH)
  ZP_CO2          (JJ)   =    PCO2          (JI)
  ZP_GC           (JJ)   =    PGC           (JI,KPATCH)
  ZP_DMAX         (JJ)   =    PDMAX         (JI,KPATCH)
END DO
IF (LHOOK) CALL DR_HOOK('PACK_CO2_INIT',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE PACK_CO2_INIT
!-------------------------------------------------------------------------------
SUBROUTINE UNPACK_CO2_INIT(KMASK,KSIZE,KPATCH)
IMPLICIT NONE
INTEGER, INTENT(IN)               :: KSIZE, KPATCH
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
!
INTEGER JJ, JI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('UNPACK_CO2_INIT',0,ZHOOK_HANDLE)
IP%XANMAX     (:,KPATCH) = XUNDEF
IP%XFZERO     (:,KPATCH) = XUNDEF
IP%XEPSO      (:,KPATCH) = XUNDEF
IP%XGAMM      (:,KPATCH) = XUNDEF
IP%XQDGAMM    (:,KPATCH) = XUNDEF
IP%XQDGMES    (:,KPATCH) = XUNDEF
IP%XT1GMES    (:,KPATCH) = XUNDEF
IP%XT2GMES    (:,KPATCH) = XUNDEF
IP%XAMAX      (:,KPATCH) = XUNDEF
IP%XQDAMAX    (:,KPATCH) = XUNDEF
IP%XT1AMAX    (:,KPATCH) = XUNDEF
IP%XT2AMAX    (:,KPATCH) = XUNDEF
IP%XAH        (:,KPATCH) = XUNDEF
IP%XBH        (:,KPATCH) = XUNDEF
IP%XTAU_WOOD  (:,KPATCH) = XUNDEF
IP%XINCREASE  (:,:,KPATCH) = XUNDEF
IP%XTURNOVER  (:,:,KPATCH) = XUNDEF
!
DO JJ=1,KSIZE
   JI                              = KMASK         (JJ)
   IP%XANMAX          (JI, KPATCH)    = YP_IP%XANMAX      (JJ,1)
   IP%XFZERO          (JI, KPATCH)    = YP_IP%XFZERO      (JJ,1)
   IP%XEPSO           (JI, KPATCH)    = YP_IP%XEPSO       (JJ,1)
   IP%XGAMM           (JI, KPATCH)    = YP_IP%XGAMM       (JJ,1)
   IP%XQDGAMM         (JI, KPATCH)    = YP_IP%XQDGAMM     (JJ,1)
   IP%XQDGMES         (JI, KPATCH)    = YP_IP%XQDGMES     (JJ,1)
   IP%XT1GMES         (JI, KPATCH)    = YP_IP%XT1GMES     (JJ,1)
   IP%XT2GMES         (JI, KPATCH)    = YP_IP%XT2GMES     (JJ,1)
   IP%XAMAX           (JI, KPATCH)    = YP_IP%XAMAX       (JJ,1)
   IP%XQDAMAX         (JI, KPATCH)    = YP_IP%XQDAMAX     (JJ,1)
   IP%XT1AMAX         (JI, KPATCH)    = YP_IP%XT1AMAX     (JJ,1)
   IP%XT2AMAX         (JI, KPATCH)    = YP_IP%XT2AMAX     (JJ,1)
   IP%XAH             (JI, KPATCH)    = YP_IP%XAH         (JJ,1)
   IP%XBH             (JI, KPATCH)    = YP_IP%XBH         (JJ,1)
   IP%XTAU_WOOD       (JI, KPATCH)    = YP_IP%XTAU_WOOD   (JJ,1)
   IP%XINCREASE       (JI, :, KPATCH) = YP_IP%XINCREASE   (JJ, :,1)
   IP%XTURNOVER       (JI, :, KPATCH) = YP_IP%XTURNOVER   (JJ, :,1)
END DO
! 
DO JJ=1,SIZE(IP%XABC)
   IP%XABC(JJ)=YP_IP%XABC(JJ)
   IP%XPOI(JJ)=YP_IP%XPOI(JJ)
END DO 
!
DEALLOCATE(ZP_GMES         )
DEALLOCATE(ZP_CO2          )
DEALLOCATE(ZP_GC           )
DEALLOCATE(ZP_DMAX         )
!
CALL ISBA_INIT_PGD_INIT(YP_IP)
!
IF (LHOOK) CALL DR_HOOK('UNPACK_CO2_INIT',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE UNPACK_CO2_INIT
!-------------------------------------------------------------------------------
!
END SUBROUTINE CO2_INIT_n
