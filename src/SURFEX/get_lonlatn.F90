!     #########
      SUBROUTINE GET_LONLAT_n (YSC, &
                               HPROGRAM)
!     ####################################
!
!!****  *GET_LONLAT_n* - routine to get some surface fields
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
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2008
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_SURFEX_n, ONLY : SURFEX_t
!
USE MODI_GET_LUOUT
USE MODI_GET_COORD_n
USE MODI_GET_SURF_SIZE_n
USE MODI_WRITE_SURF
USE MODI_INIT_IO_SURF_n
USE MODI_END_IO_SURF_n
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_CLEAN_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(SURFEX_t), INTENT(INOUT) :: YSC
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: ILUOUT
!
INTEGER            :: IRET      
 CHARACTER(LEN=100) :: YCOMMENT
!
INTEGER            :: INI      
REAL, DIMENSION(:), ALLOCATABLE :: ZLON, ZLAT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_LONLAT_N',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!-------------------------------------------------------------------------------
!
 CALL GET_SURF_SIZE_n(YSC%DTCO, YSC%U, &
                      'FULL', INI)
!
ALLOCATE(ZLON(INI))
ALLOCATE(ZLAT(INI))
!
 CALL GET_COORD_n(YSC%UG, &
                  HPROGRAM,INI,ZLON,ZLAT)      
!
 CALL IO_BUFF_CLEAN_n(YSC%IOB)
 CALL INIT_IO_SURF_n(YSC%BOP, YSC%BDD, YSC%CHE, YSC%CHI, YSC%CHS, YSC%CHN, YSC%CHU, YSC%CHT, YSC%CHW, &
                     YSC%DTCO, YSC%DTS, YSC%DTT, YSC%DTZ, YSC%DGEI, YSC%DGF, YSC%DGI, YSC%DGMI, YSC%DGMTO, &
                     YSC%DGO, YSC%DGS, YSC%DGSI, YSC%DGU, YSC%DGT, YSC%DGUT, YSC%DGW, YSC%F, YSC%FSB, &
                     YSC%GB, YSC%IOB, YSC%ICP, YSC%I, YSC%O, YSC%S, YSC%SSB, YSC%UG, YSC%U, YSC%SV, &
                     YSC%TCP, YSC%TGD, YSC%TGDO, YSC%TGR, YSC%TGRO, YSC%T, YSC%TOP, YSC%TVG, YSC%W, &
                     YSC%WSB, &
                     HPROGRAM,'FULL  ','SURF  ','WRITE')
!
YCOMMENT='XLON'
 CALL WRITE_SURF(YSC%DGU, YSC%IOB, YSC%U, &
                 HPROGRAM,'XLON',ZLON(:),IRET,HCOMMENT=YCOMMENT,HDIR='A')
!
YCOMMENT='XLAT'
 CALL WRITE_SURF(YSC%DGU, YSC%IOB, YSC%U, &
                 HPROGRAM,'XLAT',ZLAT(:),IRET,HCOMMENT=YCOMMENT,HDIR='A')
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('GET_LONLAT_N',1,ZHOOK_HANDLE)
!
!==============================================================================
!
END SUBROUTINE GET_LONLAT_n
