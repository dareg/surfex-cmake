!     #########
      SUBROUTINE PACK_PGD_ISBA(HPROGRAM,                                    &
                                 PAOSIP, PAOSIM, PAOSJP, PAOSJM,              &
                                 PHO2IP, PHO2IM, PHO2JP, PHO2JM,              &
                                 PSSO_SLOPE                                   )  
!     ##############################################################
!
!!**** *PACK_PGD_ISBA* packs ISBA physiographic fields from all surface points to ISBA points
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    03/2004
!!    Escobar J.  08/02/2005 : bug declare ILU local variable
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_ISBA_n, ONLY : I => ISBA
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
!
USE MODI_PACK_SAME_RANK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_SURF_MASK_n
!
USE MODI_GET_TYPE_DIM_n
!
USE MODI_GET_LUOUT
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
 CHARACTER(LEN=6),        INTENT(IN) :: HPROGRAM  ! Type of program
REAL,    DIMENSION(:),   INTENT(IN) :: PAOSIP    ! A/S i+ on all surface points
REAL,    DIMENSION(:),   INTENT(IN) :: PAOSIM    ! A/S i- on all surface points
REAL,    DIMENSION(:),   INTENT(IN) :: PAOSJP    ! A/S j+ on all surface points
REAL,    DIMENSION(:),   INTENT(IN) :: PAOSJM    ! A/S j- on all surface points
REAL,    DIMENSION(:),   INTENT(IN) :: PHO2IP    ! h/2 i+ on all surface points
REAL,    DIMENSION(:),   INTENT(IN) :: PHO2IM    ! h/2 i- on all surface points
REAL,    DIMENSION(:),   INTENT(IN) :: PHO2JP    ! h/2 j+ on all surface points
REAL,    DIMENSION(:),   INTENT(IN) :: PHO2JM    ! h/2 j- on all surface points
REAL,    DIMENSION(:),   INTENT(IN) :: PSSO_SLOPE! subgrid slope on all surface points
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER                        :: ILU    ! expected physical size of full surface array
INTEGER                        :: ILUOUT ! output listing logical unit
INTEGER, DIMENSION(:), POINTER :: IMASK  ! mask for packing from complete field to nature field
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PACK_PGD_ISBA',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!*    1.      Number of points and packing
!             ----------------------------
!
 CALL GET_TYPE_DIM_n(DTCO, U, &
                     'NATURE',IG%NDIM)
ALLOCATE(IMASK(IG%NDIM))
ILU=0
 CALL GET_SURF_MASK_n(DTCO, U, &
                      'NATURE',IG%NDIM,IMASK,ILU,ILUOUT)
!
!
!-------------------------------------------------------------------------------
!
!*    2.      Packing of fields
!             -----------------
!
ALLOCATE(I%XAOSIP(IG%NDIM))
ALLOCATE(I%XAOSIM(IG%NDIM))
ALLOCATE(I%XAOSJP(IG%NDIM))
ALLOCATE(I%XAOSJM(IG%NDIM))
ALLOCATE(I%XHO2IP(IG%NDIM))
ALLOCATE(I%XHO2IM(IG%NDIM))
ALLOCATE(I%XHO2JP(IG%NDIM))
ALLOCATE(I%XHO2JM(IG%NDIM))
ALLOCATE(I%XSSO_SLOPE(IG%NDIM))
 CALL PACK_SAME_RANK(IMASK,PAOSIP(:),I%XAOSIP(:))
 CALL PACK_SAME_RANK(IMASK,PAOSIM(:),I%XAOSIM(:))
 CALL PACK_SAME_RANK(IMASK,PAOSJP(:),I%XAOSJP(:))
 CALL PACK_SAME_RANK(IMASK,PAOSJM(:),I%XAOSJM(:))
 CALL PACK_SAME_RANK(IMASK,PHO2IP(:),I%XHO2IP(:))
 CALL PACK_SAME_RANK(IMASK,PHO2IM(:),I%XHO2IM(:))
 CALL PACK_SAME_RANK(IMASK,PHO2JP(:),I%XHO2JP(:))
 CALL PACK_SAME_RANK(IMASK,PHO2JM(:),I%XHO2JM(:))
 CALL PACK_SAME_RANK(IMASK,PSSO_SLOPE(:),I%XSSO_SLOPE(:))
IF (LHOOK) CALL DR_HOOK('PACK_PGD_ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PACK_PGD_ISBA
