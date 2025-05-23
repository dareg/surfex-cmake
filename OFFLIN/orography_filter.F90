!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1
!     #########
      SUBROUTINE OROGRAPHY_FILTER(HPROGRAM,HGRID,PGRID_PAR,PSEA,KOPTFILTER,KZSFILTER,PCOFILTER,PTHFILTER,PZS)
!     ##############################################################
!
!!**** *OROGRAPHY_FILTER* filters the orography
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
!!    Original    06/2004
!!    Y. Seity : 09-2018 new options for filtering 
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODI_GET_GRID_DIM
USE MODI_ZSFILTER
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=6),     INTENT(IN)    :: HPROGRAM ! name of calling program
CHARACTER(LEN=10),    INTENT(IN)    :: HGRID    ! type of grid
REAL, DIMENSION(:),   INTENT(IN)    :: PGRID_PAR! list of parameters used to define the grid
REAL, DIMENSION(:),   INTENT(IN)    :: PSEA     ! sea  fraction
INTEGER,              INTENT(IN)    :: KOPTFILTER! filtering option
INTEGER,              INTENT(IN)    :: KZSFILTER! number of filter iteration
REAL,                 INTENT(IN)    :: PCOFILTER! filtering coefficient
REAL,                 INTENT(IN)    :: PTHFILTER! filtering threshold
REAL, DIMENSION(:),   INTENT(INOUT) :: PZS      ! orography
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
LOGICAL :: GRECT    ! true when grid is rectangular
INTEGER :: IX       ! number of points in X direction
INTEGER :: IY       ! number of points in Y direction
INTEGER :: JX       ! loop counter
INTEGER :: JY       ! loop counter
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZS ! orography in a 2D array
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSEA! sea fraction in a 2D array
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*    1.     Gets the geometry of the grid
!            -----------------------------
!
IF (LHOOK) CALL DR_HOOK('OROGRAPHY_FILTER',0,ZHOOK_HANDLE)
 CALL GET_GRID_DIM(HGRID,SIZE(PGRID_PAR),PGRID_PAR,GRECT,IX,IY)
!
!-------------------------------------------------------------------------------
!
!*    2.     If grid is not rectangular, nothing is done
!            -------------------------------------------
!
IF (.NOT. GRECT .AND. LHOOK) CALL DR_HOOK('OROGRAPHY_FILTER',1,ZHOOK_HANDLE)
IF (.NOT. GRECT) RETURN
!
IF (SIZE(PZS) /= IX * IY .AND. LHOOK) CALL DR_HOOK('OROGRAPHY_FILTER',1,ZHOOK_HANDLE)
IF (SIZE(PZS) /= IX * IY) RETURN
!
!-------------------------------------------------------------------------------
!
!*    3.     Grid rectangular: orography is put in a 2D array
!            ------------------------------------------------
!
ALLOCATE(ZZS (IX,IY))
ALLOCATE(ZSEA(IX,IY))
!
DO JY=1,IY
  DO JX=1,IX
    ZZS (JX,JY) = PZS ( JX + (JY-1)*IX ) 
    ZSEA(JX,JY) = PSEA( JX + (JY-1)*IX ) 
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*    4.     Filtering in x and Y directions
!            -------------------------------
!
IF (KZSFILTER>0) CALL ZSFILTER(HPROGRAM,ZZS,(1.-ZSEA),KOPTFILTER,KZSFILTER,PCOFILTER,PTHFILTER)
!
!-------------------------------------------------------------------------------
!
!*    5.     Output field comes back into 1D vector
!            --------------------------------------
!
DO JY=1,IY
  DO JX=1,IX
    PZS ( JX + (JY-1)*IX ) = ZZS(JX,JY)
  END DO
END DO
!
DEALLOCATE(ZZS )
DEALLOCATE(ZSEA)
IF (LHOOK) CALL DR_HOOK('OROGRAPHY_FILTER',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE OROGRAPHY_FILTER
