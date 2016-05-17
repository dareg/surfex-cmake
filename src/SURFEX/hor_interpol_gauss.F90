!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE HOR_INTERPOL_GAUSS(KLUOUT,PFIELDIN,PFIELDOUT)
!     #################################################################################
!
!!****  *HOR_INTERPOL_GAUSS* - Interpolation from a gaussian grid
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      M. Jidane   Dec 2013 : initialize NNI if not already done
!!------------------------------------------------------------------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO
USE MODD_HORIBL, ONLY : LGLOBLON, LGLOBS, LGLOBN, XILO1H, XILO2H, NINLOH, &
                        XLA, XOLA, XOLO, NP, XLOPH
USE MODD_PREP,       ONLY : XLAT_OUT, XLON_OUT, LINTERP
USE MODD_GRID_GAUSS, ONLY : XILA1, XILO1, XILA2, XILO2, NINLA, NINLO, NILEN, LROTPOLE, &
                              XLAP, XLOP, XCOEF, XLAT, XLON
USE MODD_GRID_GRIB,  ONLY : NNI
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODI_HORIBL_SURF_GRIDIN
USE MODI_HORIBL_SURF_VALUE
USE MODI_HORIBL_SURF_EXTRAP
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER,            INTENT(IN)  :: KLUOUT    ! logical unit of output listing
REAL, DIMENSION(:,:), INTENT(IN)    :: PFIELDIN  ! field to interpolate horizontally
REAL, DIMENSION(:,:), INTENT(OUT)   :: PFIELDOUT ! interpolated field
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(:,:), POINTER :: ZFIELDIN0=>NULL()
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZFIELDIN 
!
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ILSMIN  
INTEGER, DIMENSION(:,:), ALLOCATABLE :: IMASKIN  ! input mask
INTEGER, DIMENSION(:), ALLOCATABLE :: IMASKOUT ! output mask
INTEGER                         :: INO, INL     ! output number of points
INTEGER                         :: JL       ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!*      1.    Allocations
!
IF (LHOOK) CALL DR_HOOK('HOR_INTERPOL_GAUSS',0,ZHOOK_HANDLE)
INO = SIZE(XLAT_OUT)
INL = SIZE(PFIELDOUT,2)
!
ALLOCATE(IMASKOUT(INO))
IMASKOUT = 1
!
IF (NNI==0) NNI=NILEN
ALLOCATE(IMASKIN (NNI,INL))
!
IF (NRANK==NPIO) THEN
  IMASKIN(:,:) = 1.
  WHERE(PFIELDIN(:,:)==XUNDEF) IMASKIN(:,:) = 0.
ENDIF

ALLOCATE(ZFIELDIN(INO,INL,12))
ALLOCATE(ILSMIN(INO,INL,12))
!
CALL HORIBL_SURF_GRIDIN(NINLA,NINLO,NNI,PFIELDIN(:,:),INO, &
                        .FALSE.,KLUOUT,LGLOBS,LGLOBN,LGLOBLON,NP, &
                        ZFIELDIN0,ZFIELDIN,ILSMIN,IMASKIN,IMASKOUT)
!
DO JL=1,SIZE(PFIELDOUT,2)
  !
  CALL HORIBL_SURF_VALUE(NNI,INO,PFIELDOUT(:,JL),LINTERP,ZFIELDIN(:,JL,:),ILSMIN(:,JL,:),&
                         XOLO,XOLA,XLA,XLOPH,IMASKIN(:,JL),IMASKOUT)
  !
ENDDO
!
!IF (LGLOBLON.OR.LGLOBN.OR.LGLOBS) THEN
!  CALL HORIBL_SURF_EXTRAP(XILA1,XILO1H,XILA2,XILO2H,NINLA,NINLO,NNI,ZFIELDIN0,&
!                         INO,NP,XLON,XLAT,PFIELDOUT,KLUOUT,LINTERP)
!ELSE
  CALL HORIBL_SURF_EXTRAP(XILA1,XILO1H,XILA2,XILO2H,NINLA,NINLO,NNI,PFIELDIN,&
                         INO,NP,XLON,XLAT,PFIELDOUT,KLUOUT,LINTERP)
!ENDIF
!
!*      5.    Deallocations
!
DEALLOCATE(IMASKIN )
DEALLOCATE(IMASKOUT)
DEALLOCATE(ZFIELDIN0)
DEALLOCATE(ZFIELDIN)
DEALLOCATE(ILSMIN)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HOR_INTERPOL_GAUSS',1,ZHOOK_HANDLE)
!
END SUBROUTINE HOR_INTERPOL_GAUSS
