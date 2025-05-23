!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########################
      SUBROUTINE WRITE_BLD_DESCRIPTION_n (HSELECT, BDD, HPROGRAM)
!     #########################
!
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!
!!
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
!!    Original    05/2012 
!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
USE MODD_BLD_DESCRIPTION_n, ONLY : BLD_DESC_t
!
USE MODI_WRITE_SURF
USE MODI_ABOR1_SFX
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
 CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSELECT
!
TYPE(BLD_DESC_t), INTENT(INOUT) :: BDD
!
 CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM
!
!
!*    0.2    Declaration of local variables
!      ------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
REAL, DIMENSION(:), ALLOCATABLE :: ZWORK
INTEGER                         :: IRESP
INTEGER                         :: I1, I2
INTEGER                         :: JL
INTEGER                         :: ITOT
 CHARACTER(LEN=100)              :: YCOMMENT
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_BLD_DESCRIPTION_n',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*    1.   Writes configuration variables of the descriptive data
!          ------------------------------------------------------
!
ALLOCATE(ZWORK(7))
!
ZWORK(1) = FLOAT(BDD%NDESC_BLD)
ZWORK(2) = FLOAT(BDD%NDESC_AGE)
ZWORK(3) = FLOAT(BDD%NDESC_USE)
ZWORK(4) = FLOAT(BDD%NDESC_WALL_LAYER)
ZWORK(5) = FLOAT(BDD%NDESC_ROOF_LAYER)
ZWORK(6) = FLOAT(BDD%NDESC_ROAD_LAYER)
ZWORK(7) = FLOAT(BDD%NDESC_FLOOR_LAYER)
!
YCOMMENT='Configuration numbers for descriptive building data'
 CALL WRITE_SURF(HSELECT, HPROGRAM,'BLD_DESC_CNF',ZWORK,IRESP,YCOMMENT,'-','Bld_dimensions  ')
DEALLOCATE(ZWORK)
!
!-------------------------------------------------------------------------------
!
!*    3.   Writes descriptive data
!          -----------------------
!
ITOT=(21+3*BDD%NDESC_ROOF_LAYER+3*BDD%NDESC_ROAD_LAYER+3*BDD%NDESC_WALL_LAYER+3*BDD%NDESC_FLOOR_LAYER)*BDD%NDESC_CODE &
      + 9*BDD%NDESC_USE+2*BDD%NDESC_AGE+BDD%NDESC_BLD
ALLOCATE(ZWORK(ITOT))
!
!
I1=0 ; I2=0
 CALL UP_DESC_IND_W(BDD%NDESC_BLD)  ; ZWORK(I1:I2) = FLOAT(BDD%NDESC_BLD_LIST(:))
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = FLOAT(BDD%NDESC_CODE_LIST(:))
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_ALB_ROOF(:)
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_ALB_ROAD(:)
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_ALB_WALL(:)
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_EMIS_ROOF(:)
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_EMIS_ROAD(:)
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_EMIS_WALL(:)
DO JL=1,BDD%NDESC_ROOF_LAYER
  CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_HC_ROOF(:,JL)
END DO
DO JL=1,BDD%NDESC_ROOF_LAYER
  CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_TC_ROOF(:,JL)
END DO
DO JL=1,BDD%NDESC_ROOF_LAYER
  CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_D_ROOF (:,JL) 
END DO
DO JL=1,BDD%NDESC_ROAD_LAYER
  CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_HC_ROAD(:,JL)
END DO
DO JL=1,BDD%NDESC_ROAD_LAYER
  CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_TC_ROAD(:,JL) 
END DO
DO JL=1,BDD%NDESC_ROAD_LAYER
  CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_D_ROAD (:,JL)
END DO
DO JL=1,BDD%NDESC_WALL_LAYER
  CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_HC_WALL(:,JL)
END DO
DO JL=1,BDD%NDESC_WALL_LAYER
  CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_TC_WALL(:,JL) 
END DO
DO JL=1,BDD%NDESC_WALL_LAYER
  CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_D_WALL (:,JL)
END DO
DO JL=1,BDD%NDESC_FLOOR_LAYER
  CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_HC_FLOOR(:,JL)
END DO
DO JL=1,BDD%NDESC_FLOOR_LAYER
  CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_TC_FLOOR(:,JL) 
END DO
DO JL=1,BDD%NDESC_FLOOR_LAYER
  CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_D_FLOOR (:,JL)
END DO
!
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_SHGC(:)
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_U_WIN(:)
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_GR(:) 
!
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_F_WASTE_CAN(:)
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_F_WATER_COND(:)
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_COP_RAT(:)
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_EFF_HEAT(:)
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_INF(:)
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_V_VENT(:) 
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_GREENROOF(:) 
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_EMIS_PANEL(:) 
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_ALB_PANEL(:) 
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_EFF_PANEL(:) 
 CALL UP_DESC_IND_W(BDD%NDESC_CODE) ; ZWORK(I1:I2) = BDD%XDESC_FRAC_PANEL(:) 
!
 CALL UP_DESC_IND_W(BDD%NDESC_USE) ; ZWORK(I1:I2) = FLOAT(BDD%NDESC_USE_LIST(:))
 CALL UP_DESC_IND_W(BDD%NDESC_USE) ; ZWORK(I1:I2) = BDD%XDESC_TCOOL_TARGET(:)
 CALL UP_DESC_IND_W(BDD%NDESC_USE) ; ZWORK(I1:I2) = BDD%XDESC_THEAT_TARGET(:)
 CALL UP_DESC_IND_W(BDD%NDESC_USE) ; ZWORK(I1:I2) = BDD%XDESC_QIN(:)
 CALL UP_DESC_IND_W(BDD%NDESC_USE) ; ZWORK(I1:I2) = BDD%XDESC_QIN_FLAT(:)
 CALL UP_DESC_IND_W(BDD%NDESC_USE) ; ZWORK(I1:I2) = BDD%XDESC_SHGC_SH(:)
 CALL UP_DESC_IND_W(BDD%NDESC_USE) ; ZWORK(I1:I2) = BDD%XDESC_SHADE(:)
 CALL UP_DESC_IND_W(BDD%NDESC_USE) ; ZWORK(I1:I2) = BDD%XDESC_NATVENT(:)
 CALL UP_DESC_IND_W(BDD%NDESC_USE) ; ZWORK(I1:I2) = BDD%XDESC_RESIDENTIAL(:)
!
 CALL UP_DESC_IND_W(BDD%NDESC_AGE) ; ZWORK(I1:I2) = FLOAT(BDD%NDESC_AGE_LIST(:))
 CALL UP_DESC_IND_W(BDD%NDESC_AGE) ; ZWORK(I1:I2) = FLOAT(BDD%NDESC_AGE_DATE(:))
!
YCOMMENT='Descriptive building data'
 CALL WRITE_SURF(HSELECT, &
                 HPROGRAM,'BLD_DESC_DAT',ZWORK,IRESP,YCOMMENT,'-','Bld_parameters  ')
DEALLOCATE(ZWORK)
!
IF (LHOOK) CALL DR_HOOK('WRITE_BLD_DESCRIPTION_n',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
CONTAINS
SUBROUTINE UP_DESC_IND_W(K)
INTEGER, INTENT(IN) :: K
I1=I2+1
I2=I2+K
END SUBROUTINE UP_DESC_IND_W
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_BLD_DESCRIPTION_n
