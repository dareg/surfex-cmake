!     #########
      SUBROUTINE AVERAGE1_OROGRAPHY (UG,USS, &
                                     KLUOUT,KNBLINES,PLAT,PLON,PVALUE,PNODATA)
!     #######################################################
!
!!**** *AVERAGE1_OROGRAPHY* computes the sum of orography, squared orography
!!                              and subgrid orography characteristics
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson         Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    12/09/95
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SSO_n, ONLY : SSO_t
!
USE MODD_PGDWORK,       ONLY : X2D_ALL, XEXTVAL2, NSIZE_ALL, X3D_ALL, N3D_ALL, NSSO
!
USE MODI_GET_MESH_INDEX
USE MODD_POINT_OVERLAY, ONLY : NOVMX
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
!
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SSO_t), INTENT(INOUT) :: USS
!
INTEGER,                 INTENT(IN)    :: KLUOUT
INTEGER,                 INTENT(IN)    :: KNBLINES
REAL, DIMENSION(:),      INTENT(IN)    :: PLAT    ! latitude of the point to add
REAL, DIMENSION(:),      INTENT(IN)    :: PLON    ! longitude of the point to add
REAL, DIMENSION(:),      INTENT(IN)    :: PVALUE  ! value of the point to add
REAL, OPTIONAL, INTENT(IN) :: PNODATA
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
INTEGER, DIMENSION(NOVMX,SIZE(PLAT)) :: IINDEX ! mesh index of all input points
                                         ! 0 indicates the point is out of the domain
INTEGER, DIMENSION(NOVMX,SIZE(PLAT)) :: ISSOX  ! X submesh index in their mesh of all input points
INTEGER, DIMENSION(NOVMX,SIZE(PLAT)) :: ISSOY  ! Y submesh index in their mesh of all input points
!
INTEGER :: JLOOP, JOVER        ! loop index on input arrays
REAL, DIMENSION(SIZE(PLAT)) :: ZVALUE
REAL :: ZNODATA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
!
!*    1.     Get position
!            ------------
!     
IF (LHOOK) CALL DR_HOOK('AVERAGE1_OROGRAPHY',0,ZHOOK_HANDLE)
!      
IF (PRESENT(PNODATA)) THEN
  ZVALUE(:) = PVALUE(:)
  ZNODATA = PNODATA
  CALL GET_MESH_INDEX(UG,KLUOUT,KNBLINES,PLAT,PLON,IINDEX,ZVALUE,ZNODATA,NSSO,ISSOX,ISSOY)
ELSE
  ZVALUE(:) = 1.
  ZNODATA = 0.
  CALL GET_MESH_INDEX(UG,KLUOUT,KNBLINES,PLAT,PLON,IINDEX,KSSO=NSSO,KISSOX=ISSOX,KISSOY=ISSOY)
ENDIF
!
!*    2.     Loop on all input data points
!            -----------------------------
!    
bloop: &
DO JLOOP = 1 , SIZE(PLAT)
!
  DO JOVER = 1, NOVMX
!
!*    3.     Tests on position
!            -----------------
!     
    IF (IINDEX(JOVER,JLOOP)==0) CYCLE bloop
!
!*    4.     Summation
!            ---------
!
    NSIZE_ALL(IINDEX(JOVER,JLOOP))=NSIZE_ALL(IINDEX(JOVER,JLOOP))+1
!
!*    5.     Orography
!            ---------
!
    X2D_ALL(IINDEX(JOVER,JLOOP),1)=X2D_ALL(IINDEX(JOVER,JLOOP),1)+PVALUE(JLOOP)
!
!*    6.     Square of Orography
!            -------------------
!
    X2D_ALL(IINDEX(JOVER,JLOOP),2)=X2D_ALL(IINDEX(JOVER,JLOOP),2)+PVALUE(JLOOP)**2
!
!*    7.     Maximum orography in a subgrid square
!            -------------------------------------
!
    N3D_ALL(IINDEX(JOVER,JLOOP),ISSOX(JOVER,JLOOP),ISSOY(JOVER,JLOOP)) = 1
    X3D_ALL(IINDEX(JOVER,JLOOP),ISSOX(JOVER,JLOOP),ISSOY(JOVER,JLOOP)) = &
         MAX (  X3D_ALL(IINDEX(JOVER,JLOOP),ISSOX(JOVER,JLOOP),ISSOY(JOVER,JLOOP)) , PVALUE(JLOOP) )   
!
!
!*    8.     Maximum orography in the mesh
!            -----------------------------
!
    XEXTVAL2(IINDEX(JOVER,JLOOP),1)=MAX(XEXTVAL2(IINDEX(JOVER,JLOOP),1),PVALUE(JLOOP))
!
!
!*    9.     Minimum orography in the mesh
!            -----------------------------
!
    XEXTVAL2(IINDEX(JOVER,JLOOP),2)=MIN(XEXTVAL2(IINDEX(JOVER,JLOOP),2),PVALUE(JLOOP))
!
!
  END DO
!
ENDDO bloop
IF (LHOOK) CALL DR_HOOK('AVERAGE1_OROGRAPHY',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE1_OROGRAPHY
