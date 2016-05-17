!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE AVERAGE1_COVER(UG,U,KLUOUT,KNBLINES,PLAT,PLON,PVALUE,PNODATA)
!     #######################################################
!
!!**** *AVERAGE1_COVER* computes the sum of values of a cover fractions
!!                              and the nature of terrain on the grid
!!                              from a data in land-cover file
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
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_SURFEX_MPI, ONLY : NRANK
USE MODD_PGDWORK, ONLY : X3D_ALL, NSIZE_ALL
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
!
USE MODI_GET_MESH_INDEX
USE MODD_POINT_OVERLAY
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
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
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
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: Z3D_ALL
REAL, DIMENSION(SIZE(PLAT)) :: ZVALUE
REAL :: ZNODATA
INTEGER :: JLOOP, JOVER, JCOV, IFOUND, ICOV, IND        ! loop index on input arrays
INTEGER :: ICOVERCLASS  ! class of cover type
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
!
!*    1.     Get position
!            ------------
!     
IF (LHOOK) CALL DR_HOOK('AVERAGE1_COVER',0,ZHOOK_HANDLE)
!
ICOV = SIZE(X3D_ALL,2)
!
IF (PRESENT(PNODATA)) THEN
  ZVALUE(:) = PVALUE(:)
  ZNODATA = PNODATA
  CALL GET_MESH_INDEX(UG,KLUOUT,KNBLINES,PLAT,PLON,IINDEX,ZVALUE,ZNODATA)
ELSE
  ZVALUE(:) = 1.
  ZNODATA = 0.
  CALL GET_MESH_INDEX(UG,KLUOUT,KNBLINES,PLAT,PLON,IINDEX)
ENDIF
!
!*    2.     Loop on all input data points
!            -----------------------------
!     
bloop: &
DO JLOOP = 1 , SIZE(PLAT)
!
!*    3.     Tests on position
!            -----------------
!    
  DO JOVER = 1, NOVMX

    IF (IINDEX(JOVER,JLOOP)==0) CYCLE bloop
!
!*    4.     Test on value meaning
!            ---------------------
!
    ICOVERCLASS = NINT(PVALUE(JLOOP))
!
    U%LCOVER(ICOVERCLASS) = .TRUE.
!
    IF (ICOVERCLASS<1 .OR. ICOVERCLASS > JPCOVER )  CYCLE
!
!*    5.     Summation
!            ---------
!
    NSIZE_ALL(IINDEX(JOVER,JLOOP))=NSIZE_ALL(IINDEX(JOVER,JLOOP))+1
!
!*    6.     Fraction of cover type
!            ----------------------
!
    IFOUND = 0
    !ICOV: number of covers already found in the domain
    DO JCOV=1,ICOV
      !if the cover read is already in the array
      IF (X3D_ALL(IINDEX(JOVER,JLOOP),JCOV,1)==ICOVERCLASS*1.) THEN
        !the number of points found is increased of 1
        X3D_ALL(IINDEX(JOVER,JLOOP),JCOV,2)=X3D_ALL(IINDEX(JOVER,JLOOP),JCOV,2)+1.
        IFOUND=1
        EXIT
      ENDIF
    ENDDO
    !if the cover is not in the array 
    IF (IFOUND==0) THEN
      !if we already have some covers for this point
      IF (X3D_ALL(IINDEX(JOVER,JLOOP),ICOV,2)/=0.) THEN
        !to save the current array
        ALLOCATE(Z3D_ALL(SIZE(X3D_ALL,1),ICOV,SIZE(X3D_ALL,3)))
        Z3D_ALL(:,:,:) = X3D_ALL(:,:,:)
        DEALLOCATE(X3D_ALL)
        !we add one cover to the size of the array
        ALLOCATE(X3D_ALL(SIZE(Z3D_ALL,1),ICOV+1,SIZE(Z3D_ALL,3)))
        X3D_ALL(:,1:ICOV,:) = Z3D_ALL(:,:,:)
        DEALLOCATE(Z3D_ALL)
        X3D_ALL(:,ICOV+1,:) = 0.
        !the number of covers already found increases
        ICOV = ICOV + 1
      ENDIF
      !first index for this point where no cover is defined
      IND = MINLOC(X3D_ALL(IINDEX(JOVER,JLOOP),:,2),1,X3D_ALL(IINDEX(JOVER,JLOOP),:,2)==0.)
      !the new cover is registered
      X3D_ALL(IINDEX(JOVER,JLOOP),IND,1) = ICOVERCLASS*1.
      X3D_ALL(IINDEX(JOVER,JLOOP),IND,2) = 1.
    ENDIF
    !
  ENDDO
!
END DO bloop
!
IF (LHOOK) CALL DR_HOOK('AVERAGE1_COVER',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE1_COVER
