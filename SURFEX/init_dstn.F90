SUBROUTINE INIT_DST_n(HPROGRAM  &! Program calling unit
       )  

USE MODD_DST_n
USE MODD_DST_SURF
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
USE MODI_GET_VEGTYPE_2_PATCH_MASK
!
USE MODI_GET_LUOUT
IMPLICIT NONE

!PASSED VARIABLES
CHARACTER(LEN=6)      :: HPROGRAM              !Passing unit

!LOCAL VARIABLES
INTEGER             :: JVEG                  ! Counter for vegetation classes
INTEGER             :: JVEG_IN               ! Vegetation index
INTEGER             :: JPATCH                ! Counter for patches
INTEGER             :: JMODE                 ! Counter for dust modes
INTEGER             :: JMODE_IDX             ! Index for dust modes
INTEGER             :: ILUOUT
INTEGER             :: IMI                   ! Current model
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!get output listing unit
IF (LHOOK) CALL DR_HOOK('INIT_DST_N',0,ZHOOK_HANDLE)
CALL GET_LUOUT(HPROGRAM,ILUOUT)

!Allocate memory for the real values which will be used by the model
ALLOCATE(XMSS_FRC_SRC(NDSTMDE))
ALLOCATE(XEMISRADIUS(NDSTMDE))
ALLOCATE(XEMISSIG(NDSTMDE))

!Get initial size distributions. This is cut and pasted
!from dead routine dstpsd.F90
!Check for different source parameterizations
IF(CEMISPARAM.eq."She84")THEN
   CRGUNITD   = 'MASS'
   XMSS_FRC_SRC_INI(:)=(/2.6e-6,0.781,0.219/)             ! [frc] Mass fraction She84 p. 75 Table 1
   XEMISRADIUS_INI(:) = 0.5d6* (/ 0.0111e-6,  2.524e-6, 42.10e-6 /)  ! [um]  Mass median radius She84 p. 75 Table 1
   XEMISSIG_INI(:) = (/ 1.89     ,  2.0     ,  2.13    /)          ! [frc] Geometric standard deviation She84 p. 75 Table 1
ELSEIF(CEMISPARAM.eq."PaG77")THEN
   CRGUNITD   = 'MASS'
   XMSS_FRC_SRC_INI(:)=(/0.036,0.957,0.007/)              ! [frc] Mass fraction BSM96 p. 73 Table 2 (ad hoc)
   XEMISRADIUS_INI(:) = 0.5d6*(/ 0.27e-6  ,  5.6e-6  ,  57.6e-6 /)   ! [um] Mass median radius PaG77 p. 2080 Table 1 
   XEMISSIG_INI(:) = (/ 1.88     ,  2.2     ,  1.62    /)          ! [frc] Geometric standard deviation PaG77 p. 2080 Table 1
ELSEIF(CEMISPARAM.eq."Dal87") THEN 
! D'Almeida, 1987 as default
   CRGUNITD   = 'MASS'
   XMSS_FRC_SRC_INI(:)=(/0.036,0.957,0.007/)              ! [frc] Mass fraction BSM96 p. 73 Table 2
   XEMISRADIUS_INI(:)=0.5d6*(/ 0.832e-6 ,  4.82e-6 , 19.38e-6 /)     ! [um] Mass median radius BSM96 p. 73 Table 2
   XEMISSIG_INI(:)=(/ 2.10     ,  1.90    ,  1.60    /)            ! [frc] Geometric standard deviation BSM96 p. 73 Table 2
ELSEIF ((CEMISPARAM.eq."alf98").OR.(CEMISPARAM.eq."EXPLI")) THEN  !  Alfaro et al 1998 as default
  IF (CEMISPARAM.eq."alf98") XFLX_MSS_FDG_FCT = 6e-4
  IF (CEMISPARAM.eq."EXPLI") XFLX_MSS_FDG_FCT = 3.5e-4
  CRGUNITD   = 'MASS'
  XMSS_FRC_SRC_INI(:)=(/0.01,0.19,0.8/)  ! [frc] Mass fraction BSM96 p. 73 Table 2
  XEMISRADIUS_INI(:)=0.5*(/ 1.5, 6.7, 14.2 /)     ! [um] Mass median radius BSM96 p. 73 Table 2
  XEMISSIG_INI(:)=(/ 1.70, 1.60, 1.50 /)       ! [frc] Geometric standard deviation BSM96 p. 73 Table 2
ELSEIF (CEMISPARAM.eq."AMMA") THEN ! Default: New distribution from AMMA
  XFLX_MSS_FDG_FCT = 60.e-4
  CRGUNITD   = 'NUMB'
  XMSS_FRC_SRC_INI(:)=(/0.008,0.092,0.99/)         ! [frc] Mass fraction 
  XEMISRADIUS_INI(:)=0.5*(/ 0.078, 0.641, 5.00 /)  ! [um] Number median radius 
  XEMISSIG_INI(:)=(/ 1.75, 1.76, 1.70 /)            ! [frc] Geometric standard deviation 
ELSEIF (CEMISPARAM.eq."CRUM") THEN ! Default: New distribution from AMMA
  XFLX_MSS_FDG_FCT = 10.e-4
  CRGUNITD   = 'NUMB'
  XMSS_FRC_SRC_INI(:)=(/0.0005,0.0029,0.9966/)     ! [frc] Mass fraction 
  XEMISRADIUS_INI(:)=0.5*(/ 0.078, 0.641, 5.00 /)  ! [um] Number median radius 
  XEMISSIG_INI(:)=(/ 1.75, 1.76, 1.70 /)           ! [frc] Geometric standard deviation   
ELSE
  WRITE(*,*) " FATAL ERROR "
  WRITE(*,*) " YOU MUST DECIDE THE EMISSIUON PARAMETERIZATION, YOU USES "
  WRITE(*,*) " CEMISPARAM = ",CEMISPARAM," AND IT IS NOT DEFINED "
  WRITE(*,*) " see init_dstn.f90 to see what dust parameterization is available. "
ENDIF


DO JMODE=1,NDSTMDE
   JMODE_IDX=JORDER(JMODE)
   XMSS_FRC_SRC(JMODE) = XMSS_FRC_SRC_INI(JMODE_IDX)
   XEMISSIG(JMODE)=XEMISSIG_INI(JMODE_IDX)
   !Get emisradius, and at the same time convert to number median radius
   IF (CRGUNITD=='MASS') THEN
   XEMISRADIUS(JMODE) =              &
          XEMISRADIUS_INI(JMODE_IDX)    &
          * EXP(-3.d0 * (LOG(XEMISSIG_INI(JMODE_IDX)))**2)  
   ELSE
   XEMISRADIUS(JMODE) = XEMISRADIUS_INI(JMODE_IDX)
   ENDIF
ENDDO
!Normalize the sum of the emissions to 1 so that all dust is
!put in one mode or the other
IF(SUM(XMSS_FRC_SRC(:)).lt.1.)THEN
   XMSS_FRC_SRC(:)=XMSS_FRC_SRC(:)/SUM(XMSS_FRC_SRC(:))
ENDIF

!Allocate memory
!ALLOCATE(NVEGNO_DST)
!Set the number of classes that can emit dust (fxm: set this elsewhere)
NVEGNO_DST = 2

!Allocate memory for the vegtype-translator
ALLOCATE(NVT_DST(NVEGNO_DST))

!Set the dust/vegtype translator vector
NVT_DST(1)  = NVT_NO
NVT_DST(2)  = NVT_ROCK

!Allocate memory for roughness lengths of erodible surfaces
ALLOCATE(Z0_EROD_DST(NVEGNO_DST))

!Set the roughness lengths corresponding to erodible surfaces
!Smooth roughness length is given to 1.d-5 (dstmbl.f90)
Z0_EROD_DST(1) = 30.d-6    !m (30 um) 
Z0_EROD_DST(2) = 200.d-6   !m (200 um) 

!Allocate memory for dust emitter surface vectors in patch vectors
IF (.NOT.ASSOCIATED(NSIZE_PATCH_DST)) ALLOCATE(NSIZE_PATCH_DST(NVEGNO_DST,NPATCH))

DO JPATCH = 1,NPATCH
   DO JVEG = 1,NVEGNO_DST
      !Count all the points in the patch where you have dust emitter vegetation
      NSIZE_PATCH_DST(JVEG,JPATCH) = COUNT(XVEGTYPE_PATCH(:,NVT_DST(JVEG),JPATCH) > 0.) 
   ENDDO
ENDDO

!Find the largest dust emitter vector in any patch
!ALLOCATE (NSIZE_LARGEST_DST)
NSIZE_LARGEST_DST = 0
DO JPATCH=1,NPATCH
   DO JVEG = 1,NVEGNO_DST
      NSIZE_LARGEST_DST = max(NSIZE_LARGEST_DST,NSIZE_PATCH_DST(JVEG,JPATCH))
   ENDDO
ENDDO

!Allocate memory for NR_PATCH_DST mask translate from patch vector to dust vector
ALLOCATE(NR_PATCH_DST(NSIZE_LARGEST_DST,NVEGNO_DST,NPATCH))

!Initialize the mask array
NR_PATCH_DST(:,:,:)=0

!Get values from the dust emitter vegetation mask
DO JPATCH=1,NPATCH
   DO JVEG=1,NVEGNO_DST
      JVEG_IN = NVT_DST(JVEG)          ! Get the real vegtype index
      CALL GET_VEGTYPE_2_PATCH_MASK(ILUOUT,&
             NSIZE_PATCH_DST(JVEG,JPATCH),             &!I Size of dust emitter vector
             NSIZE_NATURE_P(JPATCH),                   &!I Size of patch vector
             NSIZE_NATURE,                             &!I Size of nature vector
             NR_NATURE_P,                              &!I Mask from patch to nature
             XVEGTYPE_PATCH,                           &!I Fraction of vegtype of nature point within jpatch 
             NR_PATCH_DST(:NSIZE_PATCH_DST(JVEG,JPATCH),JVEG,JPATCH),  &!O Part of mask array to fill with values
             NPATCH,                                   &!I Number of possible patches
             JPATCH,                                   &!I Index of patch in question
             JVEG_IN                                  &!I Index of vegtype in question
             )  
   ENDDO !Loop on patches
ENDDO    !Loop on veg-types
IF (LHOOK) CALL DR_HOOK('INIT_DST_N',1,ZHOOK_HANDLE)

END SUBROUTINE INIT_DST_n

