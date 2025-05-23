!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #############################
      SUBROUTINE AVERAGE_DIAG_MISC_ISBA_n (DM, NDM, IO, NP, NPE)
!     #############################
!
!
!!****  *AVERAGE_DIAG_MISC_ISBA_n*  
!!
!!    PURPOSE
!!    -------
!      Average the cumulated diagnostics from all ISBA tiles
!     
!!**  METHOD
!!    ------
!
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
!!      P. Le Moigne           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2004
!!      B. Decharme  2008    New diag Total albedo, Total SWI, & Flood
!!      B. Decharme 09/2009  New diag Total soil SWI
!!      B. Decharme  2012    Averaged LAI
!!      B. Decharme  2012    New diag for DIF:
!!                           F2 stress
!!                           Root zone swi, wg and wgi
!!                           swi, wg and wgi comparable to ISBA-FR-DG2 and DG3 layers
!!                           active layer thickness over permafrost
!!                           frozen layer thickness over non-permafrost
!!      F. Wang     07/2020  New diag XWR, XSNOWALB, XSNOWRHO, XWGI, XTG,
!!                           XFRD2_SWI, XFRD3_SWI, XFRD2_TG, XFRD3_TG
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t, DIAG_MISC_ISBA_NP_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_P_t, ISBA_PE_t, ISBA_NP_t, ISBA_NPE_t
!
USE MODD_SURF_PAR,         ONLY : XUNDEF, NUNDEF
!
USE MODD_CSTS,             ONLY : XRHOLW
!
USE MODI_COMPUT_COLD_LAYERS_THICK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DM
TYPE(DIAG_MISC_ISBA_NP_t), INTENT(INOUT) :: NDM
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t), INTENT(INOUT) :: NPE
!
!*      0.2    declarations of local variables
!
TYPE(DIAG_MISC_ISBA_t), POINTER :: DMK
TYPE(ISBA_P_t), POINTER :: PK
TYPE(ISBA_PE_t), POINTER :: PEK
INTEGER                         :: JI    ! grid-cell loop counter
INTEGER                         :: JP    ! tile loop counter
INTEGER                         :: JL    ! layer loop counter
REAL, DIMENSION(SIZE(DM%XHV))   :: ZSUMDG, ZSNOW, ZSUMFRD2, ZSUMFRD3
REAL, DIMENSION(SIZE(DM%XHV))   :: ZSUMFRD2_TG, ZSUMFRD3_TG  ! for soil temperature
REAL                            :: ZWORK
INTEGER                         :: INI,INP,IDEPTH,IWORK, IMASK, JB
!
REAL, DIMENSION(SIZE(DM%XHV),IO%NGROUND_LAYER) :: ZPOND, ZTG, ZDG
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!       0.     Initialization
!              --------------
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_MISC_ISBA_N',0,ZHOOK_HANDLE)
!
IF (.NOT.DM%LSURF_MISC_BUDGET) THEN
   IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_MISC_ISBA_N',1,ZHOOK_HANDLE)
   RETURN
ENDIF
!
INI=SIZE(DM%XHV)
!
!-------------------------------------------------------------------------------
!
!       1.     Surface Miscellaneous terms
!              ---------------------------
!
DM%XHV  (:)    = 0.
DM%XPSNG(:)    = 0.
DM%XPSNV(:)    = 0.
DM%XPSN (:)    = 0.
DM%XFSAT(:)    = 0.
DM%XFFG (:)    = 0.
DM%XFFV (:)    = 0.
DM%XFF  (:)    = 0.
DM%XLAI (:)    = 0.
DM%XWR  (:)    = 0.
DM%XSNOWALB(:) = 0.
DM%XSNOWRHO(:) = 0.
DM%XTWSNOW (:) = 0.
DM%XTDSNOW (:) = 0.
DM%XTTSNOW (:) = 0.
IF (DM%LPROSNOW .AND. NPE%AL(1)%TSNOW%SCHEME=="CRO") THEN
  DM%XSNDPT_1DY     (:) = 0.
  DM%XSNDPT_3DY     (:) = 0.
  DM%XSNDPT_5DY     (:) = 0.
  DM%XSNDPT_7DY     (:) = 0.
  DM%XSNSWE_1DY     (:) = 0.
  DM%XSNSWE_3DY     (:) = 0.
  DM%XSNSWE_5DY     (:) = 0.
  DM%XSNSWE_7DY     (:) = 0.
  DM%XSNRAM_SONDE   (:) = 0.
  DM%XSN_WETTHCKN   (:) = 0.
  DM%XSN_REFRZNTHCKN(:) = 0.
  DM%XDEP_HIG       (:) = 0.
  DM%XDEP_MOD       (:) = 0.
  DM%XACC_LEV       (:) = 0.
  DM%XPRO_INF_TYP   (:) = 0.
 IF (DM%LPROBANDS) THEN
  DM%XSPEC_ALB(:,:)=0.
  DM%XDIFF_RATIO(:,:)=0.
 ENDIF 
ENDIF

ZSNOW   (:)=0.0
!   
DO JP=1,IO%NPATCH
  PK => NP%AL(JP)
  DMK => NDM%AL(JP)
  PEK => NPE%AL(JP)
!   
  DO JI=1,PK%NSIZE_P
    IMASK = PK%NR_P(JI)
!
!     Halstead coefficient
    DM%XHV  (IMASK) = DM%XHV   (IMASK) + PK%XPATCH(JI) * DMK%XHV(JI)
!
!     Snow fractions
    DM%XPSNG(IMASK) = DM%XPSNG (IMASK) + PK%XPATCH(JI) * DMK%XPSNG(JI)
    DM%XPSNV(IMASK) = DM%XPSNV (IMASK) + PK%XPATCH(JI) * DMK%XPSNV(JI)
    DM%XPSN (IMASK) = DM%XPSN  (IMASK) + PK%XPATCH(JI) * DMK%XPSN (JI)
!
!     Saturated fraction
    DM%XFSAT (IMASK) = DM%XFSAT(IMASK) + PK%XPATCH(JI) * DMK%XFSAT(JI)
!
!     Flood fractions
    DM%XFFG  (IMASK) = DM%XFFG (IMASK) + PK%XPATCH(JI) * DMK%XFFG(JI)
    DM%XFFV  (IMASK) = DM%XFFV (IMASK) + PK%XPATCH(JI) * DMK%XFFV(JI)
    DM%XFF   (IMASK) = DM%XFF  (IMASK) + PK%XPATCH(JI) * DMK%XFF (JI)
!
!     Total LAI
    IF (PEK%XLAI(JI)/=XUNDEF) DM%XLAI(IMASK) = DM%XLAI(IMASK) + PK%XPATCH(JI) * PEK%XLAI(JI)
    !
    !     Total WR
    IF (PEK%XWR(JI)/=XUNDEF) DM%XWR(IMASK) = DM%XWR(IMASK) + PK%XPATCH(JI) * PEK%XWR(JI)
    !      
    !     Snow albedo
    DM%XSNOWALB(IMASK) = DM%XSNOWALB(IMASK) + PK%XPATCH(JI) * PEK%TSNOW%ALB(JI)
    !
    !     Snow total outputs
    DM%XTWSNOW(IMASK) = DM%XTWSNOW(IMASK) + PK%XPATCH(JI) * DMK%XTWSNOW(JI)
    DM%XTDSNOW(IMASK) = DM%XTDSNOW(IMASK) + PK%XPATCH(JI) * DMK%XTDSNOW(JI)
    IF (DM%XTDSNOW(IMASK) /= XUNDEF .AND. DM%XTDSNOW(IMASK) > 0.0) THEN
      DM%XSNOWRHO(IMASK) = DM%XTWSNOW(IMASK) / DM%XTDSNOW(IMASK)
    ENDIF
    !      
    IF (DMK%XTWSNOW(JI)>0.0) THEN
      !
      DM%XTTSNOW(IMASK) = DM%XTTSNOW(IMASK) + PK%XPATCH(JI) * DMK%XTTSNOW(JI)
      ZSNOW      (IMASK) = ZSNOW    (IMASK) + PK%XPATCH(JI)
      !
      IF (DM%LPROSNOW .AND. NPE%AL(1)%TSNOW%SCHEME=="CRO") THEN
        !
        DM%XSNDPT_1DY     (IMASK) = DM%XSNDPT_1DY     (IMASK) + PK%XPATCH(JI) * DMK%XSNDPT_1DY     (JI)
        DM%XSNDPT_3DY     (IMASK) = DM%XSNDPT_3DY     (IMASK) + PK%XPATCH(JI) * DMK%XSNDPT_3DY     (JI)
        DM%XSNDPT_5DY     (IMASK) = DM%XSNDPT_5DY     (IMASK) + PK%XPATCH(JI) * DMK%XSNDPT_5DY     (JI)
        DM%XSNDPT_7DY     (IMASK) = DM%XSNDPT_7DY     (IMASK) + PK%XPATCH(JI) * DMK%XSNDPT_7DY     (JI)
        DM%XSNSWE_1DY     (IMASK) = DM%XSNSWE_1DY(     IMASK) + PK%XPATCH(JI) * DMK%XSNSWE_1DY     (JI)
        DM%XSNSWE_3DY     (IMASK) = DM%XSNSWE_3DY     (IMASK) + PK%XPATCH(JI) * DMK%XSNSWE_3DY     (JI)
        DM%XSNSWE_5DY     (IMASK) = DM%XSNSWE_5DY     (IMASK) + PK%XPATCH(JI) * DMK%XSNSWE_5DY     (JI)
        DM%XSNSWE_7DY     (IMASK) = DM%XSNSWE_7DY     (IMASK) + PK%XPATCH(JI) * DMK%XSNSWE_7DY     (JI)
        DM%XSNRAM_SONDE   (IMASK) = DM%XSNRAM_SONDE   (IMASK) + PK%XPATCH(JI) * DMK%XSNRAM_SONDE   (JI)
        DM%XSN_WETTHCKN   (IMASK) = DM%XSN_WETTHCKN   (IMASK) + PK%XPATCH(JI) * DMK%XSN_WETTHCKN   (JI)
        DM%XSN_REFRZNTHCKN(IMASK) = DM%XSN_REFRZNTHCKN(IMASK) + PK%XPATCH(JI) * DMK%XSN_REFRZNTHCKN(JI)
        !The following does not make really sense
        DM%XDEP_HIG       (IMASK) = DM%XDEP_HIG       (IMASK) + PK%XPATCH(JI) * DMK%XDEP_HIG       (JI)
        DM%XDEP_MOD       (IMASK) = DM%XDEP_MOD       (IMASK) + PK%XPATCH(JI) * DMK%XDEP_MOD       (JI)
        DM%XACC_LEV       (IMASK) = DM%XACC_LEV       (IMASK) + PK%XPATCH(JI) * DMK%XACC_LEV       (JI)
        DM%XPRO_INF_TYP   (IMASK) = DM%XPRO_INF_TYP   (IMASK) + PK%XPATCH(JI) * DMK%XPRO_INF_TYP   (JI)
        !
         
      ENDIF
!
    ENDIF
!
    ! Ugly way to keep XUNDEF instead of 0 in case of no snow
    IF ((DMK%XTWSNOW(JI)<=0.0).AND.(DM%LPROSNOW .AND. NPE%AL(1)%TSNOW%SCHEME=="CRO")) THEN
      DM%XDEP_HIG       (IMASK) = XUNDEF
      DM%XDEP_MOD       (IMASK) = XUNDEF
      DM%XACC_LEV       (IMASK) = 4
      DM%XPRO_INF_TYP   (IMASK) = 6
    ENDIF
  ENDDO
!
ENDDO
!
!
!-------------------------------------------------------------------------------
!
!       2.     Specific treatement following CISBA option
!              ------------------------------------------
!
!   Soil Wetness Index profile, Total Soil Wetness Index and 
!   Total Soil Water Content (Liquid+Solid) and Total Frozen Content
!
DM%XSWI (:,:) = 0.
DM%XTSWI(:,:) = 0.
!   
DM%XWGI (:,:) = 0.
DM%XTG  (:,:) = 0.
!
DM%XSOIL_SWI  (:) = 0.
DM%XSOIL_TSWI (:) = 0.
DM%XSOIL_TWG  (:) = 0.
DM%XSOIL_TWGI (:) = 0.
DM%XSOIL_WG   (:) = 0.
DM%XSOIL_WGI  (:) = 0.
!
ZSUMDG  (:)=0.0
!
!---------------------------------------------
IF(IO%CISBA=='DIF')THEN ! DIF case
!---------------------------------------------
!
  DM%XALT   (:) = 0. 
  DM%XFLT   (:) = 0. 

! Active and Frozen layers thickness
  ZTG(:,:)=0.0
  ZDG(:,:)=0.0
  DO JP=1,IO%NPATCH
    PK => NP%AL(JP)
    PEK => NPE%AL(JP)

    DO JL=1,IO%NGROUND_LAYER
      DO JI=1,PK%NSIZE_P
        IMASK = PK%NR_P(JI)
        DM%XTG(IMASK,JL) = DM%XTG(IMASK,JL) + PK%XPATCH(JI) * PEK%XTG(JI,JL)
        ZTG(IMASK,JL) = ZTG(IMASK,JL) + PK%XPATCH(JI) * PEK%XTG(JI,JL)
        ZDG(IMASK,JL) = ZDG(IMASK,JL) + PK%XPATCH(JI) * PK%XDG (JI,JL)
        ENDDO
     ENDDO

  ENDDO
  CALL COMPUT_COLD_LAYERS_THICK(ZDG,ZTG,DM%XALT,DM%XFLT)
!    
  ZPOND(:,:)=0.0
  DO JP=1,IO%NPATCH   
    DMK => NDM%AL(JP)
    PK => NP%AL(JP)
    PEK => NPE%AL(JP)

    DO JL = 1,IO%NGROUND_LAYER

      DO JI=1,PK%NSIZE_P
        IDEPTH = PK%NWG_LAYER(JI)
        IF(JL<=IDEPTH.AND.IDEPTH/=NUNDEF)THEN

          IMASK = PK%NR_P(JI)

          ZWORK = PK%XDZG(JI,JL)
          !Soil Ice Content
          DM%XWGI (IMASK,JL) = DM%XWGI (IMASK,JL) + ZWORK*PK%XPATCH(JI) * PEK%XWGI (JI,JL)
          !Soil Wetness Index profile
          DM%XSWI (IMASK,JL) = DM%XSWI (IMASK,JL) + ZWORK*PK%XPATCH(JI) * DMK%XSWI (JI,JL) 
          DM%XTSWI(IMASK,JL) = DM%XTSWI(IMASK,JL) + ZWORK*PK%XPATCH(JI) * DMK%XTSWI(JI,JL)
          ZPOND   (IMASK,JL) = ZPOND   (IMASK,JL) + ZWORK*PK%XPATCH(JI)
               !Total soil wetness index, total water and ice contents
          DM%XSOIL_SWI (IMASK) = DM%XSOIL_SWI (IMASK) + ZWORK * PK%XPATCH(JI) * DMK%XSWI (JI,JL)
          DM%XSOIL_TSWI(IMASK) = DM%XSOIL_TSWI(IMASK) + ZWORK * PK%XPATCH(JI) * DMK%XTSWI(JI,JL)
          ZSUMDG       (IMASK) = ZSUMDG       (IMASK) + ZWORK * PK%XPATCH(JI)
          DM%XSOIL_TWG (IMASK) = DM%XSOIL_TWG (IMASK) + ZWORK * PK%XPATCH(JI) * (PEK%XWG(JI,JL) + PEK%XWGI(JI,JL))
          DM%XSOIL_TWGI(IMASK) = DM%XSOIL_TWGI(IMASK) + ZWORK * PK%XPATCH(JI) * PEK%XWGI(JI,JL)

             ENDIF

          ENDDO

       ENDDO
    !
  ENDDO
!  
  WHERE(ZPOND(:,:)> 0.)
    DM%XWGI (:,:) = DM%XWGI (:,:) / ZPOND(:,:)
    DM%XSWI (:,:) = DM%XSWI (:,:) / ZPOND(:,:)
    DM%XTSWI(:,:) = DM%XTSWI(:,:) / ZPOND(:,:)
  ELSEWHERE
    DM%XWGI (:,:) = XUNDEF
    DM%XSWI (:,:) = XUNDEF
    DM%XTSWI(:,:) = XUNDEF
  ENDWHERE
!
! ---------------------------------------------
  IF(DM%LSURF_MISC_DIF)THEN ! LSURF_MISC_DIF case
! ---------------------------------------------
!     
    ZSUMFRD2      (:) = 0.
    ZSUMFRD3      (:) = 0.
    ZSUMFRD2_TG   (:) = 0.
    ZSUMFRD3_TG   (:) = 0.
!
    DM%XFRD2_SWI  (:) = 0.
    DM%XFRD2_TG   (:) = 0.
    DM%XFRD2_TSWI (:) = 0.
    DM%XFRD2_TWG  (:) = 0.
    DM%XFRD2_TWGI (:) = 0.
!   
    DM%XFRD3_SWI  (:) = 0.
    DM%XFRD3_TG   (:) = 0.
    DM%XFRD3_TSWI (:) = 0.
    DM%XFRD3_TWG  (:) = 0.
    DM%XFRD3_TWGI (:) = 0.

    DO JP=1,IO%NPATCH
      PK => NP%AL(JP)
      PEK => NPE%AL(JP)
      DMK => NDM%AL(JP)

      DO JI=1,PK%NSIZE_P
        IMASK = PK%NR_P(JI)
     
        DO JL = 1,IO%NGROUND_LAYER
          ZWORK = MIN(PK%XDZG(JI,JL),MAX(0.0,PK%XDG2(JI)-PK%XDG(JI,JL)+PK%XDZG(JI,JL)))
          DM%XFRD2_TG   (IMASK) = DM%XFRD2_TG   (IMASK) + ZWORK * PK%XPATCH(JI) * PEK%XTG  (JI,JL)
          ZSUMFRD2_TG   (IMASK) = ZSUMFRD2_TG   (IMASK) + ZWORK * PK%XPATCH(JI)

          ZWORK = MIN(PK%XDZG(JI,JL),MAX(0.0,PK%XDG(JI,JL)-PK%XDG2(JI)))
          DM%XFRD3_TG   (IMASK) = DM%XFRD3_TG   (IMASK) + ZWORK * PK%XPATCH(JI) * PEK%XTG  (JI,JL)
          ZSUMFRD3_TG   (IMASK) = ZSUMFRD3_TG   (IMASK) + ZWORK * PK%XPATCH(JI)

          IDEPTH= PK%NWG_LAYER(JI)

          IF(JL<=IDEPTH.AND.IDEPTH/=NUNDEF)THEN
            !
            ! ISBA-FR-DG2 comparable soil wetness index, liquid water and ice contents
            ZWORK = MIN(PK%XDZG(JI,JL),MAX(0.0,PK%XDG2(JI)-PK%XDG(JI,JL)+PK%XDZG(JI,JL)))
            DM%XFRD2_SWI  (IMASK) = DM%XFRD2_SWI  (IMASK) + ZWORK * PK%XPATCH(JI) * DMK%XSWI (JI,JL)  ! Liquid
            DM%XFRD2_TSWI (IMASK) = DM%XFRD2_TSWI (IMASK) + ZWORK * PK%XPATCH(JI) * DMK%XTSWI(JI,JL)
            DM%XFRD2_TWG  (IMASK) = DM%XFRD2_TWG  (IMASK) + ZWORK * PK%XPATCH(JI) * PEK%XWG  (JI,JL)
            DM%XFRD2_TWGI (IMASK) = DM%XFRD2_TWGI (IMASK) + ZWORK * PK%XPATCH(JI) * PEK%XWGI (JI,JL)
            ZSUMFRD2      (IMASK) = ZSUMFRD2      (IMASK) + ZWORK * PK%XPATCH(JI)
            !
            ! ISBA-FR-DG3 comparable soil wetness index, liquid water and ice contents
            ZWORK  =MIN(PK%XDZG(JI,JL),MAX(0.0,PK%XDG(JI,JL)-PK%XDG2(JI)))
            DM%XFRD3_SWI  (IMASK) = DM%XFRD3_SWI  (IMASK) + ZWORK * PK%XPATCH(JI) * DMK%XSWI (JI,JL)  ! Liquid
            DM%XFRD3_TSWI (IMASK) = DM%XFRD3_TSWI (IMASK) + ZWORK * PK%XPATCH(JI) * DMK%XTSWI(JI,JL)
            DM%XFRD3_TWG  (IMASK) = DM%XFRD3_TWG  (IMASK) + ZWORK * PK%XPATCH(JI) * PEK%XWG  (JI,JL)
            DM%XFRD3_TWGI (IMASK) = DM%XFRD3_TWGI (IMASK) + ZWORK * PK%XPATCH(JI) * PEK%XWGI (JI,JL)
            ZSUMFRD3      (IMASK) = ZSUMFRD3      (IMASK) + ZWORK * PK%XPATCH(JI)
            !
          ENDIF
        ENDDO
      ENDDO
!
    ENDDO
!    
    WHERE(ZSUMFRD2(:)>0.0)
          DM%XFRD2_SWI  (:) = DM%XFRD2_SWI  (:) / ZSUMFRD2(:)
          DM%XFRD2_TSWI (:) = DM%XFRD2_TSWI (:) / ZSUMFRD2(:)
          DM%XFRD2_TWG  (:) = DM%XFRD2_TWG  (:) / ZSUMFRD2(:)
          DM%XFRD2_TWGI (:) = DM%XFRD2_TWGI (:) / ZSUMFRD2(:)          
    ELSEWHERE
          DM%XFRD2_SWI  (:) = XUNDEF
          DM%XFRD2_TSWI (:) = XUNDEF
    ENDWHERE
    WHERE(ZSUMFRD2_TG(:)>0.0)
          DM%XFRD2_TG   (:) = DM%XFRD2_TG   (:) / ZSUMFRD2_TG(:)
    ELSEWHERE
          DM%XFRD2_TG   (:) = XUNDEF
    ENDWHERE
!    
    WHERE(ZSUMFRD3(:)>0.0)
          DM%XFRD3_TG   (:) = DM%XFRD3_TG   (:) / ZSUMFRD3(:)
          DM%XFRD3_SWI  (:) = DM%XFRD3_SWI  (:) / ZSUMFRD3(:)
          DM%XFRD3_TSWI (:) = DM%XFRD3_TSWI (:) / ZSUMFRD3(:)
          DM%XFRD3_TWG  (:) = DM%XFRD3_TWG  (:) / ZSUMFRD3(:)
          DM%XFRD3_TWGI (:) = DM%XFRD3_TWGI (:) / ZSUMFRD3(:) 
    ELSEWHERE
          DM%XFRD3_SWI  (:) = XUNDEF
          DM%XFRD3_TSWI (:) = XUNDEF
    ENDWHERE
    WHERE(ZSUMFRD3_TG(:)>0.0)
          DM%XFRD3_TG   (:) = DM%XFRD3_TG   (:) / ZSUMFRD3_TG(:)
    ELSEWHERE
          DM%XFRD3_TG   (:) = XUNDEF
    ENDWHERE
!
! ---------------------------------------------
  ENDIF ! End LSURF_MISC_DIF case
! ---------------------------------------------
!
!---------------------------------------------
ELSE ! Force-restore case
!---------------------------------------------
! 
  ZPOND(:,:) = 0.
  DO JP=1,IO%NPATCH
    PK => NP%AL(JP)
    PEK => NPE%AL(JP)
    DMK => NDM%AL(JP)
!
    DO JI=1,PK%NSIZE_P
      IMASK = PK%NR_P(JI)    
!
      ZPOND   (IMASK,1) = ZPOND   (IMASK,1) + PK%XPATCH(JI) * PK%XDG   (JI,1)
      ZPOND   (IMASK,2) = ZPOND   (IMASK,2) + PK%XPATCH(JI) * PK%XDG   (JI,2)
      DM%XWGI (IMASK,1) = DM%XWGI (IMASK,1) + PK%XPATCH(JI) * PK%XDG   (JI,1) * PEK%XWGI (JI,1)
      DM%XWGI (IMASK,2) = DM%XWGI (IMASK,2) + PK%XPATCH(JI) * PK%XDG   (JI,2) * PEK%XWGI (JI,2)
      DM%XTG  (IMASK,1) = DM%XTG  (IMASK,1) + PK%XPATCH(JI) * PEK%XTG  (JI,1)
      DM%XTG  (IMASK,2) = DM%XTG  (IMASK,2) + PK%XPATCH(JI) * PEK%XTG  (JI,2)
      DM%XSWI (IMASK,1) = DM%XSWI (IMASK,1) + PK%XPATCH(JI) * DMK%XSWI (JI,1)
      DM%XSWI (IMASK,2) = DM%XSWI (IMASK,2) + PK%XPATCH(JI) * DMK%XSWI (JI,2)
      DM%XTSWI(IMASK,1) = DM%XTSWI(IMASK,1) + PK%XPATCH(JI) * DMK%XTSWI(JI,1)
      DM%XTSWI(IMASK,2) = DM%XTSWI(IMASK,2) + PK%XPATCH(JI) * DMK%XTSWI(JI,2)
! 
      DM%XSOIL_SWI (IMASK) = DM%XSOIL_SWI (IMASK) + PK%XPATCH(JI) * PK%XDG (JI,2) * DMK%XSWI (JI,2)
      DM%XSOIL_TSWI(IMASK) = DM%XSOIL_TSWI(IMASK) + PK%XPATCH(JI) * PK%XDG (JI,2) * DMK%XTSWI(JI,2)
      DM%XSOIL_TWG (IMASK) = DM%XSOIL_TWG (IMASK) + PK%XPATCH(JI) * PK%XDG (JI,2) * (PEK%XWG(JI,2) + PEK%XWGI(JI,2))
      DM%XSOIL_TWGI(IMASK) = DM%XSOIL_TWGI(IMASK) + PK%XPATCH(JI) * PK%XDG (JI,2) * PEK%XWGI(JI,2) 
!          
      ZSUMDG       (IMASK) = ZSUMDG(IMASK) + PK%XPATCH(JI) * PK%XDG(JI,IO%NGROUND_LAYER)        
!          
     ENDDO
  ENDDO     

  WHERE(ZPOND(:,:)>0.0)
    DM%XWGI (:,:) = DM%XWGI (:,:) / ZPOND (:,:)
  ELSEWHERE
    DM%XWGI (:,:) = XUNDEF
  ENDWHERE

!     
  IF(IO%CISBA=='3-L')THEN
!          
    ZPOND(:,:)=0.0
    DO JP=1,IO%NPATCH
      DMK => NDM%AL(JP)
      PK => NP%AL(JP)
      PEK => NPE%AL(JP)    
!
      DO JI=1,PK%NSIZE_P
        IMASK = PK%NR_P(JI)       
!
        ZWORK=MAX(0.0,PK%XDG(JI,3)-PK%XDG(JI,2))
!
!           Remenber: no ice in the third layer of 3-L
        ZPOND         (IMASK,3) = ZPOND        (IMASK,3) + PK%XPATCH(JI) * ZWORK
        DM%XSWI       (IMASK,3) = DM%XSWI      (IMASK,3) + PK%XPATCH(JI) * ZWORK * DMK%XSWI (JI,3)
        DM%XSOIL_SWI  (IMASK  ) = DM%XSOIL_SWI (IMASK  ) + PK%XPATCH(JI) * ZWORK * DMK%XSWI (JI,3)  
        DM%XTSWI      (IMASK,3) = DM%XTSWI     (IMASK,3) + PK%XPATCH(JI) * ZWORK * DMK%XTSWI(JI,3)
        DM%XSOIL_TSWI (IMASK  ) = DM%XSOIL_TSWI(IMASK  ) + PK%XPATCH(JI) * ZWORK * DMK%XTSWI(JI,3)  
        DM%XSOIL_TWG  (IMASK  ) = DM%XSOIL_TWG (IMASK  ) + PK%XPATCH(JI) * ZWORK * PEK%XWG  (JI,3)  
!
       ENDDO
    ENDDO
!
    WHERE(ZPOND(:,3)>0.0)
          DM%XSWI (:,3) = DM%XSWI (:,3) / ZPOND(:,3)
          DM%XTSWI(:,3) = DM%XTSWI(:,3) / ZPOND(:,3)
    ELSEWHERE
          DM%XSWI (:,3) = XUNDEF
          DM%XTSWI(:,3) = XUNDEF
    ENDWHERE
!
  ENDIF
  
!
!---------------------------------------------
ENDIF ! End ISBA soil scheme case   
!
!---------------------------------------------
!
!       3.     Specific treatement following snow scheme
!              ------------------------------------------
!
IF(NPE%AL(1)%TSNOW%SCHEME=="CRO") THEN

  IF ( IO%LSNOWSYTRON ) THEN

    DM%XSYTMASS(:) = 0.
    DM%XSYTMASSC(:) = 0.

    DO JP=1,IO%NPATCH
      PK => NP%AL(JP)
      DMK => NDM%AL(JP)
      !   
      DO JI=1,PK%NSIZE_P
        IMASK = PK%NR_P(JI)
        !
        !     SYTRON total outputs
        DM%XSYTMASS(IMASK) = DM%XSYTMASS(IMASK) + PK%XPATCH(JI) * DMK%XSYTMASS(JI)
        DM%XSYTMASSC(IMASK) = DM%XSYTMASSC(IMASK) + PK%XPATCH(JI) * DMK%XSYTMASSC(JI)
        !
      ENDDO
      !
    ENDDO

  ENDIF

  IF ( IO%LSNOWMAK_BOOL ) THEN

    DM%XPRODCOUNT(:) = 0.

    DO JP=1,IO%NPATCH
      PK => NP%AL(JP)
      DMK => NDM%AL(JP)
      !   
      DO JI=1,PK%NSIZE_P
        IMASK = PK%NR_P(JI)
        !
        !     Snow production counter output
        DM%XPRODCOUNT(IMASK) = DM%XPRODCOUNT(IMASK) + PK%XPATCH(JI) * DMK%XPRODCOUNT(JI)
        !
      ENDDO
      !
    ENDDO
    !
  ENDIF

ENDIF


!
!       4.     Final computation for grid-cell diag
!
!              ------------------------------------
!
!Total Soil Wetness Index and Soil Water Content (m3.m-3)
WHERE(ZSUMDG(:)>0.0)
  DM%XSOIL_SWI (:) = DM%XSOIL_SWI (:)/ZSUMDG(:)
  DM%XSOIL_TSWI(:) = DM%XSOIL_TSWI(:)/ZSUMDG(:)
  DM%XSOIL_WG  (:) = DM%XSOIL_TWG (:)/ZSUMDG(:)
  DM%XSOIL_WGI (:) = DM%XSOIL_TWGI(:)/ZSUMDG(:)
ENDWHERE
!       
!Total Soil Water Content (Liquid+Solid) and Total Frozen Content (kg/m2)
DM%XSOIL_TWG (:)= DM%XSOIL_TWG (:) * XRHOLW
DM%XSOIL_TWGI(:)= DM%XSOIL_TWGI(:) * XRHOLW
!
! Snow temperature  
WHERE(ZSNOW(:)>0.0)
      DM%XTTSNOW(:) = DM%XTTSNOW(:)/ZSNOW(:)
ELSEWHERE
      DM%XTTSNOW(:) = XUNDEF
ENDWHERE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_MISC_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE_DIAG_MISC_ISBA_n
