!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ###############################################################################
SUBROUTINE APPLY_SURF_PERT (ISS, NISS, IO, S, NP, NPE)
!     ###############################################################################
!
!!****  *APPLY_SURF_PERT * - Applies surface perturbation patterns generated in pertsurf.F90
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
!!
!!    AUTHOR
!!    ------
!!     F. Bouttier
!!
!!    MODIFICATIONS
!!    -------------
!!      F. Bouttier  01/2013 : Apply random perturbations for ensembles
!!      u. Andrae    06/2023 : Make it a separate subroutine and introduce clipping for all variables
!!-------------------------------------------------------------------
!
USE MODD_SSO_n, ONLY : SSO_t, SSO_NP_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_S_t, ISBA_NP_t, ISBA_NPE_t, &
                      & ISBA_P_t, ISBA_PE_t
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(SSO_t), INTENT(INOUT) :: ISS 
TYPE(SSO_NP_t), INTENT(INOUT) :: NISS
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t), INTENT(INOUT) ::NPE
!
!
!*      0.2    declarations of local variables
!
!* forcing variables
!
TYPE(ISBA_P_t), POINTER :: PK
TYPE(ISBA_PE_t), POINTER :: PEK
TYPE(SSO_t), POINTER :: ISSK
!
! dimensions and loop counters
!
INTEGER :: JP, IMASK, JI           ! Loop counters
INTEGER :: JKV, JKL, JKC, JKA, JKZ ! Perturbation field index
!
! logical units
!
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
! --------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('APPLY_SURF_PERT',0,ZHOOK_HANDLE)
! --------------------------------------------------------------------------------------

DO JP = 1,IO%NPATCH
    PK => NP%AL(JP)
    PEK => NPE%AL(JP)
    ISSK => NISS%AL(JP)

    JKV=1*IO%NPATCH-MOD(JP,IO%NPATCH)
    JKL=2*IO%NPATCH-MOD(JP,IO%NPATCH)
    JKC=3*IO%NPATCH-MOD(JP,IO%NPATCH)
    JKA=4*IO%NPATCH-MOD(JP,IO%NPATCH)
    JKZ=5*IO%NPATCH-MOD(JP,IO%NPATCH)

    DO JI = 1,PK%NSIZE_P
      IMASK = PK%NR_P(JI)
      !
      ! Random perturbation for ensembles:
      ! reapply original perturbation patterns
      IF (ASSOCIATED(S%XPERTVEG)) THEN
        IF(PEK%XVEG(JI)/=XUNDEF) PEK%XVEG(JI) = &
        MAX(MIN(PEK%XVEG(JI) *(1.+ S%XPERTVEG(IMASK,JP) ),IO%XPERT_HIGH(JKV)),IO%XPERT_LOW(JKV))
      ENDIF
      IF (ASSOCIATED(S%XPERTLAI)) THEN
        IF(PEK%XLAI(JI)/=XUNDEF) PEK%XLAI(JI) = &
        MAX(MIN(PEK%XLAI(JI) *(1.+ S%XPERTLAI(IMASK,JP) ),IO%XPERT_HIGH(JKL)),IO%XPERT_LOW(JKL))
      ENDIF
      IF (ASSOCIATED(S%XPERTCV)) THEN
        IF(PEK%XCV(JI)/=XUNDEF)  PEK%XCV(JI)  = &
        MAX(MIN(PEK%XCV(JI) *(1.+ S%XPERTCV(IMASK,JP) ),IO%XPERT_HIGH(JKC)),IO%XPERT_LOW(JKC))
      ENDIF

      IF (ASSOCIATED(S%XPERTALB)) THEN
        IF(PEK%XALBNIR_VEG(JI)/=XUNDEF)   PEK%XALBNIR_VEG(JI)   = &
        MAX(MIN(PEK%XALBNIR_VEG(JI)  *(1.+ S%XPERTALB(IMASK,JP) ),IO%XPERT_HIGH(JKA)),IO%XPERT_LOW(JKA))
        IF(PEK%XALBVIS_VEG(JI)/=XUNDEF)   PEK%XALBVIS_VEG(JI)   = &
        MAX(MIN(PEK%XALBVIS_VEG(JI)  *(1.+ S%XPERTALB(IMASK,JP) ),IO%XPERT_HIGH(JKA)),IO%XPERT_LOW(JKA))
        IF(PEK%XALBUV_VEG(JI)/=XUNDEF)    PEK%XALBUV_VEG (JI)   = &
        MAX(MIN(PEK%XALBUV_VEG(JI)   *(1.+ S%XPERTALB(IMASK,JP) ),IO%XPERT_HIGH(JKA)),IO%XPERT_LOW(JKA))
        IF(PEK%XALBNIR_SOIL(JI)/=XUNDEF)   PEK%XALBNIR_SOIL(JI) = &
        MAX(MIN(PEK%XALBNIR_SOIL(JI) *(1.+ S%XPERTALB(IMASK,JP) ),IO%XPERT_HIGH(JKA)),IO%XPERT_LOW(JKA))
        IF(PEK%XALBVIS_SOIL(JI)/=XUNDEF)   PEK%XALBVIS_SOIL(JI) = &
        MAX(MIN(PEK%XALBVIS_SOIL(JI) *(1.+ S%XPERTALB(IMASK,JP) ),IO%XPERT_HIGH(JKA)),IO%XPERT_LOW(JKA))
        IF(PEK%XALBUV_SOIL(JI)/=XUNDEF)    PEK%XALBUV_SOIL (JI) = &
        MAX(MIN(PEK%XALBUV_SOIL(JI)  *(1.+ S%XPERTALB(IMASK,JP) ),IO%XPERT_HIGH(JKA)),IO%XPERT_LOW(JKA))
      ENDIF

      IF (ASSOCIATED(S%XPERTZ0)) THEN
        IF(PEK%XZ0(JI)/=XUNDEF)       PEK%XZ0(JI)       = &
        MAX(MIN(PEK%XZ0(JI)       *(1.+ S%XPERTZ0(IMASK,JP) ),IO%XPERT_HIGH(JKZ)),IO%XPERT_LOW(JKZ))
        IF(ISSK%XZ0EFFIP(JI)/=XUNDEF) ISSK%XZ0EFFIP(JI) = &
        MAX(MIN(ISSK%XZ0EFFIP(JI) *(1.+ S%XPERTZ0(IMASK,JP) ),IO%XPERT_HIGH(JKZ)),IO%XPERT_LOW(JKZ))
        IF(ISSK%XZ0EFFIM(JI)/=XUNDEF) ISSK%XZ0EFFIM(JI) = &
        MAX(MIN(ISSK%XZ0EFFIM(JI) *(1.+ S%XPERTZ0(IMASK,JP) ),IO%XPERT_HIGH(JKZ)),IO%XPERT_LOW(JKZ))
        IF(ISSK%XZ0EFFJP(JI)/=XUNDEF) ISSK%XZ0EFFJP(JI) = &
        MAX(MIN(ISSK%XZ0EFFJP(JI) *(1.+ S%XPERTZ0(IMASK,JP) ),IO%XPERT_HIGH(JKZ)),IO%XPERT_LOW(JKZ))
        IF(ISSK%XZ0EFFJM(JI)/=XUNDEF) ISSK%XZ0EFFJM(JI) = &
        MAX(MIN(ISSK%XZ0EFFJM(JI) *(1.+ S%XPERTZ0(IMASK,JP) ),IO%XPERT_HIGH(JKZ)),IO%XPERT_LOW(JKZ))
      ENDIF
    ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('APPLY_SURF_PERT',1,ZHOOK_HANDLE)

!==========================================================================================
END SUBROUTINE APPLY_SURF_PERT
