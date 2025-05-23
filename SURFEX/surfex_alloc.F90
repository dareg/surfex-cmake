!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE SURFEX_ALLOC(YDSURFEX)
!
USE MODD_TEB_PAR, ONLY : NTEB_PATCH_MAX
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE_ECOSG, NTILESFC
!
USE MODD_SURFEX_n, ONLY : SURFEX_t
!
USE MODD_CH_EMIS_FIELD_n, ONLY : CH_EMIS_FIELD_INIT
USE MODD_CH_SNAP_n, ONLY : CH_EMIS_SNAP_INIT
USE MODD_CH_SURF_n, ONLY : CH_SURF_INIT
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_INIT
USE MODD_DUMMY_SURF_FIELDS_n, ONLY : DUMMY_SURF_FIELDS_INIT
USE MODD_EMIS_GR_FIELD_n, ONLY : EMIS_GR_FIELD_INIT
USE MODD_SFX_GRID_n, ONLY : GRID_INIT, GRID_NP_INIT
USE MODD_CANOPY_n, ONLY : CANOPY_INIT
USE MODD_GRID_CONF_PROJ_n, ONLY : GRID_CONF_PROJ_INIT
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_INIT
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_INIT
USE MODD_SSO_n, ONLY : SSO_INIT, SSO_NP_INIT
USE MODD_SV_n, ONLY : SV_INIT
!
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_INIT
!
USE MODD_IDEAL_n, ONLY : IDEAL_INIT
!
USE MODD_DST_n, ONLY : DST_NP_INIT
USE MODD_SLT_n, ONLY : SLT_INIT
!
USE MODD_DIAG_n, ONLY : DIAG_INIT, DIAG_NP_INIT, DIAG_OPTIONS_INIT
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_INIT, DIAG_EVAP_ISBA_NP_INIT
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_INIT, DIAG_MISC_ISBA_NP_INIT
USE MODD_DIAG_OCEAN_n, ONLY : DIAG_OCEAN_INIT
USE MODD_DIAG_MISC_SEAICE_n, ONLY : DIAG_MISC_SEAICE_INIT
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DIAG_MISC_FLAKE_INIT
USE MODD_DIAG_MISC_TEB_OPTIONS_n, ONLY : DIAG_MISC_TEB_OPTIONS_INIT
USE MODD_DIAG_UTCI_TEB_n, ONLY : DIAG_UTCI_TEB_INIT
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_NP_INIT
!
USE MODD_DATA_BEM_n, ONLY : DATA_BEM_INIT
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_INIT
USE MODD_BLD_DESCRIPTION_n, ONLY : BLD_DESC_INIT
USE MODD_CH_TEB_n, ONLY : CH_TEB_INIT
USE MODD_DATA_TEB_n, ONLY : DATA_TEB_INIT
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_INIT
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_INIT
USE MODD_TEB_PANEL_n, ONLY : TEB_PANEL_INIT
USE MODD_BEM_n, ONLY : BEM_NP_INIT
USE MODD_TEB_n, ONLY : TEB_NP_INIT
!
USE MODD_CH_FLAKE_n, ONLY : CH_FLAKE_INIT
USE MODD_FLAKE_n, ONLY : FLAKE_INIT
!
USE MODD_CH_WATFLUX_n, ONLY : CH_WATFLUX_INIT
USE MODD_WATFLUX_n, ONLY : WATFLUX_INIT
!
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_INIT, CH_ISBA_NP_INIT
USE MODD_AGRI_n, ONLY : AGRI_NP_INIT
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_INIT
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_INIT, GR_BIOG_NP_INIT
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_INIT
USE MODD_ISBA_n, ONLY : ISBA_S_INIT, ISBA_K_INIT, ISBA_P_INIT, &
                        ISBA_NK_INIT, ISBA_NP_INIT, ISBA_NPE_INIT
!
USE MODD_CH_SEAFLUX_n, ONLY : CH_SEAFLUX_INIT
USE MODD_DATA_SEAFLUX_n, ONLY : DATA_SEAFLUX_INIT
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_INIT
USE MODD_OCEAN_n, ONLY : OCEAN_INIT
USE MODD_OCEAN_REL_n, ONLY : OCEAN_REL_INIT
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
TYPE (SURFEX_t), INTENT (INOUT) :: YDSURFEX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("SURFEX_ALLOC",0,ZHOOK_HANDLE)
  !
  CALL DIAG_OPTIONS_INIT(YDSURFEX%FM%DFO)
  CALL DIAG_INIT(YDSURFEX%FM%DF)
  CALL DIAG_INIT(YDSURFEX%FM%DFC)
  CALL DIAG_MISC_FLAKE_INIT(YDSURFEX%FM%DMF)
  ! 
  CALL GRID_INIT(YDSURFEX%FM%G)  
  CALL CANOPY_INIT(YDSURFEX%FM%SB)  
  CALL CH_FLAKE_INIT(YDSURFEX%FM%CHF)
  CALL FLAKE_INIT(YDSURFEX%FM%F)
  !  
  !
  CALL DIAG_OPTIONS_INIT(YDSURFEX%WM%DWO)
  CALL DIAG_INIT(YDSURFEX%WM%DW)
  CALL DIAG_INIT(YDSURFEX%WM%DWC)
  !  
  CALL GRID_INIT(YDSURFEX%WM%G)
  CALL CANOPY_INIT(YDSURFEX%WM%SB)  
  CALL CH_WATFLUX_INIT(YDSURFEX%WM%CHW)
  CALL WATFLUX_INIT(YDSURFEX%WM%W)
  !
  !
  CALL DIAG_OPTIONS_INIT(YDSURFEX%SM%SD%O) 
  CALL DIAG_INIT(YDSURFEX%SM%SD%D)
  CALL DIAG_INIT(YDSURFEX%SM%SD%DC)
  CALL DIAG_INIT(YDSURFEX%SM%SD%DI)
  CALL DIAG_INIT(YDSURFEX%SM%SD%DIC)
  CALL DIAG_OCEAN_INIT(YDSURFEX%SM%SD%GO)  
  CALL DIAG_MISC_SEAICE_INIT(YDSURFEX%SM%SD%DMI)
  !  
  CALL DATA_SEAFLUX_INIT(YDSURFEX%SM%DTS)
  CALL GRID_INIT(YDSURFEX%SM%G)
  CALL CANOPY_INIT(YDSURFEX%SM%SB)  
  CALL CH_SEAFLUX_INIT(YDSURFEX%SM%CHS)
  CALL SEAFLUX_INIT(YDSURFEX%SM%S)
  CALL OCEAN_INIT(YDSURFEX%SM%O)
  CALL OCEAN_REL_INIT(YDSURFEX%SM%OR)
  !
  !
  CALL DIAG_OPTIONS_INIT(YDSURFEX%IM%ID%O)  
  CALL DIAG_INIT(YDSURFEX%IM%ID%D)
  CALL DIAG_INIT(YDSURFEX%IM%ID%DC)
  CALL DIAG_NP_INIT(YDSURFEX%IM%ID%ND,NVEGTYPE_ECOSG)
  CALL DIAG_NP_INIT(YDSURFEX%IM%ID%NDC,NVEGTYPE_ECOSG)
  CALL DIAG_EVAP_ISBA_INIT(YDSURFEX%IM%ID%DE)
  CALL DIAG_EVAP_ISBA_INIT(YDSURFEX%IM%ID%DEC)
  CALL DIAG_EVAP_ISBA_NP_INIT(YDSURFEX%IM%ID%NDE,NVEGTYPE_ECOSG)
  CALL DIAG_EVAP_ISBA_NP_INIT(YDSURFEX%IM%ID%NDEC,NVEGTYPE_ECOSG) 
  CALL DIAG_MISC_ISBA_INIT(YDSURFEX%IM%ID%DM)
  CALL DIAG_MISC_ISBA_NP_INIT(YDSURFEX%IM%ID%NDM,NVEGTYPE_ECOSG)
  !
  CALL DATA_ISBA_INIT(YDSURFEX%IM%DTV)
  CALL CANOPY_INIT(YDSURFEX%IM%SB)
  CALL ISBA_OPTIONS_INIT(YDSURFEX%IM%O)
  CALL ISBA_S_INIT(YDSURFEX%IM%S)  
  CALL CH_ISBA_INIT(YDSURFEX%IM%CHI)
  CALL CH_ISBA_NP_INIT(YDSURFEX%IM%NCHI,NVEGTYPE_ECOSG)
  CALL GR_BIOG_INIT(YDSURFEX%IM%GB)
  CALL GR_BIOG_NP_INIT(YDSURFEX%IM%NGB,NVEGTYPE_ECOSG)
  CALL SSO_INIT(YDSURFEX%IM%ISS)
  CALL SSO_NP_INIT(YDSURFEX%IM%NISS,NVEGTYPE_ECOSG)
  CALL GRID_INIT(YDSURFEX%IM%G)
  CALL GRID_NP_INIT(YDSURFEX%IM%NG,NVEGTYPE_ECOSG)
  CALL ISBA_K_INIT(YDSURFEX%IM%K)
  CALL ISBA_NK_INIT(YDSURFEX%IM%NK,NVEGTYPE_ECOSG)
  CALL ISBA_NP_INIT(YDSURFEX%IM%NP,NVEGTYPE_ECOSG)
  CALL ISBA_NPE_INIT(YDSURFEX%IM%NPE,NVEGTYPE_ECOSG)  
  CALL AGRI_NP_INIT(YDSURFEX%IM%NAG,NVEGTYPE_ECOSG)
  !
  !
  CALL DIAG_NP_INIT(YDSURFEX%GDM%VD%ND,NTEB_PATCH_MAX)  
  CALL DIAG_EVAP_ISBA_NP_INIT(YDSURFEX%GDM%VD%NDE,NTEB_PATCH_MAX)
  CALL DIAG_EVAP_ISBA_NP_INIT(YDSURFEX%GDM%VD%NDEC,NTEB_PATCH_MAX)
  CALL DIAG_MISC_ISBA_NP_INIT(YDSURFEX%GDM%VD%NDM,NTEB_PATCH_MAX)  
  !  
  CALL DATA_ISBA_INIT(YDSURFEX%GDM%DTV)
  CALL ISBA_OPTIONS_INIT(YDSURFEX%GDM%O)
  CALL ISBA_S_INIT(YDSURFEX%GDM%S)  
  CALL GR_BIOG_INIT(YDSURFEX%GDM%GB)
  CALL ISBA_K_INIT(YDSURFEX%GDM%K)
  CALL ISBA_P_INIT(YDSURFEX%GDM%P)
  CALL ISBA_NPE_INIT(YDSURFEX%GDM%NPE,NTEB_PATCH_MAX)
  !
  !
  CALL DIAG_NP_INIT(YDSURFEX%GRM%VD%ND,NTEB_PATCH_MAX)
  CALL DIAG_EVAP_ISBA_NP_INIT(YDSURFEX%GRM%VD%NDE,NTEB_PATCH_MAX)
  CALL DIAG_EVAP_ISBA_NP_INIT(YDSURFEX%GRM%VD%NDEC,NTEB_PATCH_MAX)
  CALL DIAG_MISC_ISBA_NP_INIT(YDSURFEX%GRM%VD%NDM,NTEB_PATCH_MAX)
  !
  CALL DATA_ISBA_INIT(YDSURFEX%GRM%DTV)
  CALL ISBA_OPTIONS_INIT(YDSURFEX%GRM%O)
  CALL ISBA_S_INIT(YDSURFEX%GRM%S)  
  CALL GR_BIOG_INIT(YDSURFEX%GRM%GB)
  CALL ISBA_K_INIT(YDSURFEX%GRM%K)
  CALL ISBA_P_INIT(YDSURFEX%GRM%P)
  CALL ISBA_NPE_INIT(YDSURFEX%GRM%NPE,NTEB_PATCH_MAX)
  !
  !
  CALL DIAG_OPTIONS_INIT(YDSURFEX%TM%TD%O)
  CALL DIAG_INIT(YDSURFEX%TM%TD%D)
  CALL DIAG_MISC_TEB_OPTIONS_INIT(YDSURFEX%TM%TD%MTO)
  CALL DIAG_MISC_TEB_NP_INIT(YDSURFEX%TM%TD%NDMT,NTEB_PATCH_MAX)   
  CALL DIAG_MISC_TEB_NP_INIT(YDSURFEX%TM%TD%NDMTC,NTEB_PATCH_MAX)
  CALL DIAG_UTCI_TEB_INIT(YDSURFEX%TM%TD%DUT)
  ! 
  CALL DATA_TEB_INIT(YDSURFEX%TM%DTT)
  CALL TEB_OPTIONS_INIT(YDSURFEX%TM%TOP)
  CALL CANOPY_INIT(YDSURFEX%TM%SB)
  CALL GRID_INIT(YDSURFEX%TM%G)  
  CALL CH_TEB_INIT(YDSURFEX%TM%CHT)
  CALL TEB_PANEL_INIT(YDSURFEX%TM%TPN)
  CALL TEB_IRRIG_INIT(YDSURFEX%TM%TIR)    
  CALL TEB_NP_INIT(YDSURFEX%TM%NT,NTEB_PATCH_MAX)  
  !
  CALL DATA_BEM_INIT(YDSURFEX%TM%DTB)    
  CALL BEM_OPTIONS_INIT(YDSURFEX%TM%BOP)  
  CALL BLD_DESC_INIT(YDSURFEX%TM%BDD)
  CALL BEM_NP_INIT(YDSURFEX%TM%NB,NTEB_PATCH_MAX)  
  !
  !
  CALL DATA_COVER_INIT(YDSURFEX%DTCO)
  CALL DATA_TSZ0_INIT(YDSURFEX%DTZ)
  CALL DUMMY_SURF_FIELDS_INIT(YDSURFEX%DUU)
  !
  CALL GRID_CONF_PROJ_INIT(YDSURFEX%GCP)
  CALL SURF_ATM_GRID_INIT(YDSURFEX%UG)
  CALL SURF_ATM_INIT(YDSURFEX%U)
  CALL DIAG_OPTIONS_INIT(YDSURFEX%DUO) 
  CALL DIAG_INIT(YDSURFEX%DU) 
  CALL DIAG_INIT(YDSURFEX%DUC)
  CALL DIAG_NP_INIT(YDSURFEX%DUP,NTILESFC) 
  CALL DIAG_NP_INIT(YDSURFEX%DUPC,NTILESFC)
  CALL SSO_INIT(YDSURFEX%USS)
  CALL CANOPY_INIT(YDSURFEX%SB)
  !
  CALL DIAG_INIT(YDSURFEX%DL)
  CALL DIAG_INIT(YDSURFEX%DLC)
  CALL IDEAL_INIT(YDSURFEX%L)
  !
  CALL SV_INIT(YDSURFEX%SV)
  CALL CH_SURF_INIT(YDSURFEX%CHU)  
  CALL CH_EMIS_FIELD_INIT(YDSURFEX%CHE)
  CALL CH_EMIS_SNAP_INIT(YDSURFEX%CHN)
  CALL EMIS_GR_FIELD_INIT(YDSURFEX%EGF)  
  CALL DST_NP_INIT(YDSURFEX%NDST,NVEGTYPE_ECOSG)
  CALL SLT_INIT(YDSURFEX%SLT)
  !
IF (LHOOK) CALL DR_HOOK("SURFEX_ALLOC",1,ZHOOK_HANDLE)
!
END SUBROUTINE SURFEX_ALLOC
