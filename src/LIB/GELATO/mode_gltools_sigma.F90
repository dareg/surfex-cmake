!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!GLT_LIC The GELATO model is a seaice model used in stand-alone or embedded mode. 
!GLT_LIC  It has been developed by Meteo-France. The holder of GELATO is Meteo-France.
!GLT_LIC  
!GLT_LIC  This software is governed by the CeCILL-C license under French law and biding
!GLT_LIC  by the rules of distribution of free software. See the CeCILL-C_V1-en.txt
!GLT_LIC  (English) and CeCILL-C_V1-fr.txt (French) for details. The CeCILL is a free
!GLT_LIC  software license, explicitly compatible with the GNU GPL
!GLT_LIC  (see http://www.gnu.org/licenses/license-list.en.html#CeCILL)
!GLT_LIC  
!GLT_LIC  The CeCILL-C licence agreement grants users the right to modify and re-use the
!GLT_LIC  software governed by this free software license. The exercising of this right
!GLT_LIC  is conditional upon the obligation to make available to the community the
!GLT_LIC  modifications made to the source code of the software so as to contribute to
!GLT_LIC  its evolution.
!GLT_LIC  
!GLT_LIC  In consideration of access to the source code and the rights to copy, modify
!GLT_LIC  and redistribute granted by the license, users are provided only with a limited
!GLT_LIC  warranty and the software's author, the holder of the economic rights, and the
!GLT_LIC  successive licensors only have limited liability. In this respect, the risks
!GLT_LIC  associated with loading, using, modifying and/or developing or reproducing the
!GLT_LIC  software by the user are brought to the user's attention, given its Free
!GLT_LIC  Software status, which may make it complicated to use, with the result that its
!GLT_LIC  use is reserved for developers and experienced professionals having in-depth
!GLT_LIC  computer knowledge. Users are therefore encouraged to load and test the
!GLT_LIC  suitability of the software as regards their requirements in conditions enabling
!GLT_LIC  the security of their systems and/or data to be ensured and, more generally, to
!GLT_LIC  use and operate it in the same conditions of security. 
!GLT_LIC  
!GLT_LIC  The GELATO sofware is cureently distibuted with the SURFEX software, available at 
!GLT_LIC  http://www.cnrm.meteo.fr/surfex. The fact that you download the software deemed that
!GLT_LIC  you had knowledge of the CeCILL-C license and that you accept its terms.
!GLT_LIC  Attempts to use this software in a way not complying with CeCILL-C license
!GLT_LIC  may lead to prosecution. 
!GLT_LIC 
! =======================================================================
! ======================= MODULE mode_gltools_sigma =======================
! =======================================================================
!
!   This module contains functions that allow to know sea water density
! as a function of temperature and salinity (polynomial fits)
!   Contains also
!
! Modified : (D. Salas y Melia) 
!   Add a function to compute salinity entrapment following 
!   Cox and Weeks, JGR (1988)
!
! ------------------- BEGIN MODULE mode_gltools_sigma ---------------------

MODULE mode_gltools_sigma
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK,  JPHOOK
CONTAINS

! -------------------- END MODULE mode_gltools_sigma ----------------------


! -----------------------------------------------------------------------
! ---------------------------- FUNCTION glt_sigma ---------------------------
!
! * sigma-theta as a function of temp (deg c) and salinity (mil)
! (friedrich-levitus 3rd degree polynomial fit)

FUNCTION glt_sigma(pt,ps,nx,ny)
!
  USE modd_glt_const_thm
!
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny 
  REAL, DIMENSION(nx,ny) ::                                             &
    & glt_sigma 
  REAL, DIMENSION(nx,ny), INTENT(in) ::                                 &
    & pt,ps 

  glt_sigma(:,:) =                                                          &
     & (c1 + c3 * ps(:,:) + pt(:,:) * (c2 + c5 * ps(:,:) +                &
     & pt(:,:) * (c4 + c7 * ps(:,:) + pt(:,:) * c6))) * 1.E-03  
!
END FUNCTION glt_sigma

! ------------------------- END FUNCTION glt_sigma --------------------------
! -----------------------------------------------------------------------


! -----------------------------------------------------------------------
! -------------------------- FUNCTION glt_dsigmadt --------------------------
!
! Computes d(sigma)/dt

FUNCTION glt_dsigmadt(pt,ps)
!
  USE modd_glt_const_thm
!
  IMPLICIT NONE

  REAL ::                                                               &
    & glt_dsigmadt 
  REAL, INTENT(in) ::                                                   &
    & pt,ps 

  glt_dsigmadt = (c2 + c5 * ps + 2. * pt * (c4 + c7 * ps +                  &
    & 1.5 * pt * c6)) * 1.E-3 
!
END FUNCTION glt_dsigmadt

! ------------------------ END FUNCTION glt_dsigmadt ------------------------
! -----------------------------------------------------------------------


! -----------------------------------------------------------------------
! -------------------------- FUNCTION glt_dsigmads --------------------------
! 
! Computes d(sigma)/ds

FUNCTION glt_dsigmads(pt,ps)
!
  USE modd_glt_const_thm
!
  IMPLICIT NONE

  REAL ::                                                               &
    & glt_dsigmads 
  REAL, INTENT(in) ::                                                   &
    & pt,ps 

  glt_dsigmads = (c3 + pt * (c5 + pt * c7)) * 1.E-3
!
END FUNCTION glt_dsigmads

! ------------------------ END FUNCTION glt_dsigmads ------------------------
! -----------------------------------------------------------------------

END MODULE mode_gltools_sigma
