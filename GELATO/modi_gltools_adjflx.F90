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
! ====================== MODULE modi_gltools_adjflx =======================
! =======================================================================
!
! Goal:
! -----
!   This module contains a subroutine that prints the integral of a 
! quantity separately in both hemispheres.
!
! Created : 2001/08 (D. Salas y Melia)
!           Does part of the job formerly handled by thermo_ice routine.
! Modified: 2010/06 (D. Salas y Melia)
!           Average flux of a field (with a criterion on space) 
!
! ------------------- BEGIN MODULE modi_gltools_adjflx --------------------
!
!THXS_SFX!MODULE modi_gltools_adjflx
!THXS_SFX!INTERFACE 
!THXS_SFX!!
!THXS_SFX!FUNCTION gltools_adjflx(tpdom,ocrit,pfield)
!THXS_SFX!  USE modd_types_glt
!THXS_SFX!  USE modd_glt_param
!THXS_SFX!  TYPE(t_dom), DIMENSION(nx,ny), INTENT(in) ::  &
!THXS_SFX!        tpdom
!THXS_SFX!  LOGICAL, DIMENSION(nx,ny), INTENT(in) ::  &
!THXS_SFX!        ocrit
!THXS_SFX!  REAL, DIMENSION(nx,ny), INTENT(in) ::  &
!THXS_SFX!        pfield
!THXS_SFX!  REAL, DIMENSION(nx,ny) ::  &
!THXS_SFX!        gltools_adjflx 
!THXS_SFX!END FUNCTION gltools_adjflx
!THXS_SFX!!
!THXS_SFX!END INTERFACE
!THXS_SFX!END MODULE modi_gltools_adjflx
!
! -------------------- END MODULE modi_gltools_adjflx ---------------------
!
!
! -----------------------------------------------------------------------
! -------------------------- FUNCTION gltools_adjflx ----------------------------
!
! * Subroutine used to check global sea ice extent, area and volume in 
! both hemispheres.
!
FUNCTION gltools_adjflx(tpdom,ocrit,pfield,nx,ny,dtt)
!
  USE modd_glt_const_thm
  USE modd_types_glt, only: t_dom
!
  IMPLICIT NONE
!
  INTEGER,INTENT(IN) :: nx,ny
  REAL,INTENT(IN) :: dtt
  TYPE(t_dom), DIMENSION(nx,ny), INTENT(in) ::  &
        tpdom
  LOGICAL, DIMENSION(nx,ny), INTENT(in) ::  &
        ocrit
  REAL, DIMENSION(nx,ny), INTENT(in) ::   &
        pfield
  REAL, DIMENSION(nx,ny) ::  &
        gltools_adjflx 
!
  REAL ::  &
        zint,zsrf
!
!
! * Compute field integral
  zint = SUM( tpdom(:,:)%srf*pfield(:,:), MASK=ocrit(:,:) )
!
! * Compute surface of the domain
  zsrf = SUM( tpdom(:,:)%srf, MASK=ocrit(:,:) )
!
! * Compute the correction
  gltools_adjflx(:,:) = 0.
  IF ( zsrf>epsil1 ) THEN
      WHERE( ocrit(:,:) )
        gltools_adjflx(:,:) = dtt*zint/zsrf
      ENDWHERE
  ENDIF
!
END FUNCTION gltools_adjflx
!
! ------------------------ END FUNCTION gltools_adjflx --------------------------
! -----------------------------------------------------------------------
