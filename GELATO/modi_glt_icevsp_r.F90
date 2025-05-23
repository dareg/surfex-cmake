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
! ======================== MODULE modi_glt_icevsp_r =========================
! =======================================================================
!
! Goal:
! -----
!   This module contains a subroutine that is used to update sea ice
! salinity, taking into account desalination processes.
!
! Method:
! --------
!   Vancoppenolle et al., O. Modelling (2009)
!   We assume that the sea ice gltools_enthalpy does not change. Implicitly,
! it means that the temperature of sea ice will increase if 
! desalination occurs. The temperature change will be taken into 
! account at the next time step, in glt_vhdiff_r routine.
!
! Created : 2010/03 (D. Salas y Melia)
!
! --------------------- BEGIN MODULE modi_glt_icevsp_r ----------------------
!
!THXS_SFX!MODULE modi_glt_icevsp_r
!THXS_SFX!INTERFACE 
!THXS_SFX!!
!THXS_SFX!SUBROUTINE glt_icevsp_r( tpsit,pvsp )
!THXS_SFX!  USE modd_types_glt
!THXS_SFX!  USE modd_glt_param
!THXS_SFX!  TYPE(t_sit), DIMENSION(nt,np), INTENT(inout) ::  &
!THXS_SFX!        tpsit   
!THXS_SFX!  REAL, DIMENSION(nl,nt,np), INTENT(inout) ::  &
!THXS_SFX!        pvsp
!THXS_SFX!END SUBROUTINE glt_icevsp_r
!THXS_SFX!!
!THXS_SFX!END INTERFACE
!THXS_SFX!END MODULE modi_glt_icevsp_r
!
! ---------------------- END MODULE modi_glt_icevsp_r -----------------------
!
!
!
! -----------------------------------------------------------------------
! ------------------------- SUBROUTINE glt_icevsp_r -------------------------
!
SUBROUTINE glt_icevsp_r( tpsit,pvsp,&
  nilay,nl,np,nt,height,sf3tinv )
!
  USE modd_glt_const_thm
  USE modd_types_glt, only: t_sit
!
  IMPLICIT NONE
!
  INTEGER,INTENT(IN) :: nl,nt,np,nilay
  REAL,DIMENSION(:),INTENT(IN) :: height,sf3tinv
  TYPE(t_sit), DIMENSION(nt,np), INTENT(in) ::  &
        tpsit   
  REAL, DIMENSION(nl,nt,np), INTENT(inout) ::  &
        pvsp
!
! .. Local variables
!
  REAL, DIMENSION(np) ::  &
    zqsalt
  REAL, DIMENSION(nt,np) ::  &
    zdssi,zssieq
! Coefficients of the MY ice salinity profile (Schwarzacher, JGR 64:2357-2367,
! 1959)
  REAL, PARAMETER ::  &
    ppa=0.407, ppb=0.573
! Low salinity and high salinity transition coefficients
  REAL, PARAMETER ::  &
    ppslo=3., ppshi=4.5
  INTEGER ::  &
    jl
  REAL ::  &
    zint,z 
  REAL, DIMENSION(nilay) ::  &
    zdepth,zsalt0
  REAL, DIMENSION(nt,np) ::  &
    zalf
!
!
!
! 1. Initialization
! ==================
!
! .. We need to know at what depth wrt the ice-atm or ice-snow 
! interface every point lies.
!
DO jl=1,nilay
  zdepth(jl) = 1. - 0.5*( height(jl)+height(jl+1) )
END DO
!
! .. Compute the standard salinity profile of ice with less than
! 3 psu salinity
!
DO jl=1,nilay
  z = zdepth(jl)
  zsalt0(jl) = 1.-COS( pi*EXP( ppa/(z+ppb)*LOG(z) ) )
END DO
zint = SUM( sf3tinv(:)*zsalt0(:) )
zsalt0(:) = zsalt0(:)/zint
! 
!
!
! 2. Define salinity profiles
! ============================
!
! .. Weighing parameter for intermediate salinity
!
WHERE( tpsit(:,:)%ssi >= ppslo .AND. tpsit(:,:)%ssi < ppshi )
  zalf(:,:) = ( ppshi-tpsit(:,:)%ssi ) / ( ppshi-ppslo )
ENDWHERE
!
DO jl=1,nilay
!
  WHERE( tpsit(:,:)%ssi < ppshi )
!
! .. Low salinity: multiyear ice salinity profile
!
    WHERE( tpsit(:,:)%ssi < ppslo )
      pvsp(jl,:,:) = tpsit(:,:)%ssi * zsalt0(jl)
!
! .. Intermediate salinity: combination between a constant profile and a 
! multiyear ice salinity profile
!
    ELSEWHERE
      pvsp(jl,:,:) = tpsit(:,:)%ssi *  &
        ( zalf(:,:) * zsalt0(jl) + ( 1.-zalf(:,:) ) )
    ENDWHERE
  ELSEWHERE
!
! ..High salinity: constant profile
!
    pvsp(jl,:,:) = tpsit(:,:)%ssi
  ENDWHERE
END DO
!
! .. Snow salinity is zero !
!
pvsp(nilay+1,:,:) = 0.
!
END SUBROUTINE glt_icevsp_r
!
! ---------------------- END SUBROUTINE glt_icevsp_r ------------------------
! -----------------------------------------------------------------------
