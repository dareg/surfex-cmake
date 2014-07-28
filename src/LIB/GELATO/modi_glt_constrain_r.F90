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
! ======================= MODULE modi_glt_constrain_r =======================
! =======================================================================
!
! Goal:
! ----- 
!    This module contains a subroutine that allows to constrain 
! sea ice. So far, the only option is newtonian damping. The variables
! that can be constrained are: ice surface temperature, concentration
! and thickness.
!  
! Method:
! -------
!    Newtonian damping. This method does not ensure energy is conserved.
!
! Created : 2012/03 (D. Salas y Melia) 
! Modified: No
! 
! -------------------- BEGIN MODULE modi_glt_constrain_r --------------------
!
!THXS_SFX!MODULE modi_glt_constrain_r
!THXS_SFX!INTERFACE
!THXS_SFX!!
!THXS_SFX!SUBROUTINE glt_constrain_r( tpdom,tpmxl,tpsit,tpsil,tpdia,tpsit_d )
!THXS_SFX!  USE modd_types_glt
!THXS_SFX!  USE modd_glt_param
!THXS_SFX!  TYPE(t_dom), DIMENSION(np), INTENT(in) ::  &
!THXS_SFX!        tpdom
!THXS_SFX!  TYPE(t_mxl), DIMENSION(np), INTENT(inout) ::  &
!THXS_SFX!        tpmxl
!THXS_SFX!  TYPE(t_sit), DIMENSION(nt,np), INTENT(inout) ::  &
!THXS_SFX!        tpsit
!THXS_SFX!  TYPE(t_vtp), DIMENSION(nt,np), INTENT(inout) ::  &
!THXS_SFX!        tpsil
!THXS_SFX!  TYPE(t_dia), DIMENSION(np), INTENT(inout) ::  &
!THXS_SFX!        tpdia
!THXS_SFX!  TYPE(t_sit), DIMENSION(ntd,np), INTENT(in) ::  &
!THXS_SFX!        tpsit_d
!THXS_SFX!END SUBROUTINE glt_constrain_r
!THXS_SFX!!
!THXS_SFX!END INTERFACE
!THXS_SFX!END MODULE modi_glt_constrain_r
!
! --------------------- END MODULE modi_glt_constrain_r ---------------------
!
!
! -----------------------------------------------------------------------
! ---------------------- SUBROUTINE glt_constrain_r -------------------------
!
!
SUBROUTINE glt_constrain_r( tpdom,tpmxl,tpsit,tpsil,tpdia,tpsit_d )
  USE modd_types_glt
  USE modd_glt_param
  USE modd_glt_const_thm
  USE mode_glt_stats_r
  USE modi_gltools_newice_r
  USE mode_gltools_enthalpy
  USE modi_gltools_glterr
!
  IMPLICIT NONE 
!
  TYPE(t_dom), DIMENSION(np), INTENT(in) ::  &
        tpdom
  TYPE(t_mxl), DIMENSION(np), INTENT(inout) ::  &
        tpmxl
  TYPE(t_sit), DIMENSION(nt,np), INTENT(inout) ::  &
        tpsit
  TYPE(t_vtp), DIMENSION(nt,np), INTENT(inout) ::  &
        tpsil
  TYPE(t_dia), DIMENSION(np), INTENT(inout) ::  &
        tpdia
  TYPE(t_sit), DIMENSION(ntd,np), INTENT(in) ::  &
        tpsit_d
!
  INTEGER, PARAMETER ::  &
        ppcent=0.01
  INTEGER ::  &
        jk
  CHARACTER(5) ::  &
        ydmp
  CHARACTER(300) ::  &
        ymess
  REAL, DIMENSION(np) ::  &
        zwork,zwork2,zfsit,zhsit0,zdhsit,zdhsit0,zfac,zfacfsi,  &
        zenti0,zents0,zenti,zents
  REAL, DIMENSION(nt,np) ::  &
        zfsi,zhsi,zhsinew
!
!
!
! 1. Initialization
! ==================
!
! These error checks are now probably unnecessary, since ntd is now
! defined from cnflxin, and after a check in readnam.
!
!! .. Get first dimension of the constraint (to now if we apply mono or 
!! multi-category damping). Mono-category damping (i.e. we have e.g. just 
!! one thickness data per grid point) is applied if the first dimension 
!! of tpsit_d is found to be 1. Else, if the first dimension is equal to
!! nt (number of ice thicknesses), every category is constrained to the 
!! data concerning the same category.
!!
!! .. Error on the number of categories in the damping data
!!
!  IF ( ntd /=1 .AND. ntd /= nt ) THEN
!    WRITE( ymess, FMT='("First dimension of the damping array tpsit_d  &
!      & provided to glt_gelato routine = ",I2,". Should be equal to nt =",I2,  &
!      & "or equal to 1." )' ) ntd,nt
!    CALL gltools_glterr( 'constrain_r', TRIM(ymess), 'STOP' )
!  ENDIF
!!
!! .. Conflict between parameters in gltpar and input data dimensions
!!
!  IF ( cfsidmp(1:4)=='MONO' .OR. &
!       chsidmp(1:4)=='MONO' .OR. &
!       ctsfdmp(1:4)=='MONO' ) THEN
!    IF ( ntd > 1 ) THEN
!      WRITE( ymess, FMT='("First dimension of the damping array tpsit_d&
!        & provided to glt_gelato routine = ",I2," is > 1 and suggests&
!        & multi-category restoring. Not consistent with cfsidmp = ",A, &
!        & ", chsidmp = ",A, " and ctsfdmp = ",A," . Check this." )' )  &
!        ntd, TRIM(cfsidmp), TRIM(chsidmp), TRIM(ctsfdmp)
!      CALL gltools_glterr( 'constrain_r', TRIM(ymess), 'STOP' )
!    ELSE IF ( ntd == 1 .AND. nt == 1 ) THEN
!      WRITE( ymess, FMT='("First dimension of the damping array tpsit_d&
!        & provided to glt_gelato routine = 1, is equal to the number of&
!        & ice categories nt. Multi-category restoring can be applied.&
!        & Forcing cfsidmp, chsidmp and ctsfdmp to MULTI." )' )
!      CALL gltools_glterr( 'constrain_r', TRIM(ymess), 'WARN' )
!    ENDIF
!  ENDIF
!  IF ( cfsidmp(1:5)=='MULTI' .OR. &
!       chsidmp(1:5)=='MULTI' .OR. &
!       ctsfdmp(1:5)=='MULTI' ) THEN
!    IF ( ntd == 1 .AND. nt > 1 ) THEN
!      WRITE( ymess, FMT='("First dimension of the damping array tpsit_d&
!        & provided to glt_gelato routine = 1, not equal to the number of&
!        & ice categories nt = ",I2, ". Multi-category restoring cannot&
!        & be applied. Not consistent with cfsidmp = ",A, &
!        & ", chsidmp = ",A, " and ctsfdmp = ",A," . Check this." )' )  &
!        nt, TRIM(cfsidmp), TRIM(chsidmp), TRIM(ctsfdmp)
!      CALL gltools_glterr( 'constrain_r', TRIM(ymess), 'STOP' )
!    ENDIF
!  ENDIF
!
! .. Initial snow and ice gltools_enthalpy
!
  CALL glt_aventh( tpsit,tpsil,zenti0,zents0 )
!
! .. Initialize various fields
!
  zhsinew(:,:) = 0.
!
!
! 2. Mono-category damping
! =========================
!
  IF ( ntd /= nt .OR. ( nt == 1 .AND. ntd == 1 ) ) THEN
!
! .. Total sea ice concentration
!
    zfsit(:) = SUM( tpsit(:,:)%fsi,DIM=1 )
!
!
! 2.1. Damp sea ice concentration 
! --------------------------------
!
!   Note that concentration should be damped before thickness 
! (if both are damped)
!
    IF ( cfsidmp(1:4)=='MONO' ) THEN
!
!   We assume here that if we modify sea ice concentration, we do not 
! want to change current sea ice thickness.
!  - If total sea ice concentration zfsit is more than epsil1=1.e-10, 
! we modify the thicknesses of all ice categories to conserve volume.
!  - If total sea ice concentration is less than epsil1, we consider 
! there was no sea ice at all initially, meaning that we must create
! new sea ice from nothing.
!
      IF ( (SIZE(tpsit_d(1,:)%fsi) > 0) .AND. &
           ( MINVAL( tpsit_d(1,:)%fsi ) < 0. .OR.  &
           MAXVAL( tpsit_d(1,:)%fsi ) > 1. )) THEN
        CALL gltools_glterr( 'constrain_r',  &
          'Wrong ice concentration damping data &
        & (all values should be between 0 and 1).','STOP' )
      ENDIF
!
      zwork(:) = dtt / ( xfsidmpeft*xday2sec ) *  &
        ( MIN(tpsit_d(1,:)%fsi,xfsimax) - zfsit(:) )
!
      DO jk=1,nt
        WHERE ( zfsit(:)>epsil1 )
          ! Finally, we do not want to conserve seaice volume  ....
          ! tpsit(jk,:)%hsi = tpsit(jk,:)%hsi * tpsit(jk,:)%fsi / (tpsit(jk,:)%fsi + zwork(:))
          tpsit(jk,:)%fsi = tpsit(jk,:)%fsi + zwork(:) 
        ENDWHERE
        WHERE ( tpsit(jk,:)%fsi > xfsimax )
          ! tpsit(jk,:)%hsi = tpsit(jk,:)%hsi * tpsit(jk,:)%fsi / xfsimax
          tpsit(jk,:)%fsi = xfsimax
        ENDWHERE
        if (lp4) print*,"applying SIC constraint ",tpsit_d(1,:)%fsi, ' on zfsit=',zfsit(:),&
             " delta=", zwork(:), 'output value for category ',jk,'=', tpsit(jk,:)%fsi
      END DO
!
! If zfsit<=epsil1, the concentration of every category tpsit(jk,:)<=epsil1.
! In particular, the thinnet category has a concentration <=epsil1. If new
! ice has to appear here, we decide to increase only the concentration of the 
! thinnest category, and to assign it a small thickness, in order to limit 
! the ice volume change.
!
      zfsi(:,:) = 0.
      zhsi(:,:) = 0.
      WHERE ( zfsit(:)<=epsil1 )
        zfsi(1,:) = dtt / ( xfsidmpeft*xday2sec ) *  &
        ( MIN(tpsit_d(1,:)%fsi,xfsimax) - zfsit(:) )
!          ( tpsit_d(1,:)%fsi - zfsit(:) )
        zhsi(1,:) = xhsimin
      ENDWHERE
!
! This routine will be active only where zfsi(:,:)>=epsil1 and tpsit(:,:)<epsil1
!
      CALL gltools_newice_r( zfsi,zhsi,tpmxl,tpsit,tpsil )
!
    ENDIF
!
!
! 2.2. Damp sea ice thickness
! ----------------------------
!
    IF ( TRIM(chsidmp)=='MONO_ADD' .OR. TRIM(chsidmp)=='MONO_FAC' ) THEN
!
! .. Compute average sea ice thickness over categories
!
      zhsit0(:) = glt_avhicem_r( tpsit )

!
! .. Check thickness data toward which we would like to restore
!
      IF ( (SIZE(tpsit_d(1,:)%hsi) > 0) .AND. &
           (MAXVAL( tpsit_d(1,:)%hsi ) < -1.) ) THEN 
        CALL gltools_glterr( 'constrain_r',  &
          'Wrong ice thickness damping data (all %hsi < -1).','STOP' )
      ENDIF
!
! .. Damp sea ice thickness
!
      zdhsit0(:) = dtt / ( xhsidmpeft*xday2sec ) *  &
        ( tpsit_d(1,:)%hsi - zhsit0(:) )
!      zdhsit0(:) = dtt / ( xhsidmpeft*xday2sec ) *  &
!        ( tpsit_d(1,:)%hsi - zhsit0(:) ) / AMAX1( zfsit(:),xfsic )
!
!
! 2.2.1. First method: add the same correction to all categories
! ...............................................................
!
! .. Now, modify the thickness of the ice categories to modify the mean
! ice thickness by zdhsit. Note the contribution of a dynamic bias should 
! be proportional to the thickness of every category. This won't be the 
! case at all for thermodynamic biases (due to e.g. an atmospheric radiative 
! bias or an ocean temperature bias). In more detail (examples):
!   - ocean heat flux bias: affects all ice categories in the same way
!   - air temperature bias: affects heat conduction (broadly proportional 
! to the inverse of ice thickness)
!   - SW radiative bias: tends to affect all ice categories in the same way
!
      IF ( TRIM(chsidmp)=='MONO_ADD' ) THEN
!
! .. Update ice thickness
!
        WHERE( tpsit(:,:)%fsi>epsil1 )
          zhsinew(:,:) = AMAX1(  &
            tpsit(:,:)%hsi + SPREAD( zdhsit0(:),1,nt ), 0. )
        ELSEWHERE
          zhsinew(:,:) = 0.
        ENDWHERE
!
        WHERE( tpsit(:,:)%hsi>epsil1 .AND. zhsinew(:,:)<=epsil1 )
          tpsit(:,:)%esi = .FALSE. 
        ENDWHERE
!
!        DO jk=1,nt
!
! .. Collect all negative volumes (if any)
!
!          WHERE( zhsinew(jk,:)<0. ) 
!            zdhsit(:) = zhsinew(jk,:)*tpsit(jk,:)%fsi
!            zhsinew(jk,:) = 0.
!          ELSEWHERE
!            zdhsit(:) = 0.
!          ENDWHERE
! The residual correction cannot be taken into account if the thickest 
! category has 0 thickness ! (however, apart from possible numerical 
! problems, this cannot happen)
!          IF ( jk<nt ) THEN
!            zwork(:) = SUM( tpsit(jk+1:nt,:)%fsi, DIM=1 ) 
!            WHERE( zwork(:)>epsil1 )
!              zwork2(:) = zdhsit(:) / zwork(:)
!            ELSEWHERE
!              zwork2(:) = 0.
!            ENDWHERE
!            zhsinew(jk+1:nt,:) = zhsinew(jk+1:nt,:) +  &
!              SPREAD( zwork2(:), 1, nt-jk )
!          ENDIF
!        END DO
!
! Case without initial sea ice and damping to positive ice thickness
        zfsi(:,:) = 0.
        zhsi(:,:) = 0.
        WHERE ( zfsit(:)<=epsil1 )
          zfsi(1,:) = xfsic
          zhsi(1,:) = dtt / ( xhsidmpeft*xday2sec ) * tpsit_d(1,:)%hsi
!          zhsi(1,:) = dtt / ( xhsidmpeft*xday2sec ) *  &
!            tpsit_d(1,:)%hsi / zfsi(1,:)
        ENDWHERE
!
        tpsit(:,:)%hsi = zhsinew(:,:)
!
! This routine will be active only where zfsi(:,:)>=epsil1 and tpsit(:,:)<epsil1
!
        CALL gltools_newice_r( zfsi,zhsi,tpmxl,tpsit,tpsil )
!
!
!
! 2.2.2. Correction of all ice thicknesses by the same factor
! ............................................................
!
      ELSE IF ( TRIM(chsidmp)=='MONO_FAC' ) THEN
!
!write(145)'dhsit'
!write(145)zdhsit0
!close(145)
!
        zfac(:) = 1.
        WHERE( zhsit0(:)>epsil1 )
          zfac(:) = 1. + zdhsit0(:)/zhsit0(:)
        ENDWHERE
!
! Define a multiplicative factor for sea ice concentration (to help reducing or 
! increasing sea ice thickness)
        zfacfsi(:) = 1.
        WHERE( ABS(zfac(:)-1.) > ppcent )
! Low values of the factor: decrease sea ice concentration to contribute to reducing 
! mean sea ice thickness
          WHERE( zfac(:) < 1.-ppcent )
            zfacfsi(:) = zfac(:)/(1.-ppcent)
            zfac(:) = 1.-ppcent
! High values of the factor: increase very low sea ice concentrations, but not more than 0.15
          ELSEWHERE
            WHERE( zfsit(:)>epsil1 .AND. zfsit(:)<xfsic )
              zfacfsi(:) = EXP( dtt/(3.*xday2sec)*LOG( xfsic/zfsit(:) ) )
              zfac(:) = AMIN1( zfac(:)/zfacfsi(:),1.+ppcent )
            ENDWHERE
          ENDWHERE
        ENDWHERE
!
! We do not want to modify sea ice thickness by more than 1%, in order to avoid runaway 
! thickness of categories
        DO jk=1,nt
          WHERE ( zfsit(:)>epsil1 )
            tpsit(jk,:)%fsi = tpsit(jk,:)%fsi*zfacfsi(:)
            WHERE ( zfsit(:)>epsil1 )
              zhsinew(jk,:) = zfac(:) * tpsit(jk,:)%hsi
            ENDWHERE
          ENDWHERE
        END DO
!
! Case without initial sea ice and damping to positive ice thickness
        zfsi(:,:) = 0.
        zhsi(:,:) = 0.
        WHERE ( zfsit(:)<=epsil1 )
          zfsi(1,:) = xfsic
          zhsi(1,:) = dtt / ( xhsidmpeft*xday2sec ) *  &
            tpsit_d(1,:)%hsi / zfsi(1,:)
        ENDWHERE
!
! This routine will be active only where zfsi(:,:)>=epsil1 and tpsit(:,:)<epsil1
!
        CALL gltools_newice_r( zfsi,zhsi,tpmxl,tpsit,tpsil )
!
        tpsit(:,:)%hsi = zhsinew(:,:)
!
!
     ENDIF
!
    ELSE
       IF (TRIM(chsidmp) /= 'NONE') THEN
!
! 2.2.3. Wrong chsidmp option
! ............................
!
! Attention si h_init = 0 et dhsit0>0: décider: regarder dans quelle
! cat tombe la nouvelle thickness, initialiser tout le reste (%tsf,tpsil)
       WRITE( ymess,FMT='("Wrong value for chsidmp = ",A)' ) chsidmp
       CALL gltools_glterr( 'constrain_r',TRIM(ymess), 'STOP' )
    ENDIF
  ENDIF
!
!
! 2.3. Sea ice temperature
! -------------------------
!
    IF ( ctsfdmp(1:4)=='MONO' ) THEN
      WRITE(noutlu,*) , 'Constraining seaice temperature is not yet implemented'
    ENDIF
!
!
!
! 3. Multi-category damping
! ==========================
!
! If ntd = nt, even if nt = 1 : multi-category damping
!
  ELSE       
!
!
! 3.1. Damp sea ice concentration 
! --------------------------------
!
    IF ( cfsidmp(1:5) == 'MULTI' ) THEN
! ATTENTION AUX CAS OU %FSI INITIAL = 0 => NEWICE.
      tpsit(:,:)%fsi = dtt / ( xfsidmpeft*xday2sec ) *  &
        ( tpsit_d(:,:)%fsi - tpsit(:,:)%fsi )
    ENDIF
!
!
! 3.2. Damp sea ice thickness
! ----------------------------
!
    IF ( chsidmp(1:5) == 'MULTI' ) THEN
! ATTENTION AUX CAS OU %FSI INITIAL = 0 => NEWICE.
      tpsit(:,:)%hsi = dtt / ( xhsidmpeft*xday2sec ) *  &
        ( tpsit_d(:,:)%hsi - tpsit(:,:)%hsi )
    ENDIF
!
!
! 3.3. Damp sea ice temperature
! ------------------------------
!
    IF ( ctsfdmp(1:5) == 'MULTI' ) THEN
! ATTENTION AUX CAS OU %FSI INITIAL = 0 => NEWICE.
      tpsit(:,:)%tsf = dtt / ( xtsfdmpeft*xday2sec ) *  &
        ( tpsit_d(:,:)%tsf - tpsit(:,:)%tsf )
    ENDIF
!
  ENDIF
!
!
!
! 4. Final operations
! ====================
!
!WRITE(noutlu,*) ,glt_avhicem_r( tpsit )
!WRITE(noutlu,*) ,zhsit0
!WRITE(noutlu,*) ,glt_avhicem_r( tpsit )-zhsit0(:)
!WRITE(noutlu,*) ,zdhsit0
!  IF ( SUM( glt_avhicem_r( tpsit )-zhsit0(:)-zdhsit0(:) )>epsil5 ) THEN
!      CALL gltools_glterr( 'constrain_r',  &
!        'Thickness damping not conform to specifications.', 'STOP' )
!  ENDIF
!
! .. Diagnostic for change in snow+ice gltools_enthalpy due to damping/restoring
! (there is no separation of the effects of the different operations)
!
! Final snow and ice gltools_enthalpy
!
  CALL glt_aventh( tpsit,tpsil,zenti,zents )
!
! Diagnostic
!
  tpdia(:)%dmp = ( zenti+zents-zenti0-zents0 ) / dtt
!
END SUBROUTINE glt_constrain_r
!
! --------------------- END SUBROUTINE glt_constrain_r ----------------------
! -----------------------------------------------------------------------
