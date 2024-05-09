!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######################
      MODULE MODD_XIAS
!     ######################
!
!!****  *MODD_XIAS - a shallow layer above XIOS 
!!
!!    PURPOSE
!!    -------
!!
!!    Allow to bufferize data blocks of a field before sending it to Xios, and to
!!    have a toggle (LXIAS) indicatiing if any component (arpege or SUrfex) activates Xios
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	S.Sénési   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       04/2019
!!
!
!*       0.   DECLARATIONS
!             ------------
!
!  Basic toggles
!
LOGICAL            :: LXIAS=.FALSE.                ! Do we use XIOS for outputing diags
!
!  Setup variables 
!
LOGICAL            :: LXIAS_INVERT_LEVELS=.FALSE.  ! should invert levels for axes which name includes 'klev' ?
INTEGER            :: NBLOCK                       ! Number of blocks to bufferize (NPROMA blocks when in Arpege)
!
!  Evolving variables
!
INTEGER            :: NTIMESTEP=0                  ! Last value of timestep sent to XIOS 
!                                                  ! (needed by xias_send_block); Xios first step is 1
!
END MODULE MODD_XIAS
