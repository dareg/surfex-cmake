!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE ASSIM_ISBA_UPDATE_SNOW (IO, S, NP, NPE, HPROGRAM, KI, PSWE, HTEST )

! ------------------------------------------------------------------------------------------
!  *****************************************************************************************
!
!  Routine to update snow field for ISBA
!  Trygve Aspelien, Separating IO  06/2013
!
!
! ******************************************************************************************
! ------------------------------------------------------------------------------------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_NP_t, ISBA_NPE_t, ISBA_S_t, ISBA_PE_t, ISBA_P_t
!
USE MODD_CSTS,        ONLY : XTT
USE MODD_ASSIM,       ONLY : LSWE
USE MODD_SURF_PAR,    ONLY : XUNDEF
USE MODD_SNOW_PAR,    ONLY : XANSMIN, XANSMAX, XRHOSMIN, XRHOSMAX
!
!
USE MODI_PACK_SAME_RANK
USE MODI_ABOR1_SFX
USE MODI_ASSIM_GATHER_WRITE_INCREMENTS
!
USE YOMHOOK,          ONLY : LHOOK,DR_HOOK
USE PARKIND1,         ONLY : JPRB
!
IMPLICIT NONE
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
TYPE(ISBA_NPE_t), INTENT(INOUT) :: NPE
!
CHARACTER(LEN=6),    INTENT(IN)    :: HPROGRAM  ! program calling surf. schemes
INTEGER,             INTENT(IN)    :: KI
REAL, DIMENSION(KI), INTENT(IN)    :: PSWE
CHARACTER(LEN=2),    INTENT(IN)    :: HTEST     ! must be equal to 'OK'
!
!    Declarations of local variables
!
TYPE(ISBA_P_t), POINTER :: PK
TYPE(ISBA_PE_t), POINTER :: PEK
!
REAL, DIMENSION(KI) :: ZSWEGP    ! Grid point average of 1. guess 

CHARACTER(LEN=2)                :: CP
REAL, ALLOCATABLE, DIMENSION(:) :: ZSWE       ! Patch value of updated snow
REAL, ALLOCATABLE, DIMENSION(:) :: ZSWE_ORIG  ! Patch value of 1. guess
REAL, ALLOCATABLE, DIMENSION(:) :: ZSWEINC 
REAL, ALLOCATABLE, DIMENSION(:) :: ZSWE_IN
REAL, ALLOCATABLE, DIMENSION(:) :: ZSWEGP_P
!    Addtional snow fields with D95 snow scheme 
REAL, ALLOCATABLE, DIMENSION(:) :: ZSNA
REAL, ALLOCATABLE, DIMENSION(:) :: ZSNR
INTEGER  :: JL,JP,JI,JII,IMASK
REAL :: ZSWEMIN,ZSWEMAX,ZSWEMEAN,ZINCMIN,ZINCMAX,ZINCMEAN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! ----------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_UPDATE_SNOW',0,ZHOOK_HANDLE)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_ISBA_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
IF ( NPE%AL(1)%TSNOW%SCHEME=='D95' ) THEN
  JL = 1
ELSE
  CALL ABOR1_SFX("Update of snow is only implemented for D95")
ENDIF


!
! 1. guess grid average of SWE
ZSWEGP(:)=0
DO JP=1,IO%NPATCH
  PK  => NP%AL(JP)
  PEK => NPE%AL(JP)
  DO JI=1,PK%NSIZE_P
    IF ( PEK%TSNOW%WSNOW(JI,JL) /= XUNDEF ) THEN
      ZSWEGP(PK%NR_P(JI))=ZSWEGP(PK%NR_P(JI))+(PEK%TSNOW%WSNOW(JI,JL)*S%XPATCH(PK%NR_P(JI),JP))
    ENDIF
  ENDDO
ENDDO

! Loop patches
DO JP=1,IO%NPATCH

  PK  => NP%AL(JP)
  PEK => NPE%AL(JP)

  ALLOCATE(ZSWE(PK%NSIZE_P))
  ALLOCATE(ZSWE_ORIG(PK%NSIZE_P))
  ALLOCATE(ZSWEINC(PK%NSIZE_P))
  ALLOCATE(ZSWE_IN(PK%NSIZE_P))
  ALLOCATE(ZSWEGP_P(PK%NSIZE_P))
  ALLOCATE(ZSNA(PK%NSIZE_P))
  ALLOCATE(ZSNR(PK%NSIZE_P))


  ZSWE=0
  ! Pack ISBA fields to this patch
  CALL PACK_SAME_RANK(PK%NR_P,ZSWEGP,ZSWEGP_P)
  CALL PACK_SAME_RANK(PK%NR_P,PSWE,ZSWE_IN)
  IF ( .NOT. LSWE ) THEN
    ! Convert PSWE from meter to SWE
    WHERE( PEK%TSNOW%RHO(:,JL) /= XUNDEF )
      ZSWE_IN(:)=ZSWE_IN(:)*PEK%TSNOW%RHO(:,JL)
    ELSEWHERE
      ZSWE_IN(:)=ZSWE_IN(:) * ( 0.5 * ( XRHOSMIN + XRHOSMAX ))
    ENDWHERE
  ENDIF

  ! Only modified snow on patch if it is defined
  WHERE ( PEK%TSNOW%WSNOW(:,JL) /= XUNDEF )
    ZSWE_ORIG(:)=PEK%TSNOW%WSNOW(:,JL)

    ! Analysed values weighted to keep the ratio of snow at different patches
    ! constant
    WHERE ( ZSWEGP_P(:) > 1.0E-10)
       ZSWE(:) = ZSWE_IN(:)*ZSWE_ORIG(:)/ZSWEGP_P(:)
    ELSE WHERE
       ZSWE(:) = ZSWE_IN(:)     ! If snow analysis gives initial snow, 
                                ! then same amount on all patches
    END WHERE
    ! Set snow=0 where 1. guess = 0 and Ts>0, to avoid that the snow analysis introduce snow where it is no snow.
    WHERE ( ZSWE_ORIG(:)<1.0E-10 .AND. PEK%XTG(:,1) >XTT )
       ZSWE(:)   = 0.0
    END WHERE

    ! Snow albedo and density are given initial values in points  
    ! which get initial snow in the snow analysis 
    ZSNA(:)=PEK%TSNOW%ALB(:)
    ZSNR(:)=PEK%TSNOW%RHO(:,JL)
    WHERE ( ZSWE_ORIG(:) < 1.0E-10 .AND. ZSWE(:)>= 1.0E-10 )
      ZSNA(:)    = 0.5 * ( XANSMIN + XANSMAX )
      ZSNR(:)    = 0.5 * ( XRHOSMIN + XRHOSMAX )
    END WHERE
  END WHERE

  WRITE(CP,'(I2)') JP
  CALL ASSIM_GATHER_WRITE_INCREMENTS(HPROGRAM,"WSNOW_VEG"//CP,PEK%TSNOW%WSNOW(:,JL),ZSWE,XUNDEF,LSTAT=.TRUE.)
  CALL ASSIM_GATHER_WRITE_INCREMENTS(HPROGRAM,"ASNOW_VEG"//CP,PEK%TSNOW%ALB(:),ZSNA,XUNDEF,LSTAT=.TRUE.)
  CALL ASSIM_GATHER_WRITE_INCREMENTS(HPROGRAM,"RSNOW_VEG"//CP,PEK%TSNOW%RHO(:,JL),ZSNR,XUNDEF,LSTAT=.TRUE.)

  ! Update snow/albedo/density
  WHERE ( PEK%TSNOW%WSNOW(:,JL) /= XUNDEF )
    PEK%TSNOW%ALB(:)=ZSNA(:)
    PEK%TSNOW%RHO(:,JL)=ZSNR(:)
    PEK%TSNOW%WSNOW(:,JL)=ZSWE(:)
  END WHERE

  DEALLOCATE(ZSWE)
  DEALLOCATE(ZSWE_ORIG)
  DEALLOCATE(ZSWEINC)
  DEALLOCATE(ZSWE_IN)
  DEALLOCATE(ZSWEGP_P)
  DEALLOCATE(ZSNA)
  DEALLOCATE(ZSNR)
ENDDO

!
! -------------------------------------------------------------------------------------
 IF (LHOOK) CALL DR_HOOK('ASSIM_ISBA_UPDATE_SNOW',1,ZHOOK_HANDLE)
 END SUBROUTINE ASSIM_ISBA_UPDATE_SNOW

