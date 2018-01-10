      MODULE MODI_FLXSURF3BX
      INTERFACE
      SUBROUTINE FLXSURF3BX(CMU, CTU, RIB, FTEMP, FVAP, ILMO,
     1                   UE, FCOR, TA , QA , ZU, ZT, VA,
     1                   TG , QG , H , Z0 , Z0T,
     1                   LZZ0, LZZ0T, FM, FH, N )
      INTEGER  N
      REAL  CMU(N),CTU(N),RIB(N),FCOR(N),ILMO(N)
      REAL  FTEMP(N),FVAP(N),TA(N),QA(N),ZU(N),ZT(N),VA(N)
      REAL  TG(N),QG(N),H(N),Z0(N),UE(N)
      REAL  Z0T(N),LZZ0(N),LZZ0T(N)
      REAL  fm(N),fh(N)
      END SUBROUTINE FLXSURF3BX
      END INTERFACE
      END MODULE MODI_FLXSURF3BX
