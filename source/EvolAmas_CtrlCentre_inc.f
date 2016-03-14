c======================================================================
c     controle du centre de l'amas
c======================================================================

      double precision  Fact_rang_CC ! quelle fraction des SE 
                                     ! represente-elle le "noyau" controle ?
      parameter (Fact_Rang_CC=0.001d0)

      double precision Rmin_CC
      integer iSE_Rmin_CC
      common /common_CC/ Rmin_CC, iSE_Rmin_CC
