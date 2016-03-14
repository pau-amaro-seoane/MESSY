c======================================================================
c     commons pour EvolAmas_SuiviSegr.F
c======================================================================
c
c constante :
c ^^^^^^^^^^^
      integer N_SuiviSegrMax
      parameter (N_SuiviSegrMax=100)
c
c common :
c ^^^^^^^^
      integer N_SuiviSegr
      double precision M_SuiviSegr(N_SuiviSegrMax)
      common /SuiviSegr_Common/ M_SuiviSegr,N_SuiviSegr
