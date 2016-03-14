c======================================================================
c     commons pour EvolAmas_AniLag.F
c======================================================================
c
c constante :
c ^^^^^^^^^^^
      integer N_AniLagMax
      parameter (N_AniLagMax=10)
c
c common :
c ^^^^^^^^
      integer N_AniLag
      double precision M_AniLag(N_AniLagMax),dM_AniLag(N_AniLagMax)
      common /AniLag_Common/ M_AniLag,dM_AniLag,N_AniLag
