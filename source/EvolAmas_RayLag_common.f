c======================================================================
c     commons pour EvolAmas_RayLag.F
c======================================================================
c
c constante :
c ^^^^^^^^^^^
      integer N_RayLagMax
      parameter (N_RayLagMax=100)
c
c common :
c ^^^^^^^^
      integer N_RayLag
      double precision M_RayLag(N_RayLagMax)
      common /RayFrac_Common/ M_RayLag,N_RayLag
