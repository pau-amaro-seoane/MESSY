c======================================================================
c     2-dim array containing tabulated dynamical friction coef
c     for eta model.
c
c     TabDFCoef(j,i) contains the log10 of the omega(R_i,nu_j)
c
c     where R_i = Rmin*10**((i-1)*dlR) is the radius in N-body units
c           nu_j = nu_min*10**((j-1)*dnu) is the velocity in units of
c           the local 1D velocity dispersion
c
c     omega defined through:
c     dV/dt = -4*pi*ln(Lambda)*G^2*Rho*M*(1+m/M)/V^2 * 
c              omega(R/Rb,V/sigma_1D) * 1/V^2
c======================================================================
c
      character*256 PathFile_DFCoef
      common /PathTabDFCoef_common/ PathFile_DFCoef

      integer N_DFCoef_R, N_DFCoef_nu
      parameter (N_DFCoef_R=100, N_DFCoef_nu=100)
      double precision
     $     R_DFCoef_min, R_DFCoef_max,
     $     nu_DFCoef_min, nu_DFCoef_max 
      parameter (
     $     R_DFCoef_min=1d-5, R_DFCoef_max=1d3,
     $     nu_DFCoef_min=1d-2, nu_DFCoef_max=10d0 )

      double precision dlR_DFCoef, dlnu_DFCoef,
     $     lR_DFCoef_min, lnu_DFCoef_min
      double precision TabDFCoef(N_DFCoef_nu, N_DFCoef_R),
     $     TabSig1D(N_DFCoef_R)
      logical lTabDFCoef
      common /TabDFCoef_common/
     $     dlR_DFCoef, dlnu_DFCoef,
     $     lR_DFCoef_min, lnu_DFCoef_min,
     $     TabDFCoef, TabSig1D,
     $     lTabDFCoef
