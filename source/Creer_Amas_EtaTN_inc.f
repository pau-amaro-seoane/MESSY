      implicit none

      include 'CstesMath.f'
      include 'Param_Amas.f'
      include 'VarSE.f'
      include 'EvolAmas_TN_common.f'

      double precision acc_R
      parameter (acc_R=1.0d-9)

      character*(80) NomAmas
      common /NomsFich/ NomAmas

      double precision Ecin, Epot_st, Epot_bh
      common /energies/ Ecin, Epot_st, Epot_bh

      double precision Rinfl,Rbreak,Sigma1D_break,
     $     Rmin_cutoff
      integer iseed
      common /various/ Rinfl,Rbreak,Sigma1D_break,
     $     Rmin_cutoff,
     $     iseed
