
      implicit none

      include 'Param_Amas.f'
      include 'CstesMath.f'
      include 'VarSE.f'

      include 'EvolAmas_Etoiles_common.f'
      include 'EvolAmas_Conserv_inc.f'
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_Pertes_inc.f'
c
c constante(s) :
c ^^^^^^^^^^^^^^
      double precision FMCmax
      parameter (FMCmax=1.0d0)
