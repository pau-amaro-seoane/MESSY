      implicit none

      include 'CstesMath.f'
      include 'Param_Amas.f'
      include 'VarSE.f'
      include 'EvolAmas_Pertes_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_Compteurs_common.f'
c
c---- La liste des etoiles suivies
c
      integer Dim_Liste_SEsuiv
      parameter (Dim_Liste_SEsuiv=100)
      integer N_SEsuiv, Nmax_SEsuiv, Liste_SEsuiv(Dim_Liste_SEsuiv)
      common /common_Liste_SEsuiv/ N_SEsuiv, Nmax_SEsuiv, Liste_SEsuiv
