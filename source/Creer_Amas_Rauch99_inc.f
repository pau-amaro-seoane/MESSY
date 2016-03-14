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

      double precision expo_n, p_anis, a_min,a_max
      common /param_dist/ expo_n, p_anis, a_min,a_max

      integer iCorrec_SelfGrav
      common /control/ iCorrec_SelfGrav

                                ! Parameters for a "special" particle
      double precision m_special,a_special,e_special
      common /param_special/ m_special,a_special,e_special
