c======================================================================
c     parametres pour EvolAmas_Relax.F
c======================================================================

      double precision
     $     a_spitzer,a_bibi,    ! < constante pour le tps de relax loc.
     $     b_spitzer            ! < constante pour le tps de relax a mi-masse
      parameter (
     $     a_spitzer=0.065d0, a_bibi=rPi/128.0d0,
     $     b_spitzer=0.138d0
     $     )

      double precision Coef_lambda
      parameter (Coef_lambda = 2.34d0) ! 2.34d0 pour avoir Lambda=gamma*N_etoiles dans plummer

      integer iFichRelax        ! pour le deboguage
      parameter (iFichRelax=30)
