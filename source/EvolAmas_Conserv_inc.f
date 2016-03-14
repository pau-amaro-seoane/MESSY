c======================================================================
c     blocs commons permettant le suivi des grandeurs "conservees"
c     lors de l'evolution d'un amas : Masse et energie.
c======================================================================
                                ! Energies totales
      double precision
     $     Egrav_amas,Ecin_amas,Etot_amas,Estell_grav_amas
      common /Eamas_common/
     $     Egrav_amas,Ecin_amas,Etot_amas,Estell_grav_amas

      double precision          ! Variations dues aux processus "non conservatifs"
     $     dEtot_evap,dEtot_coll,dEtot_dechir,dEtot_disp,dEtot_captGW,
     $     dEtot_evst,dEtot_bin3c,dEtot_nk,
     $     dTtot_evap,dTtot_coll,dTtot_dechir,dTtot_disp,dTtot_captGW,
     $     dTtot_evst,dTtot_nk,
     $     dMejec_evap,dMejec_coll,dMejec_dechir,dMejec_evst,
     $     dMaccr_coll,dMaccr_dechir,dMaccr_disp,dMaccr_captGW,
     $     dMaccr_evst
      double precision          ! Valeurs initiales 
     $     Mamas_ini,Mtn_ini,Etot_ini,Ttot_ini,Estell_grav_ini,
     $     MSEtot_ini,
     $     rNbEtoiles_ini
      integer NbSE_ini

      common /Conserv_common/
     $     Mamas_ini,Mtn_ini,Etot_ini,Ttot_ini,Estell_grav_ini,
     $     MSEtot_ini,
     $     rNbEtoiles_ini,
     $     dEtot_evap,dEtot_coll,dEtot_dechir,dEtot_disp,dEtot_captGW,
     $     dEtot_evst,dEtot_bin3c,dEtot_nk,
     $     dTtot_evap,dTtot_coll,dTtot_dechir,dTtot_disp,dTtot_captGW,
     $     dTtot_evst,dTtot_nk,
     $     dMejec_evap,dMejec_coll,dMejec_dechir,dMejec_evst,
     $     dMaccr_coll,dMaccr_dechir,dMaccr_disp,dMaccr_captGW,
     $     dMaccr_evst,
     $     NbSE_ini

      character*(*) EnTeteFichCons_XDR
      parameter (EnTeteFichCons_XDR=
     $     '%%% XDR Cluster Conservation File %%%')
