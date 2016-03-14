c======================================================================
c     blocs common des temps pour un amas a la Henon
c======================================================================

      character*(*) EnTeteFichTemps_XDR
      parameter (EnTeteFichTemps_XDR='%%% XDR Cluster Times File %%%')

      double precision
     $     Tps_SE(iDimSE),      ! `temps propres` de chaque SE
     $     Tps_amas,            ! tps de l'amas (temps median)
     $     dTps_amas_inf,dTps_amas_sup, ! evaluation de la dispertion des tps (1/6, 5/6 des particules)
     $     Tps_amas_moy, Sigma_Tps_amas ! tps moyen, dispertion

      common /Temps_SE/
     $     Tps_SE,
     $     Tps_amas, dTps_amas_inf,dTps_amas_sup,
     $     Tps_amas_moy, Sigma_Tps_amas
