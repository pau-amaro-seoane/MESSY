
      integer
     $     iSauv,iSauvPart,iNbRecGrille,iNbRecArbre,
     $     iNbPasRecRelax,i_NbTirVN_PosOrb,
     $     iNb_RejetPaire,iNb_AnnulRelaxEvap,iNb_AnnulRelaxCC,
     $     iNb_CalcPot, iNb_tests_LC, iNb_SuperRenc,
     $     iNb_SuperRenc_Pi16,
     $     iNb_SuperRenc_Pi8,
     $     iNb_SuperRenc_Pi4,
     $     iNb_SuperRenc_Pi2,
     $     iNb_Coll_P_025,
     $     iNb_Coll_P_05,
     $     iNb_Coll_P_1,
     $     iNb_Coll_P_2,
     $     iNb_Kick_P_025,
     $     iNb_Kick_P_05,
     $     iNb_Kick_P_1,
     $     iNb_Kick_P_2,
     $     iNb_Renc_dR_001,
     $     iNb_Renc_dR_01,
     $     iNb_Renc_dR_05,
     $     iNb_Renc_dR_1,
     $     iDuplic

      integer*8 ! integer*8 fonctionne avec pgf77 et g77 sous linux, f77 et g77 sur SUN
     $     iPas_Evol,
     $     iNb_Relax, iNb_Coll, iNb_Kick, iNb_Dechir, iNb_Disp,
     $     iNb_CaptGW

      double precision          ! pas vraiment des compteurs...
     $                          ! mais des variables de monitoring
     $     Dev_SuperRenc_moy,  Dev2_SuperRenc_moy,
     $     dR_rel_Renc_moy, dR_rel_Renc_max,
     $     P_coll_moy, P_coll_max,
     $     P_kick_moy, P_kick_max

      common /common_compteurs/
     $     iSauv,iSauvPart,iNbRecGrille,iNbRecArbre,
     $     iNbPasRecRelax,i_NbTirVN_PosOrb,
     $     iNb_RejetPaire,iNb_AnnulRelaxEvap,iNb_AnnulRelaxCC,
     $     iNb_CalcPot, iNb_tests_LC, iNb_SuperRenc,
     $     iNb_SuperRenc_Pi16,
     $     iNb_SuperRenc_Pi8,
     $     iNb_SuperRenc_Pi4,
     $     iNb_SuperRenc_Pi2,
     $     iNb_Coll_P_025,
     $     iNb_Coll_P_05,
     $     iNb_Coll_P_1,
     $     iNb_Coll_P_2,
     $     iNb_Kick_P_025,
     $     iNb_Kick_P_05,
     $     iNb_Kick_P_1,
     $     iNb_Kick_P_2,
     $     iNb_Renc_dR_001,
     $     iNb_Renc_dR_01,
     $     iNb_Renc_dR_05,
     $     iNb_Renc_dR_1,
     $     iDuplic

      common /common_compteurs8/
     $     iNb_Relax, iNb_Coll, iNb_Kick, iNb_Dechir, iNb_Disp,
     $     iNb_CaptGW,
     $     iPas_Evol

      common /common_monitor/
     $     Dev_SuperRenc_moy, Dev2_SuperRenc_moy,
     $     dR_rel_Renc_moy, dR_rel_Renc_max,
     $     P_coll_moy, P_coll_max,
     $     P_kick_moy, P_kick_max
