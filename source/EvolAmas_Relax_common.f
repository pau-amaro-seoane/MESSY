c======================================================================
c     blocs common pour EvolAmas_Relax.F
c======================================================================

      double precision
     $     Log_Coulomb,           ! < Ln(Gamma*Netoiles)
     $     Cste_Relax,            ! < 32/Pi*Ln(Gamma*Netoiles)
     $     Correc_Log_Coul_kicks, ! < est utile si on traites explicitement les deviations a grands angles
     $     Met_moy,               ! < Masse moyenne des etoiles (en unites de masse de l'amas)
     $     TSurTrel               ! < tps evol / Tps de relax a mi-masse

      common /Relaxation_common/
     $     Log_Coulomb, Cste_Relax, Correc_Log_Coul_kicks, Met_moy,
     $     TSurTrel
c
c---- bloc common de proprietes globales pour amas, liees  a la
c     relaxation (calculees par CalcGlobRelax)
c
      double precision
     $     Rh,Rc,TrelMoy,Trel_dt,Trel_Rh,Trel_Rh0,Trel_Rh_spit,Trel_c,
     $     V2c,V2moy,Rhoc,q_Viriel,M_renc_tot,Trel_CM_moy,Trel_tst_moy
      common /BC_GlobRelax/
     $     Rh,Rc,               ! rayon de demi-masse, de coeur
     $     TrelMoy,Trel_dt,     ! tps de relax moyen, trel par moyenne des Trel_SE
     $     Trel_Rh,Trel_Rh_spit, ! tps de relax a demi-rayon (2 definitions)
     $     Trel_Rh0,            ! valeur initiale de Trel_Rh
     $     Trel_c,              ! tps de relax au centre
     $     V2c,V2moy,           ! vitesse quad au centre, moyenne
     $     Rhoc,                ! Densite centrale
     $     q_Viriel,            ! Rh*V2moy/M_amas
     $     M_renc_tot,
     $     Trel_CM_moy,
     $     Trel_tst_moy
