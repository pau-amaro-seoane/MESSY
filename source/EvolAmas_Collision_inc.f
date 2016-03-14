c
c constante(s) :
c ^^^^^^^^^^^^^^
      double precision FactRmax_Coll
      parameter (FactRmax_Coll=1.0d0)
c
c bloc common :
c ^^^^^^^^^^^^^
      double precision
     $     Vesc_coll,           ! Vitesse d'evasion (unites du code)
     $     Seff_coll,           ! Section efficace pour collisions physiques (unites du code)
     $     Seff_kick,           ! Section efficace pour deviation avec grand angle (unites du code)
     $     M_accr_coll,         ! Masse accretee sur le TN de provenance collisionnelle
     $     M_ejec_coll          ! Masse ejecteee de provenance collisionnelle

      common /Collision_common/ Vesc_coll,Seff_coll,Seff_kick,
     $     M_accr_coll,M_ejec_coll
