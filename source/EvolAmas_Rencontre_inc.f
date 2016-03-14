c======================================================================
c     parametres de la rencontre a deux particules (collision/relaxation)
c     lors d'un pas de l'evolution d'un amas
c======================================================================
c
c---- le pas de temps
c
      double precision dt_renc
c
c---- les proprietes liees aux particules
c
      double precision
     $     Mu1,Mu2,             ! M1/(M1+M2),M2/(M1+M2)
     $     V1amas(3),V2amas(3), ! Vitesses dans referentiel de l'amas (modifiees lors rencontre)
     $     Vrel(3),Vcm_renc(3), ! Vitesse relative et vitesse du CM (pas modifiees)
     $     Vrel_n,              ! Norme de la vitesse relative
     $     beta,                ! orientation du parametre d'impact
     $     b0_renc              ! G(m1+m2)/vrel**2
c
c---- les proprietes locales de l'amas
c
      double precision
     $     Psi0,                ! Le potentiel (stellaire) central (positif)
     $     R_renc,              ! Le rayon de la rencontre
     $     Met_moy_renc,        ! masse stellaire moyenne
     $     Rho_renc,Dens_renc,  ! Densite de masse, numerique d'etoiles
     $     Sigma1D_renc         ! Dispersion locale des vitesses
      integer iPG_renc ! la position de la rencontre dans la grille PG
c
c---- le type de rencontre
c
      integer RencIndef, RencPassive, Relaxation, Collision, Kick,
     $     RencUnbound
      parameter (RencIndef=-1, RencPassive=0, Relaxation=1, Collision=2,
     $     Kick=3, RencUnbound=99)

      integer iTypeRenc
c
c---- Certains resultats de la rencontre
c
      double precision Theta_relax_renc,Theta_coll_renc,Theta_kick_renc


      common /Rencontre_common/
     $     dt_renc,
     $     Mu1,Mu2,V1amas,V2amas,Vrel,Vcm_renc,Vrel_n,beta,b0_renc,
     $     Psi0,R_renc,Met_moy_renc,Rho_renc,Dens_renc,Sigma1D_renc,
     $     iPG_renc,
     $     iTypeRenc
      common /Resultats_Rencontre_common/
     $     Theta_relax_renc,Theta_coll_renc,Theta_kick_renc
