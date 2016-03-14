c
c---- parametres globaux
c
      include 'Glob_Collisions_SPH_inc.f'
c
c---- La grille formee par interpolation des resultats de simulations SPH
c
      include 'Grille_Collisions_SPH_inc.f'
c
c---- Donnees necessaires a l'utilisation de la grille dans le code
c     d'evolution d'amas
c

c     dmin_limite indique au dela de quel dmin/(R1+R2) on se permet d'appliquer
c     un traitement "rencontre distante" si la grille ne specifie pas
c     de resultat (pour l'instant, ce traitement est l'absence d'effets
c     collisionnels)
c
      real dmin_limite
      parameter (dmin_limite=0.6)
c
c     Signaux indiquant les cas/problemes lors de l'interpolation dans
c     la grille
c
      integer i_interp, i_extrap, i_inconnu, i_distant
      parameter (
     $     i_interp=i_valeur_interpolee, ! interpolation "normale"
     $     i_extrap=i_valeur_extrapolee, ! la valeur est basee sur des valeurs extrapolee (hautes/basses vitesses/ petit param d'impact)
     $     i_inconnu=i_valeur_non_specifiee, ! la valeur de la grille est non specifiee
     $     i_distant=5                   ! valeur non specifiee dans la grille mais grand parametre d'impact
     $     )
c
c     Limites en masses de la grille
c
      double precision Minf_grille_SPH, Msup_grille_SPH
c
c     Compteurs de cas speciaux
c
      integer n_dbrdM_CollSPH_b, n_dbrdM_CollSPH_h, n_extrap_CollSPH,
     $     n_distant_CollSPH, n_inconnu_CollSPH

                                ! Et un bloc common pour mettre tout ca dedans !
      common /Grille_SPH_common/ Minf_grille_SPH, Msup_grille_SPH,
     $     n_dbrdM_CollSPH_b, n_dbrdM_CollSPH_h, n_extrap_CollSPH,
     $     n_distant_CollSPH, n_inconnu_CollSPH
