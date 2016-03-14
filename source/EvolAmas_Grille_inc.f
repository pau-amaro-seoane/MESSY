c======================================================================
c     includes pour EvolAmasPG.f
c======================================================================
c
      implicit none
c
c---- Constantes math.
c
      include 'CstesMath.f'
      double precision r4Pi3
      parameter (r4Pi3=1.3333333d0*rPi)

c
c---- Nb max de super-etoiles 'SE' et autres parametres globaux
c
      include 'Param_Amas.f'
c
c---- Les proprietes des SE
c
      include 'VarSE.f'
c
c---- Donnees concernant les SE perdues
c
      include 'EvolAmas_Pertes_inc.f'
c
c---- Les tableaux contenant le profil de densite
c
      include 'EvolAmas_Grille_param.f'
      include 'EvolAmas_Grille_common.f'
c
