c======================================================================
c     includes pour EvolAmas.F
c======================================================================
c
      implicit none
c
c     PARAMETRES
c     ==========
c
c---- Constantes math.
c
      include 'CstesMath.f'
c
c---- Constantes physiques
c
      include 'CstesPhys.f'
c
c---- Nb max de super-etoiles 'SE' et autres parametres globaux
c
      include 'Param_Amas.f'
c
c---- parametres propres a EvolAmas.F
c
      include 'EvolAmas_param.f'
c
c     COMMONS
c     =======
c
c---- Variables associees aux super-etoiles
c
      include 'VarSE.f'
c
c---- Comptablisation des SE perdues par evaporation et destruction
c
      include 'EvolAmas_Pertes_inc.f'
c
c---- Les temps et pas de temps
c
      include 'EvolAmas_Temps_common.f'
c
c---- Arbre binaire pour le potentiel
c
      include 'EvolAmas_Arbre_common.f'
c
c---- proprietes du TN
c
      include 'EvolAmas_TN_common.f'
c
c---- parametres de la simulation
c
      include 'EvolAmas_Param_common.f'
c
c---- Suivi de la conservation de la masse et de l'energie
c
      include 'EvolAmas_Conserv_inc.f'
c
c---- autres commons propres a EvolAmas.F
c
      include 'EvolAmas_common.f'

