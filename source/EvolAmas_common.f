c======================================================================
c     blocs common pour EvolAmas.F
c======================================================================
c
c---- compteurs
c
      include 'EvolAmas_Compteurs_common.f'
c
c---- Facteurs de conversion fixant les dimensions
c
      include 'EvolAmas_Dimensions_common.f'
c
c---- flags et autres signaux
c
      integer
     $     flag_Tri,flag_Grille,flag_Etoiles,
     $     flag_Temps,flag_Arbre,
     $     iexit_stat           ! status de sortie
      common /common_EvolAmas_flags/
     $     flag_Tri,flag_Grille,flag_Etoiles,
     $     flag_Temps,flag_Arbre,
     $     iexit_stat
