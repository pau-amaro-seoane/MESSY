c----------------------------------------------------------------------
c     Fichier de blocs common pour "AMAS/ArbreBinAmas.f"
c----------------------------------------------------------------------
c
c---- Structure "physique" de l'arbre binaire
c
      double precision
     $     dA_arbre(iDimArbre),dB_arbre(iDimArbre)
                                ! representation du potentiel
      integer 
     $     iSE_arbre(iDimArbre), ! num. de la SE associee au noeud
     $     iAdam,               ! adresse du noeud-racine
     $     iNoeud_SE(iDimArbre) ! fournit le noeud associe a une SE
      common /ArbrePhys/
     $     dA_arbre,dB_arbre,
     $     iNoeud_SE,iSE_arbre,iAdam
c
c---- Structure "logique" de l'arbre
c
      integer
     $     iPere_arbre(0:iDimArbre),
     $                          ! adresse du pere du noeud dans les
     $                          ! tableaux *_arbre
     $     iFilsG_arbre(0:iDimArbre),iFilsD_arbre(0:iDimArbre),
     $                          ! adresses des fils du noeud
     $     iNbDecG_arbre(0:iDimArbre)
                                ! nb de descendants a gauche (y compris
                                ! noeud en question)
      common /ArbreLog/
     $     iPere_arbre,iFilsD_arbre,iFilsG_arbre,iNbDecG_arbre
c
c---- Autres...
c
                                ! Tableau listant les noeuds liberes
      integer MaxNoeuds_Libres
      parameter (MaxNoeuds_Libres=iDimArbre)
      integer NbNoeuds_Libres,Noeuds_Libres(MaxNoeuds_Libres)
                                ! compteur du nb de parcours
      integer Nb_Parcours_Arbre

      common /ArbreDivers/
     $     NbNoeuds_Libres,Noeuds_Libres,
     $     Nb_Parcours_Arbre 


      character*(*) EnTeteFichArbre_XDR
      parameter (EnTeteFichArbre_XDR=
     $     '%%% XDR Cluster Binary Tree File %%%')
