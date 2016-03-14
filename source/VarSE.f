c----------------------------------------------------------------------
c     Fichier de blocs common contenant les variables associees aux
c     super-etoiles pour simul. d'amas
c----------------------------------------------------------------------
c
c---- proprietes globales de l'amas
c
      double precision
     $     M_amas, B0
      integer
     $     NbSE, NbNoeuds
      common /Var_globales_Amas/
     $     M_amas,              ! masse totale
     $     B0,                  ! somme des -Mi/Ri
     $     NbSE,                ! Nb de SE, de SE "libres" 
     $     NbNoeuds             ! Nb de noeuds dans l''arbre binaire
c
c--- proprietes individuelles des SE
c
      double precision
     $     T_SE(iDimSE), J_SE(iDimSE), 
     $     M_SE(0:iDimArbre), R_SE(0:iDimArbre)
      common /VariablesSE/   
     $     T_SE,                ! Energie cinetique specifique
     $     J_SE,                ! Moment cinetique specifique
     $     M_SE,                ! Masse totale de la SE
     $     R_SE                 ! Rayon de la SE (couche spherique d *)
c
c     Les elements d'indice superieur a Nb_SE sont associes aux noeuds
c     vides de l'arbre (astuce technique)
c
c---- liste triee en rayon
c
      integer
     $     iListeRay
      common /Liste/
     $     iListeRay(0:iDimArbre+1) ! la procedure Tri_Avec_Arbre trie tous les
                                    ! noeuds => il faut que la liste contienne 
                                    ! suffisamment d`elements
c
c---- liste triee par rayon croissant des SE 
c     Avec les valeurs "sentinelles"
c     iListeRay(0)=0 et iListeRay(NbSE+1) = iNbSE+1 par convention
