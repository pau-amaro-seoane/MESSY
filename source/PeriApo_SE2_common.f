c
c---- proprietes de la SE
c
      double precision J22_SE,Ener_SE,Masse_SE,Ray_SE
      common /Prop_SE/
     $     J22_SE,Ener_SE,Masse_SE,Ray_SE
c
c---- encadrement des distances peri- et apo-centriques
c      
      integer Noeud_inf,Noeud_sup,Rang_inf,Rang_sup
      double precision R_inf,R_sup,A_inf,A_sup,B_inf,B_sup,Q_inf,Q_sup

      common /Encadrement_PeriApo/
     $     R_inf,R_sup,A_inf,A_sup,B_inf,B_sup,Q_inf,Q_sup,
     $     Noeud_inf,Noeud_sup,Rang_inf,Rang_sup
c
c---- proprietes du noeud courant lors du parcours
c     de l'arbre binaire
c
      integer Noeud_loc,Rang_loc
      double precision R_loc,A_loc,B_loc

      common /Noeud_courant/
     $     R_loc,A_loc,B_loc,Noeud_loc,Rang_loc
