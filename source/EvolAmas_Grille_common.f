c======================================================================
c     tableaux permettant une representation grossiere de la densite
c     pour EvolAmasXXX
c======================================================================
c

      double precision
     $     M_PG(idimPG),        ! < masse dans chaque couche
     $     V_PG(idimPG),        ! < volume de chaque couche
     $     T_PG(idimPG),        ! < energie cinetique de chaque couche
     $     Tt_PG(idimPG),       ! < energie cinetique tangentielle de chaque couche
     $     MR_PG(idimPG),       ! < somme des mi*Ri dans chaque cellule
     $     MM_PG(idimPG),       ! < somme des mi*mi dans chaque cellule
     $     Rsup_PG(0:idimPG),   ! < limite sup de chaque cellule
     $     rNSEc_PG             ! < nb moyen de SE par cellule
      integer
     $     Nmin_PG,Nmax_PG,     ! < Nb extremes tolerables de SE par cellule
     $     N_PG,                ! < nb d`elements des tableaux ci-dessus
     $     NSE_PG(idimPG),      ! < nb de SE dans chaque couche du PG
     $     NSEc_PG,             ! < nb demande de SE par cellule
     $     Npart_PG             ! < total number of particles in grid

      common /Common_Grille/
     $     M_PG, V_PG, T_PG, Tt_PG, MR_PG, MM_PG, Rsup_PG, rNSEc_PG,
     $     NSE_PG,Nmin_PG,Nmax_PG,N_PG, NSEc_PG, Npart_PG
     
c$$$      double precision Y3c_PG(idimPG)
c$$$      common /Common_Grille_Binaires/ Y3c_PG 
