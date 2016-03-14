c
c---- La grille de donnees interpolee dans les resultats de simulations SPH de collisions
c     stellaires
c
      real qualite_lim
      parameter (qualite_lim=1.0e-4) ! valeur limite indiquant si on se trouve
                                     !  dans l'enveloppe convexe des donnees
                                     ! ou non (pour triangulation de Delaunay)

      real valeur_non_specifiee ! indique qu'un element de la grille est "vierge"
      parameter (valeur_non_specifiee=1.0e30)

      integer Ntab_interp_max
      parameter (Ntab_interp_max=40**Nvar)
      real
     $     tab_interp(Ntab_interp_max,Nfct), ! tableau contenant les donnees interpolees
     $     Xmin_interp(Nvar),   ! limites inferieures du tableau dans les Nvar dimensions
     $     Xmax_interp(Nvar)   ! limites superieures du tableau dans les Nvar dimensions

      integer*1 i_valeur_non_specifiee, i_valeur_specifiee,
     $     i_valeur_interpolee, i_valeur_extrapolee
      parameter (
     $     i_valeur_non_specifiee=0,
     $     i_valeur_specifiee=1,
     $     i_valeur_extrapolee=2,
     $     i_valeur_interpolee=10 )
      integer*1
     $     tab_stat_interp(Ntab_interp_max,Nfct) ! specification de la qualite de l'interpolation (en reserve)
      integer
     $     Nvois_interp,        ! Nb de voisins utilises pour l'interpolation par noyaux
     $     N_interp(Nvar)       ! Nb d'elements de la grille d'interpolation dans les Nvar dimensions
c
c     Comme le nb de variables est laisse libre, le tableau d'interpolation
c     ne les distingue pas explicitement. Pour acceder a l'element (i1,i2,i3,i4)
c     de la fct ifct, correspondant aux valeurs x1=Xmin_interp(1)+(i1-1)*dX(1)
c     avec dX(1)=(Xmax_interp(1)-Xmin_interp(1))/(N_interp(1)-1), x2=..., x3=..., x4=...,
c     il faut employer l'indice ivar=(N1*N2*N3)*(i4-1)+(N1*N2)*(i3-1)+N1*(i2-1)+(i1-1) +1
c     dans tab_interp(ivar,ifct) avec N1=N_interp(1), etc...
c
      common tab_interp,Xmin_interp,Xmax_interp,
     $     Nvois_interp,N_interp,
     $     tab_stat_interp
      
