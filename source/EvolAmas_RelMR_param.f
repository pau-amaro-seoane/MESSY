c======================================================================
c     parametres pour Relation Masse-Rayon
c======================================================================
c
c---- Parametres physiques (ou presque)
c
      double precision M_chandra  ! Pour la relation M-R des Naines blanches (masse de Chandrasekhar)
      parameter (M_chandra=1.44d0)
                                ! Pour la relation M-R des Etoiles de neutrons
      double precision Mmin_NS, Mmax_NS, R_NS
      parameter (
     $     Mmin_NS=0.1d0,       ! 
     $     Mmax_NS=3.0d0,       ! < Selon Fryer et Kalogera 01  
     $     R_NS=1.5d-5          ! < entre 10 et 11 km; en fait, en dessous de 0.2 Msol, le rayon augmente tres fortement
     $     )
                                ! Pour la "relation M-R" des TN stellaire
      double precision R_TNstell
      parameter (R_TNstell=1.0d-20) ! Traite en point-masse

                                ! M-R exponent for very massive stars from Goodman & Tan 03, Eq. 8 : 0.475d0
      double precision MR_exp_VMS
      parameter (MR_exp_VMS = 0.475d0) 
c
c---- Parametres numeriques
c
      integer NtabMR_Max
      parameter (NtabMR_Max=100) ! Dimension (max) du tableau de la relation MR
                                 ! pour la sequence principale
      integer Next
      parameter (Next=3)        ! nb de points pour calculer l'extrapolation
                                ! de la relation M-R poure les petites et
                                ! grandes masses
      integer N_RhR_max
      parameter (N_RhR_max=64)  ! Nb max de pts dans la relation M --> Rhalf/R
