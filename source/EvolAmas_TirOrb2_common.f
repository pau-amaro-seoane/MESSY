c
c---- proprietes utiles de la particule
c
      double precision
     $     M_TO,E_TO,J_TO,Rp_TO,Ra_TO,ecc_TO,UnSurDR_TO
c
c---- tableaux permettant la representation d'une borne sup pour dP/dR
c
      integer Nmax_borne_TO
      parameter (Nmax_borne_TO=100)

      integer N_borne_TO
      double precision          ! borne_P_TH_TO contient les valeurs des probabilite
     $                          ! de selection (P_TH) a la limite inf de chaque intervalle
     $                          ! Bint_TO contient les valeur de la fonction-borne integree jusqu'a
     $                          ! la limite sup de chaque intervalle
     $     borne_P_TH_TO(Nmax_borne_TO),Bint_TO(Nmax_borne_TO),
     $     dlR_TO
c
c---- Autres
c
      integer iCompt_essais_TO ! compteur du nb d'essais

      common /TirOrb_common/
     $     M_TO,E_TO,J_TO,Rp_TO,Ra_TO,ecc_TO,UnSurDR_TO,
     $     borne_P_TH_TO,Bint_TO,dlR_TO,N_borne_TO,
     $     iCompt_essais_TO
