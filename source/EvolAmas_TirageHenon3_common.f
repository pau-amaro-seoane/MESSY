      integer Nmax_TH
      parameter (Nmax_TH=50)    ! Nb de tranches radiales employees pour la
                                ! representation de la probabilite de tirage
      integer N_TH,NbSE_par_tranche_TH
      double precision dt_TH,Pr_TH(Nmax_TH),Pr_int_TH(Nmax_TH)
      common /common_Tirage_Henon/ dt_TH,Pr_TH,Pr_int_TH,N_TH,
     $     NbSE_par_tranche_TH
                                ! dt_TH est la constante de prop entre pas de tps Dt 
                                ! et proba de tirage : P(j)=dt_TH/Dt(j)
                                ! On montre qu'il s'agit du pas de temps moyen d'avancee
                                ! de l'amas par pas (un pas=une evolution de SE)
