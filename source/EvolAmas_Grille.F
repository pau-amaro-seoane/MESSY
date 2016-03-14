c======================================================================
c     routines de calcul des profils "grossiers" pour EvolAmasXXXX
c     Version a grille "lagrangienne"
c======================================================================
c
#define _ECRIRE_PG_ 0
#define _TEST_PG_ 0
#define _VERB_PG_ 0
#define _EXTRAP_SIG_ 0 /* 1 */
#define _EXTRAP_RHO_ 0 /* 1 */
c
c----------------------------------------------------------------------
      subroutine Construire_Grille_PG (NSE_couche)
c----------------------------------------------------------------------
c     Construction des "profils grossiers" en placant NSE_couche Se par 
c     cellule de la grille
c     Applique egalement les correction de vitesses avant de construire 
c     la grille
c     LES SE DOIVENT ETRES TRIEES EN R_SE (tableau iListeRay)
c
      include 'EvolAmas_Grille_inc.f' 
      include 'EvolAmas_Param_common.f' ! pour FacMinPG et FacMaxPG
     
#if (_ECRIRE_PG_>0)
      integer iFich_PG
      parameter (iFich_PG=65)
      character*80 NomFich_PG
      integer Num_Ecrit_PG
      data Num_Ecrit_PG /0/
#endif
c
c argument :
c ^^^^^^^^^^
      integer NSE_couche
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer iSE, lSE, lMax, n_couche
      double precision
     $     Ri,Ri1,Ri_3,Ri1_3,
     $     m_couche, t_couche, tt_couche, mr_couche, mm_couche,
     $     b_couche, y_couche,
     $     MSE,RSE
c
c fonction(s) :
c ^^^^^^^^^^^^^
      logical tester_PG
c......................................................................
c
#if (_VERB_PG_>=5)
      write(0,'("|| Construction grille PG ...",$)')
      call flush_MD(0)
#endif

      N_PG = MIN(idimPG,(NbSE-NbSE_perte)/NSE_couche) ! nb de cellules
      NSEc_PG = NSE_couche
      rNSEc_PG = dble(NbSE-NbSE_perte)/dble(N_PG)   ! nb moyen de SE par cellule
      Nmin_PG = NINT(FacMinPG*rNSEc_PG)
      Nmax_PG = NINT(FacMaxPG*rNSEc_PG)
c
c---- initialisations
c     
      RSup_PG(0) = 0.5d0*R_SE(iListeRay(1))
      N_PG = 1
      n_couche  = 0
      m_couche  = 0.0d0
      t_couche  = 0.0d0
      tt_couche  = 0.0d0
      mr_couche = 0.0d0
      mm_couche = 0.0d0
      b_couche  = 0.0d0
      y_couche  = 0.0d0
      Ri1   = 0.0d0
      Ri1_3 = 0.0d0
      lMax = NINT(rNSEc_PG)
      RSE = 0.0d0

      do lSE=1,NbSE-NbSE_perte  ! on passe en revue les SE dans l`ordre
         iSE = iListeRay(lSE)   ! de leur distance au centre de l`amas

         if (lSE.GT.lMax)
     $        then              ! on vient de depasser le seuil de remplissage
                                ! de la cellule courante
            NSE_PG(N_PG) = n_couche
            M_PG  (N_PG) = m_couche
            T_PG  (N_PG) = t_couche
            Tt_PG (N_PG) = tt_couche
            MR_PG (N_PG) = mr_couche
            MM_PG (N_PG) = mm_couche
c$$$            B_PG (N_PG)  = b_couche
c$$$            Y3c_PG(N_PG) = y_couche
            Ri = 0.5d0*(RSE+R_SE(iSE)) ! Ri  : limite sup de la couche
                                       ! Ri1 : limite inf de la couche
            Rsup_PG(N_PG) = Ri
            Ri_3 = Ri**3
            V_PG(N_PG) = r4Pi3*(Ri_3-Ri1_3) 
            Ri1 = Ri
            Ri1_3 = Ri_3
            N_PG = N_PG+1       ! N_PG est le numero de la couche en "construction"
            n_couche  = 0
            m_couche  = 0.0d0
            t_couche  = 0.0d0
            tt_couche  = 0.0d0
            mr_couche = 0.0d0
            mm_couche = 0.0d0
            b_couche  = 0.0d0
            y_couche  = 0.0d0
            lMax = NINT(N_PG*rNSEc_PG)
         end if

         RSE = R_SE(iSE)
         MSE = M_SE(iSE) 
         n_couche  = n_couche  + 1
         m_couche  = m_couche  + MSE
         t_couche  = t_couche  + MSE*T_SE(iSE)
         tt_couche = tt_couche + 0.5d0*MSE*(J_SE(iSE)/R_SE(iSE))**2
         mr_couche = mr_couche + MSE*RSE
         mm_couche = mm_couche + MSE*MSE
         b_couche  = b_couche  - MSE/RSE
c$$$#if (_CHAUFF_BIN_3C_>0)
c$$$         y_couche  = y_couche + Met_SE(iSE)**2/sqrt(2.0d0*T_SE(iSE))**3
c$$$#endif
      end do
                                ! derniere cellule
      NSE_PG(N_PG) = n_couche
      M_PG  (N_PG) = m_couche
      T_PG  (N_PG) = t_couche
      Tt_PG (N_PG) = tt_couche
      MR_PG (N_PG) = mr_couche
      MM_PG (N_PG) = mm_couche
c$$$      B_PG (N_PG)  = b_couche
c$$$      Y3c_PG(N_PG) = y_couche
      Ri = 1.5d0*RSE -
     $     0.5d0*R_SE(iListeRay(NbSE-NbSE_perte-1))
      Rsup_PG(N_PG) = Ri
      Ri_3 = Ri**3
      V_PG(N_PG) = r4Pi3*(Ri_3-Ri1_3) 

#if (_VERB_PG_>=5)
      write(0,200) N_PG
 200  format(' Termine ! (',I4,' elements)')
      call flush_MD(0)
#endif
#if (_ECRIRE_PG_>0)
      Num_Ecrit_PG = Num_Ecrit_PG + 1
      write (NomFich_PG,'("Grille_PG_",I4.4,)') Num_Ecrit_PG
      open (iFich_PG,file=NomFich_PG)
      call EcrireGrilleASC(iFich_PG)
      close (iFich_PG)
#endif
#if (_TEST_PG_>0)
      if (.NOT.Tester_PG(iPG)) then
         write (0,*) '!!! PG incorrect :',iPG,' !!!'
         call Meurt('Construire_Grille_PG')
      end if
#endif
c
      end
c
c----------------------------------------------------------------------
      subroutine Construire_Grille_PG2 (Npc_opt,Npc_min,Npc_max)
c----------------------------------------------------------------------
c     Construction des "profils grossiers" en placant Npc_opt Se par 
c     cellule de la grille
c     Applique egalement les correction de vitesses avant de construire 
c     la grille
c     LES SE DOIVENT ETRES TRIEES EN R_SE (tableau iListeRay)
c
      include 'EvolAmas_Grille_inc.f' 
c
c argument :
c ^^^^^^^^^^
      integer Npc_opt,Npc_min,Npc_max
c
c constantes :
c ^^^^^^^^^^^^
      double precision factRmax,factRmin
      parameter (factRmax=2.0d0,factRmin=0.5d0)
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer iSE, lSE, n_couche,
     $     Npc_opt_,Npc_min_,Npc_max_
      double precision
     $     Rmin,Rmax,Rlim,fact_R_opt,
     $     Ri,Ri1,Ri_3,Ri1_3,
     $     m_couche, t_couche,  tt_couche, mr_couche, mm_couche,
     $     MSE,RSE
c
c fonction(s) :
c ^^^^^^^^^^^^^
      logical tester_PG
c......................................................................
c
#if (_VERB_PG_>=5)
      write(0,'("|| Construction grille PG ...",$)')
      call flush_MD(0)
#endif
c
c---- calcul des parametres
c
      Npc_opt_ = Npc_opt
      Npc_min_ = min(Npc_min,Npc_opt)
      Npc_max_ = max(Npc_max,Npc_opt)
      Rmin = factRmin*R_SE(iListeRay(1))
      Rmax = factRmax*R_SE(iListeRay(NbSE-NbSE_perte))
      fact_R_opt = exp(log(Rmax/Rmin)*Npc_opt_/(NbSE-NbSE_perte))
c
c---- initialisations
c     
      Npart_PG = 0
      n_couche  = 0
      m_couche  = 0.0d0
      t_couche  = 0.0d0
      tt_couche = 0.0d0
      mr_couche = 0.0d0
      mm_couche = 0.0d0
      Ri1   = 0.0d0
      Ri1_3 = 0.0d0
      RSE = 0.0d0
      Rsup_PG(0) = Rmin
      Rlim = Rsup_PG(0)*fact_R_opt
      N_PG = 1                  ! Numero de la cellule en cours de remplissage
c
c---- "remplissage" successif des elements de la grille
c
      do lSE=1,NbSE-NbSE_perte  ! on passe en revue les SE dans l`ordre
         iSE = iListeRay(lSE)   ! de leur distance au centre de l`amas

                                ! En a-t-on termine avec la cellule courante ?
         if (n_couche.GE.Npc_min)
     $        then              ! Il faut que la cellule contienne au moins Npc_min particules
            Ri = 0.5d0*(RSE+R_SE(iSE)) ! Ri  : limite sup de la couche
                                       ! Ri1 : limite inf de la couche
            if ( (Ri.GE.Rlim) .OR. (n_couche.GE.Npc_max) )
     $           then           ! La cellule est complete si elle a la bonne taille ou 
                                ! si elle contient le nb max de particules
   
               NSE_PG(N_PG) = n_couche
               M_PG  (N_PG) = m_couche
               T_PG  (N_PG) = t_couche
               Tt_PG (N_PG) = tt_couche
               MR_PG (N_PG) = mr_couche
               MM_PG (N_PG) = mm_couche
               Npart_PG = Npart_PG+n_couche
               Rsup_PG(N_PG) = Ri
               Ri_3 = Ri**3
               V_PG(N_PG) = r4Pi3*(Ri_3-Ri1_3) 
               Ri1 = Ri
               Ri1_3 = Ri_3
               Rlim = Ri*fact_R_opt
               n_couche  = 0
               m_couche  = 0.0d0
               t_couche  = 0.0d0
               tt_couche = 0.0d0
               mr_couche = 0.0d0
               mm_couche = 0.0d0
   
               N_PG = N_PG + 1
               if (N_PG.GT.idimPG) then
                  write(0,*) '!!! Nb d''elements max de la grille (',
     $                 idimPG,') PG atteint !!!'
                  call Meurt('Construire_Grille_PG2')
               end if

            end if
         end if

         RSE = R_SE(iSE)
         MSE = M_SE(iSE) 
         n_couche  = n_couche  + 1
         m_couche  = m_couche  + MSE
         t_couche  = t_couche  + MSE*T_SE(iSE)
         tt_couche = tt_couche + 0.5d0*MSE*(J_SE(iSE)/R_SE(iSE))**2
         mr_couche = mr_couche + MSE*RSE
         mm_couche = mm_couche + MSE*MSE
      end do

      if (n_couche.GT.0) then   ! derniere cellule
         NSE_PG(N_PG) = n_couche
         M_PG  (N_PG) = m_couche
         T_PG  (N_PG) = t_couche
         Tt_PG (N_PG) = tt_couche
         MR_PG (N_PG) = mr_couche
         MM_PG (N_PG) = mm_couche
         Ri = 1.5d0*RSE -
     $        0.5d0*R_SE(iListeRay(NbSE-NbSE_perte-1))
         Rsup_PG(N_PG) = Ri
         Ri_3 = Ri**3
         V_PG(N_PG) = r4Pi3*(Ri_3-Ri1_3) 
      end if

#if (_VERB_PG_>=5)
      write(0,200) N_PG
 200  format(' Termine ! (',I4,' elements)')
      call flush_MD(0)
#endif
#if (_TEST_PG_>0)
      if (.NOT.Tester_PG(iPG)) then
         write (0,*) '!!! PG incorrect :',iPG,' !!!'
         call Meurt('Construire_Grille_PG')
      end if
#endif
c
      end
c
c----------------------------------------------------------------------
      subroutine Ajouter_a_Grille_PG2
c----------------------------------------------------------------------
c     Add the current cluster data into the current radial grid 
c     (built with Construire_Grille_PG2). To be used to average profiles
c     over successive snapshots. But no fancy mass/energy correction.
c
      include 'EvolAmas_Grille_inc.f' 
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer iSE, lSE, n_couche, N_grid, lSE_min
      double precision
     $     Rlim,RSE,MSE,
     $     m_couche, t_couche,  tt_couche, mr_couche, mm_couche
c
c fonction(s) :
c ^^^^^^^^^^^^^
      logical tester_PG
c......................................................................
c
#if (_VERB_PG_>=5)
      write(0,'("|| Add cluster data to PG grid ...",$)')
      call flush_MD(0)
#endif

c
c---- initialisations
c     
      n_couche  = 0
      m_couche  = 0.0d0
      t_couche  = 0.0d0
      tt_couche = 0.0d0
      mr_couche = 0.0d0
      mm_couche = 0.0d0

      N_grid = 1                  ! Numero de la cellule en cours de remplissage

      Rlim = Rsup_PG(0)
      do lSE_min=1,NbSE-NbSE_perte
         iSE = iListeRay(lSE_min)
         if (R_SE(iSE).GE.Rlim) goto 98
      end do
 98   continue
c
c---- "remplissage" successif des elements de la grille
c
      Rlim = Rsup_PG(1)
      N_grid = 1
      do lSE=lSE_min,NbSE-NbSE_perte  ! on passe en revue les SE dans l`ordre
         iSE = iListeRay(lSE)         ! de leur distance au centre de l`amas

                                ! En a-t-on termine avec la cellule courante ?
         if (R_SE(iSE).GT.Rlim .OR. lSE.EQ.NbSE-NbSE_perte) then   
            NSE_PG(N_grid) = NSE_PG(N_grid) + n_couche
            M_PG  (N_grid) = M_PG  (N_grid) + m_couche
            T_PG  (N_grid) = T_PG  (N_grid) + t_couche
            Tt_PG (N_grid) = Tt_PG (N_grid) + tt_couche
            MR_PG (N_grid) = MR_PG (N_grid) + mr_couche
            MM_PG (N_grid) = MM_PG (N_grid) + mm_couche
            Npart_PG = Npart_PG + n_couche
            n_couche  = 0
            m_couche  = 0.0d0
            t_couche  = 0.0d0
            tt_couche = 0.0d0
            mr_couche = 0.0d0
            mm_couche = 0.0d0
            
            N_grid = N_grid + 1
            if (N_grid.GT.N_PG) goto 99
            Rlim = Rsup_PG(N_grid)

         end if

         RSE = R_SE(iSE)
         MSE = M_SE(iSE) 
         n_couche  = n_couche  + 1
         m_couche  = m_couche  + MSE
         t_couche  = t_couche  + MSE*T_SE(iSE)
         tt_couche = tt_couche + 0.5d0*MSE*(J_SE(iSE)/R_SE(iSE))**2
         mr_couche = mr_couche + MSE*RSE
         mm_couche = mm_couche + MSE*MSE
      end do

 99   continue

#if (_VERB_PG_>=5)
      write(0,200)
 200  format(' Done !')
      call flush_MD(0)
#endif
#if (_TEST_PG_>0)
      if (.NOT.Tester_PG(iPG)) then
         write (0,*) '!!! PG incorrect :',iPG,' !!!'
         call Meurt('Construire_Grille_PG')
      end if
#endif
c
      end
c
c----------------------------------------------------------------------
      subroutine Normaliser_Grille_PG(Nsnapshots)
c----------------------------------------------------------------------
c     Normalise grid quantities to account for the fact that Nsnapshots
c     Add the current cluster data into the current radial grid 
c     (built with Construire_Grille_PG2). To be used to average profiles
c     over successive snapshots. But no fancy mass/energy correction.
c
      include 'EvolAmas_Grille_inc.f' 
c
c argument :
c ^^^^^^^^^^
      integer Nsnapshots
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer i
      double precision renorm
c
c fonction(s) :
c ^^^^^^^^^^^^^
      logical tester_PG
c......................................................................
c
#if (_VERB_PG_>=5)
      write(0,'("|| Normalise PG grid for ",I3," snapshots ...",$)')
     $     Nsnapshots
      call flush_MD(0)
#endif

c
c---- initialisations
c     
      renorm = 1.0d0/Nsnapshots
c
c---- Renormalisation of grid quantities
c
      do i=1,N_PG
         NSE_PG(i) = NSE_PG(i) ! Should this normalised too??
         M_PG  (i) = renorm*M_PG(i)
         T_PG  (i) = renorm*T_PG(i)
         Tt_PG (i) = renorm*Tt_PG(i)
         MR_PG (i) = renorm*MR_PG(i)
         MM_PG (i) = renorm*MM_PG(i)
      end do

#if (_VERB_PG_>=5)
      write(0,200)
 200  format(' Done !')
      call flush_MD(0)
#endif
#if (_TEST_PG_>0)
      if (.NOT.Tester_PG(iPG)) then
         write (0,*) '!!! PG incorrect :',iPG,' !!!'
         call Meurt('Renormaliser_Grille_PG')
      end if
#endif
c
      end
c
c----------------------------------------------------------------------
      subroutine Construire_Grille (Npc_opt,Npc_min,Npc_max,
     $     iDimGrille,
     $     Rsup_Grille,M_Grille,MR_Grille,T_Grille,Tt_Grille,NSE_Grille,
     $     N_Grille)
c----------------------------------------------------------------------
c     Construction d'une grille des proprietes radiales de l'amas
c     en employant environ Npc_opt SE/cellule (mais au moins Npc_min
c     et au plus Npc_max)
c     LES SE DOIVENT ETRES TRIEES EN R_SE (tableau iListeRay)
c
      include 'EvolAmas_Grille_inc.f' 
c
c arguments :
c ^^^^^^^^^^^
      integer Npc_opt,Npc_min,Npc_max ! Doivent etre des variables car ces valeurs
                                      ! sont eventuellement ajustees par la routine
      integer iDimGrille
      double precision
     $     Rsup_Grille(0:iDimGrille), ! La limite radiale superieure de chaque cellule
     $     M_Grille(iDimGrille),      ! La masse dans chaque cellule
     $     MR_Grille(iDimGrille),     ! La somme des M_i*R_i dans chaque cellule
     $     T_Grille(iDimGrille),      ! L'energie cinetique dans chaque cellule
     $     Tt_Grille(iDimGrille)      ! La composante tg de l'energie cinetique dans chaque cellule
      integer
     $     NSE_Grille(iDimGrille)     ! Le nb de particule dans chaque cellule
      integer N_Grille                ! Le nb d'elements de grille cree
  
c
c constantes :
c ^^^^^^^^^^^^
      double precision factRmax,factRmin
      parameter (factRmax=2.0d0,factRmin=0.5d0)
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer iSE, lSE, n_couche
      double precision
     $     Rmin,Rmax,Rlim,fact_R_opt,
     $     Ri,Ri1,Ri_3,Ri1_3,
     $     m_couche, t_couche,  tt_couche, mr_couche,
     $     MSE,RSE
c......................................................................
c
#if (_VERB_PG_>=5)
      write(0,'("|| Construction grille ...",$)')
      call flush_MD(0)
#endif
c
c---- calcul des parametres
c
      Npc_max = min(Npc_max,NbSE-NbSE_perte)
      Npc_opt = min(Npc_opt,Npc_max)
      Npc_min = min(Npc_min,Npc_opt)
      Rmin = factRmin*R_SE(iListeRay(1))
      Rmax = factRmax*R_SE(iListeRay(NbSE-NbSE_perte))
      fact_R_opt = exp(log(Rmax/Rmin)*Npc_opt/(NbSE-NbSE_perte))
c
c---- initialisations
c     
      n_couche  = 0
      m_couche  = 0.0d0
      t_couche  = 0.0d0
      tt_couche = 0.0d0
      mr_couche = 0.0d0
      Ri1   = 0.0d0
      Ri1_3 = 0.0d0
      RSE = 0.0d0
      Rsup_Grille(0) = Rmin
      Rlim = Rsup_Grille(0)*fact_R_opt
      N_Grille = 1                  ! Numero de la cellule en cours de remplissage
c
c---- "remplissage" successif des elements de la grille
c
      do lSE=1,NbSE-NbSE_perte  ! on passe en revue les SE dans l`ordre
         iSE = iListeRay(lSE)   ! de leur distance au centre de l`amas

                                ! En a-t-on termine avec la cellule courante ?
         if (n_couche.GE.Npc_min)
     $        then              ! Il faut que la cellule contienne au moins Npc_min particules
            Ri = 0.5d0*(RSE+R_SE(iSE)) ! Ri  : limite sup de la couche
                                       ! Ri1 : limite inf de la couche
            if ( (Ri.GE.Rlim) .OR. (n_couche.GE.Npc_max) )
     $           then           ! La cellule est complete si elle a la bonne taille ou 
                                ! si elle contient le nb max de particules
   
               NSE_Grille(N_Grille) = n_couche
               M_Grille  (N_Grille) = m_couche
               T_Grille  (N_Grille) = t_couche
               Tt_Grille (N_Grille) = tt_couche
               MR_Grille (N_Grille) = mr_couche
               Rsup_Grille(N_Grille) = Ri
               Ri_3 = Ri**3
               Ri1 = Ri
               Ri1_3 = Ri_3
               Rlim = Ri*fact_R_opt
               n_couche  = 0
               m_couche  = 0.0d0
               t_couche  = 0.0d0
               tt_couche = 0.0d0
               mr_couche = 0.0d0
   
               N_Grille = N_Grille + 1
               if (N_Grille.GT.idimGrille) then
                  write(0,*) '!!! Nb d''elements max de la grille',
     $                 ' atteint !!!'
                  call Meurt('Construire_Grille')
               end if

            end if
         end if

         RSE = R_SE(iSE)
         MSE = M_SE(iSE) 
         n_couche  = n_couche  + 1
         m_couche  = m_couche  + MSE
         t_couche  = t_couche  + MSE*T_SE(iSE)
         tt_couche = tt_couche + 0.5d0*MSE*(J_SE(iSE)/R_SE(iSE))**2
         mr_couche = mr_couche + MSE*RSE
      end do

      if (n_couche.GT.0) then   ! derniere cellule
         NSE_Grille(N_Grille) = n_couche
         M_Grille  (N_Grille) = m_couche
         T_Grille  (N_Grille) = t_couche
         Tt_Grille (N_Grille) = tt_couche
         MR_Grille (N_Grille) = mr_couche
         Ri = 1.5d0*RSE -
     $        0.5d0*R_SE(iListeRay(NbSE-NbSE_perte-1))
         Rsup_Grille(N_Grille) = Ri
         Ri_3 = Ri**3
      end if

#if (_VERB_PG_>=5)
      write(0,200) N_Grille
 200  format(' Termine ! (',I4,' elements)')
      call flush_MD(0)
#endif
c
      end
c
c----------------------------------------------------------------------
      integer function index_PG (iSE)
c----------------------------------------------------------------------
c     fournit le numero de l'element de la grille PG dans lequel la SE
c     de numero iSE se trouve
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iSE
c
c fonction :
c ^^^^^^^^^^
      integer iCouche_PG
c......................................................................
c
      index_PG = iCouche_PG(R_SE(iSE))
c
      end
c
c----------------------------------------------------------------------
      integer function iCouche_PG (Ray)
c----------------------------------------------------------------------
c     fournit le numero de la couche contenant le rayon Ray
c     renvoie N_PG+1 si Ray>Rsup_PG(N_PG) et 0 si Ray<Rsup_PG(0)
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      double precision Ray
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer jbas,jhaut,jmil
c......................................................................
c
      if (Ray.GE.Revap) then
         iCouche_PG = iPG_evap
         return
      end if
c
c---- recherche par bissection selon "locate" des "Numerical Recipes"
c
      jbas = -1
      jhaut = N_PG+1
      do while (jhaut-jbas.GT.1)
         jmil=(jhaut+jbas)/2
         if (Ray.GT.Rsup_PG(jmil)) then
            jbas = jmil
         else
            jhaut = jmil
         end if
      end do
      iCouche_PG = jhaut
c
      end
c
c----------------------------------------------------------------------
      double precision  function R_PG (icouche)
c----------------------------------------------------------------------
c     fournit le rayon de la base de la couche icouche
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer icouche
c......................................................................
c
      if (icouche.LE.N_PG+1) then
         R_PG = Rsup_PG(icouche-1)
      else
         R_PG = 1.0d20
      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function Rho_PG(iPG,R)
c----------------------------------------------------------------------
c     fournit la densite pour l'element iPG de la grille PG
c     pour un rayon R (la valeur de R est utile pour inter/extra-polation)
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iPG
      double precision R
c
c fonction :
c ^^^^^^^^^^
      integer iCouche_PG
c......................................................................
c
      if (iPG.LT.0) iPG = iCouche_PG(R)
      if (iPG.LE.N_PG) then
         Rho_PG = M_PG(iPG)/V_PG(iPG)
      else
         Rho_PG = M_PG(N_PG)/
     $        (r4Pi3*(R**3-Rsup_PG(N_PG-1)**3))
      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function Rho_PG_lin(iPG,R)
c----------------------------------------------------------------------
c     fournit la densite pour l'element iPG de la grille PG
c     Fait une interpolation lineaire
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iPG
      double precision R
c
c var locales :
c ^^^^^^^^^^^^^
      double precision lR,lR1,lR2,lRho1,lRho2,A,Rcouche
      integer i1,i2
c
c fonction :
c ^^^^^^^^^^
      integer iCouche_PG
c......................................................................
c
      if (iPG.LT.0) iPG = iCouche_PG(R)

      if (iPG.GE.N_PG) then
         i1 = N_PG-1
         i2 = N_PG
      else if (iPG.LE.0) then
         i1 = 1
         i2 = 2
      else
         Rcouche = MR_PG(iPG)/M_PG(iPG)
         if (R.GT.Rcouche) then         
            i2 = min(N_PG,iPG+1)       
            i1 = i2-1
         else
            i1 = max(1,iPG-1)
            i2 = i1+1
         end if
      end if
      lR1 = LOG(MR_PG(i1)/M_PG(i1))
      lR2 = LOG(MR_PG(i2)/M_PG(i2))
      lRho1 = LOG(M_PG(i1)/V_PG(i1))
      lRho2 = LOG(M_PG(i2)/V_PG(i2))
      LR=LOG(R)
      if (i1.EQ.1 .AND. lR.LT.lR1) then
         lR=lR1                 ! ceci cree un plateau en Rho pour R->0 !!!
      end if
      A = (lR2-LR)/(lR2-lR1)
      Rho_PG_lin = EXP(A*lRho1 + (1.0d0-A)*lRho2)
c
      end
c
c----------------------------------------------------------------------
      double precision function Rho_PG_lin2(iPG,R)
c----------------------------------------------------------------------
c     fournit la densite pour l'element iPG de la grille PG
c     Fait une interpolation lineaire
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iPG
      double precision R
c
c constantes :
c ^^^^^^^^^^^^
      integer Nmc
      parameter (Nmc=10)
c
c var locales :
c ^^^^^^^^^^^^^
      double precision lRmc(Nmc),lRhomc(Nmc)
      double precision lR,lR1,lR2,lRho1,lRho2,a,b,Rcouche
      integer i1,i2,imc,imc2
c
c fonction :
c ^^^^^^^^^^
      integer iCouche_PG
c......................................................................
c
      if (iPG.LT.0) iPG = iCouche_PG(R)

      if (iPG.GE.N_PG) then
         i1 = N_PG-1
         i2 = N_PG
      else if (iPG.LE.1) then
         i1 = 1
         i2 = 2
      else
         Rcouche = MR_PG(iPG)/M_PG(iPG)
         if (R.GT.Rcouche) then         
            i2 = min(N_PG,iPG+1)       
            i1 = i2-1
         else
            i1 = max(1,iPG-1)
            i2 = i1+1
         end if
      end if

      lR1 = LOG(MR_PG(i1)/M_PG(i1))
      lR2 = LOG(MR_PG(i2)/M_PG(i2))
      lRho1 = LOG(M_PG(i1)/V_PG(i1))
      lRho2 = LOG(M_PG(i2)/V_PG(i2))
      lR=LOG(R)
      if (i1.EQ.1 .AND. lR.LT.lR1) then          
#if (_EXTRAP_RHO_>0)
                                ! extrapolation for very small R
         do imc=1,Nmc
            lRmc(imc) = LOG(MR_PG(imc)/M_PG(imc))
            lRhomc(imc) = LOG(M_PG(imc)/V_PG(imc))
         end do
         call MoindresCarres(lRmc,lRhomc,Nmc,a,b)
         Rho_PG_lin2 = EXP( lRhomc(1) + a*(LOG(R)-lRmc(1)) )
#else 
                                ! create constant density near the centre
         Rho_PG_lin2 = EXP(lRho1)
         
#endif
      else if (i2.EQ.N_PG .AND. lR.GT.lR2)
     $     then                 ! extrapolation for very large R
         do imc=1,Nmc
            imc2=N_PG-Nmc-1+imc
            lRmc(imc) = LOG(MR_PG(imc2)/M_PG(imc2))
            lRhomc(imc) = LOG(M_PG(imc2)/V_PG(imc2))
         end do
         call MoindresCarres(lRmc,lRhomc,Nmc,a,b)
         Rho_PG_lin2 = EXP( lRhomc(Nmc) + a*(LOG(R)-lRmc(Nmc)) )
      else                      ! linear interpolation in log-log
         a = (lR2-lR)/(lR2-lR1)
         Rho_PG_lin2 = EXP(a*lRho1 + (1.0d0-a)*lRho2)
      end if
c
      end
cc----------------------------------------------------------------------
      double precision function Rho_PG_spl(R)
c----------------------------------------------------------------------
c     fournit la densite pour l'element iPG de la grille PG
c     Interpole par spline cubique
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      double precision R
c
c parametre :
c ^^^^^^^^^^^
      double precision Un6
      parameter (Un6=1.0d0/6.0d0)
c
c var locales :
c ^^^^^^^^^^^^^
      double precision lR(idimPG), lRho(idimPG), Y2(idimPG),
     $     logR, pente
      integer i
c......................................................................
c
      do i=1,N_PG
         lR(i) = LOG(MR_PG(i)/M_PG(i))
         lRho(i) = LOG(M_PG(i)/V_PG(i))
      end do
      call spline(lR,lRho,N_PG,2.0d30,2.0d30,Y2)

      logR=LOG(R)
c
c--- faut-il EXTRApoler ??
c
      if (logR.LT.lR(1)) then
         pente = (lRho(2)-lRho(1))/(lR(2)-lR(1))
     $        -Un6*(lR(2)-lR(1))*Y2(2)
         Rho_PG_spl = lRho(1) + pente*(logR-lR(1))
      else if (logR.GT.lR(N_PG)) then
         pente = (lRho(N_PG)-lRho(N_PG-1))/(lR(N_PG)-lR(N_PG-1))
     $        +Un6*(lR(N_PG)-lR(N_PG-1))*Y2(N_PG-1)
         Rho_PG_spl = lRho(N_PG) + pente*(logR-lR(N_PG))
      else
c
c---- on peut interpoler sans crainte
c
         call splint(lR,lRho,Y2,N_PG,logR,Rho_PG_spl)
      end if
      Rho_PG_spl = EXP(Rho_PG_spl)
c
      end
c
c----------------------------------------------------------------------
      double precision function Rho_PG_spl2(iPG,R,Nspl)
c----------------------------------------------------------------------
c     fournit la densite pour l'element iPG de la grille PG
c     Interpole par spline cubique, en employant 2*Nspl pts pour construire
c     le spline
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iPG,Nspl
      double precision R
c
c parametre :
c ^^^^^^^^^^^
      double precision Un6
      parameter (Un6=1.0d0/6.0d0)
c
c var locales :
c ^^^^^^^^^^^^^
      double precision lR(idimPG), lRho(idimPG), Y2(idimPG),
     $     logR, pente, Rcouche
      integer i,i1,i2,N_
c
c fonction :
c ^^^^^^^^^^
      integer iCouche_PG
c......................................................................
c
      if (iPG.LT.0) iPG = iCouche_PG(R)

      if (iPG.GE.N_PG) then
         i1 = N_PG-Nspl+1
         i2 = N_PG
      else if (iPG.LE.0) then
         i1 = 1
         i2 = Nspl
      else
         Rcouche = MR_PG(iPG)/M_PG(iPG)
         if (R.GT.Rcouche) then         
            i1 = iPG-Nspl+1
            i2 = iPG+Nspl
         else
            i1 = iPG-Nspl
            i2 = iPG+Nspl-1
         end if
      end if
      i1=max(1,i1)
      i2=min(N_PG,i2)
      do i=i1,i2
         lR(i-i1+1) = LOG(MR_PG(i)/M_PG(i))
         lRho(i-i1+1) = LOG(M_PG(i)/V_PG(i))
      end do
      N_ = i2-i1+1
      call spline(lR,lRho,N_,2.0d30,2.0d30,Y2)
      logR=LOG(R)
c
c--- faut-il EXTRApoler ??
c
      if (logR.LT.lR(1)) then
         pente = (lRho(2)-lRho(1))/(lR(2)-lR(1))
     $        -Un6*(lR(2)-lR(1))*Y2(2)
         Rho_PG_spl2 = lRho(1) + pente*(logR-lR(1))
      else if (logR.GT.lR(N_)) then
         pente = (lRho(N_)-lRho(N_-1))/(lR(N_)-lR(N_-1))
     $        +Un6*(lR(N_)-lR(N_-1))*Y2(N_-1)
         Rho_PG_spl2 = lRho(N_) + pente*(logR-lR(N_))
      else
c
c---- on peut interpoler sans crainte
c
         call splint(lR,lRho,Y2,N_,logR,Rho_PG_spl2)
      end if
      Rho_PG_spl2 = EXP(Rho_PG_spl2)
c
      end
c
c----------------------------------------------------------------------
      double precision function DensSE_PG_lin(iPG,R)
c----------------------------------------------------------------------
c     Yields density of super-stars (SE) at grid element iPG
c     Proceeds by klinear interpolation.
c
c     Multiply by Net_par_SE to get stellar number density
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iPG
      double precision R
c
c var locales :
c ^^^^^^^^^^^^^
      double precision lR,lR1,lR2,lDens1,lDens2,A,Rcouche
      integer i1,i2
c
c fonction :
c ^^^^^^^^^^
      integer iCouche_PG
c......................................................................
c
      if (iPG.LT.0) iPG = iCouche_PG(R)

      if (iPG.GE.N_PG) then
         i1 = N_PG-1
         i2 = N_PG
      else if (iPG.LE.0) then
         i1 = 1
         i2 = 2
      else
         Rcouche = MR_PG(iPG)/M_PG(iPG)
         if (R.GT.Rcouche) then         
            i2 = min(N_PG,iPG+1)       
            i1 = i2-1
         else
            i1 = max(1,iPG-1)
            i2 = i1+1
         end if
      end if
      lR1 = LOG(MR_PG(i1)/M_PG(i1))
      lR2 = LOG(MR_PG(i2)/M_PG(i2))
      lDens1 = LOG(NSE_PG(i1)/V_PG(i1))
      lDens2 = LOG(NSE_PG(i2)/V_PG(i2))
      LR=LOG(R)
      if (i1.EQ.1 .AND. lR.LT.lR1)
     $     lR=lR1               ! ceci cree un plateau en Dens pour R->0 !!!
      A = (lR2-LR)/(lR2-lR1)
                                ! / Nb d'etoiles par SE (doit etre constant !!)
      DensSE_PG_lin = EXP(A*lDens1 + (1.0d0-A)*lDens2)
c
      end
c
c----------------------------------------------------------------------
      double precision function Sigma1D_PG(iPG,R)
c----------------------------------------------------------------------
c     fournit la dispertion des vitesses (1-D) pour l'element iPG de la grille PG
c     pour un rayon R (la valeur de R est utile pour inter/extra-polation)
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iPG
      double precision R
c
c var locale :
c ^^^^^^^^^^^^
      double precision arg_sqrt
      integer i
c
c constante(s) :
c ^^^^^^^^^^^^^^
      double precision RacDeuxSur3
      parameter (RacDeuxSur3=0.816497d0)
c
c fonction :
c ^^^^^^^^^^
      integer iCouche_PG
c......................................................................
c
      if (iPG.LT.0) iPG = iCouche_PG(R)
      if (iPG.LE.N_PG) then
         i = iPG
      else
         i = N_PG
      end if 
      arg_sqrt = T_PG(i)/M_PG(i)     
#if (_SURVEILLE_SQRT_>0)
      if (arg_sqrt.LE.0.0d0) then
         write(0,100)iPG,i,T_PG(i),M_PG(i)
 100     format('!!! arg_sqrt <= 0; iPG,i,T_PG(i),M_PG(i) : ',
     $        2(I3,' '),2(E13.4,' '), '!!!')
         call Meurt('Sigma1D_PG')
      end if
#endif
      Sigma1D_PG = RacDeuxSur3*SQRT(arg_sqrt)
c
      end
c
c----------------------------------------------------------------------
      double precision function Sigma1D_PG_lin2(iPG,R)
c----------------------------------------------------------------------
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iPG
      double precision R
c
c constantes :
c ^^^^^^^^^^^^
      integer Nmc
      parameter (Nmc=5)
      double precision lRac2Sur3
      parameter (lRac2Sur3=-0.2027326d0)
c
c var locales :
c ^^^^^^^^^^^^^
      double precision lRmc(Nmc),lSigmc(Nmc)
      double precision lR,lR1,lR2,lSig1,lSig2,a,b,Rcouche
      integer i1,i2,imc,imc2
c
c fonction :
c ^^^^^^^^^^
      integer iCouche_PG
c......................................................................
c
      if (iPG.LT.0) iPG = iCouche_PG(R)

      if (iPG.GE.N_PG) then
         i1 = N_PG-1
         i2 = N_PG
      else if (iPG.LE.0) then
         i1 = 1
         i2 = 2
      else
         Rcouche = MR_PG(iPG)/M_PG(iPG)
         if (R.GT.Rcouche) then         
            i2 = min(N_PG,iPG+1)       
            i1 = i2-1
         else
            i1 = max(1,iPG-1)
            i2 = i1+1
         end if
      end if
      lR1 = LOG(MR_PG(i1)/M_PG(i1))
      lR2 = LOG(MR_PG(i2)/M_PG(i2))
      lSig1 = LOG(T_PG(i1)/M_PG(i1))
      lSig2 = LOG(T_PG(i2)/M_PG(i2))
      lR = LOG(R)
      
      if (i1.EQ.1 .AND. lR.LT.lR1)
     $     then                 ! extrapolation en bas
#if (_EXTRAP_SIG_==0)  
                                ! on extrapole avec une constante
         Sigma1D_PG_lin2 = lSig1
#elif (_EXTRAP_SIG_==1)
                                ! on extrapole avec une droite
         do imc=1,Nmc
            lRmc(imc) = LOG(MR_PG(imc)/M_PG(imc))
            lSigmc(imc) = LOG(T_PG(imc)/M_PG(imc))
         end do
         call MoindresCarres(lRmc,lSigmc,Nmc,a,b)
         Sigma1D_PG_lin2 = lSigmc(1) + a*(LOG(R)-lRmc(1))
#endif
         
      else if (i2.EQ.N_PG .AND. lR.GT.lR2)
     $        then              ! extrapolation en haut         
         do imc=1,Nmc
            imc2=N_PG-Nmc-1+imc
            lRmc(imc) = LOG(MR_PG(imc2)/M_PG(imc2))
            lSigmc(imc) = LOG(T_PG(imc2)/M_PG(imc2))
         end do
         call MoindresCarres(lRmc,lSigmc,Nmc,a,b)
         Sigma1D_PG_lin2 = lSigmc(Nmc) + a*(lR-lRmc(Nmc))
      else                      ! interpolation lineaire
         a = (lR2-lR)/(lR2-lR1)
         Sigma1D_PG_lin2 = a*lSig1 + (1.0d0-a)*lSig2
      end if

      Sigma1D_PG_lin2 = EXP(0.5d0*Sigma1D_PG_lin2 + lRac2Sur3)
c
      end
c
c----------------------------------------------------------------------
      double precision function Rho_PG_SE (iSE)
c----------------------------------------------------------------------
c     fournit la densite pour la SE iSE , densite calculee selon 
c     la grille PG
c
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iSE
c
c variable locale :
c ^^^^^^^^^^^^^^^^^
      integer i_PG
c
c fonction(s) :
c ^^^^^^^^^^^^^
      integer index_PG
c......................................................................
c
      i_PG   = index_PG(iSE)
      if (i_PG.LT.N_PG) then
         Rho_PG_SE = M_PG(i_PG)/V_PG(i_PG)
      else
         Rho_PG_SE = M_PG(N_PG)/
     $        (r4Pi3*(R_SE(iSE)**3-Rsup_PG(N_PG-1)**3))
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Retirer_de_PG(R,T,M,J,iPG)
c----------------------------------------------------------------------
c     Retire la contribution de la SE de rayon=R Ecin/Masse=T Masse=M
c     de la grille PG
c     Fournit le num de la cellule de la grille corresp. (iPG)
c
      include 'EvolAmas_Grille_inc.f'
c
c variables d'appel :
c ^^^^^^^^^^^^^^^^^^^
      double precision R,T,M,J
      integer iPG
c
c bloc common part :            ! Donnees a fournir a Remettre_dans_PG
c ^^^^^^^^^^^^^^^^^^
      integer iPGret            ! numero de la cellule de retrait   
      logical lMin              ! mis a .TRUE. si cette cellule devient sous-occuppe
      common /Retrait_PG/ iPGret,lMin 
c
c var locale :
c ^^^^^^^^^^^^
      integer iPG_tst
      double precision Met
c
c fonction(s) :
c ^^^^^^^^^^^^^
      integer icouche_PG
      logical Tester_PG
c......................................................................
c
      if (iPG.LT.0) iPG = icouche_PG(R)    
      iPGret = iPG
      lMin = (NSE_PG(iPG).EQ.Nmin_PG)
      NSE_PG(iPG) = NSE_PG(iPG) - 1
      M_PG(iPG)   = M_PG(iPG)   - M
      T_PG(iPG)   = T_PG(iPG)   - M*T
      Tt_PG(iPG)  = Tt_PG(iPG)  - 0.5d0*M*(J/R)**2
      MR_PG(iPG)  = MR_PG(iPG)  - M*R
      MM_PG(iPG)  = MM_PG(iPG)  - M*M
c$$$#if (_CHAUFF_BIN_3C_>0)
c$$$      Met = M/(Coef_ConvM*Net_par_SE) ! masse de chaque etoile en masses solaires
c$$$      Y3c_PG(iPG) = Y3c_PG(iPG) - Met**2/sqrt(2.0d0*T)**3
c$$$#endif
#if (_TEST_PG_>=1)
      if (NSE_PG(iPG).LT.0) then
         write(0,90)iPG,NSE_PG(iPG)
 90      format('!!! NSE_PG(',I5,') = 'I5)
         call Meurt('Retirer_de_PG')
      end if
      if (T_PG(iPG).LT.0.0d0) then
         write(0,110)iPG,T_PG(iPG)
 110     format('!!! T_PG(',I5,') = ',E12.4)
         call Meurt('Retirer_de_PG')
      end if
      if (M_PG(iPG).LT.0.0d0) then
         write(0,120)iPG,M_PG(iPG)
 120     format('!!! M_PG(',I5,') = ',E12.4)
         call Meurt('Retirer_de_PG')
      end if
      if (MR_PG(iPG).LT.0.0d0) then
         write(0,130)iPG,MR_PG(iPG)
 130     format('!!! MR_PG(',I5,') = ',E12.4)
         call Meurt('Retirer_de_PG')
      end if
      if (MM_PG(iPG).LT.0.0d0) then
         write(0,140)iPG,MM_PG(iPG)
 140     format('!!! MM_PG(',I5,') = ',E12.4)
         call Meurt('Retirer_de_PG')
      end if
#endif

#if (_TEST_PG_>=10)
      if (.NOT.Tester_PG(iPG_tst)) then
         write(0,100)
     $        iPG_tst,iPG,R,T,M
 100     format(
     $        '!!! Grille PG erronnee en iPG = ',I3,' !!!',/,
     $        '    iPG,R,T,M : ',I3,' ',3(E12.4," ") )
         call Meurt('Retirer_de_PG')
      end if
#endif
c
      end
c
c----------------------------------------------------------------------
      subroutine Mettre_dans_PG(R,T,M,J,iPG)
c----------------------------------------------------------------------
c     Ajoute la contribution de la SE de rayon=R Ecin/Masse=T Masse=M
c     a la grille PG
c     Fournit le num de la cellule de la grille corresp. (iPG)
c
      include 'EvolAmas_Grille_inc.f'
c
c variables d'appel :
c ^^^^^^^^^^^^^^^^^^^
      double precision R,T,M,J
      integer iPG
c
c var locale :
c ^^^^^^^^^^^^
      integer iPG_tst,i
      double precision Ttot,Ttot2,Met
c
c fonction(s) :
c ^^^^^^^^^^^^^
      integer icouche_PG
c......................................................................
c
      if (R.LT.Revap) then ! on n'integre pas les SE evaporees a la grille

         if (iPG.LT.0) iPG = icouche_PG(R)

         if (iPG.GT.N_PG) then  ! debordement par le haut
#if (_VERB_PG_>=15)
            write(0,*)'|| Debordement PG par le haut'
#endif
            iPG = N_PG
            Rsup_PG(N_PG) = 1.5d0*R-0.5d0*Rsup_PG(N_PG)
            V_PG   (N_PG) = r4Pi3*(Rsup_PG(N_PG)**3-Rsup_PG(N_PG-1)**3)
         else if (iPG.LT.1) then ! debordement par le bas
#if (_VERB_PG_>=15)
            write(0,*)'|| Debordement PG par le bas'
#endif
            Rsup_PG(0) = 0.5d0*R
            iPG = 1
            V_PG(1) = r4Pi3*(Rsup_PG(1)**3-Rsup_PG(0)**3)
         end if
#if (_TEST_PG_>0) 
         if (Rsup_PG(iPG).LT.R .OR. Rsup_PG(iPG-1).GT.R) then
            write(0,*) '!!! R,Rsup_PG(iPG-1),Rsup_PG(iPG) = ',
     $           R,Rsup_PG(iPG-1),Rsup_PG(iPG),' !!!'
            call Meurt('Mettre_dans_PG')
         end if
#endif   
         NSE_PG(iPG) = NSE_PG(iPG) + 1
         M_PG(iPG)   = M_PG(iPG)   + M
         T_PG(iPG)   = T_PG(iPG)   + M*T
         Tt_PG(iPG)  = Tt_PG(iPG)  + 0.5d0*M*(J/R)**2
         MR_PG(iPG)  = MR_PG(iPG)  + M*R
         MM_PG(iPG)  = MM_PG(iPG)  + M*M
c$$$#if (_CHAUFF_BIN_3C_>0)
c$$$         Met = M/(Coef_ConvM*Net_par_SE) ! masse de chaque etoile en masses solaires
c$$$         Y3c_PG(iPG) = Y3c_PG(iPG) + Met**2/sqrt(2.0d0*T)**3
c$$$#endif
      else
         iPG=iPG_evap
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Remettre_dans_PG(R,T,M,J,iPG)
c----------------------------------------------------------------------
c     Ajoute la contribution de la SE de rayon=R Ecin/Masse=T Masse=M
c     a la grille PG
c     iPGret est le numero de la cellule d'ou la SE a ete precedemment
c     retiree.
c     Fournit le num de la cellule de la grille corresp. (iPG)
c     Reconstruit la grille si necessaire, au quel cas iRecons=1 signale
c     un sur-peuplement de la nouvelle cellule et iRecons=-1 un sous-
c     peuplement de l'ancienne cellule. iRecons=0 signifie qu'une
c     reconstruction n'a pas ete necessaire
c
      include 'EvolAmas_Grille_inc.f'
c
c variables d'appel :
c ^^^^^^^^^^^^^^^^^^^
      double precision R,T,M,J
      integer iPG
c
c bloc common part :
c ^^^^^^^^^^^^^^^^^^
      integer iPGret            ! Donnees provenant de Retirer_de_PG
      logical lMin
      common /Retrait_PG/ iPGret,lMin 
      integer iRecons
      common /Recons_PG/ iRecons
c
c var locale :
c ^^^^^^^^^^^^
      integer iPG_tst,i
      double precision Ttot,Met
c
c fonction(s) :
c ^^^^^^^^^^^^^
      integer icouche_PG
      logical Tester_PG
c......................................................................
c
      iRecons = 0

      if (iPG.LT.0) iPG = icouche_PG(R)

      if (iPG.GT.N_PG) then     ! etendre la cellule superieure de la grille
#if (_VERB_PG_>=15)
         write(0,*)'|| Debordement PG par le haut'
#endif
         iPG=N_PG
         Rsup_PG(N_PG) = 1.5d0*R-0.5d0*Rsup_PG(N_PG)
         V_PG   (N_PG) = r4Pi3*(Rsup_PG(N_PG)**3-Rsup_PG(N_PG-1)**3)
      else if(iPG.LT.1) then    ! etendre la cellule inferieure de la grille
         iPG=1
#if (_VERB_PG_>=15)
         write(0,*)'|| Debordement PG par le bas'
#endif
         Rsup_PG(0) = 0.5d0*R
         V_PG(1) = r4Pi3*(Rsup_PG(1)**3-Rsup_PG(0)**3)
      end if

      NSE_PG(iPG) = NSE_PG(iPG) + 1
      M_PG(iPG)   = M_PG(iPG)   + M
      T_PG(iPG)   = T_PG(iPG)   + M*T
      Tt_PG(iPG)  = Tt_PG(iPG)  + 0.5d0*M*(J/R)**2
      MR_PG(iPG)  = MR_PG(iPG)  + M*R
      MM_PG(iPG)  = MM_PG(iPG)  + M*M
c$$$#if (_CHAUFF_BIN_3C_>0)
c$$$      Met = M/(Coef_ConvM*Net_par_SE) ! masse de chaque etoile en masses solaires
c$$$      Y3c_PG(iPG) = Y3c_PG(iPG) + Met**2/sqrt(2.0d0*T)**3
c$$$#endif
      
      lMin = (lMin .AND. iPG.NE.iPGret)
      if (lMin) then            ! sous-peuplement d'un element de grille
#if (_VERB_PG_>0)
         write(0,105)iPGret
 105     format(
     $        '|| Nb min de SE par cellule atteint en iPG = ',I3,/,
     $        '   reconstruction PG' )
#endif
         iRecons = -1
         call TriAmas
         call Construire_Grille_PG(NSEc_PG)
         iPG = iCouche_PG(R)
      else if (NSE_PG(iPG).GT.Nmax_PG)
     $        then              ! sur-peuplement d'un element de grille
#if (_VERB_PG_>0)
         write(0,115)iPG
 115     format(
     $        '|| Nb max de SE par cellule atteint en iPG = ',I3,/,
     $        '   reconstruction PG' )
#endif
         iRecons =  1
         call TriAmas
         call Construire_Grille_PG(NSEc_PG)
         iPG = iCouche_PG(R)
      end if 

#if (_TEST_PG_>=10)
      if (.NOT.Tester_PG(iPG_tst)) then
         write(0,100)
     $        iPG_tst,iPG,R,T,M
 100     format(
     $        '!!! Grille PG erronnee en iPG = ',I3,' !!!',/,
     $        '    iPG,R,T,M : ',I3,' ',3(E12.4," ") )
         call Meurt('Remettre_dans_PG')
      end if
#endif
c
      end
c
c----------------------------------------------------------------------
      logical function Tester_PG(iPG)
c----------------------------------------------------------------------
c
      include 'EvolAmas_Grille_inc.f'
c
c constante :
c ^^^^^^^^^^^
      double precision Tol_min,Tol_max
      parameter (Tol_min=1.0d-10, Tol_max=1.0d-8)
c
c argument :
c ^^^^^^^^^^
      integer iPG
c......................................................................
c
      Tester_PG = .FALSE.

      do iPG=1,N_PG
         if (NSE_PG(iPG).LT.Nmin_PG) then
            write(0,*) '1'
            return
         else if (NSE_PG(iPG).GT.Nmax_PG) then
            write(0,*) '2'
            return
         else if (T_PG(iPG).LT.0.0d0 .OR.
     $           M_PG(iPG).LT.Tol_min) then
            write(0,*) '3'
            return
         else if (Rsup_PG(iPG).LT.Rsup_PG(iPG-1)) then
            write(0,*) '4'
            return
         end if
      end do

      Tester_PG = .TRUE.
      iPG = 0      
c
      end
c       
c----------------------------------------------------------------------
      subroutine Ecrire_Grille_PG (iOut)
c----------------------------------------------------------------------
      include 'EvolAmas_Grille_inc.f'
c
c variable d'appel :
c ^^^^^^^^^^^^^^^^^^
      integer iOut              ! num d`unite
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer i
c......................................................................
c
      write(iOut,100)N_PG,rNSEc_PG,Nmin_PG,Nmax_PG
 100  format(
     $     '# Profil de densitee grossier pour "EvolAmasXXX"',/,
     $     '# Nb de cellules : ',I5,/,
     $     '# Valeur moy,min,max de SE par cellule : ',
     $     F8.2,' ',I4,' ',I4,/,
     $     '# 1:Rsup_couche 2:M_couche 3:V_couche 4:T_couche',
     $     ' 5:Ttg_couche 6:MR_couche 7:MM_couche 8:NSE_couche')
      !!! write(iOut,110) 0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0
      do i=1,N_PG
         write(iOut,110) Rsup_PG(i),M_PG(i),V_PG(i),
     $        T_PG(i),Tt_PG(i),MR_PG(i),MM_PG(i),NSE_PG(i)
      end do
 110  format(7(E12.5,' '),I5)
c
      end
