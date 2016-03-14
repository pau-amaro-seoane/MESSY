c======================================================================
c     procedures d'initialisation pour EvolAmas.F
c======================================================================

#define _DIST_T_INI_ 4

c----------------------------------------------------------------------
      subroutine Init_Grille
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c......................................................................
c
      if (flag_Tri.NE.1) then
         call TriAmas
         flag_Tri = 1
      end if
c$$$#if (_CHAUFF_BIN_3C_>0)
c$$$      if (flag_Etoiles.LT.1) then
c$$$         write(0,*) '!!! Il faut initialiser les grandeurs stellaires',
c$$$     $        ' avant la grille !!!'
c$$$         call Meurt('Init_Grille')
c$$$      end if
c$$$#endif
      call Construire_Grille_PG(NbSECouchePot)
      flag_Grille = 1
c
      end
c
#if (_CHAUFF_BIN_3C_>10)
c----------------------------------------------------------------------
      subroutine Init_ChauffBin3c
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c......................................................................
c
      if (flag_Tri.NE.1) then
         call TriAmas
         flag_Tri = 1
      end if
      if (flag_Etoiles.LT.1) then
         write(0,*) '!!! Il faut initialiser les grandeurs stellaires',
     $        ' avant de calcluer les param du chauffage par',
     $        ' binaires !!!'
         call Meurt('Init_ChauffBin3c')
      end if
      call Determiner_Param_CB3c
c
      end
c
#endif
c----------------------------------------------------------------------
      subroutine Init_Etoiles
c----------------------------------------------------------------------
#if (_LIMITE_EXT_==0)
#   define _M_INI_ Mamas_ini
#   define _EST_INI_ Estell_grav_ini
#   define _N_INI_ rNbEtoiles_ini
#else
#   define _M_INI_ (Mamas_ini + Mext_LimExt)
#   define _EST_INI_ (Estell_grav_ini + Egrav_ext_LimExt)
#   define _N_INI_ (rNbEtoiles_ini/(1.0d0-fNbEt_ext_LimExt))
#endif
      include 'EvolAmas_inc.f'
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer i, Ntot
      double precision Mtot,Coef
c......................................................................
c         
      Coef_ConvR = Rsol_en_pc/TailleAmas_en_pc
      Net_par_SE = rNbEtoiles/(NbSE-NbSE_Perte)
      Mtot=0.0d0
c
c---- If the simulation is started from zero, we assume that the
c     the average stellar mass is computed over all particles, even in the unlikely 
c     case that some are already considered lost. (modified 29/2/08)
c
      if (i_ini.EQ.0) then
         do i=1,NbSE
            Mtot = Mtot+M_SE(i)
         end do
         Ntot = NbSE
      else
c
c---- If we continue an interupted simulation, we assume
c     that the average stellar mass is computed over the particles
c     that are still considered part of the cluster. (modified 29/2/08)
c
         do i=1,NbSE
            if (R_SE(i).LT.R_perte) then
               Mtot = Mtot+M_SE(i)
            end if
         end do
         Ntot = NbSE-NbSE_Perte
      end if
         
      Coef_ConvM = Mtot/(Ntot*Net_par_SE*MasseEtoileDef)
      Coef = Ntot*MasseEtoileDef/Mtot
                                ! Une annee en unites du code (pour evolution stellaire)
      UneAnnee = 6.707d-8/sqrt(TailleAmas_en_pc**3*Coef_ConvM/
     $     _M_INI_)
      do i=1,NbSE
         if (R_SE(i).LT.R_perte) call
     $        stellar_SetByAge_ms(i,Coef*M_SE(i),StellMet_def,0.0d0)
      end do
         
      flag_Etoiles = 1
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_Arbre
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c......................................................................
c
      if (flag_Tri.NE.1) then
         call TriAmas
         flag_Tri = 1
      end if
      Nb_Parcours_Arbre = 0
      call Construire_Arbre_Amas
      flag_Arbre = 1
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_Temps
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_Relax_common.f'
c
c constantes :
c ^^^^^^^^^^^^
      double precision cT_ini
      parameter (cT_ini=2.0d0) ! modification 5.10.00 : 2 au lieu de 10 
c
c var locales :
c ^^^^^^^^^^^^^
      integer iRang,iSE
      double precision Tps_Moy,Tps_med,T_ini
c
c fonctions :
c ^^^^^^^^^^^
      double precision random,selip,TrelLoc,RandomGauss,P_TH,PasTps_TH
c......................................................................
c
      if (flag_Etoiles.LE.0) then
         write(0,*) '!!! impossible d''executer Init_Temps avant ',
     $        'd''avoir defini le donnees stellaires !!!'
         call Meurt('Init_Temps')
      end if
      if (flag_Grille.LE.0) then
         write(0,*) '!!! impossible d''executer Init_Temps avant ',
     $        'd''avoir construit la grille !!!'
         call Meurt('Init_Temps')
      end if
      if (flag_Tri.NE.1) then
         call TriAmas
         flag_Tri = 1
      end if
#if (_DIST_T_INI_==4)
      T_ini = 0.0d0
      do iRang=1,NbSE-NbSE_perte
         T_ini = T_ini+P_TH(iRang)*PasTps_TH(iRang) ! calcul du pas de tps moyen
      end do

      T_ini=cT_ini*T_ini/Frac_Trelax
c
c     T_ini est le temps moyen de l'amas lorsque le nb de tirages
c     moyen de chaque SE est cT_ini (NON: modification 5.10.00 : ajoute /Frac_Trelax)
c
#endif
c
c---- "repartition" des temps initiaux
c
      do iRang=1,NbSE

         iSE=iListeRay(iRang)
         
         if (R_SE(iSE).GE.Revap) then ! SE evaporee
            Tps_SE(iSE) = Tps_evap
         else
#if (_DIST_T_INI_==1)
            Tps_SE(iSE) = 1.0d-10*random()
#elif (_DIST_T_INI_==2)
            Tps_SE(iSE) = Frac_Trelax*TrelLoc(R_SE(iSE))*random()
#elif (_DIST_T_INI_==3)
            Tps_SE(iSE) = PasTps_TH(iRang)*RandomGauss()
#elif (_DIST_T_INI_==4)
c
c     Imaginons que l'amas a deja evolue jusqu'au temps T_ini, alors le nb de tirages
c     subit par une particule de pas de temps Dt serait n=T_ini/Dt +/- sqrt(n).
c     Donc son temps propre sera T_ini +/- Dt*sqrt(n) = T_ini +/- sqrt(Dt*T_ini)
c     Doonc il semble realiste d'attribuer aux particules un dispertion en temps prop.
c     a la racine de leur pas de tps
c     
            Tps_SE(iSE) = sqrt(T_ini*PasTps_TH(iRang))*RandomGauss()
#else
            Tps_SE(iSE) = 0.0d0
#endif
         end if
      end do
      Tps_med = selip((NbSE-NbSE_perte)/2,NbSE,Tps_SE) 
                                ! ne marche plus vraiment avec des SE de masses differentes
      Tps_moy = 0.0d0
      
      do iSE=1,NbSE
         if (R_SE(iSE).LT.Revap) then
            Tps_moy = Tps_moy+M_SE(iSE)*Tps_SE(iSE)
            Tps_SE(iSE) = Tps_SE(iSE) - Tps_med
         end if
      end do
      Tps_moy = Tps_moy/M_amas
      Tps_Amas = 0.0d0
      
      flag_Temps = 1
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_TN
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_Etoiles_param.f'
c......................................................................
c
#if (_PRESENCE_TN_>0)
      M_ctr = M_TN_ini_def
      if (Type_CtrObj_def.EQ.'TN') then
         iType_CtrObj = TypeEtoile_BH
      else if (Type_CtrObj_def.EQ.'MS') then
         iType_CtrObj = TypeEtoile_MS
      else
         write(0,*) '!!! Type_CtrObj_def has wrong value ("',
     $        Type_CtrObj_def,'") !!!'
         call Meurt('Init_TN')
      end if
#else
      M_ctr = 0.0d0
      iType_CtrObj = 0
#endif
      M_TN = M_ctr
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_Maree
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_Maree_common.f'
c......................................................................
c
#if (_TRONC_MAREE_>0)
      R_Mar = R_Mar_ini_def
      C_Mar = 0.0d0
#else
      R_Mar = 1.0d30
      C_Mar = 0.0d0
#endif
      call Regulariser_CR_Maree
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_Conserv
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      integer i
c......................................................................
c
      MSEtot_ini = 0.0d0
      do i=1,NbSE               ! total mass contained in _all_ particles
                                ! does not change, even at cloning (is probably 1.0d0)
         MSEtot_ini = MSEtot_ini + M_SE(i)
      end do
      Estell_grav_ini=Estell_grav_amas ! part purement stellaire, 
                                       ! contribution du potentiel du TN negligee
      Etot_ini=Etot_amas
      Ttot_ini=Ecin_amas
      Mamas_ini=M_amas
      Mtn_ini=M_ctr
      NbSE_ini=NbSE
      rNbEtoiles_ini=rNbEtoiles
      dEtot_coll=0.0d0
      dEtot_evst=0.0d0
      dEtot_dechir=0.0d0
      dEtot_disp=0.0d0
      dEtot_captGW=0.0d0
      dEtot_evap=0.0d0
      dEtot_bin3c=0.0d0
      dEtot_nk=0.0d0
      dMaccr_coll=0.0d0
      dMaccr_evst=0.0d0
      dMaccr_dechir=0.0d0
      dMaccr_disp=0.0d0
      dMaccr_captGW=0.0d0
      dMejec_coll=0.0d0
      dMejec_evst=0.0d0
      dMejec_dechir=0.0d0
      dMejec_evap=0.0d0
      dTtot_coll=0.0d0
      dTtot_evst=0.0d0
      dTtot_dechir=0.0d0
      dTtot_disp=0.0d0
      dTtot_captGW=0.0d0
      dTtot_evap=0.0d0
      dTtot_nk=0.0d0
c
      end
c


