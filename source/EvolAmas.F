#define _DIFFERENT_M_FOR_GW_ 0 /* !!!!!!!!! DANGER  !!!!!!!!!!! */ 


#define _CONTROLE_CENTRE_ 0
#define _SURV_PCOLL_ 0
#define _SURV_PKICK_ 0
#define _SURV_ERR_E_ 1
#define _TEST_ZERO_MASS_ 1

#define _LIMITER_THETA_ 1
#define _THETA_GIERSZ_ 0
#define _INTERDIT_RENC_DIST_ 0

#define _TYPE_DENS_RENC_ 3

#define _ECRIT_STAT_RELAX_ 0
#define _ECRIT_RENC_RELAX_ 0

#define _VERIF_PERI_APO_ 0
#define _TEST_PERI_APO_ 0
#define _TEST_RENC_QUANTITIES_ 0

#define _PN_FLC_EXACT_ 1
#define _CALC_THETA_LC_EXACT_ 1
#define _CI_RW_AVANT_RENC_ 0
#define _ANNONCE_SE_DANS_LC_ 0
#define _ORBAV_TREL_ 0


#if (_LIMITE_EXT_==0)
#   define _M_INI_ Mamas_ini
#   define _EST_INI_ Estell_grav_ini
#   define _N_INI_ rNbEtoiles_ini
#else
#   define _M_INI_ (Mamas_ini + Mext_LimExt)
#   define _EST_INI_ (Estell_grav_ini + Egrav_ext_LimExt -Mamas_ini*B_LimExt)
#   define _N_INI_ (rNbEtoiles_ini/(1.0d0-fNbEt_ext_LimExt))
#endif


#define _PASSIF_ 0
#if (_PASSIF_==0)
c======================================================================
      program EvolAmas
c======================================================================
c     PROGRAMME D'EVOLUTION D'AMAS AVEC RELAXATION A DEUX CORPS
c     potentiel calcule par arbre binaire
c     VERSION REPRODUISANT LES PROCEDURES DE TIRAGES DE HENON
c======================================================================
c
c     Ref. :   Henon (1971), Astrophysics and space science, 14, 
c              pp 151-167
c
c              Henon (1973) in Dynamical Structure and Evolution
c              of stellar systems, 3rd advanced course of the SSAA,
c              pp 224-258
c======================================================================
      implicit none
c......................................................................
c
c---- Indiquer les specificites du code
c
      write(0,*) '>> Ecrire_CaracCode'
      call Ecrire_CaracCode
c
c---- Lire fichier de Parametres
c
      write(0,*) '>> LireParam_EvolAmas'
      call LireParam_EvolAmas
c
c---- Intialiser
c
      write(0,*) '>> Init_EvolAmas'
      call Init_EvolAmas
c
c---- Lire les fichiers d'entree et initialiser les structures de donnees
c     necessaires
c
      write(0,*) '>> Init_Amas'
      call Init_Amas
c
c---- Check "integrity" of cluster data
c
      write(0,*) '>> Check_Amas'
      call Check_Amas
c
c---- Faire la sauvegarde initiale de toute l'info utile 
c
      write(0,*) '>> Sauvegarde_Init'
      call Sauvegarde_Init
c
c---- Faire Evoluer l'amas
c
      write(0,*) '>> Evol_Amas'
      call Evol_Amas
c
c---- Terminer
c
      call quitter
c
      end

#endif /* ifndef _PASSIF_ */
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine Ecrire_CaracCode
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c......................................................................
c
      call WritePreprocDefines(iFich_Prov,'CaracCode.asc')
c
      end
c      
c----------------------------------------------------------------------
      subroutine LireParam_EvolAmas
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      integer iarg, narg, iargc
      character*256 Chaine
      logical l_stdin
c......................................................................
c
      l_stdin=.TRUE.
c
c---- lecture des parametres sur la ligne de commande
c
      narg=iargc()
      iarg=1
      if (narg.GT.0) then
         call getarg(iarg,Chaine)

         do while (iarg.LE.narg .AND. Chaine(1:1).EQ.'-')
            if (Chaine(2:index(Chaine,' ')-1).EQ.'f') then
               l_stdin=.FALSE.
               iarg=iarg+1
               call getarg(iarg,Chaine)
               open(iFich_Prov,file=Chaine,status='OLD')
            else
               write(0,*) '!!! Option ',Chaine(1:index(Chaine,' ')-1),
     $              ' non valide !!!'
               call Meurt('LireParam_EvolAmas')
            end if
            iarg=iarg+1
            if (iarg.LT.narg) call getarg(iarg,Chaine)
         end do
      end if
c
c---- Lecture du fichier de parametres
c
      if (l_stdin) then
         call Lire_Fich_Param_EvolAmas(5)
      else
         call Lire_Fich_Param_EvolAmas(iFich_Prov)
         close(iFich_Prov)
      end if
c
c---- Afficher les parametres lus
c
      open(iFich_Prov,file='_PARAMS_.asc')
      call Ecrire_Fich_Param_EvolAmas(iFich_Prov)
      close(iFich_Prov)

      if (iVerbeux.GT.0) then 
         call Ecrire_Fich_Param_EvolAmas(0)
      end if

      end
c
c----------------------------------------------------------------------
      subroutine Init_EvolAmas
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      integer ipid, istat
      character*80 hostname,date,hour
      character car
      logical lOK
c
c fonction :
c ^^^^^^^^^^
      integer LongChaine, getpid_MD
c......................................................................
c
      ipid = getpid_MD()
      call Hostname_MD(hostname)
      call QuelleDate_str(date)
      call QuelleHeure_str(hour)
c
c---- on met le pid du processus dans un fichier
c
      open(iFich_Prov,file=NomFichPID)
      write(iFich_Prov,*) ipid
      close(iFich_Prov)
c
c---- on met differentes infos dans un fichier
c
      call open_append(iFich_Prov,NomFichINFO,lOK)
      if (.NOT.lOK) then
         write(0,*) '!!! Could not open file ',NomFichINFO,' !!!'
         call Meurt('Init_EvolAmas')
      end if
      write(iFich_Prov,100) date(1:LongChaine(date)),
     $     hour(1:LongChaine(hour)), hostname(1:LongChaine(hostname)),
     $     ipid

 100  format(/,
     $     'date : ',A,/,
     $     'time : ',A,/,
     $     'host : ',A,/,
     $     'pid  : ',I8)

      call flush_MD(iFich_Prov)

      iPas_Evol = i_ini

      iSauv        = 0
      iDuplic      = 0
      iNbRecGrille = 0
      iNbRecArbre  = 0
      i_NbTirVN_PosOrb = 0
      iNb_RejetPaire = 0
      iNb_AnnulRelaxEvap = 0
      iNb_AnnulRelaxCC = 0
      iNb_CalcPot = 0
      iNb_SuperRenc = 0

      iNb_Relax  = 0
      iNb_Coll   = 0
      iNb_kick   = 0
      iNb_Dechir = 0
      iNb_Disp   = 0
      iNb_tests_LC = 0

      Dev_SuperRenc_moy = 0.0d0
      Dev2_SuperRenc_moy = 0.0d0
      iNb_SuperRenc_Pi16 = 0
      iNb_SuperRenc_Pi8  = 0
      iNb_SuperRenc_Pi4  = 0
      iNb_SuperRenc_Pi2  = 0

      iNb_Coll_P_025 = 0
      iNb_Coll_P_05  = 0
      iNb_Coll_P_1   = 0
      iNb_Coll_P_2   = 0
      P_coll_moy = 0.0d0
      P_coll_max = 0.0d0


      iNb_Kick_P_025 = 0
      iNb_Kick_P_05  = 0
      iNb_Kick_P_1   = 0
      iNb_Kick_P_2   = 0
      P_Kick_moy = 0.0d0
      P_Kick_max = 0.0d0

      dR_rel_Renc_moy = 0.0d0
      dR_rel_Renc_max = 0.0d0
      iNb_Renc_dR_001 = 0
      iNb_Renc_dR_01  = 0
      iNb_Renc_dR_05  = 0
      iNb_Renc_dR_1   = 0

      flag_Arbre   = 0
      flag_Etoiles = 0
      flag_Grille  = 0
      flag_Temps   = 0
      flag_Tri     = 0

      iexit_stat = 0
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_Amas
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_Relax_common.f'
      include 'EvolAmas_CtrlCentre_inc.f'
      include 'EvolAmas_Etoiles_param.f'
      include 'fxdr.inc'
c
c var locale :
c ^^^^^^^^^^^^
      integer i, iFich_XDR, ixdr_err, iOK
      logical lExiste
      double precision M_,Mr_
      integer*1 iTr_
c
c fonction(s) : 
c ^^^^^^^^^^^^^
      logical ReOuvrirFichBIN, ReOuvrirFichASC, ReOuvrirLectureFichXDR
      integer time_MD
#if (_EVOL_STELL_==1)
      double precision RelMR, T_MS
#endif
c......................................................................
c
c---- Lire le fichier d'amas (indispensable)
c
      if (ReOuvrirLectureFichXDR(iFich_XDR,TypeAMAS)) then
         call LitAmas_XDR(iFich_XDR)
         ixdr_err = ixdrclose(iFich_XDR)
      else
         write(0,*) '!!! impossible de continuer sans fichier',
     $        ' d''amas !!!'
         call Meurt('Init_Amas')
      end if
c
c---- Mass--Radius relation for MS stars
c
      inquire(file=NomFichMRrelMS,exist=lExiste)
      if (lExiste) then
         write(0,*) '> Read file for MS Mass--Radius relation'
         open(iFich_Prov,file=NomFichMRrelMS,status='old')
         call Read_RelMR(iFich_Prov)
         close(iFich_Prov)
      end if
c
c---- Lire le fichier de donnees stellaires (facultatif)
c
      if (ReOuvrirLectureFichXDR(iFich_XDR,TypeETOILES)) then
         call LitEtoiles_XDR(iFich_XDR) ! %SeBa% change to routine that "create" SeBa stars (from information in file)
         ixdr_err = ixdrclose(iFich_XDR)
         flag_Etoiles = 2
      else
         call Init_Etoiles ! %SeBa% call SeBa initialisation to "create" SeBa stars assuming all stars on ZAMS
      end if

#if (_EVOL_STELL_==10) /* SeBa stellar evolution */
/* This is just a temporary "solution" to test SeBa integration */
      do i=1,NbSE
         ages(i)=0.0d0
         metals(i)=0.04d0
      end do
      t_=TpsUtil_MD()
      write(6,*) '>> Creating SeBa stars...'
      call create_stars_from_zams(NbSE, Met_SE, ages, metals, iflag)
      t_=TpsUtil_MD()-t_
      write(6,*) '>> SeBa stars created; time needed: ',t_
      do i=1,NbSE
         r=stellar_radius(i)
         m=stellar_mass(i)
      end do
      t_=TpsUtil_MD()-t_
      write(6,*) '>> got masses and radii; time needed: ',t_
      do i=1,NbSE
         call evolve_star(i,1000.0d0)
      end do
      t_=TpsUtil_MD()-t_
      write(6,*) '>> stars evolved; time needed: ',t_
      do i=1,NbSE
         call evolve_star(i,10.0d7)
      end do
      t_=TpsUtil_MD()-t_
      write(6,*) '>> stars evolved; time needed: ',t_
      do i=1,NbSE
         call evolve_star(i,10.0001d7)
      end do
      t_=TpsUtil_MD()-t_
      write(6,*) '>> stars evolved; time needed: ',t_
      stop
#endif
c
c---- Lire le fichier de trou noir
c
#if (_PRESENCE_TN_>0)
      if (ReOuvrirLectureFichXDR(iFich_XDR,TypeTN)) then
         call LitTN_XDR(iFich_XDR)
         ixdr_err = ixdrclose(iFich_XDR)
      else
         call init_TN
      end if
#endif
c
c---- Lire le fichier de donnees du chp de maree
c
#if (_TRONC_MAREE_>0)
      if (ReOuvrirLectureFichXDR(iFich_XDR,TypeMAREE)) then
         call LitMaree_XDR(iFich_XDR)
         ixdr_err = ixdrclose(iFich_XDR)
      else
         call Init_Maree   
      end if
#endif
c
c---- Lire le fichier de limite exterieure
c
#if (_LIMITE_EXT_>0)
      if (ReOuvrirLectureFichXDR(iFich_XDR,TypeLimExt)) then
         call LitLimExt_XDR(iFich_XDR)
         ixdr_err = ixdrclose(iFich_XDR)
      end if
#endif
c
c---- calculs des energies totales
c
      call Calc_EnerGlob
c
c---- Lire le fichier de suivi de la conservation Masse/Energie (facultatif)
c
      if (ReOuvrirLectureFichXDR(iFich_XDR,TypeCONS)) then
         call LitCons_XDR(iFich_XDR)
         ixdr_err = ixdrclose(iFich_XDR)
      else
         call Init_Conserv
      end if 
c
c---- calcul des facteurs de conversions dimensionnelles
c
      call Calc_ConvUnites
c
c---- calcul des parametres lies a la relaxation
c
      call CalcParam_Relax
      call Calc_RelaxGlob
c
c---- Lire le fichier de la grille (facultatif)
c
      if (ReOuvrirFichBin(iFich_Prov,TypePG)) then
         call LireGrilleBIN(iFich_Prov)
         close(iFich_Prov)
         flag_Grille = 2
      else
         call Init_Grille
      end if
c
c---- Lire le fichier de l'arbre binaire (facultatif)
c
      if (ReOuvrirFichBin(iFich_Prov,TypeARBRE)) then
         call LireArbreBIN(iFich_Prov)
         close(iFich_Prov)
         flag_Arbre = 2
      else
         call Init_Arbre
      end if
c
c---- Lire le fichier du generateur de nb aleatoires (facultatif)
c
      if (ReOuvrirFichAsc(iFich_Prov,TypeRAND)) then
         call read_random(iFich_Prov)
         close(iFich_Prov)
      else
         if (iRand_Seed.LE.0) then
                                ! on prend "n'importe quoi" comme graine...
            iRand_Seed = time_MD()
            iRand_Seed = mod(iRand_Seed,7)*6003 +
     $           mod(iRand_Seed,9)*5752 + mod(iRand_Seed,13)*2904 +
     $           mod(iRand_Seed,23)*2000 +
     $           mod(iRand_Seed,6777) + mod(iRand_Seed,34)
            iRand_Seed = 2*iRand_Seed + 1
            open(iFich_prov,file='iRand_Seed.asc')
            write(iFich_prov,'(I9)')iRand_Seed
            close(iFich_prov)
         end if
         call init_random(iRand_Seed)
      end if
c
c---- Lire le fichier de chauffage par binaires 3-corps (facultatif)
c
#if (_CHAUFF_BIN_3C_>10)
      if (ReOuvrirLectureFichXDR(iFich_XDR,TypeCB3c)) then
         call LitCB3c_XDR(iFich_XDR)
         ixdr_err = ixdrclose(iFich_XDR)
      else
         call Init_ChauffBin3c
      end if
#endif
c
c---- Lire le fichier de parametres pour la loi de croissance forcee
c     du TN (facultatif)
c
#if (_FORCER_CROISSANCE_TN_>0) 
      inquire(file=NomFichCroissForcTN,exist=lExiste)
      if (.NOT.lExiste) then
         write(0,*)
     $        'Le Fichier ',NomFichCroissForcTN,' n''existe pas !!!'
         call Meurt('Init_Amas')
      end if
      open(iFich_Prov,file=NomFichCroissForcTN,status='old')
      call Lire_Param_TN_Croiss_Forc(iFich_Prov)
      close(iFich_Prov)
#endif
c
c---- Lire le fichier de parametrisation de la loi de probabilite
c     de tirage selon Henon (facultatif)
c
      if (ReOuvrirFichBin(iFich_Prov,TypeTH)) then
         call Lit_TH_bin(iFich_Prov)
         close(iFich_Prov)
      else
         call Calc_Param_Tirage_Henon
      end if
c
c---- Lire le fichier des temps (facultatif)
c
      if (ReOuvrirLectureFichXDR(iFich_XDR,TypeTEMPS)) then
         call LitTemps_XDR(iFich_XDR)
         ixdr_err = ixdrclose(iFich_XDR)
         flag_temps = 2
      else
         call Init_Temps
      end if
c
c---- File containing cumulative collision probabilities
c
#if (_WRITE_CUMUL_COLL_PROB_>0)
      if (i_ini>0) then
         if (ReOuvrirFichASC(iFich_Prov,TypeCumPcoll)) then
            call Read_CumulPcoll_ASC(iFich_Prov,iOK)
            close(iFich_Prov)
            if (iOK.NE.0) then
               write(0,*) '!!! Could not read cumulative collision',
     $              ' probability file !!!'
               call Meurt('Init_Amas')
            end if
         else
            write(0,*) '!!! Could not open cumulative collision',
     $           ' probability file !!!'
            call Meurt('Init_Amas')
         end if
      else
         call Init_CumulPcoll
      endif
#endif
c
c---- Lire le fichier de la grille de collisions SPH
c
#if (_COLLISIONS_>0 && _TYPE_COLL_>=10) /* collisions selon resultats de MES simulations SPH */
      inquire(file=NomFichGrilleSPH,exist=lExiste)
      if (.NOT.lExiste) then
         write(0,*) '!!! Le fichier de grille de collisions SPH',
     $        NomFichGrilleSPH,
     $        ' n''existe pas !!!'
         call Meurt('Init_Amas')
      end if
      open(iFich_Prov,file=NomFichGrilleSPH,status='old')
      call Lire_Grille_Collisions_SPH(iFich_Prov,'ASC')
      close(iFich_Prov)
#endif
c
c---- Initialiser la liste des SE particulieres a suivre
c
#if (_ECRIT_SUIVI_SE_PART_>0)
      call Init_SuiviSE
      inquire(file=NomFichListeSEPart,exist=lExiste)
      if (lExiste) then
         write(0,*) '> Lecture de la liste des SE a suivre'
         open(iFich_Prov,file=NomFichListeSEPart,status='old')
         call Lire_Liste_SE_a_suivre(iFich_Prov)
         close(iFich_Prov)
      end if 
#endif
c
c---- File listing "subpopulation" to track
c
#if (_FOLLOW_SUBPOP_>0)
      inquire(file=NomFichListSubPop,exist=lExiste)
      if (lExiste) then
         write(0,*) '> Reading list for subpopulation to track'
         open(iFich_Prov,file=NomFichListSubPop,status='old')
         call Read_ListSubPop(iFich_Prov)
         close(iFich_Prov)
      end if 
#endif
c
c---- File listing non-interacting stellar type combinations
c
      inquire(file=NomFichNonInterac,exist=lExiste)
      if (lExiste) then
         write(0,*) '> Reading prescriptions for non-interacting',
     $        ' stellar type combinations'
         open(iFich_Prov,file=NomFichNonInterac,status='old')
         call Read_NonInterac_Prescriptions(iFich_Prov)
         close(iFich_Prov)
      end if 
c
c---- On-the-fly grids
c
#if (_WRITE_OFGRID_>0)
      call Reinit_OFGrid
#endif

c
c---- Initialisations des variable de controle des SE les plus internes
c                              
c      if (flag_Tri.NE.1) then   ! Trier si necessaire
c         call TriAmas
c         flag_Tri = 1
c      end if
c      iSE_Rmin_CC=iListeRay(1)
c      Rmin_CC=R_SE(iSE_Rmin_CC)
c      Rang_CC=Fact_Rang_CC*NbSE
c
c---- initialisation de la routine de calcul de la repartition de la masse
c     accretee entre TN et "reservoir"
c
      call Repartir_Maccr_rsrv_TN(.TRUE.)
c
c---- Make sure the parameters for stellar evolution are initialized
c     and write them to disk for reference
c
      call InitEvolStell
      open(iFich_Prov,file='_test_ParamStellEvol.asc')
      call write_ParamEvolStell(iFich_Prov)
      close(iFich_Prov)
c
c---- Various tabulated relations are written to disk for testing and reference 
c     purposes
c
#if (_EVOL_STELL_==1) /* Simple stellar evol (MS-->remnant) */
      open(iFich_Prov,file='_test_MRrelMS.asc')
      do i=1,1000
         M_=10**(-3.5d0+i*0.001d0*6.0d0)
         write(iFich_Prov,'(2(E12.5," "))') M_,RelMR(M_,-1.0d0,
     $        TypeEtoile_MS)
      end do
      close(iFich_Prov)

      open(iFich_Prov,file='_test_TMS.asc')
      do i=1,1000
         M_=10**(-3.5d0+i*0.001d0*7.5d0)
         write(iFich_Prov,'(2(E12.5," "))') M_,T_MS(M_)
      end do
      close(iFich_Prov)

      open(iFich_Prov,file='_test_Remn.asc')
      do i=1,1000
         M_=10**(-3.5d0+i*0.001d0*6.0d0)
         call Remnant(M_,Mr_,iTr_)
         write(iFich_Prov,'(2(E12.5," "),I2)') M_,Mr_,iTr_
      end do
      close(iFich_Prov)
#endif 
#if (_TEST_RENC_QUANTITIES_>0)
      write(33,
     $     '("# 1:R_renc 2:Dens_renc 3:Met_moy_renc 4:Sigma1D_renc")')
#endif
c
      end
c
c----------------------------------------------------------------------
      subroutine Sauvegarde_Init
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c constantes :
c ^^^^^^^^^^^^
      logical lEcrase,lConserve_Sauv
      parameter (lEcrase=.FALSE.,lConserve_Sauv=.TRUE.)
      integer iConserve, iEcrase
      parameter (iConserve=5, iEcrase=3)
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      logical lExiste
      integer iAction
c.....................................................................
c
c---- Que faut-il faire avec les fichiers existant deja?
c
      if (i_ini.LE.0) then
         iAction = iEcrase ! on ecrase les fichiers sequentielles existants
      else 
         iAction = iConserve
      end if
c
c---- d'abord, initialiser certain fichiers "sequentiels"
c
      ! "LogBook" indiquant, entre autres, la liste des sauvegardes

      call OuvrirFichASC(iFich_Log,TypeLOG,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > Fichier "LogBook" deja present'
      end if                    ! ecrire l'en-tete du Logbook

      ! Fichier d'evolution des rayons lagrangiens
#if (_ECRIT_RAYLAG_>0)
                                ! Fichier specifiant les fraction de masse pour les rayons lagrangiens
      call OuvrirFichASC(iFich_Prov,TypeFRAYLAG,2,lExiste)
      if (lExiste) then
         write(0,*) ' > Fichier de fractions pourles rayons lag.'
         call Lire_FracM_RayLag(iFich_Prov)
      else
         write(0,*) ' > Pas de fichier de fractions pourles rayons lag.'
      end if
      close(iFich_Prov)

      call OuvrirFichASC(iFich_RayLag,TypeRAYLAG,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > Fichier des rayons lagrangiens deja present'
      else                      ! initialisation du fichier des rayons lagrangiens
         call Init_RayLag_ASC(iFich_RayLag)
      end if
#endif
#if (_ECRIT_LAGQUANT_>0)
                                ! Lagrange radii and various quantities summed inside Lag spheres
#if (_ECRIT_RAYLAG_==0)
      call OuvrirFichASC(iFich_Prov,TypeFRAYLAG,2,lExiste)
      if (lExiste) then
         write(0,*) ' > Fichier de fractions pourles rayons lag.'
         call Lire_FracM_RayLag(iFich_Prov)
      else
         write(0,*) ' > Pas de fichier de fractions pourles rayons lag.'
      end if
      close(iFich_Prov)
#endif

      call OuvrirFichASC(iFich_LagQuant,TypeLagQuant,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > File for Lagrange quantities already present'
      else                      ! initialisation du fichier des rayons lagrangiens
         call Init_LagQuant_ASC(iFich_LagQuant)
      end if
#endif
      ! Fichier d'evolution des anisotropies lagrangiens
#if (_ECRIT_ANILAG_>0)
      call OuvrirFichASC(iFich_Anilag,TypeANILAG,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*)
     $        ' > Fichier des anisotropies lagrangiennes deja present'
      else                      ! initialisation du fichier des anis.s lagrangiens
         call Init_Anilag_ASC(iFich_Anilag)
      end if
#endif
      ! Fichier d'evolution du spectre de masses
#if (_ECRIT_SPECMASSE_>0)
      call OuvrirFichASC(iFich_MSpec,TypeMSpec,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*)
     $        ' > Fichier de spectres de masse deja present'
      else                      ! initialisation du fichier du spectre de masse
         call Init_MassSpec_ASC(iFich_MSpec)
      end if
#endif
      ! Fichier de suivi de la segregation des masses
#if (_ECRIT_SEGR_>0)
                                ! Fichier specifiant les fraction de masse pour la segregation
      call OuvrirFichASC(iFich_Prov,TypeFSEGRLAG,2,lExiste)
      if (lExiste) then
         write(0,*)
     $        ' > Fichier de fractions Lagrangiennes pour la segregation'
         call Lire_FracM_Segr(iFich_Prov)
      else
         write(0,*) ' > Pas de fichier de fractions Lagrangiennes',
     $        ' pour la segregation'
      end if
      close(iFich_Prov)

      call OuvrirFichASC(iFich_Segr,TypeSEGR,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*)
     $        ' > Fichier de suivi de la segr. des masses deja present'
      else                      ! initialisation du fichier de suivi de la segregation des masses
         call Init_SuiviSegr_ASC(iFich_Segr)
      end if
#endif
      ! Fichier des donnees globales de relaxation
#if (_ECRIT_RELAX_>0)
      call OuvrirFichASC(iFich_Relax,TypeRELAX,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > Fichier des donnees de relaxation deja present'
      else
         call Init_FichRelax_ASC(iFich_Relax)
      end if
#endif
      ! Fichier des donnees de dechirement/disparitions ("Loss Cone")
#if ((_DECHIREMENTS_>0 || _DISPARITIONS_>0 || _CAPT_GW_>0) && _ECRIT_LC_>0)
      call OuvrirFichASC(iFich_LC,TypeLC,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > Fichier des donnees de dechirement/',
     $        'disparitions deja present'
      else
         call Init_FichLC_ASC(iFich_LC)
      end if
#endif
      ! Fichier des donnees complementaires sur les captages sur orbites
      ! a forte emission d'ondes grav.
#if (_CAPT_GW_ && _ECRIT_CAPT_GW_>0)
      call OuvrirFichASC(iFich_CaptGW,TypeCaptGW,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > Fichier des donnees de captage par ondes grav.',
     $        ' deja present'
      else
         call Init_FichCaptGW_ASC(iFich_CaptGW)
      end if
#endif
      ! Fichier des donnees des collisions
#if (_COLLISIONS_>0 && _ECRIT_COLL_>0)
      call OuvrirFichASC(iFich_Coll,TypeCOLL,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > Fichier des donnees des collisions deja present'
      else
         call Init_FichColl_ASC(iFich_Coll)
      end if
#endif
      ! Fichier des donnees des evaporations
#if (_ECRIT_EVAP_>0)
      call OuvrirFichASC(iFich_Evap,TypeEvap,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > Fichier des donnees d''evaporation deja present'
      else
         call Init_FichEvap_ASC(iFich_Evap)
      end if
#endif
      ! Fichier des donnees sur les echelles de temps
#if (_ECRIT_EVAP_>0)
      call OuvrirFichASC(iFich_Tscale,TypeTscale,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > Fichier des echelles de temps deja present'
      else
         call Init_FichTscale_ASC(iFich_Tscale)
      end if
#endif
      ! Data about "unusually" strong encounters (large-angle super-encounters)
#if (_MONITOR_STRONG_ENCOUNTERS_>0)
      call OuvrirFichASC(iFich_StrEnc,TypeStrEnc,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > Fichier pour rencontres "fortes"',
     $        ' deja present'
      else
         call Init_FichStrEnc_ASC(iFich_StrEnc)
      end if
#endif
      ! Data about exception cases in collision treatment
#if (_WRITE_COLL_EXCEPTIONS_>0)
      call OuvrirFichASC(iFich_ExcColl,TypeExcColl,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > File for collision "exceptions"',
     $        ' already present'
      else
         continue
      end if
#endif
      ! Fichier de suivi de SE particulieres
#if (_ECRIT_SUIVI_SE_PART_>0)
      call OuvrirFichASC(iFich_SuivSE,TypeSuivSE,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*)
     $        ' > Fichier de suivi de SE particulieres deja present'
      else
         call Init_FichSuiviSE_ASC(iFich_SuivSE)
      end if
#endif
      ! Fichier de suivi des differents types stellaires
#if (_ECRIT_SUIVI_TYPES_STELL_>0)
                                ! Fichier specifiant les fraction de masse
      call OuvrirFichASC(iFich_Prov,TypeFSTSLag,2,lExiste)
      if (lExiste) then
         write(0,*)
     $        ' > Fichier de fractions Lagrangiennes ',
     $        'pour les differents types stellaires'
         call Lire_FracM_SuiviTypesStell(iFich_Prov)
      else
         write(0,*) ' > Pas de fichier de fractions Lagrangiennes ',
     $        'pour les differents types stellaires'
      end if
      close(iFich_Prov)

      call Init_SuiviTypesStell
#endif
      ! Fichier d'evolution des rayons lagrangiens
#if (_MONITOR_MBINS_>0)
                                ! File for monitoring distribution of particles binned by stellar mass
      call OuvrirFichASC(iFich_Prov,TypeMbins,2,lExiste)
      if (lExiste) then
         write(0,*) ' > File specifying stellar mass bins'
         call Read_Mbins(iFich_Prov)
      else
         write(0,*) ' > No file specifying stellar mass bins'
      end if
      close(iFich_Prov)

      call OuvrirFichASC(iFich_MbinRad,TypeMbinRad,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > File to monitor Lagr. radii by mass bins present'
      else                      ! initialisation du fichier des rayons lagrangiens
         call Init_MbinsMonitor_ASC(iFich_MbinRad)
      end if
#endif
      ! Fichier de duplications
#if (_DUPLICATION_>0 && _ECRIT_DUPLIC_>0)
      call OuvrirFichASC(iFich_Duplic,TypeDuplic,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*)
     $        ' > Fichier de suivi de duplications deja present'
      else
         call Init_FichDuplic_ASC(iFich_Duplic)
      end if
#endif    
      ! Fichier des donnees d'evolution stellaire
#if (_EVOL_STELL_>0 && _ECRIT_STEVOL_>0)
      call OuvrirFichASC(iFich_StEvol,TypeStEvol,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > Fichier des donnees d''evolution ',
     $        'stellaire deja present'
      else
         call Init_StEvol_ASC(iFich_StEvol)
      end if
#endif
      ! File to keep track of natal kicks
#if (_NATKICKS_>0 && _ECRIT_NATKICKS_>0)
#   if (_EVOL_STELL_<=0 || _ECRIT_STEVOL_<=0)
      ERROR: INCONSISTENT 
#   endif
      call OuvrirFichASC(iFich_NatKcks,TypeNatKcks,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*) ' > Fichier des donnees de "natal kicks" ',
     $        'stellaire deja present'
      else
         call Init_NatKicks_ASC(iFich_NatKcks)
      end if
#endif
      ! Fichier-journal des sauvegardes demandees de l'exterieur
      call OuvrirFichASC(iFich_LogSauv,TypeLogSauv,iAction,lExiste)
      if (lExiste .AND. (iAction.EQ.iConserve)) then
         write(0,*)
     $        ' > Fichier-journal des sauvegardes deja present'
      else
         call Init_FichLogSauv_ASC(iFich_LogSauv)
      end if
c
c---- faire toutes les sauvegardes sans ecraser les fichiers deja existants
c
      call SauvegardePart_EvolAmas(lEcrase)
      call Sauvegarde_EvolAmas(lEcrase,lConserve_Sauv)
c
      end
c
c----------------------------------------------------------------------
      subroutine Check_Amas
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      integer iNbErrors,iNbErrors_Serious
c......................................................................
c
      iNbErrors = 0
      open(iFich_Prov,file='ErrorsCheckAmas.asc')
      call CheckAmas(iFich_Prov,iNbErrors,iNbErrors_Serious)
      if (iNbErrors.GT.0) then
         close(iFich_Prov)
         write (0,210) iNbErrors,iNbErrors_Serious
 210     format('!!! ',I6,' Error(s) found in cluster data',
     $        ' including ',I6,' serious one(s)',
     $        ' (see ErrorsCheckAmas.asc) !!!')
         if (iNbErrors_Serious.GT.0) call Meurt('Check_Amas')
      else
         close(iFich_Prov,STATUS='DELETE')
      end if
c
      end
c      
c----------------------------------------------------------------------
      subroutine Calc_EnerGlob
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_LimExt_common.f'
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer i,iSE
      double precision M_i,R_i,Mint,Mprec
c......................................................................
c
                                ! Trier si necessaire
      if (flag_Tri.NE.1) then
         call TriAmas
         flag_Tri = 1
      end if

      Ecin_amas = 0.0d0
      Etot_amas = 0.0d0
      Egrav_amas = 0.0d0
      Estell_grav_amas = 0.0d0
      Mint = 0.0d0
      Mprec = 0.0d0
      do i=1,NbSE-NbSE_perte
         iSE = iListeRay(i)
         M_i = M_SE(iSE)
         R_i = R_SE(iSE)
         Ecin_amas = Ecin_amas + M_i*T_SE(iSE)
         Mint = Mint+Mprec
         Estell_grav_amas = Estell_grav_amas -
     $        M_i/R_i*(Mint+FracMasseInt*M_i)
         Egrav_amas = Egrav_amas - M_i/R_i*(Mint+FracMasseInt*M_i+M_ctr)
         Mprec = M_i
      end do
      Mint = Mint+M_i
      ! Egrav_amas = Egrav_amas +B_LimExt*(Mint+M_ctr) ! Modified 14 nov 2005
      Etot_amas  = Ecin_amas + Egrav_amas 
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_ConvUnites
c----------------------------------------------------------------------
c     Doit etre appele apres que les grandeurs _ini aient ete lues
c     ou initialisees
c
      include 'EvolAmas_inc.f'
      include 'EvolAmas_LimExt_common.f'
      include 'PhysQuantities.f'
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      double precision Mtot,Mtot_en_Msol,rNet_tot
      integer i
c
c functions :
c ^^^^^^^^^^^
      double precision stellar_mass
c......................................................................
c
      Mtot_en_Msol = 0.0d0
      Mtot = 0.0d0
      do i=1,NbSE
         Mtot_en_Msol =  Mtot_en_Msol + stellar_mass(i)
         Mtot         =  Mtot + M_SE(i)
      end do
      if (rNbEtoiles.LE.0.0d0) then ! Le nb d'etoiles n'a pas pu etre lu ou initialise...

#if (_DUPLICATION_==0)
c     Si il n'y a pas de duplication, le nb d'etoiles suit le nb de SE
         rNbEtoiles = rNbEtoiles_ini*(NbSE-NbSE_Perte)/dble(NbSE)
#else
c     Si il n'y a pas de collisions ni d'evolution stellaire, la masse des etoiles est conservee
#if ( (_COLLISIONS_==0 || (_TYPE_COLL_==0 || _TYPE_COLL_==1 || _TYPE_COLL_==3)) &&  _EVOL_STELL_==0 && _ACCR_BY_STARS_==0)
         rNbEtoiles = rNbEtoiles_ini * Mtot/MSEtot_ini *
     $        (NbSE-NbSE_Perte)/dble(NbSE)
#else
c     si il n'y a pas eu de duplication, le nb d'etoiles suit le nb de SE
         if (NbSE.EQ.NbSE_ini) then
            rNbEtoiles = rNbEtoiles_ini*(NbSE-NbSE_Perte)/dble(NbSE)
         else            
            write(0,*) '!!! Impossible de determiner rNbEtoiles !!!'
            call Meurt('Calc_ConvUnites')
         end if
#endif /* ( (_COLLISIONS_==0 || (... */
#endif /* (_DUPLICATION_==0) */
      end if

      Net_par_SE   = rNbEtoiles/dble(NbSE-NbSE_Perte)
      Mtot_en_Msol = Net_par_SE*Mtot_en_Msol
      Mtot         = Mtot
      Coef_ConvM   = Mtot/Mtot_en_Msol 
      Coef_ConvR   = Rsol_en_pc/TailleAmas_en_pc
c
      UM_Ncorps = _M_INI_
      UL_Ncorps = - _M_INI_**2/(2.0d0*Estell_grav_ini)
      UT_FP     = sqrt(_M_INI_**5/(-2.0d0*_EST_INI_)**3) * ! Maybe wrong
     $     _N_INI_/
     $     LOG(Gamma_relax*_N_INI_)
                                ! Vitesse de la lumiere au carre dans les
                                ! unites du code (pour effets relativistes)
      VitLum2 = Light_c**2* Parsec/Sun_GM *
     $     TailleAmas_en_pc/(_M_INI_/Coef_ConvM)
                                ! Une annee en unites du code (pour evolution stellaire)
      UneAnnee = SideralYear / sqrt(Parsec**3/Sun_GM) /
     $     sqrt( TailleAmas_en_pc**3/(_M_INI_/Coef_ConvM) )
      UnKmS = 1.0d5 * sqrt( Parsec/Sun_GM *
     $     TailleAmas_en_pc/(_M_INI_/Coef_ConvM) )
c
      open(iFich_Prov,file='Units.asc')
      write(iFich_Prov,100) Mtot, _M_INI_, _N_INI_, Mtot_en_Msol,
     $     Coef_ConvM, Coef_ConvR, UM_Ncorps, UL_Ncorps, UT_FP,
     $     VitLum2, UneAnnee, UnKmS
 100  format(
     $     ' Tot mass in code units,                  Mtot         = ',
     $     E12.5,/,
     $     ' Tot init. mass in code units,            _M_INI_      = ',
     $     E12.5,/,
     $     ' Tot init. star numbers,                  _N_INI_      = ',
     $     E12.5,/,
     $     ' Tot mass in Msun,                        Mtot_en_Msol = ',
     $     E12.5,/,
     $     ' Solar mass in code units,                Coef_ConvM   = ',
     $     E12.5,/,
     $     ' Solar radius in code units,              Coef_ConvR   = ',
     $     E12.5,/,
     $     ' N-body mass unit, in code units,         UM_Ncorps    = ',
     $     E12.5,/,
     $     ' N-body length unit, in code units,       UL_Ncorps    = ',
     $     E12.5,/,
     $     ' Fokker-Pl. time unit, in code units,     UT_FP        = ',
     $     E12.5,/,
     $     ' Square of speed of light, in code units, VitLum2      = ',
     $     E12.5,/,
     $     ' One year, in code units,                 UneAnnee     = ',
     $     E12.5,/,
     $     ' One km/s, in code units,                 UnKmS        = ',
     $     E12.5
     $     )
      close(iFich_Prov)
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_RelaxGlob
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c......................................................................
c
                                ! Trier si necessaire
      if (flag_Tri.NE.1) then
         call TriAmas
         flag_Tri = 1
      end if
      call CalcGlobRelax
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_FichLogSauv_ASC(iFich)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
      write(iFich,100)
 100  format(
     $     '# Fichier-journal des sauvegardes demandees',
     $     ' de l''"exterieur"')
c
      end
c
c----------------------------------------------------------------------
      subroutine Evol_Amas
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'Flag_common.f'
      include 'EvolAmas_NbPas_common.f'
c
c constante(s) :
c ^^^^^^^^^^^^^^
      logical lEcrase
      parameter (lEcrase=.TRUE.)
      double precision fact_N_paquet
      parameter (fact_N_paquet=0.01d0)
c
c var locales :
c ^^^^^^^^^^^^^
      integer*8 i_inf,
     $     iPasDetEvap_suiv, iPasInfo_suiv, iPasRecalcTH_suiv,
     $     iPasRecArbre_suiv, iPasRecCB3c_suiv, iPasRecRelax_suiv,
     $     iPasRecTpsAmas_suiv, iPasSauv_suiv, iPasTestDemSauv_suiv,
     $     iPasTestDuplic_suiv, iPasSauvPart_suiv
      integer iSauvDemandee, i, N_paquet
      logical lConserve_Sauv
c
c fonctions :
c ^^^^^^^^^^^
      logical lTest_Duplic, lRepartir_Maccr_rsrv_TN
      double precision random
c......................................................................
c
      lFlag = .FALSE.
      lFlag2 = .FALSE.

      N_paquet = fact_N_paquet*NbSE_ini + 0.5d0

      if (FacNbPasInfo    .LT.1.0d30)  then
         NbPasInfo      = FacNbPasInfo    *NbSE_ini + 0.5d0
      else
         NbPasInfo      = i_Fin+1
      end if
      if (FacNbPasSauv .LT.1.0d30)  then
         NbPasSauv  = FacNbPasSauv*NbSE_ini + 0.5d0
      else
         NbPasSauv  = i_Fin+1
      end if
      if (FacNbPasSauvPart.LT.1.0d30)  then
         NbPasSauvPart= FacNbPasSauvPart     *NbSE_ini + 0.5d0
      else
         NbPasSauvPart= i_Fin+1
      end if
      if (FacNbPasRecArbre.LT.1.0d30)  then
         NbPasRecArbre  = FacNbPasRecArbre *NbSE_ini + 0.5d0
      else
         NbPasRecArbre  = i_Fin+1
      end if
      if (FacNbPasRecCB3c.LT.1.0d30)  then
         NbPasRecCB3c   = FacNbPasRecCB3c *NbSE_ini + 0.5d0
      else
         NbPasRecCB3c   = i_Fin+1
      end if
      if (FacNbPasRecRelax.LT.1.0d30)  then
         NbPasRecRelax  = FacNbPasRecRelax *NbSE_ini + 0.5d0
      else
         NbPasRecRelax  = i_Fin+1
      end if
      if (FacNbPasRecTpsAmas.LT.1.0d30)  then
         NbPasRecTpsAmas  = FacNbPasRecTpsAmas *NbSE_ini + 0.5d0
      else
         NbPasRecTpsAmas  = i_Fin+1
      end if
      if (FacNbPasDetTrelExtremes.LT.1.0d30) then
         NbPasDetTrelExtremes = FacNbPasDetTrelExtremes*NbSE_ini + 0.5d0
      else
         NbPasDetTrelExtremes = i_Fin+1
      end if
      if (FacNbPasRecalcTH .LT.1.0d30)  then
         NbPasRecalcTH = FacNbPasRecalcTH*NbSE_ini + 0.5d0
      else
         NbPasRecalcTH = i_Fin+1
      end if
      if (FacNbPasDetEvap .LT.1.0d30)  then
         NbPasDetEvap = FacNbPasDetEvap*NbSE_ini + 0.5d0
      else
         NbPasDetEvap = i_Fin+1
      end if
      if (FacNbPasTestDuplic .LT.1.0d30)  then
         NbPasTestDuplic = FacNbPasTestDuplic*NbSE_ini + 0.5d0
      else
         NbPasTestDuplic = i_Fin+1
      end if
      if (FacNbPasTestDemSauv .LT.1.0d30)  then
         NbPasTestDemSauv = FacNbPasTestDemSauv*NbSE_ini + 0.5d0
      else
         NbPasTestDemSauv = i_Fin+1
      end if
c
c---- Determination des numeros de pas correspondant au differentes "taches exceptionnelles" a faire
c
c      code produit par le script:
c
c #! /bin/zsh -f
c for var in $(cat EvolAmas.F | grep 'iPas.*_suiv' | sed 's/.*\(iPas[a-zA-Z0-9_]*_suiv\).*/\1/'| uniq | sort); do
c     nvar="NbPas"${${var%_suiv}#iPas}
c     cat <<EOF
c       $var = 0
c       do while (iPas_Evol.GE.$var)
c          $var = $var + 
c      $          $nvar
c       end do
c EOF
c done
c
      iPasDetEvap_suiv = 0
      do while (iPas_Evol.GE.iPasDetEvap_suiv)
         iPasDetEvap_suiv = iPasDetEvap_suiv + 
     $        NbPasDetEvap
      end do
      iPasInfo_suiv = 0
      do while (iPas_Evol.GE.iPasInfo_suiv)
         iPasInfo_suiv = iPasInfo_suiv + 
     $        NbPasInfo
      end do
      iPasRecalcTH_suiv = 0
      do while (iPas_Evol.GE.iPasRecalcTH_suiv)
         iPasRecalcTH_suiv = iPasRecalcTH_suiv + 
     $        NbPasRecalcTH
      end do
      iPasRecalcTH_suiv = 0
      do while (iPas_Evol.GE.iPasRecalcTH_suiv)
         iPasRecalcTH_suiv = iPasRecalcTH_suiv + 
     $        NbPasRecalcTH
      end do
      iPasRecArbre_suiv = 0
      do while (iPas_Evol.GE.iPasRecArbre_suiv)
         iPasRecArbre_suiv = iPasRecArbre_suiv + 
     $        NbPasRecArbre
      end do
      iPasRecCB3c_suiv = 0
      do while (iPas_Evol.GE.iPasRecCB3c_suiv)
         iPasRecCB3c_suiv = iPasRecCB3c_suiv + 
     $        NbPasRecCB3c
      end do
      iPasRecRelax_suiv = 0
      do while (iPas_Evol.GE.iPasRecRelax_suiv)
         iPasRecRelax_suiv = iPasRecRelax_suiv + 
     $        NbPasRecRelax
      end do
      iPasRecRelax_suiv = 0
      do while (iPas_Evol.GE.iPasRecRelax_suiv)
         iPasRecRelax_suiv = iPasRecRelax_suiv + 
     $        NbPasRecRelax
      end do
      iPasRecTpsAmas_suiv = 0
      do while (iPas_Evol.GE.iPasRecTpsAmas_suiv)
         iPasRecTpsAmas_suiv = iPasRecTpsAmas_suiv + 
     $        NbPasRecTpsAmas
      end do
      iPasTestDemSauv_suiv = 0
      do while (iPas_Evol.GE.iPasTestDemSauv_suiv)
         iPasTestDemSauv_suiv = iPasTestDemSauv_suiv + 
     $        NbPasTestDemSauv
      end do
      iPasTestDemSauv_suiv = 0
      do while (iPas_Evol.GE.iPasTestDemSauv_suiv)
         iPasTestDemSauv_suiv = iPasTestDemSauv_suiv + 
     $        NbPasTestDemSauv
      end do
      iPasTestDuplic_suiv = 0
      do while (iPas_Evol.GE.iPasTestDuplic_suiv)
         iPasTestDuplic_suiv = iPasTestDuplic_suiv + 
     $        NbPasTestDuplic
      end do
      iPasTestDuplic_suiv = 0
      do while (iPas_Evol.GE.iPasTestDuplic_suiv)
         iPasTestDuplic_suiv = iPasTestDuplic_suiv + 
     $        NbPasTestDuplic
      end do

      iSauv = 0
      iPasSauv_suiv = 0 
      do while (iPasSauv_suiv.LE.iPas_Evol)
         iPasSauv_suiv = iPasSauv_suiv+NbPasSauv
         iSauv=iSauv+1
      end do

      iSauvPart = 0
      iPasSauvPart_suiv = 0 
      do while (iPasSauvPart_suiv.LE.iPas_Evol)
         iPasSauvPart_suiv = iPasSauvPart_suiv+NbPasSauvPart
         iSauvPart=iSauvPart+1
      end do

      call Affiche_Info         ! infos affichees a l'ecran

      iPas_Evol=i_ini
      iSauvDemandee=0

      do while (iPas_Evol.LE.i_fin .AND. iSauvDemandee.LT.9)
         
         !lFlag = (iPas_Evol.EQ.204) !######

         Flag_Tri     = 0
         Flag_Grille  = 2
         Flag_Etoiles = 2
         Flag_Temps   = 2
         Flag_Arbre   = 2

#if (_PAR_PAQUETS_>0) 
         do i=1,N_paquet
#endif
c****************************************
c****************************************
c***      LE PAS D'EVOLUTION !!!      ***
c****************************************
c****************************************
#if (_AJUST_PASSIF_>0)
            if (random().LT.FracPas_AjustPassif) then
               call Ajuster_Amas_1pas
            else
               call Evoluer_Amas_1pas
            end if
#else
            call Evoluer_Amas_1pas
#endif
            iPas_Evol=iPas_Evol+1
c****************************************
c****************************************
#if (_PAR_PAQUETS_>0)
            if (flag_Grille.LT.0) goto 9 ! il ne faut pas laisser la grille se vider!
         end do
 9       continue
#endif

         if (iPas_Evol.GE.iPasDetEvap_suiv)
     $        then              ! on determine quelles sont les particules liees
            call Determiner_SE_liees
            do while (iPas_Evol.GE.iPasDetEvap_suiv)
               iPasDetEvap_suiv = iPasDetEvap_suiv+NbPasDetEvap
            end do
         end if
         if (lRepartir_Maccr_rsrv_TN())
     $        then              ! on determine la croissance du TN central
            call Repartir_Maccr_rsrv_TN(.FALSE.)
         end if
#if (_DUPLICATION_>0)
         if (iPas_Evol.GE.iPasTestDuplic_suiv .AND. lTest_Duplic())
     $        then              ! On duplique toutes les SE si necessaire
            call Dupliquer_les_SE
            do while (iPas_Evol.GE.iPasTestDuplic_suiv)
               iPasTestDuplic_suiv = iPasTestDuplic_suiv+NbPasTestDuplic
            end do
         end if
#endif
         if (flag_Grille.LT.0)  ! la grille a deborde, il faut
     $        then              ! la reconstruire
            call Init_Grille
            iNbRecGrille = iNbRecGrille+1 
         end if
#if (_CHAUFF_BIN_3C_>10)
         if (iPas_Evol.GE.iPasRecCB3c_suiv)
     $        then              ! on recalcule les parametres de chauffage par binaires 3-corps
            call Init_ChauffBin3c
            do while (iPas_Evol.GE.iPasRecCB3c_suiv)
               iPasRecCB3c_suiv = iPasRecCB3c_suiv+NbPasRecCB3c
            end do
         end if
#endif
         if (iPas_Evol.GE.iPasRecArbre_suiv)
     $        then              ! on reconstruit l'arbre
            call Init_Arbre
            iNbRecArbre = iNbRecArbre+1 
            do while (iPas_Evol.GE.iPasRecArbre_suiv)
               iPasRecArbre_suiv = iPasRecArbre_suiv+NbPasRecArbre
            end do
         end if
         if (iPas_Evol.GE.iPasRecRelax_suiv)
     $        then              ! calcul des parametres de relaxation
            call CalcParam_Relax
            iNbPasRecRelax = iNbPasRecRelax+1
            do while (iPas_Evol.GE.iPasRecRelax_suiv)
               iPasRecRelax_suiv = iPasRecRelax_suiv+NbPasRecRelax
            end do
         end if
                                ! on s'assure que le temps affiche soit tjs
                                ! correctement calcule
         if (iPas_Evol.GE.iPasRecTpsAmas_suiv)
     $        then              ! determination du temps median et autres
            call Det_TpsAmas
            flag_Temps   = 1
         end if

         if (iPas_Evol.GE.iPasRecalcTH_suiv)
     $        then              ! determination des parametres de tirage
                                ! selon methode de Henon
            if (flag_Tri.NE.1) then
               call TriAmas
               flag_Tri = 1
            end if
            if (flag_Temps.NE.1) then
               call Det_TpsAmas
               flag_Temps  = 1
            end if
            call Calc_Param_Tirage_Henon
            do while (iPas_Evol.GE.iPasRecalcTH_suiv)
               iPasRecalcTH_suiv = iPasRecalcTH_suiv+NbPasRecalcTH
            end do
         end if

         if (iPas_Evol.GE.iPasInfo_suiv .OR. iPas_Evol.GE.iPasSauv_suiv
     $        .OR. iPas_Evol.GE.iPasSauvPart_suiv)
     $        then              ! infos affichees a l'ecran
            if (flag_Temps.NE.1) then
               call Det_TpsAmas
               flag_Temps=1
            end if
            call Affiche_Info
            do while (iPas_Evol.GE.iPasInfo_suiv)
               iPasInfo_suiv = iPasInfo_suiv+NbPasInfo
            end do
         end if

         iSauvDemandee =  0
         if (iPas_Evol.GE.iPasTestDemSauv_suiv)
     $        then              ! Voir si une sauvegarde a ete demandee
            call Tester_Demande_Sauvegarde(iSauvDemandee)
            do while (iPas_Evol.GE.iPasTestDemSauv_suiv)
               iPasTestDemSauv_suiv = iPasTestDemSauv_suiv +
     $              NbPasTestDemSauv
            end do
         end if

         if (iPas_Evol.GE.iPasSauvPart_suiv .OR. iPas_Evol.GE.iPasSauv_suiv .OR. ! sauvegardes de routines ou ...
     $        iSauvDemandee.GT.0 ) then ! sauvegardes demandees de l'exterieur
            call SauvegardePart_EvolAmas(lEcrase)
            do while (iPasSauvPart_suiv.LE.iPas_Evol) ! prochaine sauvegarde de routine
               iPasSauvPart_suiv = iPasSauvPart_suiv+NbPasSauvPart
               iSauvPart=iSauvPart+1
            end do
         end if

         if (iPas_Evol.GE.iPasSauv_suiv .OR. ! sauvegardes de routines ou ...
     $        iSauvDemandee.GT.0 ) then ! sauvegardes demandees de l'exterieur
            lConserve_Sauv = (iSauvDemandee.GT.0 .OR. ! toutes les sauvegardes demandees sont preservees
     $           mod(iSauv,NbSauv_Conserve).EQ.0)
            call Sauvegarde_EvolAmas(lEcrase,lConserve_Sauv)
            do while (iPasSauv_suiv.LE.iPas_Evol) ! prochaine sauvegarde de routine
               iPasSauv_suiv = iPasSauv_suiv+NbPasSauv
               iSauv=iSauv+1
            end do
            call CalcParam_Relax ! pour etre sur de pouvoir reproduire la meme evolution
            do while (iPasSauvPart_suiv.LE.iPas_Evol) ! prochaine sauvegarde de routine
               iPasSauvPart_suiv = iPasSauvPart_suiv+NbPasSauvPart
               iSauvPart=iSauvPart+1
            end do
         end if

         if (flag_Temps.EQ.1 ) then
#if (_VERBOSITUDE_>=5)
            write(0,*) ' && Tps_Amas = ',Tps_Amas
#endif
            do while (iPas_Evol.GE.iPasRecTpsAmas_suiv)
               iPasRecTpsAmas_suiv = iPasRecTpsAmas_suiv+NbPasRecTpsAmas
            end do
         end if

      end do

c
c---- on fait une sauvegarde finale
c
      if (iSauvDemandee.EQ.0) then
         iPas_evol=iPas_evol-1
         lConserve_Sauv = .TRUE. ! sans importance ici....
         call SauvegardePart_EvolAmas(lEcrase)
         iSauvPart=iSauvPart+1
         call Sauvegarde_EvolAmas(lEcrase,lConserve_Sauv) 
         iSauv=iSauv+1
      end if
c
      end
c----------------------------------------------------------------------
      subroutine Tester_Demande_Sauvegarde(iStat)
c----------------------------------------------------------------------
c     Test si une sauvegarde est demandee; Si oui, ecrit dans le fichier
c     de journal des sauvegardes pour indiquer que la sauvegarde
c     demandee est realisee.
c
      implicit none
      include 'EvolAmas_param.f'
      include 'EvolAmas_Compteurs_common.f'
c
c argument :
c ^^^^^^^^^^
      integer iStat             ! est mis a 1 si il faut une sauvegarde; 
                                ! a 9 si il faut arreter; 
c                               ! a 19 en cas de demande d'arret par fichier _STOP_
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      character*256 Chaine
      logical lExiste
      integer iOK
      character*1 tabul
c
c fonction :
c ^^^^^^^^^^
      integer LongChaine
c......................................................................
c
      tabul = char(9)
      iStat=0
      inquire(file=NomFichDemSauv,exist=lExiste)
      if (lExiste) then 
         iStat=1

         open(iFich_Prov,file=NomFichDemSauv)

         write(Chaine,'(I10)') iPas_Evol             ! mumero du pas
         call EnleverBlancs(Chaine)
         write(iFich_LogSauv,'(A)')
     $        'iPas_Evol'      // tabul // Chaine(1:LongChaine(Chaine)) 
         call QuelleDate_str(Chaine)                 ! date
         write(iFich_LogSauv,'(A)')
     $        'date'          // tabul // Chaine(1:LongChaine(Chaine))
         call QuelleHeure_str(Chaine)                ! heure
         write(iFich_LogSauv,'(A)')
     $        'heure'         // tabul // Chaine(1:LongChaine(Chaine))
         read(iFich_Prov,'(A)',END=99,ERR=99) Chaine ! identificateur de la sauvegrde
         write(iFich_LogSauv,'(A)')
     $        'id_sauvegarde' // tabul // Chaine(1:LongChaine(Chaine))
         if (Chaine(1:LongChaine(Chaine)).EQ.'STOP') iStat=9 ! Demande d'arret de la simulation
         read(iFich_Prov,'(A)',END=9,ERR=9) Chaine   ! condition
         goto 10
 9       Chaine='_inconnue_'
 10      continue
         write(iFich_LogSauv,'(A)')
     $        'condition'     // tabul // Chaine(1:LongChaine(Chaine)) 
         read(iFich_Prov,'(A)',END=19,ERR=19) Chaine ! commentaire
         goto 20
 19      Chaine='_aucun_'
 20      continue
         write(iFich_LogSauv,'(A)')
     $        'commentaire'   // tabul // Chaine(1:LongChaine(Chaine)) 
         write(iFich_LogSauv,*)

         close(iFich_Prov,status='delete',IOSTAT=iOK)
         if (iOK.NE.0) then
            write(0,*) '!!! Impossible de supprimer le fichier ',
     $           NomFichDemSauv,' !!!'
            call Meurt('Tester_Demande_Sauvegarde')
         end if

         call flush_MD(iFich_LogSauv)

      end if
                                ! Y a-t-il un fichier d'arret ?
      lExiste=.FALSE.
      inquire(file=NomFichDemArret,exist=lExiste)
      if (lExiste) then 
         iStat=19
         open(iFich_Prov,file=NomFichDemArret)
         close(iFich_Prov,status='delete',IOSTAT=iOK)
         if (iOK.NE.0) then
            write(0,*) '!!! Impossible de supprimer le fichier ',
     $           NomFichDemArret,' !!!'
            call Meurt('Tester_Demande_Sauvegarde')
         end if
      end if

      !write(0,*) '@@@ Tester_Demande_Sauvegarde: iStat= ',iStat
      return
                                ! impossible de lire l'identificateur
 99   continue
      write(0,*) '!!! Impossible de lire l''identificateur ',
     $     'de sauvegarde dans ',NomFichDemSauv,' !!!'
      call Meurt('Tester_Demande_Sauvegarde')
c
c
      end 
c
c----------------------------------------------------------------------
      logical function lRepartir_Maccr_rsrv_TN()
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_Etoiles_param.f'
c
c var locales :
c ^^^^^^^^^^^^^
      double precision limite_delta_Mctr_Npart, Mctr_prec,
     $     limite_delta_Mctr
      save limite_delta_Mctr_Npart, Mctr_prec
      logical lPrem
      save lPrem
      data lPrem /.TRUE./
c......................................................................
c
      lRepartir_Maccr_rsrv_TN = .FALSE.
      if (iType_CtrObj.EQ.TypeEtoile_MS) return 
c
c----  initialisations
c
      if (lPrem) then
         limite_delta_Mctr_Npart = Fac_dt_accr_Npart*Mamas_ini/NbSE_ini
         Mctr_prec = M_ctr
         lPrem = .FALSE.
      end if

      limite_delta_Mctr = max(1.0d-8,
     $     min(Fac_dt_accr_M*M_ctr,limite_delta_Mctr_Npart))
      if ((M_ctr-Mctr_prec).GT.limite_delta_Mctr) then
         lRepartir_Maccr_rsrv_TN = .TRUE.
         Mctr_prec = M_ctr
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Repartir_Maccr_rsrv_TN(lInit)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c argument :
c ^^^^^^^^^^
      logical lInit ! si .TRUE., on initialise
c
c parametres :
c ^^^^^^^^^^^^
      double precision Delta_T_min
      parameter (Delta_T_min=1.0d-8)
c
c var locales :
c ^^^^^^^^^^^^^
      double precision Tps_prec, T_Eddington, M_TN_prec
      save Tps_prec, T_Eddington, M_TN_prec
      double precision st_sur_mp, Delta_T, Delta_Mmax, Delta_M_TN
      logical lOK
c......................................................................
c
      if (lInit) then
         st_sur_mp =            ! section Thomson / masse proton en unites du code !!!
     $        ( sigma_th_en_cgs / mp_en_cgs ) *
     $       ( (Msol_en_g/Coef_ConvM) / (TailleAmas_en_pc*pc_en_cm)**2 )
         T_Eddington =           ! Temps de croissance du TN au taux d'accretion maximal selon lum. d'Eddington
     $        effic_conv_lum*sqrt(VitLum2) *
     $        st_sur_mp / (r4Pi*Mu_mol_elec)
         if (effic_conv_lum.LE.0.0d0) T_Eddington = 0.0d0 ! croissance du TN non limitee
         call open_append(iFich_Prov,'divers.asc',lOK)
         write(iFich_Prov,100)T_Eddington
 100     format(
     $        '# Temps de croissance du TN au taux',
     $        ' d''accretion maximal selon lum. d''Eddington',
     $        ' (unites du code)',/,
     $        ' T_Eddington = ',E12.5)
         call flush_MD(IFich_Prov)
         close(IFich_Prov)
         Tps_prec = Tps_Amas
         M_TN_prec = M_TN
      else
         call Det_TpsAmas
         flag_Temps = 1 
         Delta_T = Tps_Amas-Tps_prec
#if (_SURV_CROISS_TPS_>0)
         if (Delta_T.LT.0.0d0) then
            write(0,*) '!!! Tps_Amas, Tps_prec = ',Tps_Amas,Tps_prec,
     $           ' !!!'
            call Meurt('Repartir_Maccr_rsrv_TN')
         end if
#endif
         if (Delta_T.GT.Delta_T_min) then

            if (T_Eddington.GT.0.0d0) then
               Delta_Mmax = M_TN_prec*(exp(Delta_T/T_Eddington)-1.0d0)
               Delta_M_TN = min(Delta_Mmax, M_ctr-M_TN)
c                               ^^^^\/^^^^  ^^^^\/^^^^
c                                    \           \
c                               Accroissement max \
c                               selon Eddington    A disposition 
c                                                 dans le reservoir
            else
               Delta_M_TN = M_ctr-M_TN
            end if
            M_TN = M_TN + Delta_M_TN
            Tps_prec = Tps_Amas
            M_TN_prec = M_TN

         end if

      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Evoluer_Amas_1pas
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
c
c var locales:
c ^^^^^^^^^^^^
      double precision dM
c
c fonction :
c ^^^^^^^^^^
      double precision M_TN_Croiss_Forc, random
c......................................................................
c
c---- Determiner la masse du TN (si croissance adiabatique)
c
#if (_PRESENCE_TN_>0 && _FORCER_CROISSANCE_TN_>0)
      dM = -M_TN
      M_TN = M_TN_Croiss_Forc(Tps_amas)
      dM = dM + M_TN
      M_ctr = M_ctr + dM
#endif
c
c---- Check for evolution of central object
c
#if (_PRESENCE_TN_>0 && _EVOL_STELL_>0)
      call Evolve_CentralObject
#endif
c
c---- Choisir la paire de SEs a faire interagir
c
      call Choisir_Paire_SE     ! determine iRang_SE1/2 iSE1/2 (parcours arbre)
c
c---- Calculer la rencontre (collision ou relaxation)
c
      call Calculer_Rencontre_Paire_SE
c
c---- Traiter l'evolution stellaire
c
#if (_EVOL_STELL_>0 || _ACCR_BY_STARS_>0)
                                ! Escape is tested (and should be prevented) in
                                ! Calculer_Rencontre_Paire_SE but if the masses
                                ! are modified now, it can change things because
                                ! the new mass of the "other" particle in the pair
                                ! affects the potential for a particle
      call Calculer_EvolStell_Paire_SE 
#endif

#if (_EVOL_UNE_SE_>0) /* on ne fait evoluer qu'une des SE !!! */
      Num_SE_Act = 1+int(2.0d0*random())
#else
      do Num_SE_Act=1,2
#endif
c
c---- Definir la SE active
c
         call DonneesPR_a_SEact
c
c---- Extraire la SE
c
         call Extraire_SEact
c
c---- Modifier la SE active en fonction des resultats de la rencontre
c     Cette phase inclut une possible "interaction" avec le trou noir
c     et le tirage de la position orbitale
c
         call Modifier_SEact
c
c---- Re-inserer la SE (et reconstruire la grille si necessaire)
c
         call Reinserer_SEact
c
c---- Mettre a jour les donnees de la SE dans les tableaux de la paire
c     (par precaution)
c
         call DonneesSEact_a_PR

#if (_EVOL_UNE_SE_==0)
      end do
#endif
c
      end

c----------------------------------------------------------------------
      subroutine Ajuster_Amas_1pas
c----------------------------------------------------------------------
c     'Evolution passive' de 2 SE : seulement deplacement orbital
c
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
c
c var locales:
c ^^^^^^^^^^^^
      integer i
      double precision  dM
c
c function :
c ^^^^^^^^^^
      double precision M_TN_Croiss_Forc
c......................................................................
c
c---- Determiner la masse du TN (si croissance adiabatique)
c
#if (_FORCER_CROISSANCE_TN_>0)
      dM = -M_TN
      M_TN = M_TN_Croiss_Forc(Tps_amas)
      dM = dM + M_TN
      M_ctr = M_ctr + dM
#endif
      do i=1,2            
c
c---- Choisir la SE a deplacer
c
         call Choisir_SEact_Ajust ! determine iRang_SE1/2 iSE1/2 (parcours arbre)
c
c---- Tirage d'une nouvelle position orbitale
c
         call Ajuster_SEact

      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine Choisir_SEact_Ajust
c----------------------------------------------------------------------
c
      include 'EvolAmas_inc.f' 
      include 'EvolAmas_SEact_inc.f'
c
c fonction :
c ^^^^^^^^^^
      double precision random
c......................................................................
c
c---- On choisit une SE completement au hasard
c
 1    continue
      iRang_act = int((NbSE-NbSE_Perte)*random())+1
      iSE_act = iListeRay(iRang_act)
      if (R_SE(iSE_act).GE.R_perte) goto 1

      M_act = M_SE(iSE_act)
      R_act = R_SE(iSE_act)
      T_act = T_SE(iSE_act)
      J_act = J_SE(iSE_act)
      Tps_act = Tps_SE(iSE_act)

      iPG_act = -1

      A_act = 0.0d0
      B_act = 0.0d0      
      Pot_act = 0.0d0
      Peri_act = 0.0d0
      Apo_act = 0.0d0
c
      end
c
c----------------------------------------------------------------------
      subroutine Ajuster_SEact
c----------------------------------------------------------------------
c
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_LimExt_common.f'
c
c var locales :
c ^^^^^^^^^^^^^
      double precision Energie
c
c......................................................................
c
c---- Extraction de l'arbre binaire 
c
      call Extraire_SE_Arbre(iSE_act,iRang_act,A_act,B_act)
      Pot_act = B_act +
     $     ( A_act+(1.0d0-FracMasseInt)*
     $     M_act )/R_act
c
c---- Extraction de la grille
c
      call Retirer_de_PG(R_act,T_act,M_act,J_act,iPG_act)
c
c---- Extraction de l'amas
c
      B0 = B0 - M_act/R_act ! mise a jour de B(R=0)
      M_amas = M_amas - M_act
      Tps_amas = Tps_amas - Tps_act/(NbSE-NbSE_perte)
c
c---- Nouvelle position sur l'orbite (si liee)
c
      Energie = T_act+Pot_act
      call CalcApo_SE_ext2(M_act,Energie,J_act,Apo_act)
      if (Apo_Act.GT.0.0d0) then
         call CalcPeri_SE_ext2(M_act,Energie,J_act,Peri_act)
         if (Peri_act.LT.0.0d0) then
            write(0,*) '!!! Pas d''orbite physique pour la SE ',
     $           'active !!!'
            write(0,*) '    i_SE,M,E,J = ',iSE_act,M_act,Energie,J_act
            call Meurt('Ajuster_SEact')
         end if
         call Determiner_PosOrb_SEact
      end if
c
c---- Reinserer la SE
c
      call Reinserer_SEact
c
      end
c
#if (_EVOL_STELL_>0)
c----------------------------------------------------------------------
      subroutine Evolve_CentralObject
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_Etoiles_param.f'
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      logical lOK
      double precision Age_en_yr, Mini
c......................................................................
c
c---- For the time being, one only checks whether VMS object reaches
c     the end of MS. In which case the simulation is stopped!!
c
      if (iType_CtrObj.EQ.TypeEtoile_MS) then
         Mini = M_TN/Coef_ConvM
         Age_en_yr = max( 0.0d0, Tps_Amas/UneAnnee - BirthDate_CtrObj ) 
         call evolve_star(id_CtrObj,Age_en_yr)
                                ! %SeBa% change to:
                                ! iSE=iSE_PR(Num_SE)
                                ! call evolve_star(iSE,UneAnnee*TpsAmas)
         if (iType_CtrObj.NE.TypeEtoile_MS) then
            call open_append(iFich_Prov,'divers.asc',lOK)
            write(0,100)          iPas_Evol,Tps_Amas,Mini,
     $           BirthDate_CtrObj,Age_en_yr
            write(iFich_Prov,105) iPas_Evol,Tps_Amas,Mini,
     $           BirthDate_CtrObj,Age_en_yr
            close(iFich_Prov)
            call exit(0)
 100        format(
     $           '********************************************',/,
     $           '* Central object has left MS',/,
     $           '* iPas_Evol        = ',I10,/,
     $           '* Tps_Amas         = ',E12.5,/,
     $           '* M_CtrObj_Msun    = ',E12.5,/,
     $           '* BirthDate_CtrObj = ',E12.5,/,
     $           '* Age_CtrObj_yr    = ',E12.5,/,
     $           '********************************************'
     $           )
 105        format(
     $           '# Central object has left MS',/,
     $           ' iPas_Evol        = ',I10,/,
     $           ' Tps_Amas         = ',E12.5,/,
     $           ' M_CtrObj_Msun    = ',E12.5,/,
     $           ' BirthDate_CtrObj = ',E12.5,/,
     $           ' Age_CtrObj_yr    = ',E12.5
     $           )
         else
            if (abs(M_TN-Coef_ConvM*Mini).GT.1.0d10) then
               write(0,*) '!!! The mass of the central object ',
     $              'has changed through stellar evolution !!!'
               write(0,*) M_TN,Coef_ConvM*Mini
               call Meurt('Evolve_CentralObject')
            end if
         end if
      end if
            
      end
c
#endif /* if (_EVOL_STELL_>0) */
c----------------------------------------------------------------------
      subroutine Choisir_Paire_SE
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
c
c var. locales :
c ^^^^^^^^^^^^^^
      double precision DeltaR_sur_R
      logical lAccepte
c
c constantes :
c ^^^^^^^^^^^^
      double precision DeltaR_sur_R_max
      parameter (DeltaR_sur_R_max=0.25d0)
c
c fonction(s) :
c ^^^^^^^^^^^^^
      double precision random
      integer ihinv_TH
c......................................................................
c
c---- On tire une paire de SE successives avec une proba. ajustee
c     sur 1/(Pas de tps) selon methode de Henon
c
c     Modified on 11 June 2003 to allow the interacting pair's members to
c     be separated by DeltaRang_Paire-1 particles
c 
      lAccepte = .FALSE.
      do while (.NOT.lAccepte)
                                ! Rang de la 1ere SE
         iRang_PR(1) = ihinv_TH(random()) ! returns a value between 1 and NbSE-NbSE_Perte-DeltaRang_Paire
         iRang_PR(2) = iRang_PR(1)+DeltaRang_Paire
         if (iRang_PR(1).LT.1) then
            write(0,*) '!!! iRang_PR(1) too small !!!'
            call Meurt('Choisir_Paire_SE')
         end if
         if (iRang_PR(2).GT.NbSE-NbSE_Perte) then
            write(0,*) '!!! iRang_PR(2) too large !!!'
            call Meurt('Choisir_Paire_SE')
         end if
         
         call Trouver_SE_Rang
     $        (iRang_PR(1),iSE_PR(1),A_PR(1),B_PR(1)) ! parcours de l'arbre
         call Trouver_SE_Rang
     $        (iRang_PR(2),iSE_PR(2),A_PR(2),B_PR(2)) ! parcours de l'arbre
#if (_INTERDIT_RENC_DIST_>0)
         DeltaR_sur_R=(R_SE(iSE_PR(2))-R_SE(iSE_PR(1)))/R_SE(iSE_PR(1))
         if (DeltaR_sur_R.GT.DeltaR_sur_R_max) then
            iNb_RejetPaire = iNb_RejetPaire + 1
            write(0,*) '> Choisir_Paire_SE : Paire rejetee car DeltaR/R ',
     $           'trop grand :', DeltaR_sur_R
         else
            lAccepte = .TRUE.
         end if
#else
         lAccepte = .TRUE.
#endif
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine Calculer_Rencontre_Paire_SE
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_Rencontre_inc.f'
      include 'EvolAmas_Grille_param.f'
      include 'EvolAmas_Grille_common.f'
#if (_WRITE_CUMUL_COLL_PROB_>0)
      include 'EvolAmas_CumulProb_common.f'
#endif

c     
c constantes :
c ^^^^^^^^^^^^
      double precision LimVcarre
      parameter (LimVcarre=-1.0d-8)
      double precision dt_min
      parameter (dt_min=1.0d-30)
      integer iRang_min, p_rho
      parameter (iRang_min=100 000 000, p_rho=8) 
c
c var locales :
c ^^^^^^^^^^^^^
      integer Num_SE, iSE_rho, irang_rho_min, irang_rho_max, p_Henon
      double precision Vcarre,Phi,V2_tg,Pcoll,Pkick,dR,
     $   A_rho,B_dum,M_rho,DR_rho, x_rand, Rmin_rho, Rmax_rho, Ener
      logical lUnbound
c
c fonctions :
c ^^^^^^^^^^^
      double precision random,RelMR,PasTps_TH,
     $     DensSE_PG_lin,Rho_PG_lin2,Sigma1D_PG_lin2,
     $     Somme_MsurR_arbre
      integer iCouche_PG, iSE
      double precision stellar_mass, stellar_radius
      logical lOrbitOverlap_Paire_SE
c......................................................................
c
c---- proprietes physiques initiales des SE de la rencontre
c
      do Num_SE=1,2
         iSE = iSE_PR(Num_SE)
         M_PR(Num_SE)    = M_SE(iSE)
         T_PR(Num_SE)    = T_SE(iSE)
         J_PR(Num_SE)    = J_SE(iSE)
         R_PR(Num_SE)    = R_SE(iSE)
         Tps_PR(Num_SE)  = Tps_SE(iSE)
      end do
c
c---- initialisations
c
      iTypeRenc = RencIndef
      Theta_relax_renc = 0.0d0
      Theta_coll_renc  = 0.0d0
      Theta_kick_renc  = 0.0d0
c
c---- Check whether particles are bound 
c
      do Num_SE=1,2
         Ener = (A_PR(Num_SE)+(1.0d0-FracMasseInt)*M_PR(Num_SE))/
     $        R_PR(Num_SE)+B_PR(Num_SE)+T_PR(Num_SE)
         if (Ener.GE.0.0d0) then
#if (_VERB_EMPECHE_EVAP_RELAX_>0)
            write(0,*) '=== step ',iPas_Evol,' particle ',
     $           iSE_PR(Num_SE),' found unboud with E = ',Ener
#endif
            iTypeRenc = RencUnbound
         end if
      end do
c
c==== Determination des proprietes physique de la rencontre
c     -----------------------------------------------------

      if (iTypeRenc.EQ.RencIndef) then ! Skip all the rest if a particle is not bound
c
c---- Potentiel central de l'amas
c     
         Psi0 = Somme_MsurR_arbre()
c
c---- le pas de temps
c
         dt_renc = 0.5d0 *
     $        ( PasTps_TH(iRang_PR(1)) + PasTps_TH(iRang_PR(2)) )
                                ! tentons d'eviter les underflow
         if (dt_renc.LT.dt_min) then
            write(0,*)
     $           '> Calculer_Rencontre_Paire_SE: dt_renc trop petit (',
     $           dt_renc,'), pousse a ',dt_min
            dt_renc = dt_min
         end if
c
c---- les proprietes locales de l'amas
c
         R_renc =
     $        ( M_PR(1)*R_PR(1)+M_PR(2)*R_PR(2) ) /
     $        ( M_PR(1)+M_PR(2) )
         iPG_renc = -1
#if (_TYPE_DENS_RENC_==1)
         Dens_renc = Net_par_SE*DensSE_PG_lin(iPG_renc,R_renc) ! on en profite pour determiner iPG_renc         
         Met_moy_renc = M_PG(iPG_renc) / (Net_par_SE*NSE_PG(iPG_renc))
         Rho_renc = Met_moy_renc*Dens_renc
#elif (_TYPE_DENS_RENC_==2)   
         Rho_renc =  Rho_PG_lin2(iPG_renc,R_renc) ! on en profite pour determiner iPG_renc  
         Met_moy_renc = M_PG(iPG_renc) / (Net_par_SE*NSE_PG(iPG_renc)) 
         Dens_renc = Rho_renc/Met_moy_renc
#elif (_TYPE_DENS_RENC_==3) /* on traite specialement les SE les plus centrales: comme Henon */
         if (iRang_PR(1).LE.iRang_min) then
c
c     modified on Nov 6, 2003 to use same formula as Henon 71a, Eq 19+24 (mine was wrong)
c
            irang_rho_min = iRang_PR(1)-p_rho
            irang_rho_max = irang_rho_min+2*p_rho+1
            do while (irang_rho_min.LT.1)
               irang_rho_min = irang_rho_min + 1
               irang_rho_max = irang_rho_max - 1
            end do
            do while (irang_rho_max.GT.(NbSE-NbSE_perte))
               irang_rho_min = irang_rho_min + 1
               irang_rho_max = irang_rho_max - 1
            end do
            if (irang_rho_min.EQ.1) then
               irang_rho_max = irang_rho_min + 2
            else if (irang_rho_max.EQ.(NbSE-NbSE_perte)) then
               irang_rho_min = irang_rho_max - 2
            end if

            call Trouver_SE_Rang(irang_rho_max,iSE_rho,A_rho,B_dum)
            M_rho = -A_rho
            Rmax_rho = R_SE(iSE_rho)
            call Trouver_SE_Rang(irang_rho_min,iSE_rho,A_rho,B_dum)
            M_rho = M_rho + A_rho + M_SE(iSE_rho)
            Rmin_rho = R_SE(iSE_rho)
            p_Henon = irang_rho_max-irang_rho_min
            Dens_renc = (p_Henon-1)*Net_par_SE /
     $           (r4Pi*R_renc**2*(Rmax_rho-Rmin_rho))
c     Dens_renc = (p_Henon-1)*Net_par_SE /
c     $        (r4PiSur3*(Rmax_rho**3-Rmin_rho**3))
c     Dens_renc = (p_Henon-1)*p_henon*Net_par_SE / (p_Henon+1) / ! Ato's formula
c     $        (r4PiSur3*(Rmax_rho**3-Rmin_rho**3))
            Met_moy_renc = M_rho / ((p_Henon+1)*Net_par_SE)
            iPG_renc = iCouche_PG(R_renc)
         else
            Rho_renc =  Rho_PG_lin2(iPG_renc,R_renc) ! on en profite pour determiner iPG_renc 
            Met_moy_renc = M_PG(iPG_renc) / (Net_par_SE*NSE_PG(iPG_renc)) 
            Dens_renc = Rho_renc/Met_moy_renc         
         end if

#else
         write(0,*) '!!! valeur de _TYPE_DENS_RENC_ non valide !!!'
         call Meurt('Calculer_Rencontre_Paire_SE')
#endif
         Sigma1D_renc = Sigma1D_PG_lin2(iPG_renc,R_renc)
#if (_TEST_RENC_QUANTITIES_>0)
         write(33,'(4(e12.5," "))')
     $        R_renc, Dens_renc, Met_moy_renc, Sigma1D_renc
#endif
c
c---- Monitoring des rencontres non-locales
c
         dR=abs(R_PR(1)-R_PR(2))/min(R_PR(1),R_PR(2))
         dR_rel_Renc_moy=dR_rel_Renc_moy+dR
         dR_rel_Renc_max=max(dR_rel_Renc_max,dR)
         if (dR.GE.0.01d0) then
            iNb_Renc_dR_001 = iNb_Renc_dR_001 + 1
            if (dR.GE.0.1d0) then
               iNb_Renc_dR_01 = iNb_Renc_dR_01 + 1
               if (dR.GE.0.5d0) then
                  iNb_Renc_dR_05 = iNb_Renc_dR_05 + 1
                  if (dR.GE.1.0d0) then
                     iNb_Renc_dR_1 = iNb_Renc_dR_1 + 1
                  end if
               end if
            end if
         end if
c
c---- Vitesses des etoiles lors de la rencontre
c
c***      V1amas(1) = J_PR(1)/R_SE(iSE_PR(1))
         V1amas(1) = J_PR(1)/R_PR(1)
         V1amas(2) = 0.0d0
         Vcarre = 2.0d0*T_PR(1)-V1amas(1)**2
         if (Vcarre.GE.0.0d0) then
            V1amas(3)  = sqrt(Vcarre)
         else if (Vcarre.GE.LimVcarre*T_PR(1)) then
            V1amas(3)  = 0.0d0
         else
            write(0,*) '!!! V_rad(1)**2 = ',Vcarre,' !!!'
            call Meurt('Calculer_Rencontre_Paire_SE')
         end if
c***  V2_tg = J_PR(2)/R_SE(iSE_PR(2))
         V2_tg = J_PR(2)/R_PR(2)
         Vcarre = 2.0d0*T_PR(2)-V2_tg**2
         if (Vcarre.GE.0.0d0) then
            V2amas(3)  = sqrt(Vcarre)
         else if (Vcarre.GE.LimVcarre*T_PR(2)) then
            V2amas(3)  = 0.0d0
         else
            write(0,*) '!!! V_rad(2)**2 = ',Vcarre,' !!!'
            call Meurt('Calculer_Rencontre_Paire_SE')
         end if
                                ! on oriente les vitesses des deux * "au hasard" 
         Phi   = r2Pi*random()
         V2amas(1) = V2_tg*cos(Phi) ! l'une par-rapport a l'autre (tout en respectant 
         V2amas(2) = V2_tg*sin(Phi) ! les modules de Vtg et Vrad)
         if (random().GT.0.5d0) then
            V2amas(3) = -V2amas(3)
         end if
                                ! vitesse relative
         Vrel(1) = V1amas(1)-V2amas(1)
         Vrel(2) = V1amas(2)-V2amas(2)
         Vrel(3) = V1amas(3)-V2amas(3)
         Vrel_n = sqrt(Vrel(1)**2+Vrel(2)**2+Vrel(3)**2)
                                ! vitesse du CM
         Mu1 = M_PR(1)/(M_PR(1)+M_PR(2))
         Mu2 = 1.0d0-Mu1
         Vcm_renc(1) = Mu1*V1amas(1) + Mu2*V2amas(1)
         Vcm_renc(2) = Mu1*V1amas(2) + Mu2*V2amas(2)
         Vcm_renc(3) = Mu1*V1amas(3) + Mu2*V2amas(3)
                                ! param d'impact pour deviations sous grands angles
         b0_renc = Coef_ConvM*
     $        (stellar_mass(iSE_PR(1))+stellar_mass(iSE_PR(2))) /
     $        (Vrel_n*Vrel_n)
c
c---- orientation (aleatoire) du parametre d'impact
c
         Beta = r2pi*random()
c
c==== Faut-il simuler une collision ou de la relaxation ?
c     ---------------------------------------------------
c
#if (_COLLISIONS_>0 || _KICKS_>0 || _WRITE_CUMUL_COLL_PROB_>0)
         call Calc_ProbaColl_Paire_SE(Pcoll,Pkick)
#   if (_WRITE_CUMUL_COLL_PROB_>0)
         CumulPcoll_SE(iSE_PR(1)) = CumulPcoll_SE(iSE_PR(1)) + Pcoll
         CumulPcoll_SE(iSE_PR(2)) = CumulPcoll_SE(iSE_PR(2)) + Pcoll
         CumulPkick_SE(iSE_PR(1)) = CumulPkick_SE(iSE_PR(1)) + Pkick
         CumulPkick_SE(iSE_PR(2)) = CumulPkick_SE(iSE_PR(2)) + Pkick
#   endif
#if (_COLLISIONS_>0)
c
c---- Monitoring des grandes proba de collision
c
         P_coll_moy=P_coll_moy+Pcoll
         P_coll_max=max(P_coll_max,Pcoll)
         if (Pcoll.GE.0.25d0) then
            iNb_Coll_P_025 = iNb_Coll_P_025 + 1
            if (Pcoll.GE.0.5d0) then
               iNb_Coll_P_05 = iNb_Coll_P_05 + 1
               if (Pcoll.GE.1.0d0) then
                  iNb_Coll_P_1 = iNb_Coll_P_1 + 1
                  if (Pcoll.GE.2.0d0) then
                     iNb_Coll_P_2 = iNb_Coll_P_2 + 1
                  end if
               end if
            end if
         end if
#endif
#if (_KICKS_>0)
c
c---- Monitoring des grandes proba de deviation a grand angles
c
         P_kick_moy=P_kick_moy+Pkick
         P_kick_max=max(P_kick_max,Pkick)
         if (Pkick.GE.0.25d0) then
            iNb_Kick_P_025 = iNb_Kick_P_025 + 1
            if (Pkick.GE.0.5d0) then
               iNb_Kick_P_05 = iNb_Kick_P_05 + 1
               if (Pkick.GE.1.0d0) then
                  iNb_Kick_P_1 = iNb_Kick_P_1 + 1
                  if (Pkick.GE.2.0d0) then
                     iNb_Kick_P_2 = iNb_Kick_P_2 + 1
                  end if
               end if
            end if
         end if
#endif
c
c---- Y'a-t-il une collision ou une deviation a grand angle ?
c
         x_rand=random()
#if (_COLLISIONS_>0 )
         if (x_rand.LT.Pcoll) then
#   if (_COLLISION_CHECK_ORBIT_OVERLAP_==1) /* Do not account for stellar radii */
            if (lOrbitOverlap_Paire_SE(0.0d0) then
#   elif (_COLLISION_CHECK_ORBIT_OVERLAP_==2) /* Account for stellar radii */
            if ( lOrbitOverlap_Paire_SE( Coef_ConvR*
     $              (stellar_radius(iSE_PR(1))+stellar_radius(iSE_PR(2)))
     $              ) ) then
#   endif
                                ! COLLISION
                                ! ^^^^^^^^^
               iTypeRenc = Collision
               call Calculer_Collision_Paire_SE
#   if (_COLLISION_CHECK_ORBIT_OVERLAP_==1 || _COLLISION_CHECK_ORBIT_OVERLAP_==2)
            end if
#   endif
         end if
#endif
#if (_KICKS_>0 )
         if (iTypeRenc.EQ.RencIndef .AND. x_rand.LT.Pkick) then
                                ! KICK
                                ! ^^^^
            iTypeRenc = Kick
            call Calculer_Kick_Paire_SE
         end if
#endif
#endif /* if (_COLLISIONS_>0 || _KICKS_>0) */

c
c     On ne fait de la relaxation que si il n'y a pas eu de processus "plus violent"
c     Cela n'est pas tres rigoureux! On particulier en ce qui concerne les deviations
c     a "grands" angles qui peuvent ne pas etre si dramatiques si la limite est fixee
c     a grand nombre de b0...
c
         if (iTypeRenc.EQ.RencIndef) then
#if (_RELAXATION_>0)
                                ! RELAXATION
                                ! ^^^^^^^^^^
            iTypeRenc = Relaxation
            call Calculer_Relaxation_Paire_SE
#else  /* Pas de relaxation */
                                ! RENCONTRE PASSIVE
                                ! ^^^^^^^^^^^^^^^^^
            iTypeRenc = RencPassive
            call Avancer_Paire_SE
#endif /* if (_RELAXATION_>0) */
         end if

      end if                    ! All of this skipped if one of the particles
                                ! is unbound
c
      end
c
#if (_COLLISIONS_>0 || _KICKS_>0 || _WRITE_CUMUL_COLL_PROB_>0)
c----------------------------------------------------------------------
      subroutine Calc_ProbaColl_Paire_SE(Pcoll, Pkick)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_Rencontre_inc.f'
      include 'EvolAmas_Collision_inc.f'
      include 'EvolAmas_Etoiles_param.f'
c
c argument :
c ^^^^^^^^^^
      double precision Pcoll, Pkick
c
c constante :
c ^^^^^^^^^^^
      double precision Pcoll_max
      parameter (Pcoll_max=1.0d0) ! la probabilite de collision maximale
c
c var locales :
c ^^^^^^^^^^^^^
      double precision M1M2, R1R2, Mass1, Mass2, Ray1, Ray2, nu1, nu2
c
c functions :
c ^^^^^^^^^^^
      double precision stellar_mass, stellar_radius
      integer*1 stellar_type
c......................................................................
c
c---- calcul de la probabilite de collision pendant dt_renc
c     entre les 2 etoiles
c
#if (_COLLISIONS_>0 || _WRITE_CUMUL_COLL_PROB_>0)
                                ! Collision
      Mass1=stellar_mass(iSE_PR(1))
      Mass2=stellar_mass(iSE_PR(2))
      M1M2 = Coef_ConvM*(Mass1+Mass2)

      Ray1=stellar_radius(iSE_PR(1))
      Ray2=stellar_radius(iSE_PR(2))
                                ! Account for possibility of tidal disruption (added 29.08.03)
      if (stellar_type(iSE_PR(1)).NE.TypeEtoile_MS .OR.
     $     stellar_type(iSE_PR(2)).NE.TypeEtoile_MS) then
         call Solve_Roche(Mass1/Mass2,nu1)
         call Solve_Roche(Mass2/Mass1,nu2)
         R1R2 = Coef_ConvR* max(Ray1+Ray2,Ray1*nu1,Ray2*nu2)
      else
         R1R2 = Coef_ConvR*(Ray1+Ray2)
      endif

      Vesc_coll = sqrt(2.0d0*M1M2/R1R2)
      Seff_coll = rPi*(FactRmax_Coll*R1R2)**2 *
     $     ( 1.0d0 + (Vesc_coll/Vrel_n)**2 )
      Pcoll = Dens_renc*Vrel_n*dt_renc * Seff_coll
      if (Pcoll.GT.1.0d0) Seff_coll = Seff_coll/Pcoll !##### essai 18.07.01
#else
      Pcoll = 0.0d0
#endif

#if (_KICKS_>0 || _WRITE_CUMUL_COLL_PROB_>0)
                                ! Deviation a grand angle
      Seff_kick = rPi*(fact_b0_kick*b0_renc)**2
      Pkick = Dens_renc*Vrel_n*dt_renc * Seff_kick
      if (Pkick.GT.1.0d0) Seff_kick = Seff_kick/Pkick !##### essai 18.07.01
#else
      Pkick = 0.0d0
#endif

#if (_SURV_PCOLL_>0)
      if (Pcoll.GT.Pcoll_max) then
         write(0,*) '!!! Pcoll = ',Pcoll,
     $        ' dans Calc_ProbaColl_Paire_SE !!!'
      end if
#endif
#if (_SURV_PKICK_>0)
      if (Pkick.GT.Pcoll_max) then
         write(0,*) '!!! Pkick = ',Pkick,
     $        ' dans Calc_ProbaColl_Paire_SE !!!'
      end if
#endif
c
      end
c
#endif /* (_COLLISIONS_>0 || _KICKS_>0 || _WRITE_CUMUL_COLL_PROB_>0) */
c
c----------------------------------------------------------------------
      logical function lOrbitOverlap_Paire_SE(Radius_Margin)
c----------------------------------------------------------------------
c     Check for orbit overlap of particle pair, i.e. whether 
c     Rapo(part_1) + Radius_Margin > Rperi(part_2) or not 
c
      implicit none
      include 'EvolAmas_PaireRenc_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision Radius_Margin
c
c local variables :
c ^^^^^^^^^^^^^^^^^ 
      double precision Discr, Ra, Rp, Ener
      integer i_int, i_ext
c......................................................................
c
c---- Let's be paranoid and check which particle is at the interior
c
      if (R_PR(1).LT.R_PR(2)) then
         i_int=1
         i_ext=2
      else
         i_int=2
         i_ext=1
      end if
c
c---- Compute (Keplerian) peri-center distance for exterior particle
c     and check if it is smaller than position of interior partricle
c
      Ener = A_PR(i_ext)/R_PR(i_ext)+B_PR(i_ext)+T_PR(i_ext)
      Discr = A_PR(i_ext)**2 + 2.0d0*(Ener-B_PR(i_ext))*J_PR(i_ext)**2
#if (_SURVEILLE_SQRT_ > 0)
      if (Discr.LT.0.0d0) then
         write(0,*) '!!! Discr=',Discr,' for Peri-center computation !!!'
         call Meurt('lOrbitOverlap_Paire_SE')
      end if
#endif
      Rp = J_PR(i_ext)**2/(-A_PR(i_ext)+sqrt(Discr))
      if (Rp .LT. R_PR(i_int)+Radius_Margin) then
         lOrbitOverlap_Paire_SE = .TRUE.
         return
      end if
c
c---- Compute (Keplerian) apo-center distance for interior particle
c     and check if it is larger than peri-center distance for exterior particle
c
      Ener = A_PR(i_int)/R_PR(i_int)+B_PR(i_int)+T_PR(i_int)
      Discr = A_PR(i_int)**2 + 2.0d0*(Ener-B_PR(i_int))*J_PR(i_int)**2
#if (_SURVEILLE_SQRT_ > 0)
      if (Discr.LT.0.0d0) then
         write(0,*) '!!! Discr=',Discr,' for Apo-center computation !!!'
         call Meurt('lOrbitOverlap_Paire_SE')
      end if
#endif
      Ra = -(-A_PR(i_int)+sqrt(Discr))/(2.0d0*(Ener-B_PR(i_int)))
      if (Rp .LT. Ra+Radius_Margin) then
         lOrbitOverlap_Paire_SE = .TRUE.
         return
      end if
c
c---- No overlap
c
      lOrbitOverlap_Paire_SE = .FALSE.
c
      end
c

#if (_COLLISIONS_ > 0)
c----------------------------------------------------------------------
      subroutine Calculer_Collision_Paire_SE
c----------------------------------------------------------------------
c     ATTENTION : La SE de numero 1 doit etre celle de rang inferieur !!!
c
      include 'EvolAmas_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_Rencontre_inc.f'
      include 'EvolAmas_Collision_inc.f'
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_SuiviTypesStell_common.f'
c
c var locales :
c ^^^^^^^^^^^^^
      double precision 
     $     b_norm,              ! Parametre d'impact normalise par R1+R2
     $     Vrel_norm,           ! vitesse relative normalisee a Vesc
     $     V1cm(2),V2cm(2),
     $     Perte_M1,Perte_M2,Perte_M,dM_accr,dM_ejec,Perte_T,
     $     UnSurR1,UnSurR2,M1M1,M1M2,M2M2,Vcarre,Somme_MsurR,
     $     Mini1,Mini2,Mfin1,Mfin2,Vrel_norm_p,ct,st,V1cm_n,V2cm_n,
     $     f1,f2
      double precision Mmerger, Age1, Age2, AgeMerger, AgeOther,
     $     centralheliumfrac
      integer iMerger, iOther
      integer*1 iTet1, iTet2
#if (_VERB_COLL_>=5)
      double precision Mmin,Mmax,Met_min,Met_max
      save Met_min,Met_max
      data Met_min,Met_max /1.0d30,0.0d0/
#endif
c
c constantes :
c ^^^^^^^^^^^^
      double precision UnSurPi
      parameter (UnSurPi=1.0d0/rPi)
c
c fonctions :
c ^^^^^^^^^^^
      double precision random
      integer*1 stellar_type
      double precision stellar_mass, stellar_radius, stellar_age,
     $     stellar_centralheliumfrac
c......................................................................
c
c---- initialisations
c
      iTet1 = stellar_type(iSE_PR(1))
      iTet2 = stellar_type(iSE_PR(2))

      b_norm = sqrt(UnSurPi*Seff_coll*random()) /
     $     (Coef_ConvR*
     $     (stellar_radius(iSE_PR(1))+stellar_radius(iSE_PR(2))))
                                ! vitesse relative normalisee
      Vrel_norm = Vrel_n/Vesc_coll
      Perte_M1 = M_PR(1)        ! pour le calcul de la perte de masse et d'enrgie grav.
      Perte_M2 = M_PR(2)
      Perte_T  =                ! pour le calcul de la perte d'energie cinetique
     $     M_PR(1)*T_PR(1)+M_PR(2)*T_PR(2)
      M1M1=M_PR(1)*M_PR(1)
      M1M2=M_PR(1)*M_PR(2)
      M2M2=M_PR(2)*M_PR(2)
      Somme_MsurR=Psi0
      Psi0=Psi0-M_PR(1)/R_PR(1)-M_PR(2)/R_PR(2)
      Mini1=stellar_mass(iSE_PR(1))
      Mini2=stellar_mass(iSE_PR(2))
      Age1=stellar_age(iSE_PR(1))
      Age2=stellar_age(iSE_PR(2))
c
c==== Calcul du resultat de la collision
c     ----------------------------------
c
      call Collision_Stellaire(iTet1,iTet2,
     $     Mini1,Mini2,Vrel_norm,b_norm,
     $     Mfin1,Mfin2,Vrel_norm_p,
     $     Theta_coll_renc)
#if (_ECRIT_COLL_>0)
      write(iFich_Coll,120) iPas_Evol, R_renc,
     $     iSE_PR(1),iSE_PR(2),iTet1,iTet2,
     $     Mini1,Mini2,Vrel_norm,b_norm,
     $     Mfin1,Mfin2,Vrel_norm_p,
     $     Theta_coll_renc,
     $     A_PR(1)/R_PR(1)+B_PR(1),
     $     A_PR(2)/R_PR(2)+B_PR(2),
     $     Age1,Age2
 120  format(I10,' ',E13.6,2(' ',I9),2(' ',I2),10(' ',E13.6),
     $     2(' ',E11.4),$)
#endif
#if (_COLL_NULLES_>0) /* on annulle les modifications */
      Mfin1 = Mini1
      Mfin2 = Mini2
      Vrel_norm_p = Vrel_norm
      Theta_coll_renc = 0.0d0
#endif
#if (_COLL_NULLES_>1) /* on ne fait meme pas "progresser" les particules */
      dt_renc = 0.0d0
#endif

#if (_VERB_COLL_>0)
      Mmin=min(Mfin1,Mfin2)
      Mmax=max(Mfin1,Mfin2)
      if (Mmax.GT.Met_max) then
         Met_max = Mmax
         write(0,100) iPas_Evol,Met_max
 100     format('** Pas ',I10,' Etoile de masse Met_max = ',E12.4,
     $        ' creee par collision')
      end if
      if (Mmin.LT.Met_min .AND. Mmin.GT.0.0d0) then
         Met_min = Mmin
         write(0,110) iPas_Evol,Met_min
 110     format('** Pas ',I10,' Etoile de masse Met_min = ',E12.4,
     $        ' creee par collision')
      end if
#endif
      ct = cos(Theta_coll_renc)
      st = sin(Theta_coll_renc)
                                ! Nouvelles masses
      M_PR(1) = M_PR(1) * Mfin1/Mini1
      M_PR(2) = M_PR(2) * Mfin2/Mini2
                                ! Nouvelles vitesses dans le referentiel du CM
                                ! (on fait l'hypothese que ce ref reste le meme pendant la collision)
      Vrel_n=Vrel_norm_p*Vesc_coll
      if (Vrel_n.GT.0.0d0) then
         V1cm_n = Mu2*Vrel_n
         V2cm_n = Mu1*Vrel_n
         V1cm(1) = ct*V1cm_n
         V1cm(2) = st*V1cm_n
         V2cm(1) = -ct*V2cm_n
         V2cm(2) = -st*V2cm_n
      else
         V1cm(1) = 0.0d0
         V1cm(2) = 0.0d0
         V2cm(1) = 0.0d0
         V2cm(2) = 0.0d0
      end if
                           
      call Vitesses_de_CM_a_Amas ! Changement de referentiel CM --> Amas
     $     (Vrel,Vcm_renc,beta,V1cm,V2cm,V1amas,V2amas)
c
c---- modification des proprietes de la paire
c
      Perte_M1 = Perte_M1-M_PR(1)
      Perte_M2 = Perte_M2-M_PR(2)
      Perte_M = Perte_M1+Perte_M2 ! perte de masse (>0)

      if (M_ctr.GT.0.0d0) then   ! masse accretee par le TN
         dM_accr = frac_accr_Coll*Perte_M
      else
         dM_accr = 0.0d0
      end if
      dM_ejec = Perte_M - dM_accr
      
      UnSurR1 = 1.0d0/R_PR(1)
      UnSurR2 = 1.0d0/R_PR(2)

      Vcarre = V1amas(1)**2+V1amas(2)**2
                                ! Energie cinetique specifique
      T_PR(1) = 0.5d0*(Vcarre+V1amas(3)**2)
                                ! Moment cinetique specifique
      J_PR(1) = sqrt(Vcarre)*R_SE(iSE_PR(1))
                                ! Temps individuel
      Tps_PR(1) = Tps_PR(1) + dt_renc 
                                ! -Masse a l'interieur du rayon de la SE
      A_PR(1) = A_PR(1) + Perte_M1 - dM_accr
      B_PR(1) = B_PR(1) + Perte_M2*UnSurR2

      Vcarre = V2amas(1)**2+V2amas(2)**2
                                ! Energie cinetique specifique
      T_PR(2) = 0.5d0*(Vcarre+V2amas(3)**2)
                                ! Moment cinetique specifique
      J_PR(2) = sqrt(Vcarre)*R_SE(iSE_PR(2))
                                ! Temps individuel
      Tps_PR(2) = Tps_PR(2) + dt_renc 
                                ! -Masse a l'interieur du rayon de la SE
      A_PR(2) = A_PR(2) + dM_ejec
c
c---- Collisions change the effective stellar age of MS stars
c
      if (iTet1.EQ.TypeEtoile_MS .AND. iTet2.EQ.TypeEtoile_MS) then

#if (_COLL_REJUV_==0) /* No rejuvenation, keep same absolute age */

         if ( (Mfin1.EQ.0.0d0 .OR. Mfin2.EQ.0.0d0) .AND.
     $        (max(Mfin1,Mfin2).GT.max(Mini1,Mini2)) ) then
c
c     merger: one uses the age of the most massive pre-collision star
c
            if (Mfin1.EQ.0.0d0) then
               iMerger = iSE_PR(2)
               Mmerger = Mfin2
               iOther = iSE_PR(1)
            else
               iMerger = iSE_PR(1)
               Mmerger = Mfin1
               iOther = iSE_PR(2)
            end if
            if (Mini1.GT.Mini2) then
               AgeMerger = Age1
               AgeOther  = Age2
            else
               AgeMerger = Age2
               AgeOther  = Age1
            end if
            call stellar_SetByAge_ms(
     $           iMerger,Mmerger,StellMet_def,AgeMerger)
            call stellar_SetByAge_ms(
     $           iOther, 0.0d0,  StellMet_def,AgeOther )
         else
            call stellar_SetByAge_ms(iSE_PR(1),Mfin1,StellMet_def,Age1)
            call stellar_SetByAge_ms(iSE_PR(2),Mfin2,StellMet_def,Age2)
         end if

#elif (_COLL_REJUV_==1) /* Reset age to 0 in case of MS-MS merger */

         if ( (Mfin1.EQ.0.0d0 .OR. Mfin2.EQ.0.0d0) .AND.
     $        (max(Mfin1,Mfin2).GT.max(Mini1,Mini2)) ) then
            if (Mfin1.EQ.0.0d0) then ! Merger
               call stellar_SetByAge_ms(iSE_PR(1),0.0d0,StellMet_def,
     $              Age1)
               call stellar_SetByAge_ms(iSE_PR(2),Mfin2,StellMet_def,
     $              0.0d0)
            else
               call stellar_SetByAge_ms(iSE_PR(2),0.0d0,StellMet_def,
     $              Age2)
               call stellar_SetByAge_ms(iSE_PR(1),Mfin1,StellMet_def,
     $              0.0d0)
            end if
         else
            call stellar_SetByAge_ms(iSE_PR(1),Mfin1,StellMet_def,Age1)
            call stellar_SetByAge_ms(iSE_PR(2),Mfin2,StellMet_def,Age2)
         end if

#elif (_COLL_REJUV_==2) /* Use of core He mass to compute effective age assuming M_He rises linearly with time */
                        /* One assumes that only hydrogen is lost during collisions                            */
      
         if ( Mfin1.EQ.0.0d0 .OR. Mfin2.EQ.0.0d0 ) then ! Merger
            if (Mfin1.EQ.0.0d0) then
               Mmerger = Mfin2
               iMerger = iSE_PR(2)
               iOther = iSE_PR(1)
            else
               Mmerger = Mfin1
               iMerger = iSE_PR(1)
               iOther = iSE_PR(2)
            end if
            centralheliumfrac  = min(1.0d0, (
     $           stellar_centralheliumfrac(iSE_PR(1))*Mini1 +
     $           stellar_centralheliumfrac(iSE_PR(2))*Mini2 )/
     $           Mmerger )
            
            call stellar_SetByHeFrac_ms(iMerger,Mmerger,StellMet_def,
     $           centralheliumfrac)
            call stellar_Destroy(iOther)
         else

            call stellar_SetByHeFrac_ms( iSE_PR(1),Mfin1,StellMet_def,
     $           min(1.0d0,
     $           stellar_centralheliumfrac(iSE_PR(1))*Mini1/Mfin1) )
            call stellar_SetByHeFrac_ms( iSE_PR(2),Mfin2,StellMet_def,
     $           min(1.0d0,
     $           stellar_centralheliumfrac(iSE_PR(2))*Mini2/Mfin2) )
         end if
#else
         write(0,*) '!!! Invalid value of _COLL_REJUV_ !!!'
         call Meurt('Calculer_Collision_Paire_SE')
#endif
      else                      ! One star at least is (was) non-MS 
         if (Mfin1.LE.0.0d0) then
            call stellar_Destroy(iSE_PR(1))
         end if
         if (Mfin2.LE.0.0d0) then
            call stellar_Destroy(iSE_PR(2))
         end if
         if (                   ! One could use the rejuvenation prescription for MS stars
     $                          ! unless it has merged with a compact object
     $        (iTet1.EQ.TypeEtoile_MS .AND. Mfin1.GT.0.0d0 .AND.
     $        Mfin1.NE.Mini1) .OR.
     $        (iTet2.EQ.TypeEtoile_MS .AND. Mfin2.GT.0.0d0 .AND.
     $        Mfin2.NE.Mini2)) then
            write(0,900) iTet1,iTet2,Mini1,Mini2,Mfin1,Mfin2
 900        format('!!! Non-trivial non-MS/MS collision (Types= ',I3,I3,
     $           ', Mini= ',E12.5,':',E12.5,
     $           ', Mfin= ',E12.5,':',E12.5,') !!!',/,
     $           '!!! I don''t know how to handle stellar product !!!')
            call Meurt('Calculer_Collision_Paire_SE')
         end if
                                ! The mass of compact objects is allowed to change

         if (iTet1.NE.TypeEtoile_MS) then
            call stellar_SetMassType_compact(iSE_PR(1),Mfin1,iTet1)
         end if
         if (iTet2.NE.TypeEtoile_MS) then
            call stellar_SetMassType_compact(iSE_PR(2),Mfin2,iTet2)
         end if
      end if

#if (_ECRIT_COLL_>0)
      write(iFich_Coll,125)
     $     stellar_age(iSE_PR(1)), stellar_age(iSE_PR(2))
 125  format(2(' ',E11.4))
      call flush_MD(iFich_Coll)
#endif
c
c---- bilan de masse et d'energie (La masse totale de l'amas est mise
c     a jour lors de la re-introduction des SE)
c
      Perte_T = Perte_T -(M_PR(1)*T_PR(1)+M_PR(2)*T_PR(2)) ! perte d'energie cinetique orbitale (>0)
      dTtot_coll = dTtot_coll - Perte_T
      dEtot_coll = dEtot_coll - Perte_T
     $     -dM_accr*Somme_MsurR 
     $     +(Perte_M1*UnSurR1+Perte_M2*UnSurR2)*(-A_PR(1)-M_PR(1)) ! -A_PR(1)=M_ctr'+Somme des Mi pour i<=iRang_PR(1)
     $     -Perte_M*B_PR(2)+(M1M2-M_PR(1)*M_PR(2))*UnSurR2         ! -B_PR(2)=Somme des Mi/Ri pour i>iRang_PR(2)
     $     +FracMasseInt*(      ! terme d'autogravitation des SE
     $     UnSurR1*(M1M1-M_PR(1)**2)+UnSurR2*(M2M2-M_PR(2)**2) )

#if (_STATIQUE_==0)
      M_ctr = M_ctr + dM_accr
#endif
                                ! Contribution  a la masse accretee et ejectee
      dMaccr_coll = dMaccr_coll+dM_accr
      dMejec_coll = dMejec_coll+dM_ejec
      if (iTet1.EQ.iTet2) then
         dMaccr_coll_STS(iTet1) = dMaccr_coll_STS(iTet1)+dM_accr
         dMejec_coll_STS(iTet1) = dMejec_coll_STS(iTet1)+dM_ejec
      else if (Perte_M.GT.0.0d0) then
         f1 = M_SE(iSE_PR(1))/(M_SE(iSE_PR(1))+M_SE(iSE_PR(2)))
         f2 = 1.0d0-f1
                                ! quite arbitrary, but what else could I do ?
         dMaccr_coll_STS(iTet1) = dMaccr_coll_STS(iTet1)+f1*dM_accr
         dMejec_coll_STS(iTet1) = dMejec_coll_STS(iTet1)+f1*dM_ejec
         dMaccr_coll_STS(iTet2) = dMaccr_coll_STS(iTet2)+f2*dM_accr
         dMejec_coll_STS(iTet2) = dMejec_coll_STS(iTet2)+f2*dM_ejec
      end if
                                ! Modification du pot. central
      Psi0=Psi0+M_PR(1)/R_PR(1)+M_PR(2)/R_PR(2)
c
      iNb_Coll = iNb_Coll+1
c
      end
c
#endif /* (_COLLISIONS_>0) */
#if _KICKS_>0
c----------------------------------------------------------------------
      subroutine Calculer_Kick_Paire_SE
c----------------------------------------------------------------------
c     ATTENTION : La SE de numero 1 doit etre celle de rang inferieur !!!
c
      include 'EvolAmas_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_Rencontre_inc.f'
      include 'EvolAmas_Collision_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      double precision b,ct,st,Vcarre,V1cm(2),V2cm(2),V1cm_n,V2cm_n
c
c constantes :
c ^^^^^^^^^^^^
      double precision UnSurPi
      parameter (UnSurPi=1.0d0/rPi)
c fonctions :
c ^^^^^^^^^^^
      double precision random
c......................................................................
c
c---- On tire un parametre d'impact (a l'infini)
c
      b = sqrt(UnSurPi*Seff_kick*random())
c
c---- angle de deviation :
c
      Theta_kick_renc = 2.0d0*atan(b0_renc/b)
      ct = cos(Theta_kick_renc)
      st = sin(Theta_kick_renc)
c     
c---- calcul de la rencontre hyperbolique avec angle de
c     deviation Theta_kick_renc
                                ! vitesses des deux particules apres la rencontre
                                ! dans le referentiel du CM
      V1cm_n = Mu2*Vrel_n
      V2cm_n = Mu1*Vrel_n
      V1cm(1) = ct*V1cm_n
      V1cm(2) = st*V1cm_n
      V2cm(1) = -ct*V2cm_n
      V2cm(2) = -st*V2cm_n
                                ! Changement de referentiel CM --> Amas
      call Vitesses_de_CM_a_Amas
     $     (Vrel,Vcm_renc,beta,V1cm,V2cm,V1amas,V2amas)
c
c---- modification des proprietes de la paire
c
      Vcarre = V1amas(1)**2+V1amas(2)**2
                                ! Energie cinetique specifique
      T_PR(1) = 0.5d0*(Vcarre+V1amas(3)**2)
                                ! Moment cinetique specifique
      J_PR(1) = sqrt(Vcarre)*R_SE(ISE_PR(1))

      Vcarre = V2amas(1)**2+V2amas(2)**2
                                ! Energie cinetique specifique
      T_PR(2) = 0.5d0*(Vcarre+V2amas(3)**2)
                                ! Moment cinetique specifique
      J_PR(2) = sqrt(Vcarre)*R_SE(iSE_PR(2))
                                ! Temps individuels
      Tps_PR(1) = Tps_PR(1) + dt_renc 
      Tps_PR(2) = Tps_PR(2) + dt_renc
c
      iNb_kick = iNb_kick+1
c
      end
c
#endif /* (_KICKS_>0) */
c
#if (_RELAXATION_>0)
c----------------------------------------------------------------------
      subroutine Calculer_Relaxation_Paire_SE
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_Rencontre_inc.f'
      include 'EvolAmas_Relax_param.f'
      include 'EvolAmas_Relax_common.f'
      include 'EvolAmas_StrEnc_param.f'
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_Etoiles_common.f'
      include 'EvolAmas_NonInterac_common.f'
c     
c constantes :
c ^^^^^^^^^^^^
      double precision LimVcarre
      parameter (LimVcarre=-1.0d-8)
      double precision Pi2,Pi4,Pi8,Pi16
      parameter (Pi2=0.5d0*rPi,Pi4=0.25d0*rPi,Pi8=0.125d0*rPi,
     $     Pi16=0.0625d0*rPi)
      double precision Theta_lim, Theta_max
      parameter (Theta_lim=Pi2, Theta_Max=Pi2)
      double precision RapRmax
      parameter (RapRmax=2.0d0)
c
c var locales :
c ^^^^^^^^^^^^^
      double precision trel_renc,ct,st,Vcarre,
     $     V1cm(2),V2cm(2),V1cm_n,V2cm_n,v1_avant(3),v2_avant(3),
     $     E1,E2, ParamImp, Tau, Theta_raw, dE,
     $     dummy1,dummy2
      integer i,NumSE, i_Strong_Encounter
      integer*1 iTypeEt1, iTypeEt2
      logical lAnnule, lPrem
      save lPrem
      data lPrem /.TRUE./
c
c fonctions :
c ^^^^^^^^^^^
      double precision Coef_Relax, RandomGauss
      double precision stellar_mass
c......................................................................
c
                                ! temps de relaxation local 
      trel_renc = 1.0d0 /
     $     ( 0.5d0*Coef_Relax(iRang_PR(1))*b0_renc**2*Vrel_n*Dens_renc )

#if (_INTERDIT_RELAX_SE_DIST_>0)
      if (R_PR(2).GT.RapRmax*R_PR(1) .OR. R_PR(1).GT.RapRmax*R_PR(2))
     $     then
         dt_renc = 0.0d0        ! carrement !!!
         goto 99                ! jump to end of routine
      end if
#endif
                                ! At least, update times of particles!
      Tps_PR(1) = Tps_PR(1) + dt_renc 
      Tps_PR(2) = Tps_PR(2) + dt_renc
c
c---- Check for (symmetrical) non-interaction
c
      iTypeEt1 = iTet_SE(iSE_PR(1))
      iTypeEt2 = iTet_SE(iSE_PR(2))
      if ( lNonInterac .AND. (.NOT.lRelaxesWith(iTypeEt1,iTypeEt2))
     $     .AND. (.NOT.lNonInteracAsymmetrical .OR.
     $           (iTypeEt1.EQ.iTypeEt2)             )   ) then
         goto 99                ! jump to end of routine
      end if

                                ! angle de diffusion gravitationnelle
#if (_THETA_GIERSZ_>0)
                                ! formule de Giersz
      Theta_relax_renc = Pi4*sqrt(dt_renc/trel_renc)
      Theta_raw = Theta_relax_renc
      if (Theta_relax_renc.GE.1.0d0) then
         Theta_relax_renc = Pi2
      else
         Theta_relax_renc=2.0d0*asin(Theta_relax_renc)
      end if

#else /* (_THETA_GIERSZ_>0) */
                                ! formule de Henon
      Theta_relax_renc = Pi2*sqrt(dt_renc/trel_renc)
      Theta_raw = Theta_relax_renc
c
c---- Monitoring des relaxations sous grand angles
c
#if (_MONITOR_STRONG_ENCOUNTERS_>0)
      i_Strong_Encounter = 0
      if (Theta_relax_renc.GT.Theta_Strong_Encounter)
     $     i_Strong_Encounter = i_Strong_Encounter + 10
#endif
      if (Theta_relax_renc.GT.Pi16) then
         iNb_SuperRenc_Pi16 = iNb_SuperRenc_Pi16+1
         if (Theta_relax_renc.GT.Pi8) then
            iNb_SuperRenc_Pi8 = iNb_SuperRenc_Pi8+1
            if (Theta_relax_renc.GT.Pi4) then
               iNb_SuperRenc_Pi4 = iNb_SuperRenc_Pi4+1
               if (Theta_relax_renc.GT.Pi2) then
                  iNb_SuperRenc_Pi2 = iNb_SuperRenc_Pi2+1
               end if
            end if
         end if
      end if
      Dev_SuperRenc_moy  = Dev_SuperRenc_moy + Theta_relax_renc
      Dev2_SuperRenc_moy = Dev2_SuperRenc_moy + Theta_relax_renc**2

#if (_LIMITER_THETA_>0)
      if (Theta_relax_renc.GT.Theta_lim) then
         Theta_relax_renc = Theta_max
      end if
#endif /* (_LIMITER_THETA_>0) */

#endif /* (_THETA_GIERSZ_>0) */

#if (_ECRIT_STAT_RELAX_>0)
      do i=1,3
         v1_avant(i)=V1amas(i)
         v2_avant(i)=V2amas(i)
      end do
#endif

      ct = cos(Theta_relax_renc)
      st = sin(Theta_relax_renc)
c     
c---- calcul de la rencontre hyperbolique avec angle de
c     deviation Theta_relax_renc
                                ! vitesses des deux particules apres la rencontre
                                ! dans le referentiel du CM
      V1cm_n = Mu2*Vrel_n
      V2cm_n = Mu1*Vrel_n
      V1cm(1) = ct*V1cm_n
      V1cm(2) = st*V1cm_n
      V2cm(1) = -ct*V2cm_n
      V2cm(2) = -st*V2cm_n
                                ! Changement de referentiel CM --> Amas
      call Vitesses_de_CM_a_Amas
     $     (Vrel,Vcm_renc,beta,V1cm,V2cm,V1amas,V2amas)
c
c---- modification des proprietes de la paire
c
      Vcarre = V1amas(1)**2+V1amas(2)**2
                                ! Energie cinetique specifique
      T_PR(1) = 0.5d0*(Vcarre+V1amas(3)**2) !# TWILIGHT ZONE HERE !!!
                                ! Moment cinetique specifique
      J_PR(1) = sqrt(Vcarre)*R_SE(ISE_PR(1))

      Vcarre = V2amas(1)**2+V2amas(2)**2
                                ! Energie cinetique specifique
      T_PR(2) = 0.5d0*(Vcarre+V2amas(3)**2) !# TWILIGHT ZONE HERE !!!
                                ! Moment cinetique specifique
      J_PR(2) = sqrt(Vcarre)*R_SE(iSE_PR(2))

#if (_MONITOR_STRONG_ENCOUNTERS_>0)
      dE = abs(T_PR(1)-T_SE(iSE_PR(1)))
      if (dE.GT.dErel_Strong_Encounter*T_SE(iSE_PR(1)))
     $     i_Strong_Encounter = i_Strong_Encounter + 1
      dE = abs(T_PR(2)-T_SE(iSE_PR(2)))
      if (dE.GT.dErel_Strong_Encounter*T_SE(iSE_PR(2)))
     $     i_Strong_Encounter = i_Strong_Encounter + 2
#endif

#if (_EMPECHE_EVAP_RELAX_>0 && _TRONC_MAREE_==0)
      lAnnule=.FALSE.
      do NumSE=1,2
         lAnnule = lAnnule .OR.
     $        ( (T_PR(NumSE) + B_PR(NumSE) +
     $        (A_PR(NumSE)+(1.0d0-FracMasseInt)*M_PR(NumSE) )/
     $        R_PR(NumSE)) .GT. 0.0d0 )
      end do      
c XXXX There seems to be a bug with evapoartion of stars even when _EMPECHE_EVAP_RELAX_>0
c XXXX Is the computation of the binding energy OK? How is it affected by the other star
c XXXX and particularly its mass change (in case of accretion)

      if (lAnnule) then
#if (_VERB_EMPECHE_EVAP_RELAX_>0)
         write(0,*)'> Relaxation step cancelled to prevent evaporation'
#endif
         do NumSE=1,2
            T_PR(NumSE) = T_SE(iSE_PR(NumSE))
            J_PR(NumSE) = J_SE(iSE_PR(NumSE))
            Tps_PR(NumSE) = Tps_SE(iSE_PR(NumSE))
         end do
         iNb_AnnulRelaxEvap = iNb_AnnulRelaxEvap + 1
      end if
#endif

#if (_CONTROLE_CENTRE_==3) /* on empeche une se de devenir autogravitante */
      lAnnule=.FALSE.
      do NumSE=1,2
         lAnnule = lAnnule .OR.
     $        ( (T_PR(NumSE)*R_PR(NumSE)-FracMasseInt*M_PR(NumSE)) .LT.
     $        0.0d0 )
      end do
      if (lAnnule) then
         write(0,*)'> Relaxation annullee pour empecher auto-gravitation'
         do NumSE=1,2
            T_PR(NumSE) = T_SE(iSE_PR(NumSE))
            J_PR(NumSE) = J_SE(iSE_PR(NumSE))
            Tps_PR(NumSE) = Tps_SE(iSE_PR(NumSE))
         end do
         iNb_AnnulRelaxCC = iNb_AnnulRelaxCC + 1
      end if
#endif        

#if (_ECRIT_STAT_RELAX_>0)
      if (lPrem) then
         open(80,file='stat_relax.asc')
         write(80,100)
 100     format('# Variations de vitesse par relax.',/,
     $        '# 1: Vavant 2: dV ')
         lPrem = .FALSE.
      end if
      write(80,110) sqrt(v1_avant(1)**2+v1_avant(2)**2+v1_avant(3)**2),
     $     sqrt((v1_avant(1)-V1amas(1))**2+(v1_avant(2)-V1amas(2))**2+
     $     (v1_avant(3)-V1amas(3))**2)
      write(80,110) sqrt(v2_avant(1)**2+v2_avant(2)**2+v2_avant(3)**2),
     $     sqrt((v2_avant(1)-V2amas(1))**2+(v2_avant(2)-V2amas(2))**2+
     $     (v2_avant(3)-V2amas(3))**2)
 110  format(2(E11.4,'
#endif

#if (_ECRIT_RENC_RELAX_>0)
      if (lPrem) then
         open(81,file='super_renc_relax.asc')
         write(81,100)
 100     format('# Super-rencontres de relaxation',/,
     $        '# 1: iPas_Evol 2: R_renc 3: Met_1 4: Met_2 5: Vrel_renc',
     $        ' 6: trel_renc 7: dt_renc')
         lPrem = .FALSE.
      end if
      write(81,800)iPas_Evol,R_renc, stellar_mass(iSE_PR(1)),
     $     stellar_mass(iSE_PR(2)),Vrel_n,trel_renc,dt_renc
 800  format(I10,6(' ',E8.3))
#endif

#if (_MONITOR_STRONG_ENCOUNTERS_>0)
      if (i_Strong_Encounter.GT.0) then
         write(iFich_StrEnc,200) iPas_Evol,i_Strong_Encounter,
     $        iSE_PR, R_PR, stellar_mass(iSE_PR(1)),
     $        stellar_mass(iSE_PR(2)),
     $        Vrel_n,trel_renc,dt_renc, Theta_raw, 
     $        T_SE(iSE_PR(1)), T_SE(iSE_PR(2)), 
     $        J_SE(iSE_PR(1)), J_SE(iSE_PR(2)),
     $        T_PR, J_PR
 200     format(I10,1x,I2,2(1x,I8),16(1x,E11.4))
         call flush_MD(iFich_StrEnc)
      end if
#endif
c
c---- Is this combination of types actually interacting? (modification Aug 16 2005)
c
c     Note that energy conservation _and_ monitoring is screwed up
c     at this point (I think)!
c
      if (lNonInterac) then
         if (.NOT.lRelaxesWith(iTypeEt1,iTypeEt2)) then
                                ! Revert to pre-encounter quantities
            T_PR(1) = T_SE(iSE_PR(1))
            J_PR(1) = J_SE(iSE_PR(1))
         end if
         if (.NOT.lRelaxesWith(iTypeEt2,iTypeEt1)) then
            T_PR(2) = T_SE(iSE_PR(2))
            J_PR(2) = J_SE(iSE_PR(2))
         end if
      end if

#if (_RELAX_NULLE_>0) /* on annulle la modification */
      do NumSE=1,2
         T_PR(NumSE) = T_SE(iSE_PR(NumSE))
         J_PR(NumSE) = J_SE(iSE_PR(NumSE))
#if (_RELAX_NULLE_>1) /* on ne fait meme pas "progresser" les particules */
         Tps_PR(NumSE) = Tps_SE(iSE_PR(NumSE))
#endif
      end do
#endif
c
 99   continue

      iNb_Relax = iNb_Relax + 1 ! Compte tous les pas de relax 
      iNb_SuperRenc = iNb_SuperRenc + 1 ! Est remis a 0 apres chaque ecriture dans le "log"
c
      end
c
#endif /* (_RELAXATION>0) */
c
c----------------------------------------------------------------------
      subroutine Init_FichStrEnc_ASC(iFich)
c----------------------------------------------------------------------
      implicit none
      include 'CstesMath.f'
      include 'EvolAmas_StrEnc_param.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
c......................................................................
c
      write(iFich,100) Theta_Strong_Encounter, dErel_Strong_Encounter
 100  format(
     $     '# Monitoring of particularly strong relaxational ',
     $     'super-encounters',/,
     $     '# Values of iStrEnc : ',/,
     $     '# 10: theta > Theta_Strong_Encounter',/,
     $     '# +1,2: dE/E >  dErel_Strong_Encounter for SE 1,2',/,
     $     '# Theta_Strong_Encounter = ',E12.5,/,
     $     '# dErel_Strong_Encounter = ',E12.5
     $     )
      call OutputFileHeader(iFich,
     $     'iPas_Evol iStrEnc iSE1 iSE2 ' //
     $     'R1 R2 Met1 Met2 ' //
     $     'Vrel_renc Trel_renc dt_renc Theta_renc ' //
     $     'Tbefore1 Tbefore2 Jbefore1 Jbefore2 ' //
     $     'Tafter1 Tafter2 Jafter1 Jafter2'
     $     )
c
      end
c     
c----------------------------------------------------------------------
      subroutine Vitesses_de_CM_a_Amas
     $     (Vrel,Vcm,beta,V1cm,V2cm,V1amas,V2amas)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $                           ! INPUT :
     $     Vrel(3),              ! vitesse relative (initiale) des 2 etoiles
     $     Vcm(3),               ! vitesse du CM de la paire (dans ref de l'amas)
     $     beta,                 ! angle d'orientation du parametre d'impact autour de Vrel
     $     V1cm(2),V2cm(2),      ! vitesses des 2 etoiles dans ref CM (dans plan orbital)
     $                           ! OUTPUT :
     $     V1amas(3),V2amas(3) ! vitesses des 2 etoiles dans ref de l'amas
c
c var locales :
c ^^^^^^^^^^^^^
      double precision Vrel_tg,Vrel_n,UnSurVrel_n,UnSurVrel_tg,
     $     cb,sb,cd,sd,cg,sg,
     $     e_y(3),e_x(3)        ! e_x est un vecteur unite // Vrel
                                ! e_y est un vecteur unite // param d'impact 
c......................................................................
c
      cb = cos(beta)
      sb = sin(beta)
      Vrel_tg = Vrel(1)**2+Vrel(2)**2
      Vrel_n  = sqrt(Vrel_tg+Vrel(3)**2)
      UnSurVrel_n = 1.0d0/Vrel_n
      Vrel_tg = sqrt(Vrel_tg)
      if (Vrel_tg.GT.1.0d-7*Vrel_n) then
         UnSurVrel_tg = 1.0d0/Vrel_tg
         cd = Vrel(1)*UnSurVrel_tg
         sd = Vrel(2)*UnSurVrel_tg
         cg = Vrel(3)*UnSurVrel_n
         sg = Vrel_tg*UnSurVrel_n
         e_x(1) = Vrel(1)*UnSurVrel_n
         e_x(2) = Vrel(2)*UnSurVrel_n
         e_x(3) = Vrel(3)*UnSurVrel_n
         e_y(1) = cd*cg*cb-sd*sb
         e_y(2) = sd*cg*cb+cd*sb
         e_y(3) = -sg*cb
      else
         e_x(1)=0.0d0
         e_x(2)=0.0d0
         e_x(3)=1.0d0
         e_y(1)=cb
         e_y(2)=sb
         e_y(3)=0.0d0
      end if
      V1amas(1) = V1cm(1)*e_x(1) + V1cm(2)*e_y(1) + Vcm(1)
      V1amas(2) = V1cm(1)*e_x(2) + V1cm(2)*e_y(2) + Vcm(2)
      V1amas(3) = V1cm(1)*e_x(3) + V1cm(2)*e_y(3) + Vcm(3)
      V2amas(1) = V2cm(1)*e_x(1) + V2cm(2)*e_y(1) + Vcm(1)
      V2amas(2) = V2cm(1)*e_x(2) + V2cm(2)*e_y(2) + Vcm(2)
      V2amas(3) = V2cm(1)*e_x(3) + V2cm(2)*e_y(3) + Vcm(3)
c
      end
c
#if (_EVOL_STELL_>0 || _ACCR_BY_STARS_>0)
c----------------------------------------------------------------------
      subroutine Calculer_EvolStell_Paire_SE
c----------------------------------------------------------------------
c     This routine implements change of mass due to stellar evolution 
c     or accretion of gas by the stars (both are exclusive of each
c     other).
c
c     Velocity kicks can be given when a NS or BH is born.
c
      include 'EvolAmas_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_Rencontre_inc.f'
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_SuiviTypesStell_common.f'
c
c constantes :
c ^^^^^^^^^^^^
      double precision Vkick_forced_ejection
      parameter (Vkick_forced_ejection = 1.0d10)
c
c var. locales :
c ^^^^^^^^^^^^^^
      double precision Age_en_yr, Perte_M1, Perte_M2, Perte_M, Perte_T,
     $     M1M1, M1M2, M2M2, Somme_MsurR, dM_accr, 
     $     UnSurR1, UnSurR2, Met_ini, Met_increased, Vescape
      integer Num_SE, iSE
      integer*1 iTet_ini, iTet1, iTet2
      double precision Vkick, ctheta, stheta, T_ini, J_ini, Vtg, Vrad,
     $     Vkick_tg, Vkick_rad 
      double precision tini_Myr,dt_Myr,frac_accr_by_MBH
c
c functions :
c ^^^^^^^^^^^
      integer*1 stellar_type
      double precision stellar_mass, stellar_dateofbirth,
     $     natal_kick_star
      double precision random
#if (_ACCR_BY_STARS_>0)
      double precision Maccr_star_Msun
#endif
c
c......................................................................
c
c---- Initialisations
c
      Perte_M1 = M_PR(1)        ! pour le calcul de la perte de masse et d'energie grav.
      Perte_M2 = M_PR(2)
      Perte_T  =                ! pour le calcul de la perte d'energie cinetique
     $     M_PR(1)*T_PR(1)+M_PR(2)*T_PR(2)
      M1M1=M_PR(1)*M_PR(1)
      M1M2=M_PR(1)*M_PR(2)
      M2M2=M_PR(2)*M_PR(2)
      Somme_MsurR = Psi0
      Psi0=Psi0-M_PR(1)/R_PR(1)-M_PR(2)/R_PR(2)
      iTet1 = stellar_type(iSE_PR(1))
      iTet2 = stellar_type(iSE_PR(2))
c
c---- Evolution (eventuelle) des deux etoiles
c
      do Num_SE=1,2
         if (M_PR(Num_SE).GT.0.0d0) then
            iSE = iSE_PR(Num_SE)
            iTet_ini = stellar_type(iSE)
            Met_ini  = stellar_mass(iSE)
            T_ini = T_PR(Num_SE)
            J_ini = J_PR(Num_SE)
#   if (_EVOL_STELL_>0)
#      if (_ACCR_BY_STARS_>0) /* Do not even compile !! */
            XXX ERROR: CANNOT HAVE _EVOL_STELL_>0 AND _ACCR_BY_STARS_>0 XXX
#      endif
            Age_en_yr = max( 0.0d0, Tps_Amas/UneAnnee -
     $           stellar_dateofbirth(iSE) )
c            Age_en_yr = max( 0.0d0, Tps_PR(Num_SE)/UneAnnee -
c     $           stellar_dateofbirth(iSE) )
            call evolve_star(iSE,Age_en_yr)
                                ! %SeBa% change to:
                                ! iSE=iSE_PR(Num_SE)
                                ! call evolve_star(iSE,UneAnnee*TpsAmas)
            
                                ! Modifier la masse/ le type de l'etoile
                                ! Accreter (un partie de) la masse perdue
#   endif /* (_EVOL_STELL_>0) */
#   if (_ACCR_BY_STARS_>0)
c
c     The star grows in mass by gas accretion (the gas is assumed to
c     come from a large distance and doesn't contribute to the potentiel
c     or creates a drag until it is accreted)
c
            tini_Myr = 1.0d-6*Tps_SE(iSE)/uneannee
            dt_Myr = 1.0d-6*dt_renc/uneannee
            Met_increased = Met_ini +
     $           Maccr_star_Msun(iSE,tini_Myr,dt_Myr)
            call stellar_SetMass_ms(iSE,Met_increased)
c
c     We assume the gas has zero linear momentum so
c     the star must slow down to conserve total momentum
c
            T_PR(Num_SE) = T_ini*(Met_ini/Met_increased)**2
            J_PR(Num_SE) = J_ini*(Met_ini/Met_increased)

#   endif /* (_ACCR_BY_STARS_>0) */

            M_PR(Num_SE) = M_PR(Num_SE) * stellar_mass(iSE)/Met_ini
c
c--- Give natal kick if a NS or a BH is formed
c
#if (_EVOL_STELL_>0 && _NATKICKS_>0)
            Vkick = natal_kick_star() ! Vkick in km/s
            if (Vkick.GT.0.0d0) then

               Vescape = sqrt(-A_PR(Num_SE)/R_PR(Num_SE)-B_PR(Num_SE))
               if (Vkick.GT.Vkick_forced_ejection) then
                  Vkick = 2.0d0*Vescape ! to make sure the star will get ejected but avoid problems with energy conservation
               else
                                ! convert to code units
                  Vkick = UnKmS * Vkick
               end if

               ! random orientation of the kick

               ctheta = 1.0d0-2.0d0*random()
               stheta = sqrt(1-ctheta**2)
               Vkick_rad = Vkick*ctheta
               Vkick_tg  = Vkick*stheta

               Vtg = J_PR(Num_SE)/R_PR(Num_SE)
               Vrad = sqrt(2.0d0*T_PR(Num_SE) - Vtg**2)
               Vtg = sqrt( Vtg**2 + Vkick_tg**2 +
     $              cos(r2Pi*random())*Vkick_tg*Vtg )
               if (random().GT.0.5d0) then
                  Vrad = Vrad + Vkick_rad
               else
                  Vrad = Vrad - Vkick_rad
               end if

               ! New values of T and J for the particle
               
               T_PR(Num_SE) = 0.5d0*(Vrad**2+Vtg**2)
               J_PR(Num_SE) = Vtg*R_PR(Num_SE)

               ! Keep track of change in energy

               dTtot_nk = dTtot_nk + M_PR(Num_SE)*(T_PR(Num_SE)-T_ini)
               dEtot_nk = dEtot_nk + M_PR(Num_SE)*(T_PR(Num_SE)-T_ini)

#   if (_ECRIT_NATKICKS_>0)
               write(iFich_NatKcks,110) iPas_Evol,iSE,
     $              Vkick/UnKmS, Vkick, Vkick_rad, Vkick_tg,
     $              T_ini, J_ini, T_PR(Num_SE), J_PR(Num_SE), Vescape
 110           format(I10,' ',I8,9(' ',E12.5))
               call flush_MD(iFich_NatKcks)
#   endif

            end if
#endif /* (_EVOL_STELL_>0 && _NATKICKS_>0) */

#if (_EVOL_STELL_>0 && _ECRIT_STEVOL_>0)
            if (stellar_type(iSE).NE.iTet_ini) then
               write(iFich_StEvol,100)
     $              iPas_Evol,iSE,Age_en_yr,R_PR(Num_SE),
     $              iTet_ini,          Met_ini,
     $              2.0d0*T_ini,
     $              stellar_type(iSE), stellar_mass(iSE),
     $              2.0d0*T_PR(Num_SE)
 100           format(I10,' ',I8,' ',E12.5,' ',E12.5,' ',
     $              2(I2,' ',E12.5,' ',E12.5,' '))
               call flush_MD(iFich_StEvol)
            end if
#endif


         end if
      end do
c
c---- bilan de masse et d'energie (La masse totale de l'amas est mise
c     a jour lors de la re-introduction des SE)
c
      Perte_M1 = Perte_M1-M_PR(1)
      Perte_M2 = Perte_M2-M_PR(2)
      Perte_M = Perte_M1+Perte_M2 ! mass "loss" (>0 for stellar evolution; <0 for accretion by stars)

      if (M_ctr.GT.0.0d0 .AND. Perte_M.GT.0.0d0) then ! fraction of mass accreted by central object
         frac_accr_by_MBH = frac_accr_EvolSt
      else
         frac_accr_by_MBH = 0.0d0
      end if
      dM_accr = frac_accr_by_MBH*Perte_M
 
      UnSurR1 = 1.0d0/R_PR(1)
      UnSurR2 = 1.0d0/R_PR(2)
                                ! -Masse a l'interieur du rayon de la SE
      A_PR(1) = A_PR(1) + Perte_M1 - dM_accr
      B_PR(1) = B_PR(1) + Perte_M2*UnSurR2
                                ! -Masse a l'interieur du rayon de la SE
      A_PR(2) = A_PR(2) + Perte_M - dM_accr

                                ! Perte d'energie cinetique orbitale (>0)
                                ! With natal kick(s), this may actually be positive.
      Perte_T = Perte_T -(M_PR(1)*T_PR(1)+M_PR(2)*T_PR(2))
      dTtot_evst = dTtot_evst - Perte_T
      dEtot_evst = dEtot_evst - Perte_T
     $     -dM_accr*Somme_MsurR 
     $     +(Perte_M1*UnSurR1+Perte_M2*UnSurR2)*(-A_PR(1)-M_PR(1)) ! -A_PR(1)=M_ctr'+Somme des Mi pour i<=iRang_PR(1)
     $     -Perte_M*B_PR(2)+(M1M2-M_PR(1)*M_PR(2))*UnSurR2         ! -B_PR(2)=Somme des Mi/Ri pour i>iRang_PR(2)
     $     +FracMasseInt*(      ! terme d'autogravitation des SE
     $     UnSurR1*(M1M1-M_PR(1)**2)+UnSurR2*(M2M2-M_PR(2)**2) )

#if (_STATIQUE_==0)
      M_ctr = M_ctr + dM_accr
#endif
      dMaccr_evst = dMaccr_evst+dM_accr
      dMejec_evst = dMejec_evst+Perte_M-dM_accr
      dMaccr_evst_STS(iTet1) = dMaccr_evst_STS(iTet1)+
     $     frac_accr_by_MBH* Perte_M1
      dMaccr_evst_STS(iTet2) = dMaccr_evst_STS(iTet2)+
     $     frac_accr_by_MBH* Perte_M2
      dMejec_evst_STS(iTet1) = dMejec_evst_STS(iTet1)+
     $     (1.0d0-frac_accr_by_MBH)* Perte_M1
      dMejec_evst_STS(iTet2) = dMejec_evst_STS(iTet2)+
     $     (1.0d0-frac_accr_by_MBH)* Perte_M2
                                ! Modification du pot. central
      Psi0=Psi0+M_PR(1)/R_PR(1)+M_PR(2)/R_PR(2)
c
      end
c
#endif /* if (_EVOL_STELL_>0 || _ACCR_BY_STARS_>0) */
c
c----------------------------------------------------------------------
      subroutine Avancer_Paire_SE
c----------------------------------------------------------------------
c
      include 'EvolAmas_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_Rencontre_inc.f'
c......................................................................
                                ! Temps individuel
      Tps_PR(1) = Tps_PR(1) + dt_renc 
      Tps_PR(2) = Tps_PR(2) + dt_renc 
c
      end
c
c----------------------------------------------------------------------
      subroutine DonneesPR_a_SEact
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
c
c functions :
c ^^^^^^^^^^^
      double precision stellar_mass, stellar_radius
      integer*1 stellar_type
c......................................................................
c
      iSE_act = iSE_PR(Num_SE_act)
      iRang_act = iRang_PR(Num_SE_act)
 
      M_act = M_PR(Num_SE_act)
      R_act = R_PR(Num_SE_act)
      T_act = T_PR(Num_SE_act)
      J_act = J_PR(Num_SE_act)
      Met_act = stellar_mass(iSE_act)
      iTet_act = stellar_type(iSE_act)
      Ret_act = stellar_radius(iSE_act)
      Tps_act = Tps_PR(Num_SE_act)

      A_act = A_PR(Num_SE_act)
      B_act = B_PR(Num_SE_act)

      iPG_act = -1
      
      Pot_act = 0.0d0
      Peri_act = 0.0d0
      Apo_act = 0.0d0
c
      end
c
c----------------------------------------------------------------------
      subroutine DonneesSEact_a_PR
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
c......................................................................
c
      iSE_PR(Num_SE_act) = iSE_act
      iRang_PR(Num_SE_act) = iRang_act
 
      M_PR(Num_SE_act) = M_act
      R_PR(Num_SE_act) = R_act
      T_PR(Num_SE_act) = T_act
      J_PR(Num_SE_act) = J_act
      Tps_PR(Num_SE_act) = Tps_act

      A_PR(Num_SE_act) = A_act
      B_PR(Num_SE_act) = B_act
c
      end
c
c----------------------------------------------------------------------
      subroutine Extraire_SEact
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_Rencontre_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      double precision A,B
c......................................................................
c
c---- Upgrade contribution of star to the grid      
c
#if (_WRITE_OFGRID_>0)
      call Upgrade_OFGrid(iSE_act,dt_renc)
#endif
c
c==== Extraction de l'arbre binaire 
c     -----------------------------
c
      call Extraire_SE_Arbre(iSE_act,iRang_act,A,B)
c
c---- Dans le cas de la SE 1, les valeurs de A_PR et B_PR sont toujours valables :
c     elles tiennent compte des eventuelles modifications collisionnelles des masses
c     des 2 SE. D'autre part la SE 2 n'a pas ete bougee.
c
c     Dans le cas de la SE 2, il faut tenir compte de la nouvelle position
c     (et de l'eventuelle nouvelle masse) de la SE 1. Ceci est fait correctement par
c     le parcours de l'arbre. Par contre, cet arbre ne tient pas encore compte de
c     l'eventuelle modification de masse de la SE 2
c     
      if (Num_SE_act.EQ.2) then
         A_act = A+M_SE(iSE_act)-M_act
         B_act = B
      end if
      Pot_act = B_act +
     $     ( A_act+(1.0d0-FracMasseInt)*
     $     M_act )/R_SE(iSE_act)
c
c==== Extraction de la grille
c     -----------------------
c
      iPG_act = -1
      call Retirer_de_PG(R_SE(iSE_act),T_SE(iSE_act),M_SE(iSE_act),
     $     J_SE(iSE_act),iPG_act) 
c
c==== Extraction de l'amas
c     --------------------
c
      B0 = B0 - M_SE(iSE_act)/R_SE(iSE_act) ! mise a jour de B(R=0)
      Tps_amas = Tps_amas - Tps_SE(iSE_act)/(NbSE-NbSE_perte) ! mise a jour approx du temps de l'amas
      M_amas = M_amas - M_SE(iSE_act)
c
      end
c
c----------------------------------------------------------------------
      subroutine Modifier_SEact
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_Maree_common.f'
      include 'EvolAmas_LossCone_inc.f'
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_Etoiles_common.f'
      include 'EvolAmas_NonInterac_common.f'
      include 'EvolAmas_Rencontre_inc.f'
c
c var locales :
c ^^^^^^^^^^^^
      logical lPerte,l_ModifPot_dM2,lChgt,lNonRelax
      double precision E_act,dM2
      integer iLC
#if (_VERIF_PERI_APO_>0) 
      double precision la, lp, dl,
     $     Rp_1,Vr2p_1,Rp_2,Vr2p_2, Ra_1,Vr2a_1,Ra_2,Vr2a_2
#endif
#if (_TEST_PERI_APO_>0)
      double precision A_test,B_test,Vr2_peri_test,Vr2_apo_test,
     $     r_typ_test,Vr2_typ_test
      integer iRang
      double precision fact_test_vr2_periapo
      parameter ( fact_test_vr2_periapo=1.0d-4 )
#endif
c
c fonctions :
c ^^^^^^^^^^^
      double precision R_tronc_mar, Pot_Amas_Arbre_r
      logical SEact_A_Suivre
c......................................................................
c
      lPerte = .FALSE.
      iLC = iSAFE
c
c---- Modification du potentiel due a la modification de
c     masse de la seconde particule
c
c     PAS TRES ELEGANT !!!         XXXX 
c
      dM2 = M_PR(2)-M_SE(iSE_PR(2))
      l_ModifPot_dM2 = (Num_SE_Act.EQ.1 .AND. dM2.NE.0.0d0)
      if (l_ModifPot_dM2) then
         call Modif_Pot_Arbre(iNoeud_SE(iSE_PR(2)),dM2)
                                ! M_SE de la particule 2 est aussi modifiee !!
         B0 = B0+dM2/R_SE(iSE_PR(2))
      end if
c
c---- Is the star non relaxing?
c
      lNonRelax = .FALSE.
      if (lNonInterac) then
         lNonRelax = .NOT.lRelaxesWith(
     $        iTet_act, iTet_SE(mod(Num_SE_Act,2)+1) )
      end if
c
c---- Y a-t-il eu destruction collisionnelle totale ?
c
      if (Met_act.LE.0.0d0) then
         call PerteSEact_par_Collision
         lPerte = .TRUE.
      else                      ! Pas de destruction par collision
c
c---- Eventuellement, chauffage par binaires
c
#if (_CHAUFF_BIN_3C_>0)
         call Chauffage_SEact_par_bin_3c
#endif
                                ! calcul APOCENTRE
         E_act = T_act+Pot_act
         call CalcApo_SE_ext2
     $        (M_act,E_act,J_act,Apo_act)
c
c---- Y a-t-il evaporation ?
c

#if (_TRONC_MAREE_>0)
         if (Apo_act.LE.0.0d0 .OR. Apo_act.GT.R_tronc_mar()) then
#else
         if (Apo_act.LE.0.0d0) then
            if (E_act.LT.0d0) then
               write(0,*) '!!! E_act, Apo_act = ', E_act, Apo_act, ' !!!'
               call Meurt('Modifier_SEact')
            end if
#   if(_EMPECHE_EVAP_RELAX_>0)
            if (iTypeRenc.NE.RencUnbound) then
               write(0,100) iPas_Evol, iSE_act, Num_SE_act, E_act
 100           format('=== step ',I11,' evap. of part ',
     $              I7,' ( Num_SE_cat=',I1,') with E_act = ', E12.5,
     $              ' (but _EMPECHE_EVAP_RELAX_>0)')
               !call Meurt('Modifier_SEact')
            end if
#   endif
#endif
            call PerteSEact_par_Evap
            lPerte = .TRUE.
         else                   ! Pas d'evaporation

                                ! calcul PERICENTRE
            call CalcPeri_SE_ext2
     $           (M_act,E_act,J_act,Peri_act)
#if (_TEST_PERI_APO_>0)
            call Determiner_PotRang_Arbre_r(Peri_act,A_test,B_test,iRang)
            Vr2_peri_test = 2.0d0*(E_act-B_test)-(2.0d0*(A_test-FracMasseInt*
     $           M_act)+J_act**2/Peri_act)/Peri_act
            call Determiner_PotRang_Arbre_r(Apo_act,A_test,B_test,iRang)
            Vr2_apo_test = 2.0d0*(E_act-B_test)-(2.0d0*(A_test-FracMasseInt*
     $           M_act)+J_act**2/Apo_act)/Apo_act
            r_typ_test = 0.5d0*(Apo_act+Peri_act)
            Vr2_typ_test = J_act**2/(2.0d0*Apo_act*Peri_act) *
     $           (Apo_act-r_typ_test)*(r_typ_test-Peri_act)/r_typ_test**2
            if ( abs(Vr2_peri_test).GT.fact_test_vr2_periapo*Vr2_typ_test .OR.
     $           abs(Vr2_apo_test).GT.fact_test_vr2_periapo*Vr2_typ_test ) then
c$$$         write(75,'(2(E12.5," "))')
c$$$     $        Vr2_peri_test/Vr2_typ_test,Vr2_apo_test/Vr2_typ_test
               write(0,*) '!!! Vr2_peri_test,Vr2_apo_test/Vr2_typ_test = ',
     $              Vr2_peri_test/Vr2_typ_test,Vr2_apo_test/Vr2_typ_test,
     $              ' !!!'
               call Meurt('Modifier_SEact')
            end if
#endif
            if (Peri_act.LT.0.0d0) then
               write(0,*) '!!! Pas d''orbite physique pour la SE ',
     $              'active !!!'
               write(0,*) '    i_SE,M,E,J = ', iSE_act,
     $              M_act,E_act,J_act
               call Meurt('Modifier_SEact')
            end if
#if (_VERIF_PERI_APO_>0) 
c
c---- On teste le changement de signe de Vrad**2 au pericentre et a l'apocentre
c
            lp=log10(Peri_act)
            la=log10(Apo_act)
            dl=0.001d0*(la-lp)
            Rp_1 = 10**(lp-dl)  ! position juste avant pericentre
            Rp_2 = 10**(lp+dl)  ! position juste apres pericentre
            Ra_1 = 10**(la-dl)  ! position juste avant apocentre
            Ra_2 = 10**(la+dl)  ! position juste apres apocentre
            Vr2p_1 = 2.0d0*E_act - 2.0d0*Pot_Amas_Arbre_r(Rp_1)
     $           + 2.0d0*FracMasseInt*M_act/Rp_1
     $           - (J_act/Rp_1)**2 
            Vr2p_2 = 2.0d0*E_act - 2.0d0*Pot_Amas_Arbre_r(Rp_2)
     $           + 2.0d0*FracMasseInt*M_act/Rp_2
     $           - (J_act/Rp_2)**2 
            Vr2a_1 = 2.0d0*E_act - 2.0d0*Pot_Amas_Arbre_r(Ra_1)
     $           + 2.0d0*FracMasseInt*M_act/Ra_1
     $           - (J_act/Ra_1)**2 
            Vr2a_2 = 2.0d0*E_act - 2.0d0*Pot_Amas_Arbre_r(Ra_2)
     $           + 2.0d0*FracMasseInt*M_act/Ra_2
     $           - (J_act/Ra_2)**2 
            if (Vr2p_1.GT.0.0d0 .OR. Vr2p_2.LT.0.0d0 .OR.
     $           Vr2a_2.GT.0.0d0 .OR. Vr2a_1.LT.0.0d0) then
               write(0,100) Peri_act, Rp_1,Vr2p_1,Rp_2,Vr2p_2,
     $              Apo_act, Ra_1,Vr2a_1,Ra_2,Vr2a_2
 100           format('!!! Erreur lors du calcul Peri/Apocentres !!!',/,
     $              '    Peri_act : Rp_1,Vr2p_1 | Rp_2,Vr2p_2 : ',E12.5,
     $              ' : ',E12.5,' , ',E12.5,' | ',E12.5,' , ',E12.5,/,
     $              '    Apo_act  : Ra_1,Vr2a_1 | Ra_2,Vr2a_2 : ',E12.5,
     $              ' : ',E12.5,' , ',E12.5,' | ',E12.5,' , ',E12.5)
               call Meurt('Modifier_SEact')
            end if
#endif /* if (_VERIF_PERI_APO_>0) */

#if (_TEST_MODIF_PROP_SE_ACT_>0) /* test 26.4.2001 */
            call test_modif_prop_orb_SE_act (.TRUE.,M_act,
     $           E_act,J_act,T_act,pot_act,
     $           Peri_act,Apo_act,lChgt)
#endif
c
c---- Y a-t-il entree dans le loss-cone ?
c
#if (_DECHIREMENTS_>0 || _DISPARITIONS_>0 || _CAPT_GW_>0)
            if (M_TN.GT.0.0d0 .AND. .NOT.lNonRelax) then

c
c     On determine d'abord si il y a disparition (franchissement de l'horizon)
c     ou dechirement et, dans ce cas, avec quels
c     parametres (en fait le seul parametre est Beta)
c
               call Determiner_LossCone_SEact(iLC)

               if (iLC.EQ.iSAFE) then 
                  continue
#   if (_DECHIREMENTS_>0)
               else if (mod(iLC,100).EQ.iDISRP) ! tidal disruption or collision with central star
     $                 then                     ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                  call Interac_CentrObject_SEact
c     
c---- Y a-t-il destruction totale par le TN ?
c     C'est en fait la seule possibilite traitee correctement par le code
c     dans son etat actuel. Car une destruction partielle modifie l'orbite
c     ce qui necessiterait de re-calculer Peri- et Apocentre
c     (+ autres modifications eventuelles)
c
                  if (Met_act.LE.0.0d0) then 
                     call PerteSEact_par_Dechirement
                     lPerte = .TRUE.
                  end if
#   endif
#   if (_DISPARITIONS_>0)
               else if (mod(iLC,100).EQ.iCROSS) ! direct horizon crossing
     $                 then                     ! ^^^^^^^^^^^^^^^^^^^^^^^
                  call Calculer_Disp_SEact
                  call PerteSEact_par_Disparition
                  lPerte = .TRUE.
#   endif
#   if (_CAPT_GW_>0)
               else if (mod(iLC,100).EQ.iGWCPT) ! captage par emission de GW
     $                 then                     ! ^^^^^^^^^^^^^^^^^^^^^^^^^^
                  call Calculer_CaptGW_SEact 
                  call PerteSEact_par_CaptageGW
                  lPerte = .TRUE.
#   endif
               else
                  write(0,*) '!!! Invalid outcome for ',
     $                 'Determiner_LossCone_SEact: iLC = ',iLC,' !!!'
                  call Meurt('Modifier_SEact')
               end if

            end if
#endif /* (_DECHIREMENTS_>0 || _DISPARITIONS_>0 || _CAPT_GW_>0) */
         end if ! (Apo_act.LE.0.0d0)
      end if ! (Met_act.LE.0.0d0)
c
c---- Determination d'une nouvelle position orbitale 
c
#if (_STATIQUE_==0)
#if (_DEPL_ORBITAL_>0)
      if (.NOT.lPerte) then
         call Determiner_PosOrb_SEact
      end if
#endif
#else /* (_STATIQUE_!=0) */
                                ! on annule toutes les modifications
      M_act = M_SE(iSE_act)
      T_act = T_SE(iSE_act)
      J_act = J_SE(iSE_act)
      R_act = R_SE(iSE_act)
                                !One cannot easily undo modifications to stellar properties because they are done in-place!

                                ! on fait tout de meme progresser l'etoile
                                ! pour avoir une idee du temps "equivalent" ecoule
      !Tps_act = Tps_SE(iSE_act)
                                ! on determine une nouvelle position sur la meme orbite
#if (_DEPL_ORBITAL_>0)
      E_act = T_act+Pot_act
      call CalcApo_SE_ext2
     $     (M_act,E_act,J_act,Apo_act)
      if (Apo_act.GT.0.0d0) then
         call CalcPeri_SE_ext2
     $        (M_act,E_act,J_act,Peri_act)
         call Determiner_PosOrb_SEact
      end if
#endif
#endif /* if (_STATIQUE_==0) */
c
c---- On annule la modification du potentiel due a la
c     modification de masse de la seconde SE
c     car elle sera refaite lors de la re-introduction
c     de celle-ci dans l'arbre...
c
c     PAS TRES ELEGANT !!!
c
      if (l_ModifPot_dM2) then
         call Modif_Pot_Arbre (iNoeud_SE(iSE_PR(2)),-dM2)
         B0 = B0-dM2/R_SE(iSE_PR(2))
      end if
c
c---- Eventuellement, consigner certaines proprietes de la SE
c
#if (_ECRIT_SUIVI_SE_PART_>0)
      if (SEact_A_Suivre()) call Ecrire_SuiviSE(iFich_SuivSE)
#endif
c
      end
c
#if (_CHAUFF_BIN_3C_>0)
c----------------------------------------------------------------------
      subroutine Chauffage_SEact_par_bin_3c
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_Rencontre_inc.f'
      include 'EvolAmas_ChauffBin3c_common.f'
      include 'EvolAmas_Grille_param.f'
      include 'EvolAmas_Grille_common.f'
c
c constante(s) :
c ^^^^^^^^^^^^^^
      double precision dErel_max
      parameter (dErel_max=1.0d1) ! augmentation relative maximale admissible de l'energie cinetique
c
c var locales :
c ^^^^^^^^^^^^^
      double precision T_prec,dE,V2_c,Met_moy, Rho, Sigma
c
c fonctions :
c ^^^^^^^^^^^
      double precision Somme_MsurR_Arbre, Sigma1D_PG, Sigma1D_PG_lin2
c......................................................................
c
#if (_CHAUFF_BIN_3C_==1)
                                ! dispersion centrale des vitesses (3D)
      V2_c = 2.0d0*T_PG(1)/M_PG(1) !## Somme_MsurR_Arbre() !## potentiel central
                                ! augmentation d'energie specifique due aux binaires "3-corps" 
      Met_moy_loc =             ! masse stellaire moyenne (ponderation par la masse)
     $     (MM_PG(iPG_renc)/(Net_par_SE*M_PG(iPG_renc)))

      dE =  dt_renc *
     $     V2_c * C4_1D * Rho_renc**2 *
     $     Met_moy_loc**3 / Sigma1D_PG(iPG_renc,R_renc)**9
#elif (_CHAUFF_BIN_3C_==2)
      V2_c = 0.99d0*V2_c + 0.01d0*2.0d0*T_PG(1)/M_PG(1)
                                ! version ultra-simplifiee valable pour des masses stellaires egales
                                ! "etendue" pour tenir compte d'un spectre de masse si une masse domine au centre
      Met_moy = M_PG(iPG_renc)/(Net_par_SE*NSE_PG(iPG_renc))
      Rho = Rho_renc
      Sigma = Sigma1D_PG_lin2(iPG_renc,R_renc)
      dE = dt_renc * V2_c * C4_1D *  Rho**2 * Met_moy**3 / Sigma**9
#endif
                                ! verifier que l'augmentation d'energie n'est pas trop forte
      T_prec = T_act
      if (dE.GT.-dErel_max*Etot_amas/M_amas) then
         write(0,*) '!!! dE_chauffage, Eliaison_moy = ',dE,
     $        -Etot_amas/M_amas,' !!!'
         call Meurt('Chauffage_SEact_par_bin_3c')
      end if
                                ! augmenter la vitesse en fct du chauffage (isotropie)
      T_act = T_prec + dE
      J_act = sqrt(T_act/T_prec)*J_act
                                ! suivi de l'energie fournie par les binaires
      dEtot_bin3c = dEtot_bin3c + dE*M_act
c
      end
c
#endif /* (_CHAUFF_BIN_3C_>0) */

#if (_DECHIREMENTS_>0 || _DISPARITIONS_>0 || _CAPT_GW_>0)
c----------------------------------------------------------------------
      subroutine Determiner_LossCone_SEact(iFlagLC)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_Rencontre_inc.f'
      include 'EvolAmas_LossCone_inc.f'
      include 'EvolAmas_GW_common.f'
      include 'EvolAmas_Relax_common.f'
c     
c     argument :
c     ^^^^^^^^^^
      integer iFlagLC           ! Est mis > 0 si il y a dechirement (ou collision avec VMS centrale) 
                                !      et < 0 si il y a traversee de l'horizon
                                !         = 0 sinon 
c     
c     constantes :
c     ^^^^^^^^^^^^
      double precision DemiPi
      parameter (DemiPi=0.5d0*rPi)
      double precision Coef_Secu_Diff
      parameter (Coef_Secu_Diff=1.0d30) 
      double precision Coef_Sep_FE  ! coefficient de separtion des cas empty/full loss-cone
      parameter (Coef_Sep_FE=1.0d-30)! une valeur de 0 correspond a 100% de cas traites dans la limite "empty"
                                     ! une valeur tres grande, a 100% de cas traites dans la limite "full"
c                                    ! la valeur normale est 1     
c     var locales :
c     ^^^^^^^^^^^^^
      double precision
     $     Theta_disp,Theta_dechir,Theta_lc,Vlc,
     $     Porb_harm, Porb_quad, Trel_quad,
     $     Theta_relax_loc, Theta_relax_typ, Theta_relax_test,
     $     Theta_V, Vpart_renc(3), Vtg, V, Vcarre, fN_orb,
     $     Pdechir_1, Pdechir_n, Pdechir_n_sup, delta_PA, dT_PA,
     $     Vtg_ini, R_peri, Trel_loc, Trel_renc
      integer N_PA, Nreal_PA, iIssue
      logical lInit_GW
c  
c     fonction :
c     ^^^^^^^^^^
      double precision random,TrelLoc,Somme_MsurR_Arbre
      logical lTest_GW_Capture
c......................................................................
c     
c---- Initialisations
c
      iFlagLC = iSAFE
c
c---- calcul de l'angle d'ouverture du LC (pour la traversee de l'horizon et le dechirement de maree)
c
#if (_DISPARTIONS_==2 || _DECHIREMENTS_==2) /* Remove particle if semi-major axis smaller than given (fixed) value */
      if ((Apo_act+Peri_act).LT.2.0d0*R_tid_disr_def_NB) then
#   if (_DECHIREMENTS_==2)
         iFlagLC = iRLX_DISRP
#   endif
#   if (_DISPARTIONS_==2)
         iFlagLC = iRLX_CROSS
#   endif
      return
#endif
      call Calc_LC_SEact(Theta_disp,Theta_dechir)
#if (_DECHIREMENTS_==0)
      Theta_dechir = 0.0d0
#endif
#if (_DISPARITIONS_==0)
      Theta_disp   = 0.0d0
#endif
      Theta_lc = max(Theta_disp,Theta_dechir)
c
c---- Compute orbital period (and orbit-averaged relaxation time)
c
#if (_ECRIT_LC_>0 || (_CAPT_GW_>0 && _TYPE_TREL_GW_==1))
      call Calc_PorbTrel_SEact(Porb_quad,Trel_quad)
#else
      call Calc_Porb_SEact(Porb_quad)
#endif
      fN_orb = dt_renc/Porb_quad ! nb d'orbites
#if (_CAPT_GW_>0)
#  if (_TYPE_TREL_GW_==0) /* By default, use value from encounter */
      Trel_over_Porb_GW = fN_orb*(DemiPi/Theta_relax_renc)**2
#  elif (_TYPE_TREL_GW_==1) /* Integration over orbit; scaled by 1/mass of star */
      Trel_over_Porb_GW = Met_moy/Met_act * Trel_quad/Porb_quad
#  elif (_TYPE_TREL_GW_==2) /* fixed value (!!) */
      Trel_over_Porb_GW = 1.0d9*UneAnnee*Trlx_GW_def_Gyr/Porb_quad
#  else 
      PROBLEM XXXX
#  endif
#endif
c
c---- || Si la rencontre etait une collision ou un kick, pas d'autre chose a faire que de
c     || regarder si l'etoile se trouve ou non dans le LC. Il n'y a pas de
c     || diffusion a traiter de maniere subtile !
c
      if (iTypeRenc.EQ.Collision .OR. iTypeRenc.EQ.Kick) then
#if (_CI_RW_AVANT_RENC_==0) /* on prend les condition APRES rencontre comme c.i. du parcours aleat. */
         V = sqrt(2.0d0*T_act)
         Vtg = J_act/R_act
#else /* on prend les condition AVANT rencontre comme c.i. du parcours aleat. */
         V = sqrt(2.0d0*T_SE(iSE_act))
         Vtg = J_SE(iSE_act)/R_SE(iSE_act)
#endif
                                ! Is the star tidally disrupted?
                                ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         Theta_V = asin(Vtg/V)  ! direction de la vitesse de la SE (la fct asin est sans doute inutile)
         if (Theta_V.LT.Theta_lc) then ! on est entre dans le LC
            if (iTypeRenc.EQ.Collision) then
               if (Theta_disp.LT.Theta_dechir) then
                  iFlagLC = iCOLL_DISRP
               else
                  iFlagLC = iCOLL_CROSS
               end if
            else
               if (Theta_disp.LT.Theta_dechir) then
                  iFlagLC = iKICK_DISRP
               else
                  iFlagLC = iKICK_CROSS
               end if
            end if
#if (_PRESENCE_TN_>0 && _CAPT_GW_>0)
         else                   ! Is the star captured by grav wave emission?
            lInit_GW=.TRUE.     ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            if (lTest_GW_Capture(lInit_GW,V**2,Vtg**2)) then
               if (iTypeRenc.EQ.Collision) then
                  iFlagLC = iCOLL_GWCPT
               else
                  iFlagLC = iKICK_GWCPT
               end if
            end if
#endif
         end if
c
c---- || Si la rencontre etait un "super-encounter" de relaxation, il faut tenter de
c     || tenir compte du processus diffusif sous-jacent qui va determiner le risque
c     || de dechirement.
c
      else if (iTypeRenc.EQ.Relaxation) then
c
c---- Test de degrossissage pour savoir si on est dans la limite
c     full loss-cone
c
c         Porb_harm =            ! Borne inf. sur la periode orbitale (approx. harmonique)
c     $        rPi*Peri_act*Apo_act/J_act
                                ! angle quadratique total de diffusion ...
         Theta_relax_typ  =     ! ...valeur typique 
     $        DemiPi*sqrt(Frac_Trelax)
         Theta_relax_loc  =     ! ...en fct du tps de relaxation local
     $        DemiPi*sqrt(dt_renc/TrelLoc(R_renc))
         
         Theta_relax_test = Theta_relax_loc ! l'angle quad de diffusion employe dans les tests

c$$$c!#######\\\\# essai 3.5.2001
c$$$c
c$$$c---- Calcul de la probabilite de destruction base sur le travail de Spurzem
c$$$c     avec approximation isotrope
c$$$c
c$$$            Omega = 0.25d0*Theta_lc**2
c$$$            T_in  = C_T_in*Omega*TrelLoc(R_renc)
c$$$            T_out = C_T_out*Porb_quad
c$$$            Plc = T_in/(T_in+T_out)
c$$$            Tlc = (T_out+Plc*T_in)/Omega ! Devrait y a avoir un facteur 1/Plc encore...
c$$$            P_dechir = dt_renc/Tlc
c$$$            write(88,*) P_dechir
c$$$            if (random().LT.P_dechir) then
c$$$               iFlagLC = EMPTY
c$$$               V = sqrt(2.0d0*T_act)
c$$$               Vtg = 0.5d0*V*sin(Theta_lc)
c$$$            end if
c$$$            goto 7
c$$$c!######/////#

#if (_CI_RW_AVANT_RENC_==0) /* on prend les condition APRES rencontre comme c.i. du parcours aleat. */
         V = sqrt(2.0d0*T_act)
         Vtg = J_act/R_act
#else /* on prend les condition AVANT rencontre comme c.i. du parcours aleat. */
         V = sqrt(2.0d0*T_SE(iSE_act))
         Vtg = J_SE(iSE_act)/R_SE(iSE_act)
#endif
         Theta_V = asin(Vtg/V)  ! direction de la vitesse de la SE (la fct asin est sans doute inutile)
         if ((Theta_lc+Coef_Secu_Diff*Theta_relax_test).GT.Theta_V)
     $        then
c
c---- Il y a risque de rentree dans le LC, il faut donc calculer
c     le parcours aleatoire de la direction de la vitesse
c     Remarquons que ce test n'est pas vraiment correct car il faudrait
c     le faire dans le referentiel de la rencontre... A moins qu'il ne faille
c     employer un autre angle de diffusion...
c
            iNb_tests_LC = iNb_tests_LC + 1
#if (_CI_RW_AVANT_RENC_==0) /* on prend les condition APRES rencontre comme c.i. du parcours aleat. */
            if (Num_SE_act.EQ.1) then
               Vpart_renc(1) = V1amas(1)-Vcm_renc(1)
               Vpart_renc(2) = V1amas(2)-Vcm_renc(2)
               Vpart_renc(3) = V1amas(3)-Vcm_renc(3)
            else
               Vpart_renc(1) = V2amas(1)-Vcm_renc(1)
               Vpart_renc(2) = V2amas(2)-Vcm_renc(2)
               Vpart_renc(3) = V2amas(3)-Vcm_renc(3)
            end if
#else /* on prend les condition AVANT rencontre comme c.i. du parcours aleat. */
            if (Num_SE_act.EQ.1) then
               Vpart_renc(1) = Mu2*Vrel(1)
               Vpart_renc(2) = Mu2*Vrel(2)
               Vpart_renc(3) = Mu2*Vrel(3)
            else
               Vpart_renc(1) = -Mu1*Vrel(1)
               Vpart_renc(2) = -Mu1*Vrel(2)
               Vpart_renc(3) = -Mu1*Vrel(3)
            end if
#endif
c
c---- Si possible, on emploie un pas plus grand dans le parcours aleatoire
c     afin d'avoir moins de pas a calculer. Le probleme est que les cond ini
c     de dechirements ainsi obtenues ne sont pas exactes (on penetre trop
c     profondement dans le LC...)
c
            delta_PA = Theta_relax_renc/sqrt(fN_orb)
            N_PA = int(fN_orb+0.5d0)
            Vlc = V*sin(Theta_lc)
            call ParcoursAlea_LC(Vpart_renc,Vcm_renc,
     $           Vlc,delta_PA,N_PA,Nreal_PA,Vtg,Vcarre,iIssue)
            if (iIssue.EQ.1)
     $           then           ! il y a dechirement dans regime "empty"
               if (Theta_disp.LT.Theta_dechir) then
                  iFlagLC = iRLX_DISRP
               else
                  iFlagLC = iRLX_CROSS
               end if
            else if (iIssue.EQ.2)
     $              then        ! captage par orbite a forte emission d'ondes grav.
               iFlagLC = iRLX_GWCPT
#if (_ANNONCE_SE_DANS_LC_>=1)
            else
               Vtg_ini = J_act/R_act
               if (Vtg_ini.LT.Vlc) then
                  write(0,400) iPas_Evol,iSE_act,Vlc,Vtg_ini
 400              format(
     $                 '> Determiner_LossCone_SEact: iPas_Evol = ',
     $                 I10,
     $                 ' iSE = ',I7,' laissee dans le LC (Vlc = ',
     $                 E12.5,') avec Vtg_ini = ',E12.5)
               end if
#endif
            end if
         end if

c$$$c!###############\\\\\\! essai 9.5.01
c$$$         if (iFlagLC.NE.0) then
c$$$            dT_PA = 0.5d0*M_act*(Vcarre-V*V)
c$$$            if (dT_PA.LT.0.0d0) then ! perte d'energie => on donne l'energie correspondante a l'autre SE!
c$$$               if (Num_SE_act.EQ.1) then
c$$$                  T_PR(2) = T_PR(2) - dT_PA/M_PR(2)
c$$$                  dTtot_dechir = dTtot_dechir - dT_PA
c$$$                  dEtot_dechir = dEtot_dechir - dT_PA
c$$$               else             ! l'autre SE a deje ete re-introduite
c$$$                  if (R_SE(iSE_PR(1)).LT.R_perte) then
c$$$                     T_SE(iSE_PR(1)) = T_SE(iSE_PR(1)) -
c$$$     $                    dT_PA/M_SE(iSE_PR(1))
c$$$                     dTtot_dechir = dTtot_dechir - dT_PA
c$$$                     dEtot_dechir = dEtot_dechir - dT_PA
c$$$                  end if
c$$$               end if
c$$$            end if
c$$$         end if
c$$$C!###############//////!

      end if ! if (iTypeRenc.EQ.Collision)
c$$$$ 7    continue !######
c      
      if (iFlagLC.NE.iSAFE) then
         if (mod(iFlagLC,100).NE.iGWCPT) then ! On exclue les cas de captage par GW
#if (_DECHIR_NULS_==0)
            J_act = R_act*Vtg
#endif
c
c     on ne modifie pas la norme de la vitesse, donc E, cela permet d'evacuer
c     les problemes de conservation de l'energie et est coherent avec le fait
c     que l'entree dans le LC est essentiellement reglee par la diffusion en J
c     (la limite du LC est quasi-independante de E)
c
            R_peri = J_act**2 / (2.0d0*M_ctr) ! Approximation parabolique
c !!! adapt here for VMS
c
c     Beta = R_star/R_peri * (M_BH/M_star)**(1/3)
c
            Beta_dechir = Coef_ConvR*Ret_act/R_peri
     $           * (M_TN/(Coef_ConvM*Met_act))**0.33333333d0
c
c     Il faut voir si il s'agit d'un dechirement ou d'une traversee
c     de l'horizon.
c     Si l'horizon est dans la zone de Roche, on admet
c     qu'il y a forcement dechirement, bien que l'on puisse imaginer
c     une etoile plongeant diretement dans l'horizon "en ligne droite",
c     mais cela semble peu probable!
c
         else
            Beta_dechir = -1.0d0
         end if

#if (_ECRIT_LC_>0)
         if (Theta_relax_loc.GT.0.0d0) then
            Trel_loc = dt_renc*(DemiPi/Theta_relax_loc)**2
         else
            Trel_loc = -1.0d30
         end if
         if (Theta_relax_renc.GT.0.0d0) then
            Trel_renc = dt_renc*(DemiPi/Theta_relax_renc)**2
         else
            Trel_renc = -1.0d30
         end if
         write(iFich_LC,130) iPas_Evol,iSE_act,
     $        R_act,
     $        T_act+Pot_act,J_act,
     $        iTet_act,Met_act,Theta_disp,
     $        Theta_dechir,dt_renc,Porb_quad,Beta_dechir,iFlagLC,
     $        0.5d0*(Vcarre-V*V),-Somme_MsurR_Arbre(),
     $        Trel_loc, Trel_renc,
     $        Trel_quad, Peri_act,Apo_act
      
 130     format(I10,' ',I7,' ',3(' ',E10.3),' ',I3,6(' ',E10.3),' ',I3,
     $        7(' ',E10.3))
         call flush_MD(iFich_LC)
#endif
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_LC_SEact(Theta_Disp,Theta_Dechir)
c----------------------------------------------------------------------
c     Calcul des angles sous-tendus par l'horizon et la zone de Roche
c     tels que "vus" par la SE active
c
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_LossCone_inc.f'
      include 'EvolAmas_Etoiles_param.f'
c
c argument :
c ^^^^^^^^^^
      double precision Theta_Disp, Theta_Dechir
c     
c constantes :
c ^^^^^^^^^^^^
      double precision coef_lc
      parameter (coef_lc=1.0d0) ! Facteur sur la taille du LC; 1.0d0 en principe
c
c fonctions :
c ^^^^^^^^^^^
      double precision
     $     Rmar_TN, Rhoriz_eff_TN, ThetaCrit_SEact,
     $     R_VMS
c......................................................................
c
      if (iType_CtrObj.EQ.TypeEtoile_BH) then

#if (_CONST_R_TID_DISRUP_==0)
         Theta_dechir = ThetaCrit_SEact(
     $        coef_lc*Rmar_TN(Met_act,Ret_act)/
     $        Beta_Dechir_crit )
#else /* Same tidal disruption radius for all stars */
         Theta_dechir = ThetaCrit_SEact(
     $        coef_lc*R_tid_disr_def_NB/UL_Ncorps) ! In principle UL_Ncorps=1, but this is safer
#endif
         Theta_disp   = ThetaCrit_SEact( Rhoriz_eff_TN() )

      else if (iType_CtrObj.EQ.TypeEtoile_MS) then

         Theta_dechir = ThetaCrit_SEact(coef_lc*(R_VMS()+
     $        Coef_ConvR*Ret_act))
         Theta_disp   = 0.0d0
      else
         write(0,*) '!!! Don''t know how to deal with central',
     $        ' object of type ',iType_CtrObj,' !!!'
         call Meurt('Calc_LC_SEact')
      end if

c
      end
c
c----------------------------------------------------------------------
      double precision function ThetaCrit_SEact(R)
c----------------------------------------------------------------------
c     Gives the "loss cone" aperture angle for the acrive SE to reach
c     pericenter distance R
c
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R
c
c fonctions :
c ^^^^^^^^^^^
      double precision Somme_MsurR_Arbre, Pot_Amas_Arbre_r
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision sin2
c......................................................................
c
#if (_CALC_THETA_LC_EXACT_>0)   
                                ! Modified 23 1 2007 to correct mistake [sqrt(asin()) instead of asin(sqrt())]
                                ! and include test on argument of asin()

                                ! Not really exact if there is an extebded gas "reservoir" around the MBH
c      sin2 = (R/R_act)**2 *(
c     $     1.0d0 + (
c     $     M_ctr/R*(1.0d0-R/R_act) +
c     $     Pot_act + Somme_MsurR_Arbre() ) / T_act
c     $     )                              ! Not really exact if there is an extended gas "reservoir" around the MBH

c
c---- The active SE has been extracted from the tree so it does not contribute 
c     any more to the potential Pot_Amas_Arbre_r
c
      sin2 = (R/R_act)**2 *(    ! Corrected formula (?), 2 May 2008
     $     ( Pot_act - Pot_Amas_Arbre_r(R) ) / T_act + 1.0d0
     $     )
      sin2 = min(sin2,1.0d0)    ! sin2 can formally be > 1 if all orbits enter the 
                                ! "disruption radius", even when the velocity is initially
                                ! tangential
      if (sin2.LT.0.0d0) then
         write(0,*) '!!! sin2 = ',sin2,' !!!'
         write(0,200) R,R_act,M_act,T_act,J_act,
     $        M_ctr,Pot_act,Pot_Amas_Arbre_r(R)
 200     format(
     $        '  R,R_act,M_act,T_act,J_act : ',5(E14.7,' '),/,
     $        '  M_ctr,Pot_act,Pot_Amas_Arbre_r() : ',
     $        3(E14.7,' '))
         call Meurt('ThetaCrit_SEact')
      end if
      ThetaCrit_SEact = ASIN(sqrt(sin2))

#else /* Tres bonne approximation */

#  if (_CI_RW_AVANT_RENC_==0) /* on prend les condition APRES rencontre comme c.i. du parcours aleat. */
      ThetaCrit_SEact = sqrt(
     $     R*M_ctr / (T_act*R_act**2) )
#  else /* on prend les condition AVANT rencontre comme c.i. du parcours aleat. */
      ThetaCrit_SEact = sqrt(
     $     R*M_ctr / (T_SE(iSE_act)*R_SE(iSE_act)**2) )
#  endif
#endif /* if (_CALC_THETA_LC_EXACT_>0) */
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_Porb_SEact(Porb)! Compute orbital period
                                      ! by Gauss-Chebychev quadrature
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_LossCone_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision Porb
c
c var. locales :
c ^^^^^^^^^^^^^^
      double precision Peri,Apo
c     
c constantes :
c ^^^^^^^^^^^^
      double precision eps_IntegrOrb
      parameter (eps_IntegrOrb=1.0d-4) ! precision de l'integration orbitale de Porb (2.0d-3)
c......................................................................
c
                                ! Valeur "exacte" de Porb pour la nouvelle orbite
                                ! calculee par integration orbitale de 1/Vrad
#if (_CI_RW_AVANT_RENC_==0) /* on prend les condition APRES rencontre comme c.i. du parcours aleat. */
      call Calc_Porb_SE(M_act,T_act+Pot_act,J_act,
     $     Peri_act,Apo_act,eps_IntegrOrb,Porb) 
#else /* on prend les condition AVANT rencontre comme c.i. du parcours aleat. */
      call CalcPeri_SE_ext2(M_SE(iSE_act),
     $     T_SE(iSE_act)+Pot_act,J_SE(iSE_act),Peri)
      call CalcApo_SE_ext2(M_SE(iSE_act),
     $     T_SE(iSE_act)+Pot_act,J_SE(iSE_act),Apo)
      call Calc_Porb_SE(M_SE(iSE_act),
     $     T_SE(iSE_act)+Pot_act,J_SE(iSE_act),
     $     Peri,Apo,eps_IntegrOrb,Porb)
#endif
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_PorbTrel_SEact(Porb,Trel) ! Compute orbital period and orbit-averaged relaxation time 
                                                ! by Gauss-Chebychev quadrature
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_LossCone_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision Porb,Trel
c
c var. locales :
c ^^^^^^^^^^^^^^
      double precision Peri,Apo
c     
c constantes :
c ^^^^^^^^^^^^
      double precision eps_IntegrOrb
      parameter (eps_IntegrOrb=1.0d-4) ! precision de l'integration orbitale de Porb (2.0d-3)
c......................................................................
c
                                ! Valeur "exacte" de Porb pour la nouvelle orbite
                                ! calculee par integration orbitale de 1/Vrad
#if (_CI_RW_AVANT_RENC_==0) /* on prend les condition APRES rencontre comme c.i. du parcours aleat. */
      call Calc_PorbTrel_SE(M_act,T_act+Pot_act,J_act,
     $     Peri_act,Apo_act,eps_IntegrOrb,Porb,Trel) 
#else /* on prend les condition AVANT rencontre comme c.i. du parcours aleat. */
      call CalcPeri_SE_ext2(M_SE(iSE_act),
     $     T_SE(iSE_act)+Pot_act,J_SE(iSE_act),Peri)
      call CalcApo_SE_ext2(M_SE(iSE_act),
     $     T_SE(iSE_act)+Pot_act,J_SE(iSE_act),Apo)
      call Calc_PorbTrel_SE(M_SE(iSE_act),
     $     T_SE(iSE_act)+Pot_act,J_SE(iSE_act),
     $     Peri,Apo,eps_IntegrOrb,Porb,Trel)
#endif
c
      end
c
#endif  /* if (_DECHIREMENTS_>0 || _DISPARITIONS_>0 || _CAPT_GW_>0) */

#if (_CAPT_GW_>0)
c----------------------------------------------------------------------
      logical function lTest_GW_Capture(lInit,V2,Vtg2)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_GW_common.f'
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_EvolStell_common.f'
c
c constants :
c ^^^^^^^^^^^
      double precision DemiPi2
      parameter (DemiPi2=0.25d0*rPi*rPi)
      double precision coef_T_GW
      parameter (coef_T_GW=5.0d0/(128.0d0*rPi))
      double precision cf_GW_1, cf_GW_2
      parameter (cf_GW_1=73.0d0/24.0d0,cf_GW_2=37.0d0/96.0d0)
c     
c argument :
c ^^^^^^^^^^
      logical lInit
      double precision V2,Vtg2
c
c local variables :
c ^^^^^^^^^^^^^^^^^^
      double precision Met_GW
c......................................................................
c
      if (lInit) then
         GMsurR_GW = M_TN/R_act
         GMsurR2_GW = GMsurR_GW/R_act
                                ! coefficient dans le calcul du rapport T_GW/Porb
#if (_DIFFERENT_M_FOR_GW_==0)
         CT_GW = coef_T_GW *
     $        sqrt(VitLum2**5/M_TN**3)/(Coef_ConvM*Met_act)
#else 
         if (iTet_act.EQ.TypeEtoile_WD) then
            Met_GW = M_WD
         else if (iTet_act.EQ.TypeEtoile_NS) then
            Met_GW = M_NS
         else if (iTet_act.EQ.TypeEtoile_BH) then
            Met_GW = M_BH
         else
            Met_GW = Met_act
         end if

         CT_GW = coef_T_GW *
     $        sqrt(VitLum2**5/M_TN**3)/(Coef_ConvM*Met_GW)
#endif
         lInit = .FALSE.
      end if

      El_GW = GMsurR_GW-0.5d0*V2 ! Energie de liaison au TN
      if (El_GW.LT.0.0d0) then
         lTest_GW_Capture = .FALSE.
         return
      end if

      a_GW = 0.5d0*M_TN/El_GW
      e_GW = 1.0d0-Vtg2/(GMsurR2_GW*a_GW)
      if (e_GW.LT.0.0d0) then
         write(0,*) '!!! e^2 = ',e_GW,' V,Vtg =,',
     $        sqrt(V2),sqrt(Vtg2),' !!!'
         call Meurt('lTest_GW_Capture')
      end if
      e_GW = sqrt(e_GW)
      Rperi_GW = a_GW*(1.0d0-e_GW)
      theta2_GW = 2.0d0*GMsurR2_GW*Rperi_GW/V2
      T_GW_sur_Porb = CT_GW *
     $     sqrt( a_GW**5 * ((1.0d0+e_GW)*(1.0d0-e_GW))**7 ) /
     $     (1.0d0+e_GW*e_GW*(cf_GW_1+cf_GW_2*e_GW*e_GW))

      lTest_GW_Capture = (DemiPi2*T_GW_sur_Porb .LT.
     $     Coef_Trlx_GW_Capt*theta2_GW*Trel_over_Porb_GW
     $     .AND. T_GW_sur_Porb .LT. 1.0d9*UneAnnee*Tmax_GW_Capt_Gyr)
c
      end
c
#endif /* (_CAPT_GW_>0) */
#if (_DECHIREMENTS_>0)
c----------------------------------------------------------------------
      subroutine Interac_CentrObject_SEact
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_LossCone_inc.f'
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_SuiviTypesStell_common.f'
c
c var locales :
c ^^^^^^^^^^^^^
      double precision Mfin,dEorb,T_prec,
     $     Perte_M,Perte_T,UnSurR,MM,dM_accr,somme_MsurR,
     $     Mass_VMS, HeFrac_VMS
      integer*1 iTet
c
c fonctions :
c ^^^^^^^^^^^
      double precision Somme_MsurR_Arbre,stellar_age,
     $     stellar_centralheliumfrac
c......................................................................
c
c---- initialisations
c
      iTet = iTet_act
      UnSurR = 1.0d0/R_act
      MM = M_act**2
      Perte_M = M_act
      Perte_T = M_act*T_act

#if (_COLL_REJUV_==2)
      if (iType_CtrObj.EQ.TypeEtoile_MS) then
c
c     Mass of Helium in the core of the star (contribution to core of central star after merger)
c
         HeFrac_VMS = stellar_centralheliumfrac(iSE_act)*Met_act
      end if
#endif   
                                ! on calcule la somme de M/R en tenant compte
                                ! des masses eventuellement modifiees
      somme_MSurR = Somme_MsurR_Arbre()+M_act*UnSurR
                                ! La correction suivante est en principe nulle
                                ! puisque la modification de masse de la SE 2
                                ! a ete provisoirement appliquee a l'arbre
      if (Num_SE_Act.EQ.1) somme_MSurR = somme_MSurR
     $     - M_SE(iSE_PR(2))/R_SE(iSE_PR(2))+M_PR(2)/R_PR(2)

      if (iType_CtrObj.EQ.TypeEtoile_BH) then
c
c---- Calcul du resultat du dechirement par le TN
c
         call Dechirement_par_TN(
     $        iTet_act,Met_act,beta_dechir,
     $        Mfin,dEorb)
      else if (iType_CtrObj.EQ.TypeEtoile_MS) then
c         call Collision_with_CentralMSStar(XXX) !!! Adapt for VMS
c
c---- For the time being, we assume perfect merger with the central star
c
         Mfin = 0.0d0 
      else
         write(0,*) '!!! Don''t know how to deal with central',
     $        ' object of type ',iType_CtrObj,' !!!'
         call Meurt('Interac_CentrObject_SEact')
      end if
c     
c---- modification des proprietes de la particule active
c
                                ! Diminution de la masse
      M_act = M_act*Mfin/Met_act
      Perte_M = Perte_M - M_act
      if (Mfin.GT.0.0d0) then
#if (_COLL_REJUV_==0 || _COLL_REJUV_==1) /* No rejuvenation, keep same absolute age */
         call stellar_SetByAge_ms(iSE_act,Mfin,StellMet_def,
     $        stellar_age(iSE_act))
#elif (_COLL_REJUV_==2) /* Use of core He mass to compute effective age assuming M_He rises linearly with time */
         call stellar_SetByHeFrac_ms(iSE_act,Mfin,StellMet_def,
     $        min(1.0d0,
     $        stellar_centralheliumfrac(iSE_act)*Met_act/Mfin) )
#else
         write(0,*) '!!! Invalid value of _COLL_REJUV_ !!!'
         call Meurt('Interac_CentrObject_SEact')
#endif
      end if

      Met_act = Mfin

      if (M_ctr.GT.0.0d0) then   ! masse accretee par le TN
         dM_accr = frac_accr_Dechir*Perte_M
      else
         dM_accr = 0.0d0
      end if
                                ! Diminution de la vitesse (perte d'energie orbitale)
      if (Mfin.GT.0.0d0) then
         T_prec = T_act
         T_act = T_prec - dEorb
         J_act = sqrt(T_act/T_prec)*J_act
      else
         T_act = 0.0d0
         J_act = 0.0d0
      end if
                                ! -Masse a l'interieur du rayon de la SE
      A_act = A_act + Perte_M - dM_accr
c
c---- bilan de masse et d'energie (La masse totale de l'amas est mise
c     a jour lors de la re-introduction de la SE)
c
      Perte_T = Perte_T -M_act*T_act ! perte d'energie cinetique orbitale (>0)
      dTtot_dechir = dTtot_dechir - Perte_T
      dEtot_dechir = dEtot_dechir - Perte_T
     $     -dM_accr*somme_MsurR
     $     +Perte_M*UnSurR*(-A_act-M_act) ! -A_act=M_ctr'+Somme des Mi pour i<=iRang_act
     $     -Perte_M*B_act ! -B_act=Somme des Mi/Ri pour i>iRang_act
     $     +FracMasseInt*UnSurR*(MM-M_act**2) ! terme d'autogravitation de la SE

#if (_STATIQUE_==0)
      M_ctr = M_ctr+dM_accr     ! Change of M_TN is done in Repartir_Maccr_rsrv_TN 
                                ! unless the central object is on the MS, in which case
                                ! its state is update below
#endif
      dMaccr_dechir = dMaccr_dechir+dM_accr
      dMejec_dechir = dMejec_dechir+Perte_M-dM_accr
      dMaccr_dechir_STS(iTet) = dMaccr_dechir_STS(iTet)+dM_accr
      dMejec_dechir_STS(iTet) = dMejec_dechir_STS(iTet)+Perte_M-dM_accr
c
c---- Update mass and effective age of central object if it is a MS star
c
      if (iType_CtrObj.EQ.TypeEtoile_MS) then
         Mass_VMS = (M_TN+dM_accr)/Coef_ConvM ! New mass in M_sun
#if (_COLL_REJUV_==0) /* No rejuvenation, keep same absolute age  */
         call stellar_SetByAge_ms(id_CtrObj,Mass_VMS,
     $        StellMet_def,stellar_age(id_CtrObj))
#elif (_COLL_REJUV_==1) /* Reset age to 0 in case of MS-MS merger */
         call stellar_SetByAge_ms(id_CtrObj,Mass_VMS,
     $        StellMet_def,0.0d0)
#elif (_COLL_REJUV_==2) /* Use of core He mass to compute effective age assuming M_He rises linearly with time */
         !!write(99,*) stellar_centralheliumfrac(id_CtrObj), Mass_VMS, M_TN/Coef_ConvM, HeFrac_VMS !##########
         HeFrac_VMS = min(1.0d0, (HeFrac_VMS +
     $        stellar_centralheliumfrac(id_CtrObj)*M_TN/Coef_ConvM)/
     $        Mass_VMS)
         call stellar_SetByHeFrac_ms(id_CtrObj,Mass_VMS,
     $        StellMet_def,HeFrac_VMS)
         !!write(99,*) stellar_centralheliumfrac(id_CtrObj), Mass_VMS, M_TN/Coef_ConvM, HeFrac_VMS !##########
#else
         write(0,*) '!!! Invalid value of _COLL_REJUV_ !!!'
         call Meurt('Interac_CentrObject_SEact')
#endif
      end if
c
      iNb_Dechir = iNb_Dechir+1
c
      end
c
#endif /* if (_DECHIREMENTS_>0) */ 

#if (_DISPARITIONS_>0)
c----------------------------------------------------------------------
      subroutine Calculer_Disp_SEact
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_SuiviTypesStell_common.f'
c
c var locales :
c ^^^^^^^^^^^^^
      double precision 
     $     Perte_M,Perte_T,UnSurR,MM,dM_accr,somme_MsurR
      integer*1 iTet
c
c fonctions :
c ^^^^^^^^^^^
      double precision Somme_MsurR_Arbre
c......................................................................
c
c---- initialisations
c
      UnSurR = 1.0d0/R_act
      MM = M_act**2
      Perte_M = M_act
      Perte_T = M_act*T_act
                                ! on calcule la somme de M/R en tenant compte
                                ! des masses eventuellement modifiees
      somme_MSurR = Somme_MsurR_Arbre()+M_act*UnSurR
                                ! La correction suivante est en principe nulle
                                ! puisque la modification de masse de la SE 2
                                ! a ete provisoirement appliquee a l'arbre
      if (Num_SE_Act.EQ.1) somme_MSurR = somme_MSurR
     $     - M_SE(iSE_PR(2))/R_SE(iSE_PR(2))+M_PR(2)/R_PR(2)
c
c---- modification des proprietes de la particule active
c
      M_act = 0.0d0
      T_act = 0.0d0
      J_act = 0.0d0
      dM_accr = Perte_M
c
c---- bilan de masse et d'energie (La masse totale de l'amas est mise
c     a jour lors de la re-introduction de la SE)
c
      dTtot_disp = dTtot_disp - Perte_T
      dEtot_disp = dEtot_disp - Perte_T
     $     -dM_accr*somme_MsurR
     $     +Perte_M*UnSurR*(-A_act) ! -A_act=M_ctr'+Somme des Mi pour i<=iRang_act
     $     -Perte_M*B_act ! -B_act=Somme des Mi/Ri pour i>iRang_act
     $     +FracMasseInt*UnSurR*MM ! terme d'autogravitation de la SE

#if (_STATIQUE_==0)
      M_ctr = M_ctr+dM_accr
      M_TN = M_TN+dM_accr
#endif
      dMaccr_disp = dMaccr_disp+dM_accr
      iTet = iTet_act
      dMaccr_disp_STS(iTet) = dMaccr_disp_STS(iTet)+dM_accr
c
      iNb_Disp = iNb_Disp+1
c
      end
c
#endif /* if (_DISPARITIONS_>0) */ 
#if (_CAPT_GW_>0)
c----------------------------------------------------------------------
      subroutine Calculer_CaptGW_SEact
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_PaireRenc_inc.f'
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_SuiviTypesStell_common.f'
c
c var locales :
c ^^^^^^^^^^^^^
      double precision 
     $     Perte_M,Perte_T,UnSurR,MM,dM_accr,somme_MsurR
      integer*1 iTet
c
c fonctions :
c ^^^^^^^^^^^
      double precision Somme_MsurR_Arbre
c......................................................................
c
c---- initialisations
c
      iTet = iTet_act
      UnSurR = 1.0d0/R_act
      MM = M_act**2
      Perte_M = M_act
      Perte_T = M_act*T_act
                                ! on calcule la somme de M/R en tenant compte
                                ! des masses eventuellement modifiees
      somme_MSurR = Somme_MsurR_Arbre()+M_act*UnSurR
                                ! La correction suivante est en principe nulle
                                ! puisque la modification de masse de la SE 2
                                ! a ete provisoirement appliquee a l'arbre
      if (Num_SE_Act.EQ.1) somme_MSurR = somme_MSurR
     $     - M_SE(iSE_PR(2))/R_SE(iSE_PR(2))+M_PR(2)/R_PR(2)
c
c---- modification des proprietes de la particule active
c
      M_act = 0.0d0
      T_act = 0.0d0
      J_act = 0.0d0
      dM_accr = Perte_M
c
c---- bilan de masse et d'energie (La masse totale de l'amas est mise
c     a jour lors de la re-introduction de la SE)
c
      dTtot_captGW = dTtot_captGW - Perte_T
      dEtot_captGW = dEtot_captGW - Perte_T
     $     -dM_accr*somme_MsurR
     $     +Perte_M*UnSurR*(-A_act) ! -A_act=M_ctr'+Somme des Mi pour i<=iRang_act
     $     -Perte_M*B_act ! -B_act=Somme des Mi/Ri pour i>iRang_act
     $     +FracMasseInt*UnSurR*MM ! terme d'autogravitation de la SE

#if (_STATIQUE_==0)
      M_ctr = M_ctr+dM_accr
      M_TN = M_TN+dM_accr
#endif
      dMaccr_captGW = dMaccr_captGW+dM_accr
      dMaccr_captGW_STS(iTet) = dMaccr_captGW_STS(iTet)+dM_accr
c
      iNb_CaptGW = iNb_CaptGW+1
c
      end
c
#endif /* if (_CAPT_GW_>0) */ 
c----------------------------------------------------------------------
      subroutine PerteSEact_par_Evap
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_SuiviTypesStell_common.f'
      include 'EvolAmas_Rencontre_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      integer*1 iTet
      integer iTypeEvap
c......................................................................
c
c---- bilan de masse et d'energie (La masse totale de l'amas est mise
c     a jour lors de la re-introduction des SE)
c
      iTypeEvap = 1
#if (_VERB_PERTES_>0)
      write(0,*) '~~~ Perte de la SE ',iSE_act,' (Num_SE_act=',
     $     Num_SE_act,') par Evaporation'
#endif
      iTet = iTet_act
      dMejec_evap = dMejec_evap + M_act
      dMejec_evap_STS(iTet) = dMejec_evap_STS(iTet) + M_act
      dTtot_evap = dTtot_evap - M_act*T_act
      dEtot_evap = dEtot_evap - M_act*
     $     (T_act+Pot_act)

#if (_ECRIT_EVAP_>0)
      if (iTypeRenc.EQ.RencUnbound) iTypeEvap = 3
      write(iFich_Evap,100) iPas_Evol,iTypeEvap,iSE_act,M_act,
     $     R_act,T_act,Pot_act,J_act
 100  format(I10,' ',I2,' ',I8,5(' ',E11.4))
      call flush_MD(iFich_Evap)
#endif

      NbSE_perte = NbSE_perte+1
      NbSE_evap = NbSE_evap+1
      rNbEtoiles = rNbEtoiles-Net_par_SE
      R_act = R_evap
      Tps_act = Tps_perte
      call stellar_Destroy(iSE_act)
c
      end
c
c----------------------------------------------------------------------
      subroutine PerteSEact_par_Collision
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
c......................................................................
c
#if (_VERB_PERTES_>0)
      write(0,*) '~~~ Perte de la SE ',iSE_act,' (Num_SE_act=',
     $     Num_SE_act,') par Collision'
#endif
      NbSE_perte = NbSE_perte+1
      NbSE_dest_coll = NbSE_dest_coll+1
      rNbEtoiles = rNbEtoiles-Net_par_SE
      R_act = R_dest_coll
      Tps_act = Tps_perte
      call stellar_Destroy(iSE_act)
c
      end
c
c----------------------------------------------------------------------
      subroutine PerteSEact_par_Dechirement
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
c......................................................................
c
#if (_VERB_PERTES_>0)
      write(0,*) '~~~ Perte de la SE ',iSE_act,' (Num_SE_act=',
     $     Num_SE_act,') par Dechirement'
#endif
      NbSE_perte = NbSE_perte+1
      NbSE_dest_dechir = NbSE_dest_dechir+1
      rNbEtoiles = rNbEtoiles-Net_par_SE
      R_act = R_dechir
      Tps_act = Tps_perte
      call stellar_Destroy(iSE_act)
c
      end
c
c----------------------------------------------------------------------
      subroutine PerteSEact_par_CaptageGW
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
c......................................................................
c
#if (_VERB_PERTES_>0)
      write(0,*) '~~~ Perte de la SE ',iSE_act,' (Num_SE_act=',
     $     Num_SE_act,') par Captage par ondes. grav.'
#endif
      NbSE_perte = NbSE_perte+1
      NbSE_dest_captGW = NbSE_dest_captGW+1
      rNbEtoiles = rNbEtoiles-Net_par_SE
      R_act = R_captGW
      Tps_act = Tps_perte
      call stellar_Destroy(iSE_act)
c
      end
c
c----------------------------------------------------------------------
      subroutine PerteSEact_par_Disparition
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
c......................................................................
c
#if (_VERB_PERTES_>0)
      write(0,*) '~~~ Perte de la SE ',iSE_act,' (Num_SE_act=',
     $     Num_SE_act,') par Disparition'
#endif
      NbSE_perte = NbSE_perte+1
      NbSE_dest_disp = NbSE_dest_disp+1
      rNbEtoiles = rNbEtoiles-Net_par_SE
      R_act = R_disp
      Tps_act = Tps_perte
      call stellar_Destroy(iSE_act)
c
      end
c
c----------------------------------------------------------------------
      subroutine Determiner_PosOrb_SEact
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_PeriApo_inc.f'
      include 'EvolAmas_CtrlCentre_inc.f'
      include 'EvolAmas_TirOrb2_common.f'
      include 'Flag_common.f'
c
c constantes :
c ^^^^^^^^^^^^
      integer Nmax_tir_pos
      parameter (Nmax_tir_pos=200)
c
c var locales :
c ^^^^^^^^^^^^^
      double precision Rnouveau,Vrad2nouveau,Energie,Eself
      integer itir_pos
c......................................................................
c
      Energie = T_act+Pot_act
      itir_pos = 0

 1    continue
      call TirerRorb_SE (
     $     M_act,Energie,J_act,
     $     iRang_act,
     $     Peri_act,Apo_act,A_peri,A_apo,iRang_peri,
     $     Rnouveau,Vrad2nouveau )
      i_NbTirVN_PosOrb = i_NbTirVN_PosOrb + iCompt_essais_TO
      itir_pos = itir_pos + 1
      
#if (_CONTROLE_CENTRE_==1)
      if (iRang_act.LT.Rang_CC .AND. Rnouveau.LT.FactR_CC*R_act) then
                                ! on laisse la SE a son ancienne position
         if (itir_pos.EQ.Nmax_tir_pos) then
            write(0,*) '> Nb max d''essais de tirage ',
     $           'd''une position orbitale (CONTROLE_CENTRE) atteint'
            Rnouveau = R_act
            Vrad2nouveau = 2.0d0*T_act
     $        - (J_act/R_act)**2
         else
            goto 1
         end if
      end if
#elif (_CONTROLE_CENTRE_==2)
      if (iRang_act.LT.Rang_CC .AND. Rnouveau.LT.Rmin_CC) then
         Rnouveau = R_act       ! on laisse la SE a son ancienne position
         Vrad2nouveau = 2.0d0*T_act
     $        - (J_act/R_act)**2
      end if
#elif (_CONTROLE_CENTRE_==3) /* on empeche la SE de devenir autogravitante */
      Eself = 0.5d0*( Vrad2nouveau + (J_act/Rnouveau)**2 ) -
     $     FracMasseInt*M_act/Rnouveau
      if (Eself .LT. 0.0d0) then
         if (itir_pos.EQ.Nmax_tir_pos) then ! on laisse la SE a son ancienne position
            write(0,*) '> Nb max d''essais de tirage ',
     $           'd''une position orbitale (CONTROLE_CENTRE) atteint'
            Rnouveau = R_act
            Vrad2nouveau = 2.0d0*T_act
     $           - (J_act/R_act)**2
         else
            goto 1
         end if
      end if
#endif

      if (Rnouveau.LT.Rmin_CC) then
         Rmin_CC = Rnouveau
         iSE_Rmin_CC = iSE_act
      end if

      R_act = Rnouveau
      T_act = 0.5d0*( Vrad2nouveau + (J_act/Rnouveau)**2 )

      Pot_act = Energie - T_act
c
      end
c
c----------------------------------------------------------------------
      subroutine Reinserer_SEact
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_SEact_inc.f'
      include 'EvolAmas_Grille_param.f'
      include 'EvolAmas_Grille_common.f'
c
c var. locales :
c ^^^^^^^^^^^^^^
      double precision A,B
      integer iPG_prec
c......................................................................
c
c---- Mise a jour des prop. de la SE active
c
      M_SE(iSE_act) = M_act
      T_SE(iSE_act) = T_act
      J_SE(iSE_act) = J_act
      R_SE(iSE_act) = R_act
      Tps_SE(iSE_act) = Tps_act
#if (_TEST_ZERO_MASS_==1)
      if (M_SE(iSE_act).LE.0.0d0 .AND. R_SE(iSE_act).LT.R_perte) then
         write(0,*) '!!! How funny, active SE ',iSE_act,
     $        ' has zero mass !!!'
         call Meurt('Reinserer_SEact')
      end if
#endif
c
c==== Reinsertion dans l'arbre
c     ------------------------
c
      call Inserer_SE_Arbre(iSE_act,iRang_act,A,B)
c
c==== Reinsertion dans la grille
c     --------------------------
c
      iPG_prec = iPG_act
      iPG_act = -1
      call Mettre_dans_PG(R_act,T_act,M_act,J_act,iPG_act)
c
c---- on determine si il y a lieu de reconstruire la grille
c
      if  (NSE_PG(iPG_prec).LT.Nmin_PG) then
         flag_Grille = -1
#if (_VERBOSITUDE_>5)
         write(0,100)iPG_prec,NSE_PG(iPG_prec),Nmin_PG
 100     format('|| Nb min de SE par cellule atteint NSE_PG(',
     $        I4,') = ',I3,'< ',I3)
#endif    
      end if
      if ( iPG_act.GE.1 .AND. iPG_act.LE.N_PG
     $     .AND. (NSE_PG(iPG_act).GT.Nmax_PG) ) then
         flag_Grille = -2
#if (_VERBOSITUDE_>5)
         write(0,110)iPG_act,NSE_PG(iPG_act),Nmax_PG
 110     format('|| Nb max de SE par cellule atteint NSE_PG(',
     $        I4,') = ',I3,'> ',I3)
#endif    
      endif
c
c==== Reinsertion dans l'amas
c     -----------------------
c
      if (R_act.LT.R_perte) then
         B0 = B0 + M_act/R_act  ! mise a jour de B(R=0)
         M_amas = M_amas + M_act
         Tps_amas = Tps_amas + Tps_act/(NbSE-NbSE_perte)
      end if
c
      end
c
c@@@@@@@@@@@@@@@ FIN REFONTE @@@@@@@@@@@@@@@@@@@@@@@@
c----------------------------------------------------------------------
      subroutine Affiche_Info
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_Glob_common.f'
      include 'EvolAmas_GlobTech_common.f'
      include 'EvolAmas_NbPas_common.f'
      include 'EvolAmas_LimExt_common.f'
      include 'EvolAmas_Etoiles_param.f'
      include 'EvolAmas_NonInterac_common.f'
c
c constantes :
c ^^^^^^^^^^^^
      double precision Max_err_E
      parameter (Max_err_E=1.0d-9)
      double precision prec_Rc
      parameter (prec_Rc=0.01d0)
c
c fonction(s) :
c ^^^^^^^^^^^^^
      real TpsUtil_MD
      double precision Somme_MsurR_Arbre
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      double precision UnSurN,UnSurN_rel,Rc, etot_, ilast_step_info
      save ilast_step_info
      data ilast_step_info /0/
c......................................................................
c
      if (ilast_step_info.EQ.0) ilast_step_info = i_ini
      call Calc_EnerGlob
      call Calc_Rcoeur(Rc,prec_Rc) ! calcul la densite et la dispersion centrale
                                   ! (dans common de 'EvolAmas_Glob_common.f')
      Pot_ctr = -Somme_MsurR_Arbre() ! potentiel central (contribution stellaire

      call QuelleDate_MD(jour,mois,annee)
      call QuelleHeure_MD(heure,minute,sec)

      if (ilast_step_info.LT.iPas_Evol) then
         UnSurN = 1.0d0/dble(iPas_Evol-ilast_step_info)
      else
         UnSurN = 1.0d30
      end if
      ilast_step_info = iPas_Evol

      Frac_Pas_Relax = iNb_SuperRenc*UnSurN
      if (iNb_SuperRenc.EQ.0) iNb_SuperRenc=-1
      UnSurN_rel = 1.0d0/dble(iNb_SuperRenc)

      Tps_CPU = TpsUtil_MD()
#if (_LIMITE_EXT_==0)
      Viriel = -(Ecin_amas+Etot_amas)/Etot_amas
#else
      etot_ = Ecin_amas + Ecin_ext_LimExt +
     $     Egrav_amas + Egrav_ext_LimExt - B_LimExt*(M_amas+M_ctr)
      Viriel = -(Ecin_amas+Ecin_ext_LimExt+etot_)/etot_
#endif
                                ! Note that the energy change due to natal kicks
                                ! is included in dEtot_evst (and dTtot_evst)
      Erreur_rel_E = ( Etot_ini-Etot_amas
     $     +dEtot_evap+dEtot_coll+dEtot_evst
     $     +dEtot_dechir+dEtot_disp
     $     +dEtot_captGW+dEtot_bin3c)
      Erreur_rel_M = ( Mamas_ini-M_amas
     $     -dMejec_evap-dMejec_coll-dMejec_evst-dMejec_dechir
     $     -dMaccr_coll-dMaccr_evst-dMaccr_dechir-dMaccr_disp
     $     -dMaccr_captGW )
#if ( _ACCR_BY_STARS_>0)
      Erreur_rel_E = Erreur_rel_E / min(Etot_amas,Etot_ini)
      Erreur_rel_M = Erreur_rel_M / max(M_amas,Mamas_ini)
#else
      Erreur_rel_E = Erreur_rel_E / Etot_ini
      Erreur_rel_M = Erreur_rel_M / Mamas_ini
#endif
      Nb_moy_VN_pos_orb_par_pas = i_NbTirVN_PosOrb*UnSurN
      Nb_moy_parcrs_arbre_par_pas = Nb_Parcours_Arbre*UnSurN
      if (i_ini.LT.iPas_Evol) then
         Nb_moy_test_LC_par_pas = dble(iNb_tests_LC)/
     $        dble(iPas_Evol-i_ini)
      else
         Nb_moy_test_LC_par_pas = 0.0d0
      end if
c
c     Monitoring des super-rencontres pour simuler la relaxation
c
      Dev_SuperRenc_moy   = Dev_SuperRenc_moy*UnSurN_rel
      Dev2_SuperRenc_moy  = sqrt(Dev2_SuperRenc_moy*UnSurN_rel)
      Frac_SuperRenc_Pi16 = iNb_SuperRenc_Pi16*UnSurN_rel
      Frac_SuperRenc_Pi8  = iNb_SuperRenc_Pi8 *UnSurN_rel
      Frac_SuperRenc_Pi4  = iNb_SuperRenc_Pi4 *UnSurN_rel
      Frac_SuperRenc_Pi2  = iNb_SuperRenc_Pi2 *UnSurN_rel
c
c     Monitoring des cas a grande probabilite de collisions
c
#if (_COLLISIONS_>0)
      P_coll_moy          = P_coll_moy*UnSurN
      Frac_Coll_P_025     = iNb_Coll_P_025*UnSurN
      Frac_Coll_P_05      = iNb_Coll_P_05 *UnSurN
      Frac_Coll_P_1       = iNb_Coll_P_1  *UnSurN
      Frac_Coll_P_2       = iNb_Coll_P_2  *UnSurN
#endif
c
c     Monitoring des cas a grande probabilite de dif. gr. angles
c
#if (_KICKS_>0)
      P_kick_moy          = P_kick_moy*UnSurN
      Frac_Kick_P_025     = iNb_Kick_P_025*UnSurN
      Frac_Kick_P_05      = iNb_Kick_P_05 *UnSurN
      Frac_Kick_P_1       = iNb_Kick_P_1  *UnSurN
      Frac_Kick_P_2       = iNb_Kick_P_2  *UnSurN
#endif
c
c     Monitoring des rencontres "non-locales"
c
      dR_rel_Renc_moy  = dR_rel_Renc_moy*UnSurN
      Frac_Renc_dR_001 = iNb_Renc_dR_001*UnSurN
      Frac_Renc_dR_01  = iNb_Renc_dR_01 *UnSurN
      Frac_Renc_dR_05  = iNb_Renc_dR_05 *UnSurN
      Frac_Renc_dR_1   = iNb_Renc_dR_1  *UnSurN

c
c---- ecriture a l'ecran
c
      write (0,100) iPas_Evol
 100  format(/,
     $     '|=======================================================|',/,
     $     '|=>>>>>>>>>>>>>>> iPas_Evol : ',I10.10,' <<<<<<<<<<<<<<=|',/,
     $     '|=======================================================|' )
      call Ecrire_Info_Phys(0)
      call Ecrire_Info_Tech(0)
      call flush_MD(0)
c
c---- ecriture dans le fichier "Log"
c
      write (iFich_Log,100) iPas_Evol
      call Ecrire_Info_Phys(iFich_Log)
      call Ecrire_Info_Tech(iFich_Log)
      call flush_MD(iFich_Log)
c
c---- Re-initialisation des compteurs
c
      i_NbTirVN_PosOrb = 0
      Nb_Parcours_Arbre = 0
      iNb_SuperRenc = 0
      iNb_RejetPaire = 0
      iNb_AnnulRelaxEvap = 0
      iNb_AnnulRelaxCC = 0

      Dev_SuperRenc_moy = 0.0d0
      Dev2_SuperRenc_moy = 0.0d0
      iNb_SuperRenc_Pi16 = 0
      iNb_SuperRenc_Pi8  = 0
      iNb_SuperRenc_Pi4  = 0
      iNb_SuperRenc_Pi2  = 0
      iNb_Coll_P_025 = 0
      iNb_Coll_P_05  = 0
      iNb_Coll_P_1   = 0
      iNb_Coll_P_2   = 0
      P_coll_moy = 0.0d0
      P_coll_max = 0.0d0
      iNb_Kick_P_025 = 0
      iNb_Kick_P_05  = 0
      iNb_Kick_P_1   = 0
      iNb_Kick_P_2   = 0
      P_kick_moy = 0.0d0
      P_kick_max = 0.0d0
      dR_rel_Renc_moy = 0.0d0
      dR_rel_Renc_max = 0.0d0
      iNb_Renc_dR_001 = 0
      iNb_Renc_dR_01  = 0
      iNb_Renc_dR_05  = 0
      iNb_Renc_dR_1   = 0

#if (_SURV_ERR_E_>0 && _FORCER_CROISSANCE_TN_==0 && _EVOL_UNE_SE_==0 && _STATIQUE_==0)
      if (.NOT.lNonInteracAsymmetrical .AND.
     $     abs(Erreur_rel_E).GT.Max_err_E) then
         write(0,*) '!!! Energie non conservee !!!'
         call Meurt('Affiche_Info')
      end if
#endif
c
      end
c
c----------------------------------------------------------------------
      subroutine Quitter
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c......................................................................
c
c---- Fermeture des fichiers
c
      close(iFich_Log)
#if (_ECRIT_RAYLAG_>0)
      close(iFich_RayLag)
#endif
#if (_ECRIT_LAGQUANT_>0)
      close(iFich_LagQuant)
#endif
#if (_ECRIT_ANILAG_>0)
      close(iFich_AniLag)
#endif
#if (_ECRIT_SPECMASSE_>0)
      close(iFich_MSpec)
#endif
#if (_ECRIT_SEGR_>0)
      close(iFich_Segr)
#endif
#if (_ECRIT_RELAX_>0)
      close(iFich_Relax)
#endif
#if (_ECRIT_COLL_>0)
      close(iFich_Coll)
#endif
#if (_ECRIT_STEVOL_>0)
      close(iFich_StEvol)
#endif

      write(0,*) ' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      write(0,*) ' @@@@ SIMULATION TERMINEE @@@@'
      write(0,*) ' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
      call exit(iexit_stat)
c
      end
c
