#define _RECALC_ENERGLOB_ 0
c======================================================================
c     routine d'entree-sortie specifiques a EvolAmas.F
c======================================================================

c
c     dans les routines d'ouverture de fichier qui suivent,
c     iAction indique l'action a realiser en fonction de l'etat du
c     fichier, selon la convention :
c
c     iAction | Fichier existe | Fichier n'existe pas
c     -----------------------------------------------
c         0   | rien            | rien
c         1   | rien            | creer
c         2   | ouvrir au debut | rien < These is an input file; the name may or may not start
c         3   | ouvrir au debut | creer  with the name (suffixe, such as '_EvolAmas%%') 
c         4   | ouvrir au bout  | rien   of the simulation
c         5   | ouvrir au bout  | creer
c
c----------------------------------------------------------------------
      subroutine OuvrirFichASC(iFich,TypeFich,iAction,lExiste)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
      character*(*) TypeFich
      integer iAction
      logical lExiste
c
c var locales :
c ^^^^^^^^^^^^^
      character*32 ChaineProv
      character*160 NomFich
      common /common_NomFich_EvolAmas/ NomFich
      integer iOK
      logical lOK
c......................................................................
c
      lOK = .TRUE.
c
c---- former le nom du fichier et voir si il existe
c
      if ( iFich.GE.50 .AND. iFich.NE.iFich_Prov) then ! fichiers sequentiels
         write(ChaineProv,100)TypeFich
 100     format('%%',A,'.asc')
      else                      ! fichiers lies au pas
         write(ChaineProv,110)iPas_Evol,TypeFich
 110     format('%',I10.10,'%',A,'.asc')
      end if
      call ConcatChaines(NomSimul,ChaineProv,NomFich)
      inquire(file=NomFich,exist=lExiste)

      if      (iAction.EQ.0) then
         continue
      else if (iAction.EQ.1) then
         if (.NOT.lExiste) then
            open(iFich,file=NomFich,form='formatted',IOSTAT=iOK)
            lOK = (iOK.EQ.0)
         end if
      else if (iAction.EQ.2) then
         if (.NOT.lExiste) then ! Maybe a the file exits with a name 
                                ! that does not include the infamous suffix
            write(NomFich,115)TypeFich
 115        format(A,'.asc')
            inquire(file=NomFich,exist=lExiste)
         end if
         if (lExiste) then
            open(iFich,file=NomFich,form='formatted',IOSTAT=iOK)
            lOK = (iOK.EQ.0)
         end if
      else if (iAction.EQ.3) then
         open(iFich,file=NomFich,form='formatted',IOSTAT=iOK)
         lOK = (iOK.EQ.0)
      else if (iAction.EQ.4) then
         if (lExiste) then
            call open_append(iFich,NomFich,lOK)
         end if
      else if (iAction.EQ.5) then
         call open_append(iFich,NomFich,lOK)
      else 
         write(0,*) '!!! mauvaise valeur de iAction : ',iAction,' !!!'
         call Meurt('OuvrirFichASC')
      end if
      if (.NOT.lOK) then
         write(0,*) '!!! Cannot open file as requested : !!!'
         call Meurt('OuvrirFichASC')
      end if
         
c
      end
c
c----------------------------------------------------------------------
      subroutine OuvrirFichBIN(iFich,TypeFich,iAction,lExiste)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
      character*(*) TypeFich
      integer iAction
      logical lExiste
c
c var locales :
c ^^^^^^^^^^^^^
      character*32 ChaineProv
      character*160 NomFich
      common /common_NomFich_EvolAmas/ NomFich
      integer iOK
      logical lOK
c......................................................................
c
c---- former le nom du fichier et voir si il existe
c
      lOK=.TRUE.
      if ( iFich.GE.50 .AND. iFich.NE.iFich_Prov) then ! fichiers sequentiels
         write(ChaineProv,100)TypeFich
 100     format('%%',A,'.bin')
      else                      ! fichiers lies au pas
         write(ChaineProv,110)iPas_Evol,TypeFich
 110     format('%',I10.10,'%',A,'.bin')
      end if
      call ConcatChaines(NomSimul,ChaineProv,NomFich)
      inquire(file=NomFich,exist=lExiste)

      if      (iAction.EQ.0) then
         continue
      else if (iAction.EQ.1) then
         if (.NOT.lExiste) then
            open(iFich,file=NomFich,form='unformatted',IOSTAT=iOK)
            lOK = (iOK.EQ.0)
         end if
      else if (iAction.EQ.2) then
         if (lExiste) then
            open(iFich,file=NomFich,form='unformatted',IOSTAT=iOK)
            lOK = (iOK.EQ.0)
         end if
      else if (iAction.EQ.3) then
         open(iFich,file=NomFich,form='unformatted',IOSTAT=iOK)
         lOK = (iOK.EQ.0)
      else if (iAction.EQ.4) then
         if (lExiste) then
            call open_append_BIN(iFich,NomFich,lOK)
         end if
      else if (iAction.EQ.5) then
         call open_append_BIN(iFich,NomFich,lOK)
      else 
         write(0,*) '!!! mauvaise valeur de iAction !!!'
         call Meurt('OuvrirFichBIN')
      end if
      if (.NOT.lOK) then
         write(0,*) '!!! Cannot open file as requested :',
     $            TypeFich,iAction,' !!!'
         call Meurt('OuvrirFichBIN')
      end if
c
      end
c
c----------------------------------------------------------------------
      logical function ReOuvrirFichBIN(iFich,TypeFich)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
      character*(*) TypeFich
      character*160 NomFich
      common /common_NomFich_EvolAmas/ NomFich
c......................................................................
c
      call OuvrirFichBIN(iFich,TypeFich,2,ReOuvrirFichBin)
      if (ReOuvrirFichBin) then
         write(0,100)TypeFich,NomFich
 100     format('> fichier de Type ',A,' present :',/,'   ',A40)
      else
         write(0,110)TypeFich,NomFich
 110     format('> fichier de Type ',A,' ABSENT :',/,'   ',A40)
      end if
c
      end
c
c----------------------------------------------------------------------
      logical function ReOuvrirFichASC(iFich,TypeFich)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
      character*(*) TypeFich
      character*160 NomFich
      common /common_NomFich_EvolAmas/ NomFich
c......................................................................
c
      call OuvrirFichASC(iFich,TypeFich,2,ReOuvrirFichAsc)
      if (ReOuvrirFichAsc) then
         write(0,100)TypeFich,NomFich
 100     format('> fichier de Type ',A,' present :',/,'   ',A40)
      else
         write(0,110)TypeFich,NomFich
 110     format('> fichier de Type ',A,' ABSENT :',/,'   ',A40)
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine OuvrirEcritureFichXDR(iFich_XDR,TypeFich,
     $     lExiste,lEcrase)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich_XDR
      character*(*) TypeFich
      logical lExiste, lEcrase
c
c var locales :
c ^^^^^^^^^^^^^
      character*32 ChaineProv
      character*160 NomFich
c......................................................................
c
c---- former le nom du fichier
c
      write(ChaineProv,110)iPas_Evol,TypeFich
 110  format('%',I10.10,'%',A,'.xdr')
      call ConcatChaines(NomSimul,ChaineProv,NomFich)
      inquire(file=NomFich,exist=lExiste)
c
c---- ouvrir le fichier
c
      if (.NOT.lExiste .OR. lEcrase)
     $     iFich_XDR=initxdr(NomFich,'w',.TRUE.)
c
      end
c
c----------------------------------------------------------------------
      subroutine OuvrirLectureFichXDR(iFich_XDR,TypeFich,lExiste)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'fxdr.inc'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich_XDR
      character*(*) TypeFich
      logical lExiste
c
c var locales :
c ^^^^^^^^^^^^^
      character*32 ChaineProv
      character*160 NomFich
      common /common_NomFich_EvolAmas/ NomFich
c......................................................................
c
c---- former le nom du fichier
c
      write(ChaineProv,110)iPas_Evol,TypeFich
 110  format('%',I10.10,'%',A,'.xdr')
      call ConcatChaines(NomSimul,ChaineProv,NomFich)
      inquire(file=NomFich,exist=lExiste)
c
c---- ouvrir le fichier
c
      if (lExiste) then
         iFich_XDR = initxdr(NomFich,'r',.TRUE.)
      else
         iFich_XDR = -1
      end if
c
      end
c
c----------------------------------------------------------------------
      logical function ReOuvrirLectureFichXDR(iFich_XDR,TypeFich)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich_XDR
      character*(*) TypeFich

      character*160 NomFich
      common /common_NomFich_EvolAmas/ NomFich
c......................................................................
c
      call OuvrirLectureFichXDR(iFich_XDR,TypeFich,
     $     ReOuvrirLectureFichXDR)
      if (ReOuvrirLectureFichXDR) then
         write(0,100)TypeFich,NomFich
 100     format('> fichier de Type ',A,' present :',/,'   ',A40)
      else
         write(0,110)TypeFich,NomFich
 110     format('> fichier de Type ',A,' ABSENT :',/,'   ',A40)
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine SupprimerFichXDR(iPas,TypeFich,iComport)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer*8 iPas            ! numero du pas correspondant au fichier
                                ! a effacer
      character*(*) TypeFich
      integer iComport          ! Determine le comportement si le fichier n'existe pas
                                ! iComport=0  => on continue en faisant semblant de rien et que la vie est belle
                                ! iComport=1  => on rale quand meme un peu
                                ! iComport>=2 => on boude definitivement en plaquant tout
                                ! iComport<0  => on ne dit rien, mais n'en pensant pas moins, on modifie iComport,
                                !                le rendant positif pour signaler le probleme
c
c var locales :
c ^^^^^^^^^^^^^
      logical lExiste
      character*160 NomFich
      common /common_NomFich_EvolAmas/ NomFich
c
c fonctions :
c ^^^^^^^^^^^
      integer LongChaine
c......................................................................
c
c---- former le nom du fichier
c
      write(NomFich,110)NomSimul(1:LongChaine(NomSimul)),iPas,TypeFich
 110  format(A,'%',I10.10,'%',A,'.xdr')
c
c---- traitement si le fichier n'existe pas
c
      inquire(file=NomFich,exist=lExiste)

      if (.NOT.lExiste) then
         if (iComport.LT.0) then
            iComport = -iComport
         else if (iComport.GT.0) then
            write(0,*) '!!! Le fichier a effacer ',
     $           NomFich(1:LongChaine(NomFich)),' n''existe pas !!!'
            if (iComport.GT.1) call Meurt('SupprimerFichXDR')
         end if
      else
c
c---- Effacement si le fichier existe
c
         open(iFich_Prov,file=NomFich,status='old')
         close(iFich_Prov,status='delete')
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine SupprimerFichASC(iPas,TypeFich,iComport)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer*8 iPas            ! numero du pas correspondant au fichier
                                ! a effacer
      character*(*) TypeFich
      integer iComport          ! Determine le comportement si le fichier n'existe pas
                                ! iComport=0  => on continue en faisant semblant de rien et que la vie est belle
                                ! iComport=1  => on rale quand meme un peu
                                ! iComport>=2 => on boude definitivement en plaquant tout
                                ! iComport<0  => on ne dit rien, mais n'en pensant pas moins, on modifie iComport,
                                !                le rendant positif pour signaler le probleme
c
c var locales :
c ^^^^^^^^^^^^^
      logical lExiste
      character*160 NomFich
      common /common_NomFich_EvolAmas/ NomFich
c
c fonctions :
c ^^^^^^^^^^^
      integer LongChaine
c......................................................................
c
c---- former le nom du fichier
c
      write(NomFich,110)NomSimul(1:LongChaine(NomSimul)),iPas,TypeFich
 110  format(A,'%',I10.10,'%',A,'.asc')
c
c---- traitement si le fichier n'existe pas
c
      inquire(file=NomFich,exist=lExiste)

      if (.NOT.lExiste) then
         if (iComport.LT.0) then
            iComport = -iComport
         else if (iComport.GT.0) then
            write(0,*) '!!! Le fichier a effacer ',
     $           NomFich(1:LongChaine(NomFich)),' n''existe pas !!!'
            if (iComport.GT.1) call Meurt('SupprimerFichASC')
         end if
      else
c
c---- Effacement si le fichier existe
c
         open(iFich_Prov,file=NomFich,status='old')
         close(iFich_Prov,status='delete')
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_Info_Phys(NumFich)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_Glob_common.f'
c
c argument :
c ^^^^^^^^^^
      integer NumFich
c
c fonctions :
c ^^^^^^^^^^^
      double precision R_tronc_mar, stellar_age
#if (_FORCER_CROISSANCE_TN_>0)
      double precision M_TN_Croiss_Forc
#endif
#if (_ACCR_BY_STARS_>0)
      double precision dMaccrdt_clust_MsunMyr
#endif
c......................................................................
c
#if (_FORCER_CROISSANCE_TN_>0)
      M_TN = M_TN_Croiss_Forc(Tps_amas)
      M_ctr = M_TN
#endif
      write (NumFich,100)
     $     iPas_Evol,
     $     Tps_amas, Tps_amas/UT_FP, Tps_amas/UneAnnee,
     $     dTps_amas_inf,dTps_amas_sup,
     $     Tps_amas_moy,Sigma_Tps_amas,
     $     NbSE-NbSE_perte,
     $     (NbSE-NbSE_perte)*Net_par_SE,
     $     M_amas,
     $     M_amas/Coef_ConvM,
     $     M_TN, M_ctr-M_TN,
     $     iType_CtrObj, stellar_age(id_CtrObj),
     $     dMejec_evap, dMejec_coll, dMejec_evst, dMejec_dechir,
     $     dMaccr_coll, dMaccr_evst, dMaccr_dechir, dMaccr_disp,
     $     dMaccr_captGW, 
     $     Etot_amas,  Ecin_amas,
     $     Estell_grav_amas,
     $     dEtot_evap, dEtot_coll, dEtot_evst,
     $     dEtot_dechir,dEtot_disp,dEtot_captGW,
     $     dEtot_bin3c,dEtot_nk,
     $     R_tronc_mar(),
     $     Pot_ctr,
     $     Rho_ctr,
     $     Sigma3D_ctr
#if (_ACCR_BY_STARS_>0)
      write (NumFich,110)
     $     tAccrByStar_1Msun_Myr,
     $     dMaccrdt_clust_MsunMyr(1.0d-6*Tps_Amas/UneAnnee)
#endif
      write (NumFich,150)
     $     iNb_Relax, iNb_Coll, iNb_Kick, iNb_Dechir, iNb_Disp,
     $     iNb_CaptGW
 100  format(
     $     '|======================================================|',/,
     $     '|=                 DONNEES PHYSIQUES                  =|',/,
     $     '|======================================================|',/,
     $     '|  iPas_Evol    : ',I10,/,
     $     '|  Tps_Ucode    : ',G20.10,' Tps_Ufokpl : ',E14.7,
     $     ' Tps_en_yr : ',E11.4,/,
     $     '|  dTps_inf     : ',E14.7,' dTps_sup    : ',E14.7,/,
     $     '|  Tps_moy      : ',G20.10,' Sigma_Tps : ',E14.7,/,
     $     '|  NbSE_subsist : ',I9,/,
     $     '|  Net_subsist  : ',E14.7,/,
     $     '|  M_amas       : ',E11.4,/,
     $     '|  M_amas_en_Msol : ',E11.4,/,
     $     '|  M_TN         : ',E11.4,' M_rsrv      : ',E11.4,/,
     $     '|  iType_CtrObj : ',I2,' Age_CtrObj_yr : ',E11.4,/,
     $     '|  dMejec_evap  : ',E11.4,' dMejec_coll : ',E11.4,
     $     ' dMejec_evst : ',E11.4,' dMejec_dechir : ',E11.4,/,
     $     '|                             dMaccr_coll : ',E11.4,
     $     ' dMaccr_evst : ',E11.4,' dMaccr_dechir : ',E11.4,
     $     ' dMaccr_disp : ',E11.4,' dMaccr_captGW : ',E11.4,/,
     $     '|  Etot_amas    : ',E14.7,' Ecin_amas   : ',E14.7,/,
     $     '|  Estell_grav  : ',E14.7,/,
     $     '|  dEtot_evap   : ',E11.4,' dEtot_coll  : ',E11.4,
     $     ' dEtot_evst  : ',E11.4,
     $     ' dEtot_dechir  : ',E11.4,' dEtot_disp  : ',E11.4,
     $     '  dEtot_captGW : ',E11.4,' dEtot_bin3c : ',E11.4,
     $     ' dEtot_nk  : ',E11.4,/,
     $     '|  R_maree      : ',E14.7,/,
     $     '|  Pot_ctr      : ',E14.7,/,
     $     '|  Rho_ctr      : ',E14.7,/,
     $     '|  Sigma3D_ctr  : ',E14.7
     $     )
 110  format(
     $     '|  tAccrByStar_1Msun_Myr  : ',E12.5,/,
     $     '|  dMaccrdt_clust_MsunMyr : ',E12.5
     $     )
 150  format(
     $     '|  Nb_relax     : ',I10,' Nb_coll : ',I10,' Nb_kick : ',I10,
     $     ' Nb_dechir : ',I8,' Nb_disp : ',I6,' Nb_captGW : ',I6
     $     )
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_Info_Tech(NumFich)
c----------------------------------------------------------------------
      include 'EvolAmas_inc.f'
      include 'EvolAmas_GlobTech_common.f'
#if (_COLLISIONS_>0 && _TYPE_COLL_>=10)
      include 'EvolAmas_Collisions_SPH_inc.f'
#endif
c
c argument :
c ^^^^^^^^^^
      integer NumFich
c......................................................................
c      
      write (NumFich,100)
     $     iPas_Evol,
     $     jour,mois,annee,
     $     heure,minute,sec,
     $     Tps_CPU,
     $     Viriel,
     $     Erreur_rel_E,
     $     Erreur_rel_M,
     $     Nb_moy_VN_pos_orb_par_pas,
     $     Nb_moy_parcrs_arbre_par_pas,
     $     Nb_moy_test_LC_par_pas,
     $     iNb_RejetPaire,
     $     iNb_AnnulRelaxEvap,
     $     iNb_AnnulRelaxCC
 100  format(
     $     '|======================================================|',/,
     $     '|=                 DONNEES TECHNIQUES                 =|',/,
     $     '|======================================================|',/,
     $     '|  iPasEvol : ',I10,/,
     $     '|  Date  : ',I2.2,'.',I2.2,'.',I4.4,/,
     $     '|  Heure : ',I2.2,'h',I2.2,'m',I2.2,'s',/,
     $     '|  Tcpu  : ',G10.4,/,
     $     '|--- Conservations :',/,
     $     '|  Viriel       : ',E9.2,/,
     $     '|  Erreur_rel_E : ',E9.2,/,
     $     '|  Erreur_rel_M : ',E9.2,/,
     $     '|--- Efficacite du code :',/,
     $     '|  Nb_moy_VN_pos_orb_par_pas   : ',F6.1,/,  
     $     '|  Nb_moy_parcrs_arbre_par_pas : ',F6.1,/,
     $     '|  Nb_moy_test_LC_par_pas      : ',F6.3,/,
     $     '|--- Nb de cas exceptionnels depuis dernier affichage :',/,
     $     '|  Nb_RejetPaire      : ',I5,/,
     $     '|  Nb_AnnulRelaxEvap  : ',I5,/,
     $     '|  Nb_AnnulRelaxCC    : ',I5
     $     )
c
c     Monitoring des super-rencontres pour simuler la relaxation
c
      write (NumFich,105)
     $     Frac_Pas_relax,
     $     Dev_SuperRenc_moy,
     $     Dev2_SuperRenc_moy,
     $     Frac_SuperRenc_Pi16,
     $     Frac_SuperRenc_Pi8,
     $     Frac_SuperRenc_Pi4,
     $     Frac_SuperRenc_Pi2
 105  format(
     $     '|--- Frac des pas pour relax :',/,
     $     '|  Frac_Pas_Relax          : ',E10.2,/,
     $     '|--- Deviation moy/quad moy lors des super-rencontres :',/,
     $     '|  Dev_SuperRenc_moy       : ',E10.2,/,
     $     '|  Dev2_SuperRenc_moy      : ',E10.2,/,
     $     '|--- Fraction de super-rencontres avec deviation > ',
     $     'Pi/16, /8, /4, /2 :',/,
     $     '|  Frac_SuperRenc_Pi16     : ',F6.4,/,
     $     '|  Frac_SuperRenc_Pi8      : ',F6.4,/,
     $     '|  Frac_SuperRenc_Pi4      : ',F6.4,/,
     $     '|  Frac_SuperRenc_Pi2      : ',F6.4
     $     )
c
c     Monitoring des cas a grande probabilite de collisions
c
#if (_COLLISIONS_>0)
      write (NumFich,106)
     $     P_coll_moy, P_coll_max,
     $     Frac_Coll_P_025,
     $     Frac_Coll_P_05,
     $     Frac_Coll_P_1,
     $     Frac_Coll_P_2
 106  format(
     $     '|--- Moyenne/max des proba de collision :',/,
     $     '|  P_coll_moy               : ',E10.2,/,
     $     '|  P_coll_max               : ',E10.2,/,
     $     '|--- Fraction des rencontres avec Pcoll > ',
     $     '0.25, 0.5, 1, 2 :',/,
     $     '|  Frac_Coll_P_025         : ',F6.4,/,
     $     '|  Frac_Coll_P_05          : ',F6.4,/,
     $     '|  Frac_Coll_P_1           : ',F6.4,/,
     $     '|  Frac_Coll_P_2           : ',F6.4
     $     )
#endif
c
c     Monitoring des cas a grande probabilite de diffusion grand angle
c
#if (_KICKS_>0)
      write (NumFich,107)
     $     P_kick_moy, P_kick_max,
     $     Frac_Kick_P_025,
     $     Frac_Kick_P_05,
     $     Frac_Kick_P_1,
     $     Frac_Kick_P_2
 107  format(
     $     '|--- Moyenne/max des proba de dif. gr. angles (kicks) :',/,
     $     '|  P_kick_moy               : ',E10.2,/,
     $     '|  P_kick_max               : ',E10.2,/,
     $     '|--- Fraction des rencontres avec Pkick > ',
     $     '0.25, 0.5, 1, 2 :',/,
     $     '|  Frac_Kick_P_025         : ',F6.4,/,
     $     '|  Frac_Kick_P_05          : ',F6.4,/,
     $     '|  Frac_Kick_P_1           : ',F6.4,/,
     $     '|  Frac_Kick_P_2           : ',F6.4
     $     )
#endif
c
c     Monitoring des rencontres "non-locales"
c
      write (NumFich,109)
     $     dR_rel_Renc_moy,
     $     dR_rel_Renc_max,
     $     Frac_Renc_dR_001,
     $     Frac_Renc_dR_01,
     $     Frac_Renc_dR_05,
     $     Frac_Renc_dR_1
 109  format(
     $     '|--- Moyenne/max des dR/Rmin des rencontes :',/,
     $     '|  dR_rel_Renc_moy         : ',E10.2,/,
     $     '|  dR_rel_Renc_max         : ',E10.2,/,
     $     '|--- Fraction des rencontres avec dR/Rmin > ',
     $     '0.01, 0.1, 0.5, 1 :',/,
     $     '|  Frac_Renc_dR_001        : ',F7.5,/,
     $     '|  Frac_Renc_dR_01         : ',F7.5,/,
     $     '|  Frac_Renc_dR_05         : ',F7.5,/,
     $     '|  Frac_Renc_dR_1          : ',F7.5
     $     )
      
#if (_COLLISIONS_>0 && _TYPE_COLL_>=10)  /* collisions selon resultats de MES simultaions SPH */
c
c     compteurs de cas speciaux lors de l'emploi des resultats SPH
c     pour les collisions
c
      write (NumFich,110) n_dbrdM_CollSPH_b, n_dbrdM_CollSPH_h,
     $     n_extrap_CollSPH, n_distant_CollSPH, n_inconnu_CollSPH
 110  format(
     $     '|  n_debordM_CollSPH_bas       : ',I7,/,
     $     '|  n_debordM_CollSPH_haut      : ',I7,/,
     $     '|  n_extrap_CollSPH            : ',I7,/,
     $     '|  n_distant_CollSPH           : ',I7,/,
     $     '|  n_inconnu_CollSPH           : ',I7
     $     )
#endif
c
      end
c
c----------------------------------------------------------------------
      subroutine Sauvegarde_EvolAmas(lEcrase,lConserver_Sauv)
c----------------------------------------------------------------------
c     sauvegarde generale de l'etat du systeme
c
c     Toutes les donnees essentielles sont sauvegardees
c     au format XDR afin de pouvoir etre traitees sur n'importe quel
c     ordinateur, independemment de la plateforme sur laquelle la
c     simulation a ete faite.
c
      include 'EvolAmas_inc.f'
      include 'fxdr.inc'
c
c argument(s) :
c ^^^^^^^^^^^^^
      logical lEcrase           ! faut-il ecraser les fichiers deja existants ?
      logical lConserver_Sauv   ! faudra-t-il conserver cette sauvegarde (fichiers XDR volumineux) ?
c
c constantes :
c ^^^^^^^^^^^^
      integer iFich
      parameter (iFich=11)
#if (_ECRIT_GRILLE_>0)
      integer Npc_opt_def,Npc_min_def,Npc_max_def ! Nb de part/cellule dans la grille
      parameter (Npc_opt_def=50,Npc_min_def=10,Npc_max_def=1000)
      integer iDimGrille
      parameter (iDimGrille=iDimSE/Npc_min_def+1)
#endif
      integer iUsePercent_max
      parameter (iUsePercent_max=98)
c
c var locales :
c ^^^^^^^^^^^^^
      integer iFich_XDR, ixdr_err
      integer iAction
      logical lExiste
      integer annee,mois,jour,heure,minute,sec
#if (_ECRIT_PROFILPOT_>0)
      integer i,j
      double precision Rpp(iDimSE),Upp(iDimSE),dUpp(iDimSE)
#endif
#if (_ECRIT_GRILLE_>0)
      integer Npc_opt,Npc_min,Npc_max,N_Grille,NSE_Grille(iDimGrille)
      double precision
     $     Rsup_Grille(0:iDimGrille),M_Grille(iDimGrille),
     $     MR_Grille(iDimGrille),T_Grille(iDimGrille),
     $     Tt_Grille(iDimGrille)
#endif
      logical lConserver_Sauv_prec
      save lConserver_Sauv_prec
      data lConserver_Sauv_prec /.TRUE./
      integer*8 iPas_Sauv_Prec
      save iPas_Sauv_Prec
      data iPas_Sauv_Prec /0/
      integer df_MD,iOK,iUsePercent
      integer*8 iOneK_Blocks,iUsed,iAvailable
      character*256 Filesystem,Mounted_on
c......................................................................
c
c---- One tests wether there is some disk space available
c     If not, it's better to stop the simulation.
c
#if (_CONDOR_==0 && _NO_DF_==0)
      iOK = df_MD('.',Filesystem,iOneK_Blocks,iUsed,iAvailable,
     $     iUsePercent,Mounted_on)
#else
      iOK = 0
      iUsePercent=-1
#endif
      if (iOK.NE.0) then
         write(0,*) '!!! Error with df_MD !!!'
         call Meurt('Sauvegarde_EvolAmas')
      end if
      if (iUsePercent.GE.iUsePercent_max) then
         write(0,*) '!!! Disk nearly full (file system ',
     $        Filesystem(1:index(FileSystem,' ')-1),') !!!'
         call exit(1)
      end if
      
      call QuelleDate_MD(jour,mois,annee)
      call QuelleHeure_MD(heure,minute,sec)

      write(0,120)iSauv,iPas_Evol,
     $     jour,mois,annee,heure,minute,sec
 120  format(/,'iSauv = ',I4,' iPas_Evol =',I10,/,
     $     '< Date : ',I2.2,'.',I2.2,'.',I4.4,
     $     ' Heure : ',I2.2,':',I2.2,':',I2.2,' >',/,
     $     '  > Ecriture des donnees ...')
      if (lEcrase) then
         iAction = 3
      else
         iAction = 1
      end if
c
c     -------------------------------------------------------
c==== Donnees physiques indispensables (decrivent le systeme)
c     -------------------------------------------------------
c
c---- Fichier XDR d'amas
c
      call OuvrirEcritureFichXDR(iFich_XDR,TypeAMAS,lExiste,lEcrase)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)' for star cluster'
         call EcritAmas_XDR(iFich_XDR)
         write(0,110)
         ixdr_err = ixdrclose(iFich_XDR)
      end if
      if (.NOT.lConserver_Sauv_Prec) call
     $     SupprimerFichXDR(iPas_Sauv_Prec,TypeAMAS,1)
c
c---- Fichier XDR des temps et pas de temps
c
      call OuvrirEcritureFichXDR(iFich_XDR,TypeTEMPS,lExiste,lEcrase)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for times'
         call EcritTemps_XDR(iFich_XDR)
         write(0,110)
         ixdr_err = ixdrclose(iFich_XDR)
      end if
      if (.NOT.lConserver_Sauv_Prec) call
     $     SupprimerFichXDR(iPas_Sauv_Prec,TypeTEMPS,1)
c
c---- ASCII file for cumulative collision probabilites
c
#if (_WRITE_CUMUL_COLL_PROB_>0)
      call OuvrirFichASC(iFich_Prov,TypeCumPcoll,iAction,lExiste)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for cumulative collision probabilites'
         call Write_CumulPcoll_ASC(iFich_Prov)
         write(0,110)
         call flush_MD(iFich_Prov)
         close(iFich_Prov)
         if (.NOT.lConserver_Sauv_Prec) call
     $        SupprimerFichASC(iPas_Sauv_Prec,TypeCumPcoll,1)
      end if
#endif
c
c---- Fichier XDR de prop. stellaires
c
#if (_COLLISIONS_>0 || _ECRIT_ETOILES_>0)
      call OuvrirEcritureFichXDR(iFich_XDR,TypeEtoiles,lExiste,lEcrase)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for stellar data'
         call EcritEtoiles_XDR(iFich_XDR)
         write(0,110)
         ixdr_err = ixdrclose(iFich_XDR)
      if (.NOT.lConserver_Sauv_Prec) call
     $     SupprimerFichXDR(iPas_Sauv_Prec,TypeEtoiles,1)
      end if
#endif
c
c---- Fichier XDR de donnees sur le trou noir
c
#if (_PRESENCE_TN_>0)
      call OuvrirEcritureFichXDR(iFich_XDR,TypeTN,lExiste,lEcrase)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for central massive object'
         call EcritTN_XDR(iFich_XDR)
         write(0,110)
         ixdr_err = ixdrclose(iFich_XDR)
      end if
      if (.NOT.lConserver_Sauv_Prec) call
     $     SupprimerFichXDR(iPas_Sauv_Prec,TypeTN,1)
#endif
c
c---- Fichier XDR de donnees sur le chp de maree
c
#if (_TRONC_MAREE_>0)
      call OuvrirEcritureFichXDR(iFich_XDR,TypeMAREE,lExiste,lEcrase)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for tidal truncation'
         call EcritMaree_XDR(iFich_XDR)
         write(0,110)
         ixdr_err = ixdrclose(iFich_XDR)
      end if
      if (.NOT.lConserver_Sauv_Prec) call
     $     SupprimerFichXDR(iPas_Sauv_Prec,TypeMAREE,1)
#endif
c
c---- Fichier XDR de donnees sur la limite exterieure (troncation)
c
#if (_LIMITE_EXT_>0)
      call OuvrirEcritureFichXDR(iFich_XDR,TypeLIMEXT,lExiste,lEcrase)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for outer wall'
         call EcritLimExt_XDR(iFich_XDR)
         write(0,110)
         ixdr_err = ixdrclose(iFich_XDR)
      end if
      if (.NOT.lConserver_Sauv_Prec) call
     $     SupprimerFichXDR(iPas_Sauv_Prec,TypeLIMEXT,1)
#endif
c
c---- Fichier XDR de donnees le chauffage par binaires "3-corps"
c
#if (_CHAUFF_BIN_3C_>0 && _ECRIT_CB3C_>0 )
      call OuvrirEcritureFichXDR(iFich_XDR,TypeCB3CE,lExiste,lEcrase)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for heating by 3-body binaries'
         call EcritCB3c_XDR(iFich_XDR)
         write(0,110)
         ixdr_err = ixdrclose(iFich_XDR)
      end if
      if (.NOT.lConserver_Sauv_Prec) call
     $     SupprimerFichXDR(iPas_Sauv_Prec,TypeCB3CE,1)
#endif
c
c---- Fichier XDR de donnees de suivi de la conservation
c     masse/energie
c
      call OuvrirEcritureFichXDR(iFich_XDR,TypeCONS,lExiste,lEcrase)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for conserved quantities'
         call EcritCons_XDR(iFich_XDR)
         write(0,110)
         ixdr_err = ixdrclose(iFich_XDR)
      end if
      if (.NOT.lConserver_Sauv_Prec) call
     $     SupprimerFichXDR(iPas_Sauv_Prec,TypeCONS,1)
c
c     ---------------------------
c==== Autres donnees facultatives
c     ---------------------------
c
c---- Fichier d'etat du generateur de nb aleatoires
c
#if (_ECRIT_RAND_>0)
      call OuvrirFichASC(iFich_Prov,TypeRAND,iAction,lExiste)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for state of random generator'
         call save_random(iFich_Prov)
         call flush_MD(iFich_Prov)
         close(iFich_Prov)
         write(0,110)
      end if
#endif
c
c---- Fichier xdr de grille PG
c
#if (_ECRIT_PG_>0)
      call OuvrirEcritureFichXDR(iFich_XDR,TypePG,lExiste,lEcrase)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for PG grid'
         call EcritGrille_XDR(iFich_XDR)
         write(0,110)
         ixdr_err = ixdrclose(iFich_XDR)
      end if
#endif
c
c---- Profil du potentiel (debuggage)
c
#if (_ECRIT_PROFILPOT_>0)
      call OuvrirFichBIN(iFich_Prov,'ProfPot',iAction,lExiste)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for potential profile'
         j=0
         do i=1,NbSE
            if (R_SE(i).LT.Revap) then
               j=j+1
               Rpp(j) = R_SE(i)
               call Pot_Arbre(Rpp(j),Upp(j),dUpp(j))
            end if
         end do
         write(iFich)j,(Rpp(i),i=1,j),(Upp(i),i=1,j),(dUpp(i),i=1,j)
         call flush_MD(iFich_Prov)
         close(iFich_Prov)
         write(0,110)
      end if
#endif
c
c---- Fichier XDR de l'arbre
c
#if (_ECRIT_ARBRE_>0)
      call OuvrirEcritureFichXDR(iFich_XDR,TypeARBRE,lExiste,lEcrase)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for binary tree'
         call EcritArbre_XDR(iFich_XDR)
         call flush_MD(iFich)
         write(0,110)
         ixdr_err = ixdrclose(iFich_XDR)
      end if
#endif
c
c---- Parametres determinant la loi de tirage des SE
c
#if (_ECRIT_TH_>0)
      call OuvrirFichBIN(iFich_Prov,TypeTH,iAction,lExiste)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'of parameters for particle selection (TH)'
         call Ecrit_TH_bin(iFich_Prov)
         call flush_MD(iFich_Prov)
         close(iFich_Prov)
         write(0,110)
      end if
#endif
c
 100  format('/ Writing file ',A,$)
 110  format(' /')
c
c---- On se souvient du numero de la sauvegarde et
c     de la necessite (ou non) de la conserver
c
      iPas_Sauv_Prec = iPas_evol
      lConserver_Sauv_Prec = lConserver_Sauv
c
      end
c
c----------------------------------------------------------------------
      subroutine SauvegardePart_EvolAmas(lEcrase)
c----------------------------------------------------------------------
c     sauvegarde partielle
c
      include 'EvolAmas_inc.f'
c
c argument(s) :
c ^^^^^^^^^^^^^
      logical lEcrase           ! faut-il ecraser les fichiers deja existants ?
c
c constantes :
c ^^^^^^^^^^^^
      integer iFich
      parameter (iFich=11)
#if (_ECRIT_GRILLE_>0)
      integer Npc_opt_def,Npc_min_def,Npc_max_def ! Nb de part/cellule dans la grille
      parameter (Npc_opt_def=50,Npc_min_def=10,Npc_max_def=1000)
      integer iDimGrille
      parameter (iDimGrille=iDimSE/Npc_min_def+1)
#endif
      integer iUsePercent_max
      parameter (iUsePercent_max=99)
c
c var locales :
c ^^^^^^^^^^^^^
      integer annee,mois,jour,heure,minute,sec
#if (_ECRIT_GRILLE_>0)
      integer Npc_opt,Npc_min,Npc_max,N_Grille,NSE_Grille(iDimGrille)
      double precision
     $     Rsup_Grille(0:iDimGrille),M_Grille(iDimGrille),
     $     MR_Grille(iDimGrille),T_Grille(iDimGrille),
     $     Tt_Grille(iDimGrille)
#endif
      integer df_MD,iOK,iUsePercent
      integer*8 iOneK_Blocks,iUsed,iAvailable
      character*256 Filesystem,Mounted_on
c......................................................................
c
c---- One tests wether there is some disk space available
c     If not, it's better to stop the simulation.
c
#if (_CONDOR_==0 && _NO_DF_==0)
      iOK = df_MD('.',Filesystem,iOneK_Blocks,iUsed,iAvailable,
     $     iUsePercent,Mounted_on)
#else
      iOK = 0
      iUsePercent=-1
#endif
      if (iOK.NE.0) then
         write(0,*) '!!! Error with df_MD !!!'
         call Meurt('SauvegardePart_EvolAmas')
      end if
      if (iUsePercent.GE.iUsePercent_max) then
         write(0,*) '!!! Disk nearly full (file system ',
     $        Filesystem(1:index(FileSystem,' ')-1),') !!!'
         call exit(1)
      end if
      
      call QuelleDate_MD(jour,mois,annee)
      call QuelleHeure_MD(heure,minute,sec)

      write(0,120)iSauvPart,iPas_Evol,
     $     jour,mois,annee,heure,minute,sec
 120  format(/,'iSauvPart = ',I4,' iPas_Evol =',I10,/,
     $     '< Date : ',I2.2,'.',I2.2,'.',I4.4,
     $     ' Heure : ',I2.2,':',I2.2,':',I2.2,' >',/,
     $     '  > Ecriture des donnees ...')

      call SauvegardePart_(lEcrase)
c
 100  format('/ Ecriture fichier ',A,$)
 110  format(' /')
c
      end
c
c----------------------------------------------------------------------
      subroutine SauvegardePart_(lEcrase)
c----------------------------------------------------------------------
c     sauvegarde partielle
c
      include 'EvolAmas_inc.f'
c
c argument(s) :
c ^^^^^^^^^^^^^
      logical lEcrase 
c
c constantes :
c ^^^^^^^^^^^^
      integer iFich
      parameter (iFich=11)
#if (_ECRIT_GRILLE_>0)
      integer Npc_opt_def,Npc_min_def,Npc_max_def ! Nb de part/cellule dans la grille
      parameter (Npc_opt_def=50,Npc_min_def=10,Npc_max_def=1000)
      integer iDimGrille
      parameter (iDimGrille=iDimSE/Npc_min_def+1)
#endif
c
c var locales :
c ^^^^^^^^^^^^^
      integer iFich_XDR, ixdr_err, iAction
      logical lExiste
#if (_ECRIT_GRILLE_>0)
      integer Npc_opt,Npc_min,Npc_max,N_Grille,NSE_Grille(iDimGrille)
      double precision
     $     Rsup_Grille(0:iDimGrille),M_Grille(iDimGrille),
     $     MR_Grille(iDimGrille),T_Grille(iDimGrille),
     $     Tt_Grille(iDimGrille)
#endif
c......................................................................
c
      if (lEcrase) then
         iAction = 3
      else
         iAction = 1
      end if
c
c---- Fichier xdr de grille "optimisee" (pour representation graphique)
c
#if (_ECRIT_GRILLE_>0)
      Npc_opt=Npc_opt_def
      Npc_min=Npc_min_def
      Npc_max=Npc_max_def
      call Construire_Grille (Npc_opt,Npc_min,Npc_max,
     $     iDimGrille,
     $     Rsup_Grille,M_Grille,MR_Grille,T_Grille,Tt_Grille,NSE_Grille,
     $     N_Grille)
      call OuvrirEcritureFichXDR(iFich_XDR,TypeGRILLE,lExiste,lEcrase)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for radial grids'
         call EcritGrille_XDR2(iFich_XDR,
     $        N_Grille,
     $        Rsup_Grille,M_Grille,MR_Grille,T_Grille,Tt_Grille,
     $        NSE_Grille)
         write(0,110)
         ixdr_err = ixdrclose(iFich_XDR)
      end if
#endif
c
c
c---- "On-the-fly" grid: attempt to have continuously updated "oversampled" grid for various profiles
c
#if (_WRITE_OFGRID_>0)
      call OuvrirFichASC(iFich,TypeOFGrid,iAction,lExiste)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'for ''on-the-fly'' radial grids'
         call write_OFGrid(iFich)
         call flush_MD(iFich)
         close(iFich)
         call Reinit_OFGrid
         write(0,110)
      end if   
#endif
c
c     --------------------------
c==== fichiers ASCII sequentiels
c     --------------------------
c
c---- Fichier ascii des rayons lagrangiens
c
#if (_ECRIT_RAYLAG_>0)
      write(0,100)'for Lagrange radii'
      if (flag_Tri.NE.1) then
         call TriAmas
         flag_Tri = 1
      end if
      call Ecrire_RayLag_ASC(iFich_RayLag)
      call flush_MD(iFich_RayLag)
      write(0,110)
#endif
c
c---- Lagrange radii and various quantities summed inside Lag spheres
c
#if (_ECRIT_LAGQUANT_>0)
      write(0,100)'for Lagrange quantities'
      if (flag_Tri.NE.1) then
         call TriAmas
         flag_Tri = 1
      end if
      call Ecrire_LagQuant_ASC(iFich_LagQuant)
      call flush_MD(iFich_LagQuant)
      write(0,110)
#endif
c
c---- Fichier ascii des anisotropies lagrangiens
c
#if (_ECRIT_ANILAG_>0)
      write(0,100)'for anisotropy'
      if (flag_Tri.NE.1) then
         call TriAmas
         flag_Tri = 1
      end if
      call Ecrire_AniLag_ASC(iFich_AniLag)
      call flush_MD(iFich_AniLag)
      write(0,110)
#endif
c
c---- Fichier ascii du spectre de masse
c
#if (_ECRIT_SPECMASSE_>0)
      write(0,100)'for mass spectrum'
      call Write_MassSpec_ASC(iFich_MSpec)
      call flush_MD(iFich_MSpec)
      write(0,110)
#endif
c
c---- Fichier ascii de segregation
c
#if (_ECRIT_SEGR_>0)
      write(0,100)'to monitor mass segregation'
      if (flag_Tri.NE.1) then
         call TriAmas
         flag_Tri = 1
      end if
      call Ecrire_SuiviSegr_ASC(iFich_Segr)
      call flush_MD(iFich_Segr)
      write(0,110)
#endif
c
c---- Fichier ascii des grandeurs globales de relaxation
c
#if (_ECRIT_RELAX_>0)
      write(0,100)'for global realaxtion quantities'
      call Calc_RelaxGlob
      call Ecrire_Relax_ASC(iFich_Relax)
      call flush_MD(iFich_Relax)
      write(0,110)
#endif
c
c---- Fichier de suivi des differents types stellaires
c
#if (_ECRIT_SUIVI_TYPES_STELL_>0)
      write(0,100)'to follow stellar types'
      call SuiviTypesStell
      write(0,110)
#endif
c
c---- Fichier de suivi de la distrib en R par bin de masse
c
#if (_MONITOR_MBINS_>0)
      write(0,100)'to Lagrange radii per mass bins'
      call Write_MbinsMonitor_ASC(iFich_MbinRad)
      call flush_MD(iFich_MbinRad)
      write(0,110)
#endif
c
c---- File to follow a (small) sub-population of particles
c
#if (_FOLLOW_SUBPOP_>0)
      call OuvrirFichASC(iFich,TypeSubPop,iAction,lExiste)
      if (.NOT.lExiste .OR. lEcrase) then
         write(0,100)'to follow subpopulation'
         call write_SubPop(iFich)
         call flush_MD(iFich)
         close(iFich)
         write(0,110)
      end if
#endif
 100  format('/ Writing file ',A,$)
 110  format(' /')
c
      end
c
c-------------------------------------------------------------------------------
      subroutine open_append(iFile,Name,lOK)
c-------------------------------------------------------------------------------
      implicit none
c
c arguments:
c """"""""""
      integer iFile
      character*(*) Name
      logical lOK
c
c local vars :
c """"""""""""
      logical lExist
      integer iOK
c......................................................................
c
      lOK = .TRUE.
      inquire(file=Name,exist=lExist)
      if (lExist) then
         open(iFile,file=Name,IOSTAT=iOK,form='formatted',
     $        access='append')
      else
         open(iFile,file=Name,IOSTAT=iOK,form='formatted',
     $        STATUS='new')
      end if
      if (iOK.NE.0) then
         close(iFile)
         lOK = .FALSE.
      end if
c
      end
c
c-------------------------------------------------------------------------------
      subroutine open_append_BIN(iFile,Name,lOK)
c-------------------------------------------------------------------------------
      implicit none
c
c arguments:
c """"""""""
      integer iFile
      character*(*) Name
      logical lOK
c
c local vars :
c """"""""""""
      logical lExist
      integer iOK
c......................................................................
c
      lOK = .TRUE.
      inquire(file=Name,exist=lExist)
      if (lExist) then
         open(iFile,file=Name,IOSTAT=iOK,form='unformatted',
     $        access='append')
      else
         open(iFile,file=Name,IOSTAT=iOK,form='unformatted',
     $        STATUS='new')
      end if
      if (iOK.NE.0) then
         close(iFile)
         lOK = .FALSE.
      end if
c
      end
c$$$
c$$$c-------------------------------------------------------------------------------
c$$$      subroutine open_append(iFile,Name,lOK)
c$$$c-------------------------------------------------------------------------------
c$$$      implicit none
c$$$c
c$$$c arguments:
c$$$c """"""""""
c$$$      integer iFile
c$$$      character*(*) Name
c$$$      logical lOK
c$$$c
c$$$c local vars :
c$$$c """"""""""""
c$$$      logical lExist
c$$$      character one_car
c$$$      integer iOK
c$$$c......................................................................
c$$$c
c$$$      lOK = .TRUE.
c$$$      inquire(file=Name,exist=lExist)
c$$$      open(iFile,file=Name,IOSTAT=iOK)
c$$$      if (iOK.NE.0) goto 999
c$$$      if (lExist) then
c$$$         do while (.TRUE.)
c$$$            read(ifile,'(A)',END=9,ERR=999) one_car
c$$$         end do
c$$$      end if
c$$$ 9    continue
c$$$      return
c$$$ 999  continue
c$$$      close(iFile)
c$$$      lOK = .FALSE.
c$$$c
c$$$      end
