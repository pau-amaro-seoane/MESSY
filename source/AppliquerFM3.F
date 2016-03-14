      program AppliquerFM3
c====================================================================================================
c     Give stellar masses to a stellar cluster according to a piecewise power-law IMF.
c
c     usage is:
c
c     appliquerfm3 [-m MlimList] [-a ExponentList] [-age age_in_yrs] [-Tsf SF_time_scale_in_yrs] [-iseed random_seed] FichIn FichOut
c
c     format of MlimList is like (with quotes)  "0.01,0.08,0.5,200"
c     format of ExponentList is like            "0.3,1.3,2.3" (one value less than MlimList)
c
c     Default values (Heggie's collaborative experiment) :
c         limit masses "0.1,1.5"
c         exponent "2.35"
c
c     One can produce a stellar population of age age_in_yrs with switch -age
c     (a very simple stellar evolution model is used with direct MS --> compact remnant transition and no MS evolution)
c     In this case, one may adopt a parametrised SFR expressed by SF(t) ~ exp(-t/SF_time_scale_in_yrs)
c     with switch -Tsf.
c     
c     -iseed switch allows one to impose a given seed for the random-number generator.
c
c     FichIn and FichOut are files of stellar data 
c     of format "_ETOILES.xdr"
c====================================================================================================

c====================================================================================================
c     Applique une fct des masses a un amas d'etoiles
c
c     appliquerfm3 [-Tsf SF_time_scale_in_yrs] [-age age_in_yrs] [-m ListeMlim] [-a ListeExposants] FichIn FichOut
c
c     Applique une fonction des masses qui est une loi de puissance par morceaux
c     Le format de ListeMlim      est du type "0.01,0.08,0.5,200"
c     Le format de ListeExposants est du type "0.3,1.3,2.3" (une valeur de moins que ListeExposants
c
c     Valeurs par defaut (collaborative experiment de D.Heggie) :
c         ListeMlim 0.1,1.5
c         Exposant 2.35
c
c     On peut vieillir la population jusqu'a un age donne et utiliser une SF-history du type SF(t) ~ exp(-t/Tsf)
c     Par defaut, l'age est 0 et Tsf->0 (burst instantane)
c     
c     La masse totale est donnees en masses solaires; cette options prend le dessus sur
c     le nb total d'etoiles qui est ajuste en consequence
c
c     FichIn et FichOut sont des fichiers de donnees stellaires
c     du type "_ETOILES.xdr"
c====================================================================================================
      include 'AppliquerFM3_inc.f'
      include 'AppliquerFM3_Param_FM_inc.f'
      include 'EvolAmas_Maree_common.f'
      include 'fxdr.inc'
c
c constante(s) :
c ^^^^^^^^^^^^^^
      integer iIn
      parameter (iIn=10)
c
c var locales :
c ^^^^^^^^^^^^^
      character*80 NomFich,NomFich2
      integer iSE,i_arg,iXDR, ixdr_err
      double precision Ecin_0,Egrav_0,M_tot_,
     $     renR,renT,renJ,M_ZAMS,coef_SFH
      logical lEtoiles,lMaree
c
c fonction(s) :
c ^^^^^^^^^^^^^
      integer LongChaine
      double precision random
c......................................................................
c
c---- Lecture des parametres de la FM specifies sur la ligne de commande
c
      i_arg = 1                 ! numero du prochain parametre a lire
      call Lect_Param_FM(i_arg)

      call getarg(i_arg,NomFich)
      i_arg=i_arg+1

      write(0,*) ' >> Lecture de donnees au format XDR, fichier ',
     $     NomFich(1:LongChaine(NomFich))
      
      iXDR = initxdr(NomFich,'r',.TRUE.)
      call LitAmas_XDR(iXDR)
      ixdr_err = ixdrclose(iXDR)
      NomFich = NomFich(1:LongChaine(NomFich)-4)
      
      call ConcatChaines(NomFich,'_ETOILES.xdr',NomFich2)
      inquire(file=NomFich2,exist=lEtoiles)
      if (lEtoiles) then
         write(0,*) ' + Fichier "ETOILES" present'
         iXDR = initxdr(NomFich2,'r',.TRUE.)
         call LitEtoiles_XDR(iXDR)
         ixdr_err = ixdrclose(iXDR)
      else
         write(0,*) ' + Fichier "ETOILES" absent'
      end if
      
      call ConcatChaines(NomFich,'_MAREE.xdr',NomFich2)
      inquire(file=NomFich2,exist=lMaree)
      if (lMaree) then
         write(0,*) ' + Fichier "MAREE" present'
         iXDR = initxdr(NomFich2,'r',.TRUE.)
         call LitMaree_XDR(iXDR)
         ixdr_err = ixdrclose(iXDR)
      else
         write(0,*) ' + Fichier "MAREE" absent'
      end if
c
c---- determination des energies
c
      call TriAmas
      call Calc_EnerGlob
      Ecin_0  = Ecin_amas
      Egrav_0 = Egrav_amas
c
c---- tirage des masses des etoiles
c
      call Tirage_IMF_LPPM(N_morc, Mlim, alpha, NbSE, Met_SE)

      M_tot_ = 0.0d0
      if (Tsf.GT.0.0d0) then
         coef_SFH = 1.0d0-exp(-Age_Amas/Tsf)
      end if
      do iSE=1,NbSE
         M_ZAMS=Met_SE(iSE)
c
c---- On tire la date de naissance (sur la ZAMS) en fct de la SF-history (exponentielle)
c
         if (Tsf.GT.0.0d0) then
            DNet_SE(iSE) = -Tsf*log(1.0d0-random()*coef_SFH) ! date de naissance...
         else
            DNet_SE(iSE) = 0.0d0
         end if
c
c---- On "vieilli" l'etoile jusqu'a l'age specifie
c
         call Evoluer_Etoile(M_ZAMS,Age_Amas-DNet_SE(iSE),
     $        Met_SE(iSE),iTet_SE(iSE))
         M_tot_ = M_tot_+Met_SE(iSE)
      end do
c
c---- masse des SE pour obtenir la bonne masse totale
c
      do iSE=1,NbSE
         M_SE(iSE) = 1.0d0/M_tot_*Met_SE(iSE)
      end do
c
c---- informer sur les differents types stellaires
c
      call Statistiques
c
c---- on decale les dates de naissance pour mettre l'age a 0
c
      do iSE=1,NbSE
         DNet_SE(iSE) = DNet_SE(iSE)-Age_amas
      end do
c
c---- calcul energie et viriel
c
      call Calc_EnerGlob
      call Affiche_Ener
c
c---- renormalisation rayons/vitesses pour retrouver
c     les valeurs initiales des energies
c
      if (lRenorm) then
         renR = Egrav_amas/Egrav_0
         renT = Ecin_0/Ecin_amas
         renJ = renR*sqrt(RenT)
         write(0,120) renR,renT,renJ
 120     format(' > coefficient de renormalisation de R,T,J : ',
     $        3(F6.4,' '))
         do iSE=1,NbSE
            R_SE(iSE)=renR*R_SE(iSE)
            J_SE(iSE)=renJ*J_SE(iSE)
            T_SE(iSE)=renT*T_SE(iSE)
         end do
         if (lMaree) R_mar=min(1.0d20,renR*R_mar)
      end if
c
c---- calcul energie et viriel
c
      call Calc_EnerGlob
      call Affiche_Ener
c
c---- ecriture des fichiers
c
      call getarg(i_arg,NomFich)
      i_arg=i_arg+1

      write(0,*) ' >> Ecriture de donnees au format XDR, fichier ',
     $     NomFich(1:LongChaine(NomFich))
      
      iXDR = initxdr(NomFich,'w',.TRUE.)
      call EcritAmas_XDR(iXDR)
      ixdr_err = ixdrclose(iXDR)
      NomFich = NomFich(1:LongChaine(NomFich)-4)
      
      call ConcatChaines(NomFich,'_ETOILES.xdr',NomFich2)
      iXDR = initxdr(NomFich2,'w',.TRUE.)
      call EcritEtoiles_XDR(iXDR)
      ixdr_err = ixdrclose(iXDR)
      
      if (lMaree) then
         call ConcatChaines(NomFich,'_MAREE.xdr',NomFich2)
         iXDR = initxdr(NomFich2,'w',.TRUE.)
         call EcritMaree_XDR(iXDR)
         ixdr_err = ixdrclose(iXDR)
      end if
c     
      end
c
c----------------------------------------------------------------------
      subroutine Lect_Param_FM(i_arg)
c----------------------------------------------------------------------
      include 'AppliquerFM3_inc.f'
      include 'AppliquerFM3_Param_FM_inc.f'
      include 'EvolAmas_EvolStell_common.f'
c
c parametre :
c ^^^^^^^^^^^
      integer i_arg
c
c var locale :
c ^^^^^^^^^^^^
      character*64 Chaine, form_str
      integer l,i1,i2,i,iseed
      character*256 ListeMlim, ListeExposants
c......................................................................
c
      Age_Amas=0.0d0 
      Tsf=0.0d0
      ListeExposants='2.35'
      ListeMlim='0.1,1.5'
      iseed = 676583
      lRenorm = .TRUE.

      call InitEvolStell

      call getarg(i_arg,Chaine)

      do while (Chaine(1:1).EQ.'-')
        
         if (Chaine.EQ.'-m') then
            i_arg=i_arg+1
            call getarg(i_arg,ListeMlim)
         else if (Chaine.EQ.'-a') then
            i_arg=i_arg+1
            call getarg(i_arg,ListeExposants)
         else if (Chaine.EQ.'-age') then
            i_arg=i_arg+1
            call getarg(i_arg,Chaine)
            read(Chaine,*) Age_Amas
         else if (Chaine.EQ.'-Tsf') then
            i_arg=i_arg+1
            call getarg(i_arg,Chaine)
            read(Chaine,*) Tsf
         else if (Chaine.EQ.'-iseed') then
            i_arg = i_arg+1
            call getarg(i_arg,Chaine)
            read(Chaine,*) iseed
         else if (Chaine.EQ.'-Mzams_max_WD') then
            i_arg=i_arg+1
            call getarg(i_arg,Chaine)
            read(Chaine,*) Mzams_max_WD
         else if (Chaine.EQ.'-Mzams_max_NS') then
            i_arg=i_arg+1
            call getarg(i_arg,Chaine)
            read(Chaine,*) Mzams_max_NS
         else if (Chaine.EQ.'-M_max_WD') then
            i_arg=i_arg+1
            call getarg(i_arg,Chaine)
            read(Chaine,*) M_max_WD
         else if (Chaine.EQ.'-M_max_NS') then
            i_arg=i_arg+1
            call getarg(i_arg,Chaine)
            read(Chaine,*) M_max_NS
         else if (Chaine.EQ.'-M_WD') then
            i_arg=i_arg+1
            call getarg(i_arg,Chaine)
            read(Chaine,*) M_WD
         else if (Chaine.EQ.'-M_NS') then
            i_arg=i_arg+1
            call getarg(i_arg,Chaine)
            read(Chaine,*) M_NS
         else if (Chaine.EQ.'-M_BH') then
            i_arg=i_arg+1
            call getarg(i_arg,Chaine)
            read(Chaine,*) M_BH
         else if (Chaine.EQ.'-M_WD') then
            i_arg=i_arg+1
            call getarg(i_arg,Chaine)
            read(Chaine,*) M_WD
         else if (Chaine.EQ.'-Mremn_min_BH') then
            i_arg=i_arg+1
            call getarg(i_arg,Chaine)
            read(Chaine,*) Mremn_min_BH
         else if (Chaine.EQ.'-NoRenorm') then
            lRenorm = .FALSE.
         else
            write(0,110) Chaine
 110        format('!!! Option ',A8,' inconnue !!!')
            call exit(1)
         end if
         i_arg = i_arg+1
         call getarg(i_arg,Chaine)

      end do
c
c---- Determination des limites en masse et exposants
c
      N_morc = 0
      l=index(ListeExposants,' ')-1
      i2=1
      do while (.TRUE.)
         i1=i2
         do while (i2.LT.l .AND. ListeExposants(i2+1:i2+1).NE.',')
            i2=i2+1
         end do
         N_morc=N_morc+1
         read(ListeExposants(i1:i2),*) alpha(N_morc)
         i2=i2+2
         if (i2.GT.l) goto 9
      end do
 9    continue

      N_morc = -1
      l=index(ListeMlim,' ')-1
      i2=1
      do while (.TRUE.)
         i1=i2
         do while (i2.LT.l .AND. ListeMlim(i2+1:i2+1).NE.',')
            i2=i2+1
         end do
         N_morc=N_morc+1
         read(ListeMlim(i1:i2),*) Mlim(N_morc)
         i2=i2+2
         if (i2.GT.l) goto 19
      end do
 19   continue

      write(0,100)
 100  format('> Parametres de la fct des masses :')
      write(form_str,130) N_morc+1
 130  format('("   - Mlim       = ",',I2,'(E11.3," "))')
      write(0,form_str) (Mlim(i),i=0,N_morc)
      write(form_str,140) N_morc
 140  format('("   - Exposants  = ",',I2,'(E11.3," "))')
      write(0,form_str) (alpha(i),i=1,N_morc)
      write(0,150) Age_Amas
 150  format('   - Age_Amas   = ',E12.4,' ans')
      write(0,160) Tsf
 160  format('   - Tsf        = ',E12.4,' ans')

      call Init_random(iseed)
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_EnerGlob
c----------------------------------------------------------------------
      include 'AppliquerFM3_inc.f'
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer i,iSE
      double precision M_i,R_i,Mint,Mprec
c......................................................................
c
      Ecin_amas = 0.0d0
      Etot_amas = 0.0d0
      Egrav_amas = 0.0d0
      Mint = 0.0d0
      Mprec = 0.0d0
      do i=1,NbSE-NbSE_perte
         iSE = iListeRay(i)
         M_i = M_SE(iSE)
         R_i = R_SE(iSE)
         Ecin_amas = Ecin_amas + M_i*T_SE(iSE)
         Mint = Mint+Mprec
         Egrav_amas = Egrav_amas - M_i/R_i*(Mint+FracMasseInt*M_i)
         Mprec = M_i
      end do

      Etot_amas  = Ecin_amas + Egrav_amas
c
      end
c
c----------------------------------------------------------------------
      subroutine Affiche_Ener
c----------------------------------------------------------------------
      include 'AppliquerFM3_inc.f'
c
c var locale :
c ^^^^^^^^^^^^
      double precision Viriel
c......................................................................
c
      write(0,90) Egrav_amas,Ecin_amas,Etot_amas
 90   format(' Energie grav, cin, tot      : ',3(E12.4,' | '))
      Viriel = (2.0d0*Ecin_amas+Egrav_amas)/ABS(Etot_amas)
      write(0,100) Viriel
 100  format(' Viriel (2Ecin+Egrav)/|Etot| : ',E12.4)
c
      end 
c  
c----------------------------------------------------------------------
      subroutine Statistiques
c----------------------------------------------------------------------
      include 'AppliquerFM3_inc.f'
      include 'AppliquerFM3_Param_FM_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      integer NTypes
      parameter (Ntypes=5)
      integer NSE(NTypes)
      double precision Mtot(Ntypes),Mmin(Ntypes),Mmax(Ntypes),
     $     Met, M_tot, AgeMoy, coef_convM
      integer i,iSE
      integer iType
      character*128 Designation, Design_TypeEtoile(NTypes)
c
c fonction :
c ^^^^^^^^^^
      integer LongChaine
c......................................................................
c
      Design_TypeEtoile(TypeEtoile_MS) = 'Etoiles de sequence principale'
      Design_TypeEtoile(TypeEtoile_WD) = 'Naines blanches'
      Design_TypeEtoile(TypeEtoile_NS) = 'Etoiles de neutrons'
      Design_TypeEtoile(TypeEtoile_BH) = 'Trous noirs'
      Design_TypeEtoile(TypeEtoile_RG) = 'Geantes'

      do iType=1,Ntypes
         NSE(iType)=0
         Mtot(iType)=0.0d0
         Mmin(iType)=1.0d30
         Mmax(iType)=-1.0d30
      end do
      AgeMoy = 0.0d0

      coef_convM=0.0d0
      M_tot=0.0d0
      do i=1,NbSE-NbSE_perte
         iSE = iListeRay(i)
         AgeMoy = AgeMoy + Age_Amas-DNet_SE(iSE)
         Met = Met_SE(iSE)
         iType = iTet_SE(iSE)
         NSE(iType) = NSE(iType)+1
         Mtot(iType) = Mtot(iType)+M_SE(iSE)
         Mmin(iType) = min(Mmin(iType),Met)
         Mmax(iType) = max(Mmax(iType),Met)
         coef_convM=coef_convM+Met_SE(iSE)
         M_tot = M_tot+M_SE(iSE)
      end do
      AgeMoy = AgeMoy/(NbSE-NbSE_perte)
      coef_convM=coef_convM/M_tot

      write(0,90) AgeMoy, coef_convM*M_tot/(NbSE-NbSE_perte)
 90   format('> Statistique sur les etoiles :',/,
     $     '      Age   moyen    (yrs) : ',E12.5,/,
     $     '      Masse moyenne (Msol) : ',E12.5)
      do iType=1,NTypes
         Designation = Design_TypeEtoile(iType)
         if (NSE(iType).GT.0) then
            write(0,100) Designation(1:LongChaine(Designation)),
     $        dble(NSE(iType))/dble(NbSE-NbSE_perte),
     $        Mtot(iType)/M_tot,
     $        coef_convM*Mtot(iType)/NSE(iType),
     $        Mmin(iType),Mmax(iType)
         else
            write(0,100) Designation(1:LongChaine(Designation)),
     $           0.0d0,
     $           0.0d0,
     $           0.0d0,
     $           0.0d0,0.0d0
         end if
      end do

 100  format(
     $     '   Type Stellaire : ',A,/,
     $     '        Fraction en nombre   : ',F7.5,/,
     $     '        Fraction en masse    : ',F7.5,/,
     $     '        Masse moyenne (Msol) : ',E11.4,/,
     $     '        Masses extr. (Msol)  : ',E11.4,' ',E11.4 )
c
      end
c
c----------------------------------------------------------------------
      subroutine Meurt(Mess)
c----------------------------------------------------------------------
      character*(*) Mess
     
      write(0,*) '!!! Erreur : ',Mess,' !!!'
      call exit(1)

      end
c
