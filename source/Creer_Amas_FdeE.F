      program Creer_Amas_FdeE

#define _AJUSTE_VIR_  /* Si definie, "corrige" les vitesses pour forcer la virielisation */
#define _LAXISTE_     /* Si definie, laisse passer un amas mal virielise */
#define _MODE_TIRAGE_POS_ 2 /* 1 : SE regulierement espacees en M(r) */
                            /* 2 : tirage R selon rho(R)             */
#define _INTEGRER_POISSON_ 1 /* ne pas toucher ! */

c======================================================================
c     Cree un amas a la Henon pour une fonction de distribution f(E)
c     donnee
c     Version "b" employant les unites "standard" selon Heggie & Mathieu
c     ou M=1 G=1 Etot=-0.25
c======================================================================
c
c     Type de f(E) disponibles :
c
c     1 : King
c     2 : Wooley & Dickens
c     3 : Polytrope
c     4 : Plummer
c     5 : Dehnen gamma=0
c     6 : Hernquist (Dehnen gamma=1)
c     7 : Jaffe (Dehnen gamma=2)
c     8 : Isochrone
c

#define _TYPE_ 4 

/* Malheureusement, il faut definire les noms a part car la variable _TYPE_ */
/* doit etre de 1 caractere au plus pour permettre les tests (merci cpp !) */
#if (_TYPE_==1)
#   define _NOM_ 'King'
#elif (_TYPE_==2)
#   define _NOM_ 'WooleyDickens'
#elif (_TYPE_==3)
#   define _NOM_ 'Polytrope'
#elif (_TYPE_==4)
#   define _NOM_ 'Plummer'
#elif (_TYPE_==5)
#   define _NOM_ 'Dehnen_Gamma0'
#   define _INTEGRER_POISSON_ 0
#elif (_TYPE_==6)
#   define _NOM_ 'Hernquist'
#   define _INTEGRER_POISSON_ 0
#elif (_TYPE_==7)
#   define _NOM_ 'Jaffe'
#   define _INTEGRER_POISSON_ 0
#elif (_TYPE_==8)
#   define _NOM_ 'Isochrone'
#   define _INTEGRER_POISSON_ 0
#else
      ERREUR !!!!!!
#endif

#if (_TYPE_==5)
#define _GAMMA_ 0.0d0
#define _ETA_ 3.0d0
#endif
#if (_TYPE_==6)
#define _GAMMA_ 1.0d0
#define _ETA_ 2.0d0
#endif
#if (_TYPE_==7)
#define _GAMMA_ 2.0d0
#define _ETA_ 1.0d0
#endif
c
      call Initialiser
c
c---- on entre les parametres determinant l'amas
c
      call Saisir_Parametres
c
c---- on calcule les grandeurs des SE de l'amas
c
      call Calc_Amas
c
c---- on ecrit l'amas
c
      call Ecrire_Amas
c
      end
c
c======================================================================
c
c----------------------------------------------------------------------
      subroutine Initialiser
c----------------------------------------------------------------------
c
      include 'Creer_Amas_FdeE_inc.f'
      include 'EvolAmas_TN_common.f'
c......................................................................
c
      M_amas = 1.0d0
#if (_TYPE_==1)
      NbParamFD=1
#elif (_TYPE_==2)
      NbParamFD=1
#elif (_TYPE_==3)
      NbParamFD=1
#elif (_TYPE_==4)
      NbParamFD=0
#elif (_TYPE_==5)
      NbParamFD=0
#elif (_TYPE_==6)
      NbParamFD=0
#elif (_TYPE_==7)
      NbParamFD=0
#elif (_TYPE_==8)
      NbParamFD=0
#else
      ERREUR !!!!!!!!!!!!!
#endif
c
      end
c
c----------------------------------------------------------------------
      subroutine Saisir_Parametres
c----------------------------------------------------------------------
c
      include 'Creer_Amas_FdeE_inc.f'
c
      integer i,narg,iargc,iarg
      character*256 arg
      integer LongChaine
c......................................................................
c
      iRand_Seed = 79347
      NbSE = -1
      NomOut = ''
      ParamFD(1)=1.0d40
c 
c---- lecture des options sur la ligne de commande
c
      iarg = 0
      narg = iargc()
      do while(iarg.LT.narg)
         iarg = iarg+1
         call getarg(iarg,arg)
         if (arg(1:LongChaine(arg)).EQ.'-iseed') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) iRand_Seed
         else if (arg(1:LongChaine(arg)).EQ.'-N') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) NbSE
         else if (arg(1:LongChaine(arg)).EQ.'-Name') then
            iarg = iarg+1
            call getarg(iarg,arg)
            NomOut=arg
#if (_TYPE_==1)
         else if (arg(1:LongChaine(arg)).EQ.'-W') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) ParamFD(1)
#endif
         else
            write(0,*) '!!! option ',arg(1:LongChaine(arg)),
     $           ' non comprise !!!'
            call exit(1)
         end if
      end do

      write(0,*)'>> Creation d''un amas de type ',_NOM_
      write(0,*)
      if (NbSE.LE.0) then
         write(0,'("--> Nb de SE                      : ",$)') 
         read(*,*) NbSE
      end if
      if (NbSE.GT.iDimSE) then
         write(0,*) '!!! Le nb max de particules est ',iDimSE,' !!!'
         call exit(1)
      end if
      if (NomOut.EQ.'') then
         write(0,'("--> Nom du fichier d''amas a creer : ",$)') 
         read(*,'(A)') NomOut
      end if
      i=index(NomOut,' ')-1
      if (NomOut(i-3:i).EQ.'.xdr') NomOut=NomOut(1:i-4)

      if (ParamFD(1).EQ.1.0d40) then
#if (_TYPE_==1)
         write(0,'("--> Valeur de Wo                  : ",$)')
         read(*,*)ParamFD(1)
#else
         if (NbParamFD.GT.0) then
            write(0,'("--> Parametres de la FD (",I2,")      : ",$)')
     $           NbParamFD
            read(*,*)(ParamFD(i),i=1,NbParamFD)
         end if
#endif
      end if
      write(0,*)

      write(0,100) NomOut(1:LongChaine(NomOut)), NbSE
 100  format(
     $     ' Nom_Fichier : ',A,/,
     $     ' NbSE        : ',I8)
      if (NbParamFD.GT.0) then
         write(0,110)
 110     format(' ParamFD     :',$)
         write(0,120)(ParamFD(i),i=1,NbParamFD)
 120     format(' ',E12.5,$)
         write(0,*) ''
      end if
      write(0,130) iRand_Seed
 130  format(
     $     ' iRand_Seed  : ',I10)
      write(0,*) ''
c
      call Init_random(iRand_Seed)
      NbNoeuds = NbSE
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_Amas
c----------------------------------------------------------------------
c
      include 'Creer_Amas_FdeE_inc.f'
c......................................................................
c
      call Calc_Tab
      call Normaliser_Tab
      call Tirage_SE
#ifdef _ECRIRE_TAB_
      call Ecrire_tab
#endif
c
      end
c
c--------------------------------------------------------
      subroutine Calc_Tab       
c--------------------------------------------------------
c     Cette routine calcul les tableaux R,Rho,Psi
c     R et Rho sont non normalises
c
      include 'Creer_Amas_FdeE_inc.f'
c
c constantes :
c ^^^^^^^^^^^
      double precision tiny_mass
      parameter (tiny_mass=1.0d-10)
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      double precision
     $     Rayon,
     $     Vect_PsiM(2),        ! Vect_PsiM(1) : potentiel Psi
     $                          ! Vect_PsiM(2) : M(u) masse entre 0 et r
     $                          !
     $     drRK                 ! pas d integration pour RK
      double precision
     $     Ri1,int_Omega_i,int_Omega_i1 ! afin d'integrer l energie pot totale
      double precision m_, r_, CoefR, old_R, old_Rho, old_M, old_Psi
      integer iPlat             ! permet de compter le nb de pas d integration
                                ! sans changement notable de M(r)
c
c blocs common particuliers :
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      double precision Val_FctRho
      common /ValeurDensitee/ Val_FctRho ! pour echange avec FctRho
c
c routines et fonctions :
c ^^^^^^^^^^^^^^^^^^^^^^^
      double precision FctRho, Fonc_Dist
      external RoutDerivPsiM
#if (_INTEGRER_POISSON_==0)
      double precision Rho_de_R, Psi_de_R, M_de_R
#endif
c.......................................................
c
c---- initialisations :
c
      iPlat     = 0
      NbElemTab = 0
      TabPsi(0) = 1.0d0
      TabM(0)   = 0.0d0
      TabR(0)   = 0.0d0
      Omega       = 0.0d0
      int_Omega_i = 0.0d0
      lRayonFini = .FALSE.
      Dist_typ = 1.0d0
      Dist_min = 1.0d-6*Dist_typ

#if (_TYPE_==1)
      TabPsi(0)=ParamFD(1)
#elif (_TYPE_==2)
      TabPsi(0)=1.0d0 ! a verifier absolument !!!!
#elif (_TYPE_==3)
      TabPsi(0)=1.0d0 ! a verifier absolument !!!!
#elif (_TYPE_==4)
      TabPsi(0)=1.0d0
#elif (_TYPE_==5)
      m_ = tiny_mass
      r_ = m_**0.3333333333d0
      TabM(0) = m_
      TabR(0) = r_
      TabPsi(0) = 0.5d0*(1.0d0-r_**2/(1.0d0+r_)**2) ! Cf Tremaine et al. 1994, AJ 107, 634
      Dist_typ = 1.0d0
      Dist_min = 1.0d-2*r_
#elif (_TYPE_==6)
      TabM(0) = tiny_mass
      m_ = TabM(0)
      r_ = (m_+sqrt(m_))/(1.0d0-m_)
      TabR(0) = r_
      TabPsi(0)=1.0d0/(1.0d0+r_) ! Cf Tremaine et al. 1994, AJ 107, 634
      Dist_typ = 1.0d0
      Dist_min = 1.0d-2*r_
#elif (_TYPE_==7)
      TabM(0) = tiny_mass
      m_ = TabM(0)
      r_ = m_/(1.0d0-m_)
      TabR(0) = r_
      TabPsi(0)= log(1.0d0+1.0d0/r_) ! Cf Tremaine et al. 1994, AJ 107, 634
      Dist_typ = 1.0d0
      Dist_min = 1.0d-2*r_
#elif (_TYPE_==8)
      m_ = tiny_mass
      r_ = (4.0d0*m_)**0.3333333333d0
      TabM(0) = m_
      TabR(0) = r_
      TabPsi(0) = 1.0d0/(1.0d0+sqrt(1.0d0+r_*r_)) ! Cf BT87, p38
      Dist_typ = 1.0d0
      Dist_min = 1.0d-2*r_
#else 
      ERREUR !!!!!!!!!!!!!
#endif
      TabRho(0) = FctRho(TabPsi(0))
      Rho0      = TabRho(0)
      Vect_PsiM(1) = TabPsi(0)
      Vect_PsiM(2) = TabM(0)
      Rayon        = TabR(0)

      call init_RKDOPRI54_1pas(PrecRK)
      drRK = 10.0d0*Dist_min      ! valeur initiale du pas d integration RK
c
c---- boucle de calcul des tableaux
c
      CoefR=1.2d0
 1    continue                  ! WHILE ... 

c
c---- on teste si on a atteint la limite exterieure de l'amas
c     2 criteres : - soit la masse atteint un plateau
c                  - soit le pot. reduit atteint la valeur 0
c
                                ! 1er critere
      if (
     $     (NbElemTab.GT.10) .AND.
     $     ((tabM(NbElemTab)-tabM(NbElemTab-1)) .LT.
     $     1d-8*tabM(NbElemTab-1)) ) then
         iPlat = iPlat+1
         if (iPlat.GE.10) then
            write(*,*) 'Plateau atteint pour M(r)'
            goto 99
         end if
      else
         iPlat = 0
      end if
                                ! 2eme critere

      if (tabPsi(NbElemTab).LE.0.0d0) then
         write(*,*) 'Zero de psi atteint'
         lRayonFini=.TRUE.
         goto 99
      end if
      
      if (NbElemTab.EQ.DimTab) then
         write(0,*) '!!! Tableaux ''tabXXX'' trop petits !!!'
         stop
      end if
c     
c---- pas d'integration selon Runge-Kutta
c
      Ri1 = Rayon
#if (_INTEGRER_POISSON_>0)
      call RKDOPRI54_1pas (Vect_PsiM,2,Rayon,drRK,RoutDerivPsiM)
#else /* on connait explicitement Psi(R), M(R) et Rho(R) ! */
      old_R = tabR(NbElemTab)
      old_Rho = tabRho(NbElemTab)
      old_Psi = tabPsi(NbElemTab)
      old_M = tabM(NbElemTab)
      coefR=min(1.8d0,coefR**2)
 7    continue
      if (coefR.LT.1.0001d0) then
         write(0,*) '!!! Impossible d''incrementer R !!!'
         call Meurt('Calc_Tab')
      end if
      Rayon = coefR*old_R
      Val_FctRho = Rho_de_R(Rayon)
      Vect_PsiM(1) = Psi_de_R(Rayon)
      Vect_PsiM(2) = M_de_R(Rayon)
      if (abs(old_Rho-Val_FctRho).GT.0.05d0*max(1.0d-20,old_Rho)) then
         coefR = sqrt(coefR)
         goto 7
      end if
      if (abs(old_Psi-Vect_PsiM(1)).GT.0.05d0*max(1.0d-20,old_Psi)) then
         coefR = sqrt(coefR)
         goto 7
      end if
      if (abs(old_M-Vect_PsiM(2)).GT.0.05d0*max(1.0d-20,old_M)) then
         coefR = sqrt(coefR)
         goto 7
      end if
#endif
      NbElemTab = NbElemTab + 1
      tabR  (NbElemTab) = Rayon
      tabPsi(NbElemTab) = Vect_PsiM(1)
      tabM  (NbElemTab) = Vect_PsiM(2)
      tabRho(NbElemTab) = Val_FctRho
      int_Omega_i1 = int_Omega_i
      int_Omega_i  = Rayon**2*Val_FctRho*Vect_PsiM(1)
                                ! on integre l energie potentielle selon methode
                                ! des trapezes
      Omega = Omega + 0.5d0*(Ri1-Rayon)*(int_Omega_i+int_Omega_i1)
      
      goto 1
                                ! END WHILE
 99   continue

      Omega = 2.0d0*rPi*Omega
c
      end
c
#if (_INTEGRER_POISSON_==0)
#if (_TYPE_==5) /* Modele de Dehnen avec gamma=0 */
c----------------------------------------------------------------------
      double precision function Rho_de_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c......................................................................
c
      Rho_de_R = 0.2387324146d0/(1+R)**4
c
      end
c
c----------------------------------------------------------------------
      double precision function M_de_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c......................................................................
c
      M_de_R = R*R*R/(1.0d0+R)**3
c
      end
c
c----------------------------------------------------------------------
      double precision function Psi_de_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c......................................................................
c
      Psi_de_R = 0.5d0*(1.0d0-R*R/(1.0d0+R)**2)
c
      end
c
#elif (_TYPE_==6) /* Modele de Hernquist */
c----------------------------------------------------------------------
      double precision function Rho_de_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c......................................................................
c
      Rho_de_R = 0.1591549431d0/(R*(1.0d0+R)**3)
c
      end
c
c----------------------------------------------------------------------
      double precision function M_de_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c......................................................................
c
      M_de_R = R*R/(1.0d0+R)**2
c
      end
c
c----------------------------------------------------------------------
      double precision function Psi_de_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c......................................................................
c
      Psi_de_R = 1.0d0/(1.0d0+R)
c
      end
c
#elif (_TYPE_==7) /* Modele de Jaffe */
c----------------------------------------------------------------------
      double precision function Rho_de_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c......................................................................
c
      Rho_de_R = 0.07957747155d0/(R*R*(1.0d0+R)**2)
c
      end
c
c----------------------------------------------------------------------
      double precision function M_de_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c......................................................................
c
      M_de_R = R/(1.0d0+R)
c
      end
c
c----------------------------------------------------------------------
      double precision function Psi_de_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c......................................................................
c
      Psi_de_R = log(1.0d0+1.0d0/R)
c
      end
c
#elif (_TYPE_==8) /* Modele isochrone de Henon */
c----------------------------------------------------------------------
      double precision function Rho_de_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c
c var locales :
c ^^^^^^^^^^^^^
      double precision a, a2
c......................................................................
c
      a2=1.0d0+R*R
      a=sqrt(a2)
      Rho_de_R = 0.07957747155d0*
     $     (3.0d0*(1.0d0+a)*a2-R*R*(1.0d0+3.0d0*a)) / (((1.0d0+a)*a)**3)
c
      end
c
c----------------------------------------------------------------------
      double precision function M_de_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c
c constantes :
c ^^^^^^^^^^^^
      double precision eps
      parameter (eps=1.0d-9)
c
c var locales :
c ^^^^^^^^^^^^^
      integer Nmax_tab
      parameter (Nmax_tab=10000)
      integer Ntab
      double precision x_tab(Nmax_tab), y_tab(Nmax_tab),
     $     y2_tab(Nmax_tab)
      logical lPrem
      save Ntab, x_tab, y_tab, y2_tab, lPrem
      data lPrem /.TRUE./
      double precision Rini,Rfin,Mini,Mfin,hini,hmin
      integer nok,nbad,i
c
c fonctions :
c ^^^^^^^^^^^
      external Deriv_Masse
c......................................................................
c
      if (lPrem) then
         Rini = 1.0d-7
         Mini =  0.25d0*Rini**3
         Rfin = 1.0d6
         hini = 1.0d-6
         hmin = 1.0d-12
         call odeint_rk(
     $     Mini,1,Rini,Rfin,eps,hini,hmin,nok,nbad,Deriv_Masse,
     $     Nmax_tab,Ntab,hini,x_tab,y_tab)
         call spline(x_tab,y_tab,Ntab,1.0d30,1.0d30,y2_tab)
         Mfin=y_tab(Ntab)
         do i=1,Ntab
            y_tab(i)=y_tab(i)/Mfin
         end do
         lPrem = .FALSE.
      end if

      if (R.LT.x_tab(1)) then 
         M_de_R = 0.0d0
      else if (R.GT.x_tab(Ntab)) then 
         M_de_R = 1.0d0
      else
         call splint(x_tab,y_tab,y2_tab,Ntab,R,M_de_R)
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Deriv_Masse(R,M,dMdR)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R,M,dMdR
c
c constante :
c ^^^^^^^^^^^
      double precision QuatrePi
      parameter (QuatrePi=.12566370614359172464D+02)
c
c fonctions :
c ^^^^^^^^^^^
      double precision Rho_de_R
c......................................................................
c
      dMdR=QuatrePi*R*R*Rho_de_R(R)
c
      end
c
c----------------------------------------------------------------------
      double precision function Psi_de_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c......................................................................
c
      Psi_de_R = 1.0d0/(1.0d0+sqrt(1.0d0+R*R))
c
      end
c
#endif 
#endif /* if (_INTEGRER_POISSON_==0) */
#ifdef _ECRIRE_TAB_
c
c-----------------------------------------------------------------
      subroutine Ecrire_Tab
c-----------------------------------------------------------------
c
      include 'Creer_Amas_FdeE_inc.f'
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      character*80 NomFich
      integer i
      double precision RenRho
c.................................................................
c
      RenRho=RenR**(-3)
      call ConcatChaines(NomOut,'.tab',NomFich)
      open(iFich,file=NomFich,form='formatted')
      write(iFich,100)NbElemTab
 100  format (
     $     '# Tableaux pour amas cree par Creer_Amas_FdeE.F',/,
     $     '# Nb d''elements : ',I6,/,
     $     '# WARNING: NOT ACCURATE (use Integrer_Amas.F)',/,
     $     '# N-body normalisation is done for particles,',
     $     ' not for the arrays !',/,
     $     '# 1:R 2:Rho 3:M 4:Psi'
     $     )
      do i=1,NbElemTab
         write(iFich,110) RenR*tabR(i),RenRho*tabRho(i),tabM(i),
     $        tabPsi(i)
 110     format(4(' ',E14.6))
      end do
      close(iFich)
c
      end
c
#endif
c
c-------------------------------------------------------
      subroutine RoutDerivPsiM(PsiM,r,der_PsiM)
c-------------------------------------------------------
      include 'Creer_Amas_FdeE_inc.f'
c
c variables d'appel :
c ^^^^^^^^^^^^^^^^^^^
      double precision PsiM(2),r,der_PsiM(2)
c
c fonction(s) :
c ^^^^^^^^^^^^^
      double precision FctRho
c.......................................................
c
      if (r.GT.1.0d-15) then
         der_PsiM(1)=-PsiM(2)/(r*r)
      else
         der_PsiM(1)=-0.33333333d0*r4Pi*Rho0
      end if
      der_PsiM(2)=r4Pi*FctRho(PsiM(1))*(r*r)
c
      end
c
c------------------------------------------------------------
      double precision function FctRho(Psi)
c------------------------------------------------------------
c
      include 'Creer_Amas_FdeE_inc.f'
c
c variables d'appel :
c ^^^^^^^^^^^^^^^^^^^
      double precision Psi
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      double precision Vmax
c
c blocs common particuliers :
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      double precision Psi_Int
      common /Psi_Cour/ Psi_Int
      double precision Val_FctRho
      common /ValeurDensitee/ Val_FctRho 
c
c fonctions :
c ^^^^^^^^^^^
      double precision Integrand
      external Integrand
      double precision Int_RKDOPRI54
c
c...........................................................
c
      Psi_Int = Psi
      if (Psi.GT.0.0d0) then
         Vmax = sqrt(2.0d0*Psi)
c
c---- Calcul de la densite pour le potentiel Psi
c     par integration num de la FD sur les vitesses
c     => il s'agit d'integrer :
c        4*Pi*Integrale( V**2*FD(-0.5*V**2+Psi) ; V : 0->Vmax )
c
         call qromb(Integrand,0.0d0,Vmax,1.0d-6,Val_FctRho)
         Val_FctRho = r4Pi * Val_FctRho
c$$$         Val_FctRho = r4Pi *
c$$$     $        Int_RKDOPRI54(Integrand,0.0d0,Vmax,0.0001d0*Vmax,PrecRK)
      else
         Val_FctRho = 0.0d0
      end if

      FctRho = Val_FctRho
c
      end
c
c--------------------------------------------------
      double precision function Integrand(V)
c--------------------------------------------------
c
      implicit none
c
c variables d'appel :
c ^^^^^^^^^^^^^^^^^^^
      double precision V
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      double precision Valeur
c
c blocs common particuliers :
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      double precision Psi_Int
      common /Psi_Cour/ Psi_Int
c
c fonction de distribution F(e) :
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      double precision Fonc_Dist
c..................................................
c
      Valeur = V**2*Fonc_Dist(-0.5d0*V**2+Psi_Int)
      Integrand = Valeur

#ifdef _SURVEILLER_INT_
      write(66,*)V,Valeur
#endif
      end
c      
c------------------------------------------------------------
      subroutine Normaliser_Tab
c------------------------------------------------------------
c
      include 'Creer_Amas_FdeE_inc.f'
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      double precision Eta,Eta2,UnSurEta    ! facteur de normalisation
      integer i
c............................................................
c
      Eta = tabM(NbElemTab)/M_amas
      Eta2 = Eta**2
      UnSurEta = 1.0/Eta
c$$$      Omega=UnSurEta*Omega 
c$$$      write(0,600) UnSurEta,UnSurEta,Eta2
c$$$ 600  format(
c$$$     $     '> Renormalisation factors: ',/,
c$$$     $     '   Renorm_M  : ',E12.5,/,
c$$$     $     '   Renorm_R  : ',E12.5,/,
c$$$     $     '   Renorm_Rho: ',E12.5,/)
      do i=1,NbElemTab
         tabR(i)   = UnSurEta*tabR(i)
         tabRho(i) = Eta2*tabRho(i)
         tabM(i)   = UnSurEta*tabM(i)
      end do
      if (lRayonFini) then
                                ! Cas ou Eps_max < 0
                                ! Rayon fini
         Rayon_Amas = (
     $        tabR(NbElemTab)*tabPsi(NbElemTab-1) -
     $        tabR(NbElemTab-1)*tabPsi(NbElemTab) 
     $        ) / ( tabPsi(NbElemTab-1)-tabPsi(NbElemTab) )
         tabR(NbElemTab) = Rayon_Amas
         Omega = Omega-0.5d0*M_amas**2/Rayon_Amas

      else
                                ! Cas ou Eps_max = 0
                                ! Rayon infini ( en principe)
         Rayon_Amas = tabR(NbElemTab) 
      end if
      tabM(NbElemTab)   = M_amas 
      tabPsi(NbElemTab) = 0.0d0
      tabRho(NbElemTab) = 0.0d0
c
      end
c    
c------------------------------------------------------------
      subroutine Tirage_SE
c------------------------------------------------------------
c
      include 'Creer_Amas_FdeE_inc.f'
      include 'EvolAmas_Etoiles_common.f'
      include 'EvolAmas_Temps_common.f'
      include 'EvolAmas_Maree_common.f'
      include 'EvolAmas_Pertes_inc.f'
      include 'EvolAmas_Conserv_inc.f'
c
c constantes :
c ^^^^^^^^^^^^
      double precision UnSurRac2
      parameter ( UnSurRac2 = .70710678118654757274D+00 )
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      double precision
     $     R_i, Mr_i, dMSE_dMTab, V_i, Psi_i,
     $     M__SE, F_VN, V2, Vmax, Fmax, f_V, VradSE, VtgSE
      integer iSE,iTab
      integer itick, inext_tick

c
c fonction(s) :
c ^^^^^^^^^^^^^
      double precision  Fonc_Dist, random
c............................................................
c
c---- initialisations
c
      Tps_amas   = 0.0d0
      Egrav_amas = 0.0d0
      Ecin_amas  = 0.0d0
      Etot_amas  = 0.0d0
      NbSE_perte  = 0
      M__SE = M_amas/NbSE
      iTab=0
      Mr_i=-0.5d0*M__SE
      
      write(0,90)
 90   format(/,
     $     '###########################',/,
     $     '## Creating cluster..... ##',/,
     $     '###########################',//,
     $     '0    10   20   30   40   50   60   70   80   90  100%',/,
     $     '|----|----|----|----|----|----|----|----|----|----|',/,
     $     $)
      inext_tick=0
      itick=0

      do iSE=1,NbSE
         if (iSE.GE.inext_tick) then
            write(0,95)
 95         format('^',$)
            itick=itick+1
            inext_tick = 0.02d0*NbSE*itick
         end if
#if (_MODE_TIRAGE_POS_==1) /* SE regulierement espacees en M(r) */
         Mr_i=Mr_i+M__SE
c
c---- on cherche l'emplacement adequat dans le tableau TabM tq
c     TabM(iTab) <= Mr(i) < TabM(iTab+1)
c
 1       continue
         if (TabM(iTab+1).LT.Mr_i) then
            iTab=iTab+1
            if (iTab.GT.NbElemTab) then
               write(0,*) ' !!! Indice iTab > NbElemTab dans routine ',
     $              'Tirage_SE   !!!'
               stop
            end if
            goto 1
         end if
#elif (_MODE_TIRAGE_POS_==2) /* R_SE tire aleatoirement selon rho(R) */
         Mr_i=M_amas*random()
c
c---- on determine R_SE par inversion de la relation Mint(r)
c
         call locate(TabM,NbElemTab,Mr_i,iTab)
#endif /* if (_MODE_TIRAGE_POS_==1) */
         dMSE_dMtab  = (Mr_i-TabM(iTab))/(TabM(iTab+1)-TabM(iTab))
         R_i   = TabR(iTab)   + (TabR(iTab+1)-TabR(iTab))    *dMSE_dMTab
         Psi_i = TabPsi(iTab) + (TabPsi(iTab+1)-TabPsi(iTab))*dMSE_dMTab
c
c---- tirage de la vitesse selon methode de Von Neumann
c     
c$$$                                ! calculons la valeur max de distribution
c$$$                                ! de vitesse locale
c$$$                                ! max de 0.5*v**2*FD(0.5*v**2+Psi(r))
c$$$         E_max = Fct_Emax(v)
c$$$         Fmax = (-E_max+Psi_i)*Fonc_Dist(E_max) 

                                ! plus rustique mais plus simple :
                                ! borne superieure pour la fct de dist de V <
                                ! Vmax**2 * FDmax(r) !!
         Fmax = Psi_i * Fonc_Dist(Psi_i) 
c                \         \______________valeur max de la FD en R_i
c                 \
c                  \______________________valeur max de 0.5*V**2 (Eps=0)

         if (Fmax.LE.0.0d0) then
            write(0,*) '!!! Fmax = ',Fmax,' !!!'
            write(0,*) R_i, Psi_i, Fonc_Dist(Psi_i)
         end if

         Vmax = sqrt(2.0*Psi_i)

 2       continue
         V_i  = Vmax*random()
         F_VN = Fmax*random()
         V2 = 0.5d0*V_i*V_i
         f_V = V2*Fonc_Dist(-V2+Psi_i)

         if (f_V.GT.Fmax) then
            write(0,*) '!!! f_V(=',f_V,') > Fmax(=',Fmax,') !!!'
            stop
         end if

         if (F_VN.GT.f_V) then
            goto 2
         end if
c
c---- on attribue ses proprietes a la SE :
c
         M_SE   (iSE) = M__SE
         R_SE   (iSE) = R_i
 3       continue
         VradSE = V_i*(1.0-2.0*random())
	 VtgSE = sqrt(V2-0.5d0*VradSE**2)

                                ! faisons en sorte d`eviter les orbites purement radiales
                                ! ou purement circulaires afin de ne pas avoir
                                ! de complication numeriques
         if ( (VtgSE .LT.1.0d-6*V_i) .OR. (VradSE.LT.1.0d-6*V_i) ) then
            goto 3
         endif
         J_SE(iSE) = R_i*VtgSE*rRacDeux
         T_SE(iSE) = V2
         Met_SE(iSE) = 1.0d0    ! a default d'autre chose
#ifdef _CALC_DIR_ECIN_
         Ecin_amas = Ecin_amas + M__SE*V2
#endif
                                ! REMARQUE :
                                ! La valeur du moment cinetique des SE est
                                ! fixee dans la routine "CalculerE"
      end do
      write(0,*) ''
c
c---- Modification des vitesses pour satisfaire le viriel
c
      call TriAmas
      call Calc_EnerGlob
      call Affiche_Ener
      call CorrigeViriel
      call Calc_EnerGlob
      call Affiche_Ener
c
c---- renormalisation afin de correspondre aux unites "standard"
c
      renR=-2.0d0*Egrav_amas
      renT=1.0d0/renR
      renJ=renR*sqrt(renT)
      do iSE=1,NbSE
         R_SE(iSE)=renR*R_SE(iSE)
         J_SE(iSE)=renJ*J_SE(iSE)
         T_SE(iSE)=renT*T_SE(iSE)
      end do
      R_mar = renR*tabR(NbElemTab)
      C_mar = R_Mar/M_amas**0.3333333333d0
      call Calc_EnerGlob
      call Affiche_Ener
      call Affiche_Taille
c
      end
c 
c----------------------------------------------------------------------
      subroutine Calc_EnerGlob
c----------------------------------------------------------------------
      include 'Creer_Amas_FdeE_inc.f'
      include 'EvolAmas_Conserv_inc.f'
      include 'EvolAmas_Pertes_inc.f'
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
      dEtot_evap = 0.0d0
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
                                ! Energie des SE libres
      do i=NbSE-NbSE_perte+1,NbSE
         iSE = iListeRay(i)
         dEtot_evap = dEtot_evap + M_SE(iSE)*T_SE(iSE)
      end do

      Etot_amas  = Ecin_amas + Egrav_amas
c
      end
c  
c------------------------------------------------------------
      subroutine Affiche_Ener
c------------------------------------------------------------
c
      include 'Creer_Amas_FdeE_inc.f'
      include 'EvolAmas_Conserv_inc.f'
c............................................................
c
      write(*,90) Egrav_amas,Ecin_amas,Etot_amas
 90   format(' Energie grav, cin, tot      : ',3(E12.4,' | '))
      Viriel = (2.0d0*Ecin_amas+Egrav_amas)/ABS(Etot_amas)
      write(*,100) Viriel
 100  format(' Viriel (2Ecin+Egrav)/|Etot| : ',E12.4)
c
      end 
c  
c------------------------------------------------------------
      subroutine Affiche_Taille
c------------------------------------------------------------
c
      include 'Creer_Amas_FdeE_inc.f'
      include 'EvolAmas_Maree_common.f'
      include 'EvolAmas_Pertes_inc.f'
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer i,iSE
      double precision Rmin,Rmax,Rhalf,Mint
c............................................................
c
      Rmin=rEnorme
      Rmax=-rEnorme
      Mint=0.0d0
      Rhalf = 0.0d0
      do i=1,NbSE-NbSE_perte
         iSE=iListeRay(i)
         Mint=Mint+M_SE(iSE)
         if (Rhalf.EQ.0.0d0  .AND. Mint.GT.0.5d0*M_amas) then
            Rhalf = R_SE(iListeRay(i-1)) +
     $           (R_SE(iSE)-R_SE(iListeRay(i-1)))/M_SE(iSE)*
     $           (0.5d0*M_amas-Mint+M_SE(iSE))
         end if
         Rmin=min(Rmin,R_SE(iSE))
         Rmax=max(Rmax,R_SE(iSE))
      end do
      write(0,100) Rmin,Rmax,Rhalf,R_mar
 100  format(' > Rmin,Rmax,Rhalf,R_mar = ',4(E12.4," "))
c
      end
c
c------------------------------------------------------------
      subroutine CorrigeViriel
c------------------------------------------------------------
c
      include 'Creer_Amas_FdeE_inc.f'
      include 'EvolAmas_Conserv_inc.f'
c   
c variables locales :
c ^^^^^^^^^^^^^^^^^^^    
      double precision Fact_Corr_Vit 
      integer iSE
c............................................................
c
c---- testons si l'amas est virialise
c 
#ifndef _LAXISTE_
      if (ABS(Viriel).GT.tol_Vir/sqrt(0.0001*NbSE)) then
         write(0,*) '!!! La valeur du viriel est innacceptable, on stoppe tout !!!'
         stop
      end if
#endif	
#ifdef _AJUSTE_VIR_
                                ! on modifie les vitesses afin d`obtenir un amas
                                ! parfaitement virialise (en principe l`energie
                                ! gravitationnelle est plus 'fiable')
      
      write(*,*) ' Modification des vitesses des ',
     $     'SE afin de satisfaire le viriel ...'
      Fact_Corr_Vit = sqrt(-Egrav_amas/(2.0d0*Ecin_amas))
      do iSE=1,NbSE
         T_SE(iSE) = Fact_Corr_Vit**2*T_SE(iSE) 
         J_SE(iSE) = Fact_Corr_Vit*J_SE(iSE)
      end do
#endif
c
      end 
c   
c-----------------------------------------------------------------
      double precision function Fonc_Dist(Eps)
c-----------------------------------------------------------------
c     La fonction de distribution exprimee en fct de 
c     l'energie reduite F(Epsilon)
c
      include 'Creer_Amas_FdeE_inc.f'
c
c variables d'appel :
c ^^^^^^^^^^^^^^^^^^^
      double precision Eps
c
c fonctions :
c ^^^^^^^^^^^
      double precision dawson, erfcc
c
c var locales :
c ^^^^^^^^^^^^^
      double precision sq2e
c.................................................................
c
      if (Eps.GT.0.0d0 .AND. Eps.LE.1.0001d0*TabPsi(0)) then
#if (_TYPE_==1)
         Fonc_Dist = exp(Eps)-1.0d0
#elif (_TYPE_==2)
         Fonc_Dist = exp(ParamFD(1)*Eps)
#elif (_TYPE_==3)
         Fonc_Dist = Eps**ParamFD(1) ! exposant de 3.5 correspond au modele de Plummer
#elif (_TYPE_==4)
         Fonc_Dist = Eps**3.50d0
#elif (_TYPE_==5) 
         sq2e=sqrt(2.0d0*Eps)
         if (Eps.GT.1.0d-5) then
            Fonc_Dist =         ! Cf Tremaine et al. 1994, AJ 107, 634, Eq. 21
     $           0.02418865082d0*( 2.0d0*sq2e*
     $           (3.0d0-4.0d0*Eps)/(1.0d0-2.0d0*Eps) +
     $           3.0d0*log((1.0d0-sq2e)/(1.0d0+sq2e)) )
         else
            Fonc_Dist = 0.02418865082d0*( 3.2d0*sq2e*Eps**2 ) ! From Series development of Fonc_Dist/sq2e near 0
         end if
#elif (_TYPE_==6)
         if (Eps.GT.1.0d0) then
            write(0,*) '!!! Appel de FD(Eps) Hernquist avec Eps = ',Eps,
     $           ' !!!'
            call Meurt('Fonc_Dist')
         end if
         Fonc_Dist =      ! Cf Tremaine et al. 1994, AJ 107, 634, Eq. 21
     $        0.002850659838d0*( 3.0d0*asin(sqrt(Eps)) -
     $        sqrt(Eps*(1.0d0-Eps))*
     $        (3.0d0+2.0d0*Eps*(1.0d0+Eps*(-12.0d0+8.0d0*Eps))) ) /
     $        (1.0d0-Eps)**2.5d0 
#elif (_TYPE_==7) 
         Fonc_Dist =      ! Cf Tremaine et al. 1994, AJ 107, 634, Eq. 21
     $        0.02021065203d0*(
     $                        -exp(Eps)*erf(sqrt(Eps)) +
     $        0.7071067812d0 * exp(Eps)**2*erf(sqrt(2.0d0*Eps)) -
     $        1.128379167d0  * dawson(sqrt(Eps)) +
     $        0.7978845608d0 * dawson(sqrt(2.0d0*Eps)) )
#elif (_TYPE_==8) 
         if (Eps.GT.1.0d0) then
            write(0,*) '!!! Appel de FD(Eps) isochrone avec Eps = ',Eps,
     $           ' !!!'
            call Meurt('Fonc_Dist')
         end if
         Fonc_Dist =      ! Cf BT87 p239
     $        0.002850659838d0*sqrt(Eps)/(2.0d0*(1.0d0-Eps))**4 * (
     $        27.0d0+Eps*(-66.0d0+Eps*(320.0d0+Eps*(-240+Eps*64.0d0))) +
     $        3.0d0*(-9.0d0+Eps*(28.0d0+16.0d0*Eps)) *
     $        asin(sqrt(Eps))/sqrt(Eps*(1.0d0-Eps)) )
#else
         ERREUR !!!!!!!!!
#endif
      else
         Fonc_Dist = 0.0d0
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_Amas
c----------------------------------------------------------------------
c
      include 'Creer_Amas_FdeE_inc.f'
      include 'EvolAmas_Conserv_inc.f'
      include 'EvolAmas_Maree_common.f'
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      character*80 NomFich
      character*80 Form_
      integer i, iXDR, ixdr_err
      logical lexiste
c
c fonction :
c ^^^^^^^^^^
      integer initxdr, ixdrclose
c......................................................................
c 
c---- fichier d'amas (XDR)
c
      call ConcatChaines(NomOut,'.xdr',NomFich)
      inquire(file=NomFich,exist=lExiste)
      if (lExiste) then
         write(0,*) '!!! Le fichier ',NomFich(1:index(NomFich,' ')-1),
     $        ' existe deja !!!'
         call exit(1)
      end if
      !open(iFich,file=Nomfich,form='unformatted',status='new')
      iXDR = initxdr(NomFich,'w',.TRUE.)
      call EcritAmas_XDR(iXDR)
      ixdr_err = ixdrclose(iXDR)

#if (_TYPE_==1 || _TYPE_==2 || _TYPE_==3 )
c
c---- fichier de maree
c
      call ConcatChaines(NomOut,'_MAREE.xdr',NomFich)
      inquire(file=NomFich,exist=lExiste)
      if (lExiste) then
         write(0,*) '!!! Le fichier ',NomFich(1:index(NomFich,' ')-1),
     $        ' existe deja !!!'
         call exit(1)
      end if
      !open(iFich,file=NomFich,form='unformatted',status='new')
      iXDR = initxdr(NomFich,'w',.TRUE.)
      call EcritMaree_XDR(iXDR)
      ixdr_err = ixdrclose(iXDR)
#endif
c
c---- fichier de documentation
c
      call ConcatChaines(NomOut,'.info',NomFich)
      open(iFich,file=NomFich,form='formatted')
      write(iFich,100)
     $     NomOut,
     $     NbSE,iRand_Seed,
     $     M_amas,Psi0,Rayon_Amas,
     $     Egrav_amas,Ecin_amas,Etot_amas,
     $     Viriel,
     $     Omega,
     $     R_mar
 100  format(
     $     ' Parametres de fichier d''amas cree par Creer_Amas_FdeE',/,
     $     '=======================================================',/,
     $     ' Nom du fichier d''amas       : ',A,/,
     $     ' Nom du type d''amas          : ',_NOM_,/,
     $     ' Nb de super-etoiles         : ',I8,/,
     $     ' Graine du gen. aleatoire    : ',I10,/,
     $     ' Masse de l''amas             : ',E12.4,/,
     $     ' Potentiel red. au centre    : ',E12.4,/,
     $     ' Rayon de l''amas             : ',E12.4,/,
     $     ' Energie grav, cin, tot      : ',3(E12.4,' | '),/,
     $     ' Viriel (2Ecin+Egrav)/|Etot| : ',E12.4,/,
     $     ' Energie pot. integree       : ',E12.4,/,
     $     ' Rayon de maree              : ',E12.4
     $     )
      if (NbParamFD.GT.0) then
         write(Form_,110)NbParamFD
 110     format(
     $        '('' Parametres du potentiel     :'',/,',I2,'('' '',E12.4))'
     $        )
         write(iFich,Form_) (ParamFD(i),i=1,NbParamFD)
      end if
      
      close(iFich)  
c
      end
c
c----------------------------------------------------------------------
      SUBROUTINE locate(xx,n,x,j)
c----------------------------------------------------------------------
      INTEGER j,n
      DOUBLE PRECISION x,xx(0:n)
      INTEGER jl,jm,ju
      jl=-1
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END
c----------------------------------------------------------------------
      subroutine Erreur
c----------------------------------------------------------------------
c     
      write(0,*)' UNE ERREUR S''EST PRODUITE !!!'
      stop
c
      end
c
c----------------------------------------------------------------------
      block data init
      include 'Creer_Amas_FdeE_inc.f'
      data renR,renJ,renT /1.0d0,1.0d0,1.0d0/
      end
c----------------------------------------------------------------------
      
