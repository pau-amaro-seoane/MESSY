c======================================================================
c     routines pour suivi de l'evolution de la segregation des masses
c     pour EvolAmas.F
c======================================================================
c----------------------------------------------------------------------
      subroutine Lire_FracM_Segr(iIn)
c----------------------------------------------------------------------
c     read the fractional masses for thelagragian radii, in file iIn
c
      include 'EvolAmas_SuiviSegr_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iIn
c
c var locales :
c ^^^^^^^^^^^^^
      integer i
      double precision frac, frac_prec
c......................................................................
c
      N_SuiviSegr = 0
      frac_prec = 0.0d0
      do while(.TRUE.)
         read(iIn,*,ERR=99,END=99) frac
         if (frac.GT.1.0d0 .OR. frac.LE.0.0d0) then
            write(0,*) '!!! Erreur dans le fichier de fractions ',
     $           'lagrangiennes : les fractions doivent etre ',
     $           'dans ]0;1] !!!',frac
            call Meurt('Lire_FracM_Segr')
         end if
         if (frac.LE.frac_prec) then
            write(0,*) '!!! Erreur dans le fichier de fractions ',
     $           'lagrangiennes : les fractions doivent etre ',
     $           'en ordre croissant !!!'
            call Meurt('Lire_FracM_Segr')
         end if
         if (N_SuiviSegr.EQ.N_SuiviSegrMax) then
            write(0,*) '!!! Erreur dans le fichier de fractions ',
     $           'lagrangiennes : ',N_SuiviSegrMax,
     $           ' fractions au maximum !!!'
            call Meurt('Lire_FracM_Segr')
         end if
         N_SuiviSegr=N_SuiviSegr+1
         M_SuiviSegr(N_SuiviSegr)=frac
         frac_prec = frac
      end do
         
 99   continue
      if (N_SuiviSegr.LT.1) then
         write(0,*) '!!! Erreur dans le fichier de fractions ',
     $        'lagrangiennes : aucune fraction specifiee !!!'
         call Meurt('Lire_FracM_Segr')
      end if
         
      do i=N_SuiviSegr+1,N_SuiviSegrMax
         M_SuiviSegr(i) = -1.0d0
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_SuiviSegr_ASC(iOut)
c----------------------------------------------------------------------
c
      include 'EvolAmas_SuiviSegr_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iOut
c
c var locales :
c ^^^^^^^^^^^^^
      integer i
c......................................................................
c
c---- ecriture en-tete
c
      write(iOut,100)N_SuiviSegr
 100  format(
     $     '# Segregation des masses pour amas',/,
     $     '# NbCouches : ',I3 )
      write(iOut,105)
 105  format('# 1: iPas_Evol 2: Tps_amas || Masse stellaire moyenne :',$)
      do i=1,N_SuiviSegr
         write(iOut,110)i+2,
     $        NINT(10000.0d0*M_SuiviSegr(i)),M_SuiviSegr(i)
 110     format(I3,' : m',I5.5,' (M<=',F6.4,') ',$)
      end do
      write(iOut,*)
      call flush_MD(iOut)
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_SuiviSegr_ASC(iOut)
c----------------------------------------------------------------------
c     L'amas doit etre trie dans iListeRay.
c
      include 'EvolAmas_SuiviSegr_inc.f'
      include 'EvolAmas_Compteurs_common.f'
      include 'EvolAmas_Temps_common.f'
c
c argument :
c ^^^^^^^^^^
      integer iOut
c
c var locales :
c ^^^^^^^^^^^^^
      integer lSE, iSE, icouche
      double precision Msomme, Msup, Nsomme
c
c functions :
c ^^^^^^^^^^^
      double precision stellar_mass
c......................................................................
c
      write(iOut,100)iPas_Evol, Tps_amas
 100  format(I10,' ',E14.6,' ',$)

      icouche = 0
      Msomme = 0.0d0
      Nsomme = 0.0d0
      lSE = 0

      do icouche=1,N_SuiviSegr
         Msup = M_SuiviSegr(icouche)*M_amas
         do while (Msomme.LE.Msup .AND.
     $        Msomme.LT.0.999999d0*M_amas)
            lSE = lSE + 1
            iSE = iListeRay(lSE)
            Msomme = Msomme + M_SE(iSE)
            Nsomme = Nsomme + M_SE(iSE)/stellar_mass(iSE)
         end do
         write(iOut,110) Msomme/Nsomme
 110     format(E14.5,' ',$)
      end do

      write(iOut,*)
      call flush_MD(iOut)
c
      end
c
c----------------------------------------------------------------------
      block data BD_Init_SuiviSegr
c----------------------------------------------------------------------
      include 'EvolAmas_SuiviSegr_common.f'
#if (_RAY_LAG_SPECIAL_==0)
      ! valeurs standard
      data N_SuiviSegr /8/
      data M_SuiviSegr /
     $     0.005,0.01,0.05,0.10,0.25,0.50,0.90,1.00,
     $     92*-1.0d0 /
#elif (_RAY_LAG_SPECIAL_==1)
      ! valeurs speciales pour comparaison dans modeles avec masse totale
      ! variee d'un facteur 2**n
      data N_SuiviSegr /11/
      data M_SuiviSegr /
     $     0.0025,0.005,0.01,0.02,0.04,0.08,0.16,0.25,0.32,0.64,1.00,
     $     89*-1.0d0 /
#endif
      end
c----------------------------------------------------------------------
