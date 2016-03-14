c======================================================================
c     routines pour suivi de l'evolution des rayons lagrangiens
c     pour EvolAmas.F
c======================================================================
c----------------------------------------------------------------------
      subroutine Lire_FracM_RayLag(iIn)
c----------------------------------------------------------------------
c     read the fractional masses for thelagragian radii, in file iIn
c
      include 'EvolAmas_RayLag_inc.f'
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
      N_RayLag = 0
      frac_prec = 0.0d0
      do while(.TRUE.)
         !write(0,*) 'xxxx' !#########
         read(iIn,*,ERR=99,END=99) frac
         !write(0,*) 'yyyyy' !#########
         if (frac.GT.1.0d0 .OR. frac.LE.0.0d0) then
            write(0,*) '!!! Erreur dans le fichier de fractions ',
     $           'lagrangiennes : les fractions doivent etre ',
     $           'dans ]0;1] !!!',frac
            call Meurt('Lire_FracM_RayLag')
         end if
         if (frac.LE.frac_prec) then
            write(0,*) '!!! Erreur dans le fichier de fractions ',
     $           'lagrangiennes : les fractions doivent etre ',
     $           'en ordre croissant !!!'
            call Meurt('Lire_FracM_RayLag')
         end if
         if (N_RayLag.EQ.N_RayLagMax) then
            write(0,*) '!!! Erreur dans le fichier de fractions ',
     $           'lagrangiennes : ',N_RayLagMax,
     $           ' fractions au maximum !!!'
            call Meurt('Lire_FracM_RayLag')
         end if
         N_RayLag=N_RayLag+1
         M_RayLag(N_RayLag)=frac
         frac_prec = frac
      end do
         
 99   continue
      if (N_RayLag.LT.1) then
         write(0,*) '!!! Erreur dans le fichier de fractions ',
     $        'lagrangiennes : aucune fraction specifiee !!!'
         call Meurt('Lire_FracM_RayLag')
      end if
         
      do i=N_RayLag+1,N_RayLagMax
         M_RayLag(i) = -1.0d0
      end do
c
      end
c----------------------------------------------------------------------
      subroutine Init_RayLag_ASC(iOut)
c----------------------------------------------------------------------
c     Initialise l'ecriture pour un fichier indiquant le rayon des
c     couches contenant des fractions sucessives de la masse de l'amas. 
c
      include 'EvolAmas_RayLag_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iOut
c
c var locales :
c ^^^^^^^^^^^^^
      integer i, expo_coef10
      double precision coef10
      character*256 frm
c......................................................................
c
c---- ecriture en-tete
c
      expo_coef10 = -( int(log10(M_RayLag(1))+100.0d0)-100 )
      expo_coef10 = max(expo_coef10,3)
      coef10 = 10.0d0**expo_coef10
      write(frm,200) expo_coef10+1,expo_coef10+1
 200  format('(I3,": R",I',I1,'.',I1,'," (",F8.6,") ",$)')
      write(iOut,100)N_RayLag
 100  format(
     $     '# Rayons massiques pour amas',/,
     $     '# NbCouches : ',I3 )
      write(iOut,105)
 105  format('# 1: iPas_Evol 2: NbSE_liees 3: M_liee 4: Tps_amas',
     $     ' || FracMasse :',$)
      do i=1,N_RayLag
         write(iOut,frm)i+4,
     $        NINT(coef10*M_RayLag(i)),M_RayLag(i)
      end do
      write(iOut,*)
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_RayLag_ASC(iOut)
c----------------------------------------------------------------------
c     L'amas doit etre trie dans iListeRay.
c
      include 'EvolAmas_RayLag_inc.f'
      include 'EvolAmas_Temps_common.f'
      include 'EvolAmas_Compteurs_common.f'
c
c argument :
c ^^^^^^^^^^
      integer iOut
c
c var locales :
c ^^^^^^^^^^^^^
      integer lSE, iSE, icouche
      double precision Msomme, Mlim
c......................................................................
c
      write(iOut,100)iPas_evol, NbSE-NbSE_perte, M_amas, Tps_amas
 100  format(I10,' ',I7,' ',2(E14.5,' '),$)

      Mlim = 0.0d0
      icouche = 0
      Msomme = 0.0d0
      lSE = 0 
      iSE = 0 ! useless

      do icouche=1,N_RayLag
         Mlim = min(M_RayLag(icouche),0.9999d0)*M_amas
         do while (Msomme.LT.Mlim .AND. Msomme.LT.M_amas)
            lSE = lSE + 1
            iSE = iListeRay(lSE)
            Msomme = Msomme + M_SE(iSE)
         end do
         write(iOut,110) R_SE(iSE)
 110     format(E14.5,' ',$)
      end do

      write(iOut,*)
c
      end
c
c----------------------------------------------------------------------
      block data BD_Init_RayLag
c----------------------------------------------------------------------
      include 'EvolAmas_RayLag_common.f'
#if (_RAY_LAG_SPECIAL_==0)
      ! valeurs standard
      data N_RayLag /17/
      data M_RayLag /
     $     0.001,0.005,0.01,0.02,0.05,0.10,0.15,0.20,0.30,
     $     0.40, 0.50,0.60,0.70,0.80,0.90,0.95,0.99,
     $     83*-1.0d0 /
#elif (_RAY_LAG_SPECIAL_==1)
      ! valeurs speciales pour comparaison dans modeles avec masse totale
      ! variee d'un facteur 2**n
      data N_RayLag /14/
      data M_RayLag /
     $     0.00125,0.0025,0.005,0.01,0.02,0.04,0.08,0.16,0.25,0.32,0.64,
     $     0.75,0.90,0.95,
     $     86*-1.0d0 /
#endif
      end
c----------------------------------------------------------------------
