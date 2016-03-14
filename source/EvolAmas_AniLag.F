c======================================================================
c     routines pour suivi de l'evolution de l'anisotropie
c     pour EvolAmas.F
c======================================================================

c----------------------------------------------------------------------
      subroutine Init_AniLag_ASC(iOut)
c----------------------------------------------------------------------
c     Initialise l'ecriture pour un fichier indiquant le rayon des
c     couches contenant des fractions sucessives de la masse de l'amas. 
c
      include 'EvolAmas_AniLag_inc.f'
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
      write(iOut,100)N_AniLag
 100  format(
     $     '# Anisotropie pour amas',/,
     $     '# anisotropie  B = 2-<Vtan**2>/<Vrad**2>',/,
     $     '# NbCouches : ',I3 )
      write(iOut,105)
 105  format('# 1: Tps_amas || Anisotropie :',$)
      do i=1,N_AniLag
         write(iOut,110)i+1,
     $        NINT(1000.0d0*M_AniLag(i)),M_AniLag(i)-dM_AniLag(i),
     $        M_AniLag(i)
 110     format(I3,' : B',I4.4,' (M=',F5.3,'-',F5.3,') ',$)
      end do
      write(iOut,*)
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_AniLag_ASC(iOut)
c----------------------------------------------------------------------
c     L'amas doit etre trie dans iListeRay.
c
      include 'EvolAmas_AniLag_inc.f'
      include 'EvolAmas_Temps_common.f'
c
c argument :
c ^^^^^^^^^^
      integer iOut
c
c var locales :
c ^^^^^^^^^^^^^
      integer lSE, iSE, icouche
      double precision Msomme, Minf,Msup, Vtan2,Vrad2, SVtan2,SVrad2
c......................................................................
c
      write(iOut,100)Tps_amas
 100  format(E14.6,' ',$)

      icouche = 0
      Msomme = 0.0d0
      lSE = 0
      iSE = 0

      do icouche=1,N_AniLag
         Msup = min(M_AniLag(icouche),0.9999d0)*M_amas
         Minf = max(1d-10*M_amas,Msup-dM_AniLag(icouche)*M_amas)
         do while (Msomme.LT.Minf .AND. Msomme.LT.M_amas)
            lSE = lSE + 1
            iSE = iListeRay(lSE)
            Msomme = Msomme + M_SE(iSE)
         end do
         SVtan2 = 0.0d0
         SVrad2 = 0.0d0
         do while (Msomme.LE.Msup .AND. Msomme.LT.M_amas)
            Vtan2 = (J_SE(iSE)/R_SE(iSE))**2
            Vrad2 = 2.0d0*T_SE(iSE)-Vtan2
            SVtan2 = SVtan2+M_SE(iSE)*Vtan2
            SVrad2 = SVrad2+M_SE(iSE)*Vrad2
            lSE = lSE + 1
            iSE = iListeRay(lSE)
            Msomme = Msomme + M_SE(iSE)
         end do
         write(iOut,110) 2.0d0-SVtan2/SVrad2
 110     format(E14.5,' ',$)
      end do

      write(iOut,*)
c
      end
c
c----------------------------------------------------------------------
      block data BD_Init_AniLag
c----------------------------------------------------------------------
      include 'EvolAmas_AniLag_common.f'
      data N_AniLag /6/
      data M_AniLag /
     $     0.05d0,0.1d0,0.2d0,0.5d0,0.75d0,0.95d0,
     $     4*-1.0d0 /
      data dM_AniLag /
     $     6*0.05d0,
     $     4*-1.0d0 /
      end
c----------------------------------------------------------------------
