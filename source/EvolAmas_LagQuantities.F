c======================================================================
c     Routines to write Lagrange radii and various quantities summed 
c     in Lagrange spheres
c======================================================================
c
c----------------------------------------------------------------------
      subroutine Init_LagQuant_ASC(iOut)
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
c constants :
c ^^^^^^^^^^^
      integer N_quantities
      parameter (N_quantities=7)
c
c var locales :
c ^^^^^^^^^^^^^
      integer i,iR,iQ
      integer expo_coef10
      double precision coef10
      character*2000 header
      character*3 prefix, prefixs(N_quantities)
      character*3 num_str
      character*20 num_str2, frm
      data prefixs /'R','N','M','Tr','Tt','Vr2','Vt2'/
c
c functions :
c ^^^^^^^^^^^
      integer LongChaine
c......................................................................
c
c---- ecriture en-tete
c
      expo_coef10 = -( int(log10(M_RayLag(1))+100.0d0)-100 )
      expo_coef10 = max(expo_coef10,3)
      coef10 = 10.0d0**expo_coef10

      write(frm,200) expo_coef10+1,expo_coef10+1
 200  format('(I',I1,'.',I1,')')

      write(iOut,100)N_RayLag
 100  format(
     $     '# Lagrangian quantities',/,
     $     '# =====================',/,
     $     '# N_Lag_Rad : ',I3,/,
     $     '#',/,
     $     '# R = Radius of Lagrange sphere (containing approximately',
     $     ' a prescribed fraction of the remaining total mass)',/,
     $     '# N = Number of particles with R_i<= R (i.e. "within R")',/,
     $     '# M = Mass within R',/,
     $     '# Tr = 0.5*sum m_i*Vrad_i**2 for all particles ',
     $     'within the Lagrange radius (Vrad_i is the component',
     $     ' of velocity in the radial direction, i.e. away or',
     $     ' towards the centre of the cluster)',/,
     $     '# Tt = 0.5*sum m_i*Vtan_i**2 (Vtan_i is the component of',
     $     ' velocity in the tangential direction)',/,
     $     '# For an isotropic cluster, one has',
     $     ' Tr = 0.5*Tt (approximately)',/,
     $     '# Vr_2= sum Vrad_i**2',
     $     ' for all particles within the Lagrange radius',/,
     $     '# Vt_2 = sum Vtan_i**2',
     $     ' for all particles within the Lagrange radius',/,
     $     '#',/,
     $     '# N-body units are used with',
     $     ' Mclust(t=0)=1, G=1 and Eclust(t=0)=-1/4',/,
     $     '#')
      write(iOut,105)
 105  format('# 1: iPas_Evol ',$)
      
      i=2
      do iR=1,N_RayLag
         do iQ=1,N_quantities
            prefix=prefixs(iQ)
            write(num_str,'(I3)') i
            write(num_str2,frm) NINT(coef10*M_RayLag(iR))
            write(iOut,115) i, prefix(1:LongChaine(prefix)),
     $           num_str2(1:LongChaine(num_str2))
            i=i+1
         end do
      end do
 115  format(' ',I3,': ',A,'_',A,$)
      write(iOut,*)
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_LagQuant_ASC(iOut)
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
      double precision Mlim, Mclust
      integer N
      double precision M, Tr, Tt, Vr2, Vt2
      double precision Mi, Tri, Tti, Vr2i, Vt2i
c......................................................................
c
      write(iOut,100)iPas_evol
 100  format(I10,' ',$)

      Mclust = 0.0d0
      do lSE=1,NbSE-NbSE_perte
         Mclust = Mclust + M_SE(iListeRay(lSE))
      end do
         
      Mlim = 0.0d0
      icouche = 0

      N = 0
      M = 0.0d0
      Tr = 0.0d0
      Tt = 0.0d0
      Vr2 = 0.0d0
      Vt2 = 0.0d0

      lSE = 0

      do icouche=1,N_RayLag
         Mlim = min(M_RayLag(icouche),0.9999d0)*Mclust
         do while (M.LT.Mlim .AND. M.LT.Mclust)
            lSE = lSE + 1
            iSE = iListeRay(lSE)
            N = N + 1
            M = M + M_SE(iSE)
            Vt2i = (J_SE(iSE)/R_SE(iSE))**2
            Vr2i = 2.0d0*T_SE(iSE)-Vt2i
            Tri = 0.5d0*M_SE(iSE)*Vr2i
            Tti = 0.5d0*M_SE(iSE)*Vt2i
            Tr = Tr + Tri
            Tt = Tt + Tti
            Vr2 = Vr2 + Vr2i
            Vt2 = Vt2 + Vt2i
         end do
         write(iOut,110) R_SE(iSE),N,M,Tr,Tt,Vr2,Vt2
 110     format(E16.7,' ',I8,5(' ',E16.7),$)
      end do

      write(iOut,*)
c
      end
c
