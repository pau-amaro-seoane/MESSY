c
#if (_SANS_PG_==0)
c----------------------------------------------------------------------
      subroutine Calc_PorbTrel_SE
     $     (Mpart,Epart,Jpart,Peri_part,Apo_part,eps, Porb,Trel)
c----------------------------------------------------------------------
c     calcule la periode orbitale et le temps de relaxation moyen 
c     pour une particule extraite de l'amas
c
      include 'EvolAmas_IntegrOrb_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $     Mpart,Epart,Jpart,Peri_part,Apo_part, ! Proprietes physique de la part.
     $     eps,                 ! Precision de l'integration
     $     Porb,Trel            ! Periode orbitale, moyenne orbitale du "tps de relax"
c
c var locales :
c ^^^^^^^^^^^^^
      integer j,n,kbase,iPG,iRang
      double precision V2,Rho,sigma,Porb_prec,Trel_prec,x,R,
     $     R_p,R_m,Vr,dPorb
c
c fonctions :
c ^^^^^^^^^^^
      double precision Rho_PG_lin, Sigma1D_PG_lin2, Coef_Relax
      integer iCouche_PG
c......................................................................
c
c---- au premier appel, on calcule les abscisses pour G-C
c
      if (x_GC(1).LE.-2.0d0) call Init_GC
         
      Porb = 0.0d0
      Trel = 0.0d0
c
c---- on integre dt/dr de Rperi a Rapo pour calculer la periode orbitale
c     et 1/trel pour calculer le tps de relax
c 
      R_p = 0.5d0*(Apo_part+Peri_part)
      R_m = 0.5d0*(Apo_part-Peri_part)
      Trel_prec = 1.0d20
      Porb_prec = 1.0d20
      n = 0                     ! ordre de l'integration
      do while ( (abs(Porb-Porb_prec).GT.0.5d0*eps*abs(Porb)) .AND.
     $     (abs(Trel-Trel_prec).GT.0.5d0*eps*abs(Trel)) )
         if (n.EQ.n_GC_max) then
            write(0,*)'!!! n_GC_max atteint !!!'
            call Meurt('Calc_PorbTrel_SE')
         end if
         n = n + 1              ! on augmente l'ordre
         kbase = ((n-1)*n)/2
         Trel_prec = Trel
         Porb_prec = Porb
         Trel = 0.0d0
         Porb = 0.0d0
         do j = 1,n
            x = x_GC(kbase+j) 
            R = R_p+R_m*x       ! evaluer l'integrand a l'abscisse de G-C
            call Calc_RangVitesse(Mpart,Epart,Jpart,R,iRang,V2)
            Vr = V2-(Jpart/R)**2
#if (_SURVEILLE_SQRT_>0)
            if (Vr.LT.0.0d0) then
               write(0,*) '!!! Vr**2 = ',Vr,' !!!'
               call Meurt('Calc_PorbTrel_SE')
            end if
#endif
            Vr=sqrt(Vr)
            iPG = iCouche_PG(R)
            Rho = Rho_PG_lin(iPG,R)
            sigma = Sigma1D_PG_lin2(iPG,R)
            dPorb = sqrt(1.0d0-x*x)/Vr ! c'est l'integrand divise par la fct poids
            Porb = Porb + dPorb        ! W(x)=1/sqrt(1-x**2)
            Trel = Trel + Coef_Relax(iRang)*dPorb*Rho/
     $           (sqrt(V2+3.0d0*sigma**2))**3
         end do
         Porb = Porb/dble(n)    ! on multiplie par le poids constant
         Trel = Trel/dble(n)
      end do

                                ! Normalisations
      Trel = Porb/(Met_moy*Trel)
      Porb = r2Pi*R_m*Porb
c
      end
c
#endif /* (_SANS_PG_==0) */
c
c----------------------------------------------------------------------
      subroutine Calc_Porb_SE
     $     (Mpart,Epart,Jpart,Peri_part,Apo_part,eps, Porb)
c----------------------------------------------------------------------
c     calcule la periode orbitale 
c     pour une particule extraite de l'amas
c
      include 'EvolAmas_IntegrOrb_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $     Mpart,Epart,Jpart,Peri_part,Apo_part, ! Proprietes physique de la part.
     $     eps,                 ! Precision de l'integration
     $     Porb                 ! Periode orbitale
c     
c var locales :
c ^^^^^^^^^^^^^
      integer j,n,kbase,iRang
      double precision V2,Porb_prec,x,R,R_p,R_m,Vr,dPorb
c......................................................................
c
c---- au premier appel, on calcule les abscisses pour G-C
c
      if (x_GC(1).LE.-2.0d0) call Init_GC
         
      Porb = 0.0d0
c
c---- on integre dt/dr de Rperi a Rapo pour calculer la periode orbitale
c     et 1/trel pour calculer le tps de relax
c 
      R_p = 0.5d0*(Apo_part+Peri_part)
      R_m = 0.5d0*(Apo_part-Peri_part)
      Porb_prec = 1.0d20
      n = 0                     ! ordre de l'integration
      do while ( (abs(Porb-Porb_prec).GT.0.5d0*eps*abs(Porb)) )
         if (n.EQ.n_GC_max) then
            write(0,*)'!!! n_GC_max atteint !!!'
            call Meurt('Calc_Porb_SE')
         end if
         n = n + 1              ! on augmente l'ordre
         kbase = ((n-1)*n)/2
         Porb_prec = Porb
         Porb = 0.0d0
         do j = 1,n
            x = x_GC(kbase+j) 
            R = R_p+R_m*x       ! evaluer l'integrand a l'abscisse de G-C
            call Calc_RangVitesse(Mpart,Epart,Jpart,R,iRang,V2)
            Vr = V2-(Jpart/R)**2
#if (_SURVEILLE_SQRT_>0)
            if (Vr.LT.0.0d0) then
               write(0,*) '!!! Vr**2 = ',Vr,' !!!'
               call Meurt('Calc_Porb_SE')
            end if
#endif
            Vr=sqrt(Vr)
            dPorb = sqrt(1.0d0-x*x)/Vr ! c'est l'integrand divise par la fct poids
            Porb = Porb + dPorb        ! W(x)=1/sqrt(1-x**2)
         end do
         Porb = Porb/dble(n)    ! on multiplie par le poids constant
      end do
                                ! Normalisations
      Porb = r2Pi*R_m*Porb
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_GC
c----------------------------------------------------------------------
      include 'EvolAmas_IntegrOrb_inc.f'
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer j,k,n
c......................................................................
c
      k = 0
      do n=1,n_GC_max
         do j=1,n
            k=k+1
            x_GC(k) = cos(rPi*(dble(j)-0.5d0)/dble(n))
         end do
      end do
c
      end
c    
c----------------------------------------------------------------------
      subroutine Calc_RangVitesse(Mpart,Epart,Jpart,R, iRang,V2)
c----------------------------------------------------------------------
      include 'EvolAmas_IntegrOrb_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $     Mpart,Epart,Jpart,R  ! En entree
      integer
     $     iRang                ! En sortie
      double precision 
     $     V2          
c
c var locales :
c ^^^^^^^^^^^^^
      double precision A,B
c......................................................................
c
      call Determiner_PotRang_Arbre_r(R,A,B,iRang) ! parcours de l'arbre
      V2  = 2.0d0*(Epart-B)-2.0d0*(A-FracMasseInt*Mpart)/R
c
      end
c
c----------------------------------------------------------------------
      block data BD_Init_IntegrOrb
      include 'EvolAmas_IntegrOrb_GC_inc.f'
      data x_GC /k_GC_max*-2.0d0/
      end
c----------------------------------------------------------------------
