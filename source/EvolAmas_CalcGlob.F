#define _VERB_RCOEUR_ 0
#define _TYPE_CALC_RCOEUR_ 3

#if (_TYPE_CALC_RCOEUR_==1)
c----------------------------------------------------------------------
      subroutine Calc_Rcoeur(Rc,prec)
c----------------------------------------------------------------------
c     l'amas doit etre trie !!
c
      include 'EvolAmas_CalcGlob_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision Rc,prec  ! prec est la predcision relative a atteindre
c
c constantes :
c ^^^^^^^^^^^^
      integer Niter_max
      parameter (Niter_max=20)
      double precision fMtot_Mc_ini
      parameter (
     $     fMtot_Mc_ini=0.1d0  ! la fraction de la masse totale definissant initialement le coeur
     $     )
c
c var locales :
c ^^^^^^^^^^^^^
      double precision Rc_prec,Mint,Mc,Mc_prec
      integer iter,iRang,iSE
c......................................................................

c
c---- Valeur initiale du rayon de coeur
c
      if (Rc.LE.0.0d0) then     ! si la routine est appelee avec une valeur pour Rc,
                                ! on prend celle-ci comme valeur initiale
         Mint = 0.0d0
         Mc = fMtot_Mc_ini*M_amas
         iRang = 0
         do while (Mint.LT.Mc) 
            iRang = iRang+1
            iSE = iListeRay(iRang)
            Mint = Mint+M_Se(iSE)
         end do
         Rc = 0.5d0*(R_SE(iSE)+R_SE(iListeRay(iRang-1)))
      end if
c
c---- Calcul iteratif du rayon de coeur
c
      Rc_prec = 0.0d0
      iter = 0
#if (_VERB_RCOEUR_>0)
         write(0,100) iter,Rc
 100     format(' > Calc_Rcoeur : iter = ',I2,' R_coeur = ',E12.4)
#endif
      do while (abs(Rc-Rc_prec).GT.prec*0.5d0*(Rc+Rc_prec))
         iter = iter+1
         if (iter.GT.Niter_max) then
            write(0,*) '!!! Nb maximal d''iterations atteint dans',
     $           ' Calc_Rcoeur !!!'
            return
         end if
         Mc_prec = Mc
         Rc_prec = Rc
         call Raffine_Rcoeur(Mc_prec,Rc,Mc)
#if (_VERB_RCOEUR_>0)
         write(0,100) iter,Rc
#endif
      end do
c
      end
c
c------------------------------------------------------------
      subroutine Raffine_Rcoeur(Mc_prec,Rc,Mc)
c------------------------------------------------------------
      include 'EvolAmas_CalcGlob_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision Mc_prec,Rc,Mc
c
c constantes :
c ^^^^^^^^^^^^
      double precision TroisSur4Pi,fMc_centre
      parameter (
     $     TroisSur4Pi = 3.0d0/(4.0d0*rPi),
     $     fMc_centre=0.01d0     ! la fraction de la masse du coeur definissant le "centre"
     $     )
c
c var locales :
c ^^^^^^^^^^^^^
      double precision Tc,M_centre,Sigma2_0,Rho4Pi3_0
      integer iRang,iSE
c......................................................................
c
      M_centre = fMc_centre*Mc_prec
      Tc = 0.0d0
      Mc = 0.0d0
      iRang = 0
      do while (Mc.LT.M_centre)
         iRang = iRang+1
         iSE =iListeRay(iRang)
         Mc = Mc+M_SE(iSE) 
         Tc = Tc+M_SE(iSE)*T_SE(iSE)
      end do
      Sigma2_0 = 2.0d0*Tc/Mc    ! dispertion centrale des vitesses
      Rho4Pi3_0 = Mc/(R_SE(iSE)**3) ! densite*4*Pi/3
      Rc = sqrt(Sigma2_0/Rho4Pi3_0) ! Nouvelle valeur du rayon de coeur

      if (Rc.GT.R_SE(iListeRay(NbSE-NbSE_perte))) then
         write(0,*) '!!! Rayon de coeur (',Rc,') > Taille amas (',
     $        R_SE(iListeRay(NbSE-NbSE_perte)),') !!!'
         Mc = M_amas
      else
                                ! Determination de la masse du coeur, c'est a dire
                                ! la masse contenue dans le rayon de coeur
         do while (R_SE(iSE).LT.Rc)
            iRang = iRang+1
            iSE =iListeRay(iRang)
            Mc = Mc+M_SE(iSE) 
         end do
      end if

      Rho_ctr = TroisSur4Pi*Rho4Pi3_0
      Sigma3D_ctr = sqrt(Sigma2_0)
#if (_VERB_RCOEUR_>1)
      write(0,100) Rc,Mc,Sigma3D_ctr,Rho_ctr
 100  format('Rc,Mc,sigma_3D(0),Rho(0) ='4(' ',E12.4))
#endif
         
c
      end
c
#elif (_TYPE_CALC_RCOEUR_==2)
c----------------------------------------------------------------------
      subroutine Calc_Rcoeur(Rc,prec)
c----------------------------------------------------------------------
c     l'amas doit etre trie !!
c
      include 'EvolAmas_CalcGlob_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision Rc,prec  ! prec est la precision relative a atteindre
c
c constantes :
c ^^^^^^^^^^^^
      integer NSE_min
      parameter (NSE_min=20)
      double precision TroisSur4Pi
      parameter (TroisSur4Pi = 3.0d0/(4.0d0*rPi))
      double precision fRi
      parameter ( fRi = 5.0d0 )
c
c var locales :
c ^^^^^^^^^^^^^
      double precision Mint,Tint,Rc_prec
      integer iRang,iSE
c......................................................................
c
      Mint = 0.0d0
      Tint = 0.0d0
      do iRang=1,NSE_min
         iSE=iListeRay(iRang)
         Mint=Mint+M_SE(iSE)
         Tint=Tint+M_SE(iSE)*T_SE(iSE)
      end do

      Rc_prec = 0.0d0
      do iRang=NSE_min+1,NbSE-NbSE_perte
         Rc = sqrt(2.0d0*Tint*R_SE(iSE)**3/Mint**2)
         if ( Rc.LT.fRi*R_SE(iSE) .AND.
     $        (abs(Rc-Rc_prec).LE.prec*0.5d0*(Rc+Rc_prec)) ) goto 9
         Rc_prec = Rc
         iSE=iListeRay(iRang)
         Mint=Mint+M_SE(iSE)
         Tint=Tint+M_SE(iSE)*T_SE(iSE)
      end do

      write(0,*) '!!! Impossible de calculer R_coeur !!!'
      call Meurt('Calc_Rcoeur')
      return

 9    continue                  ! sortie normale
      Rho_ctr = TroisSur4Pi*Mint/R_SE(iSE)**3
      Sigma3D_ctr = sqrt(2.0d0*Tint/Mint)
c
      end
#elif (_TYPE_CALC_RCOEUR_==3)
c
c----------------------------------------------------------------------
      subroutine Calc_Rcoeur(Rc,prec)
c----------------------------------------------------------------------
c     l'amas doit etre trie !!
c
      include 'EvolAmas_CalcGlob_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision Rc,prec  ! prec est la predcision relative a atteindre
c
c constantes :
c ^^^^^^^^^^^^
      integer NSE_min_def, NSE_max_def
      parameter (NSE_min_def=50, NSE_max_def=500) 
                                ! nb min/max de particules intervenant dans la determination
                                ! des proprietes centrales
      double precision fMtot_Mc
      parameter (fMtot_Mc=0.001d0) ! la fraction de la masse totale definissant la region centrale
      double precision TroisSur4Pi
      parameter (TroisSur4Pi = 3.0d0/(4.0d0*rPi))
c
c var locales :
c ^^^^^^^^^^^^^
      double precision Tint,Mint,Mlim
      integer iRang,iSE,NSE_min
c......................................................................
c
      iSE=0 ! useless
      Mint = 0.0d0
      Tint = 0.0d0
      iRang = 0
      Mlim = fMtot_Mc*M_amas
      NSE_min = min(NSE_min_def,(NbSE-NbSE_perte)/2+1)
      do while (iRang.LT.NSE_min .OR.
     $     (Mint.LT.Mlim .AND. iRang.LT.NSE_max_def))
         iRang = iRang+1
         iSE = iListeRay(iRang)
         Mint=Mint+M_SE(iSE)
         Tint=Tint+M_SE(iSE)*T_SE(iSE)
      end do
      Rho_ctr = TroisSur4Pi*Mint/R_SE(iSE)**3
      Sigma3D_ctr = sqrt(2.0d0*Tint/Mint)
      Rc = sqrt(2.0d0*Tint*R_SE(iSE)**3/Mint**2)
c
      end
c
#endif /* if (_TYPE_CALC_RCOEUR_==1) */
c
c----------------------------------------------------------------------
      subroutine Calc_Rhalf(Rh)
c----------------------------------------------------------------------
c     l'amas doit etre trie !!
c
      include 'EvolAmas_CalcGlob_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision Rh
c
c var locales :
c ^^^^^^^^^^^^^
      integer i,iSE
      double precision Mint,Mdemi
c......................................................................
c
      Mdemi=0.5d0*M_amas
      Mint=0.0d0
      i=0
      iSE=0 ! useless
      do while (Mint.LT.Mdemi)
         i=i+1
         iSE=iListeRay(i)
         Mint=Mint+M_SE(iSE)
      end do
      Rh = 0.5d0*(R_SE(iSE)+R_SE(iListeRay(i-1)))
c
      end
c
