      program Creer_Amas_EtaTN

#define _AJUSTE_VIR_  /* Si definie, "corrige" les vitesses pour forcer la virielisation */
#define _LAXISTE_     /* Si definie, laisse passer un amas mal virielise */
#define _MODE_TIRAGE_POS_ 1 /* 1 : SE regulierement espacees en M(r) */
                            /* 2 : tirage R selon rho(R)             */
#define _INTEGRER_POISSON_ 1 /* ne pas toucher ! */

c======================================================================
c     Create initial cluster file for MC simulation according to 
c     eta-model with central BH. Cf Tremaine et al. (94), AJ 107, 634 
c======================================================================
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
      call Ecrire_Fichiers
c
      end
c
c======================================================================
c
c----------------------------------------------------------------------
      subroutine Saisir_Parametres
c----------------------------------------------------------------------
c
      include 'Creer_Amas_EtaTN_inc.f'
      include 'EvolAmas_LimExt_common.f'
c
c var locales :
c ^^^^^^^^^^^^^
      character*20 arg
      integer iarg,narg
      double precision my_eta,my_mu
c
c fonctions :
c ^^^^^^^^^^^
      integer iargc, LongChaine
      double precision Sigma1D_of_R
c......................................................................
c
      my_eta = -1.0d0
      NbSE = -1
      NomAmas = ''
      iseed = 67654
      Rmin_cutoff = 0.0d0       ! No minimum radius per default
      Mext_LimExt = 0.0d0       ! No radial truncation (reflective wall or such) by default
c
c---- lecture des options sur la ligne de commande
c
      iarg = 0
      narg = iargc()
      do while(iarg.LT.narg)
         iarg = iarg+1
         call getarg(iarg,arg)
         if (arg(1:LongChaine(arg)).EQ.'-M_TN' .OR.
     $        arg(1:LongChaine(arg)).EQ.'-mu') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) my_mu
         else if (arg(1:LongChaine(arg)).EQ.'-eta') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) my_eta
         else if (arg(1:LongChaine(arg)).EQ.'-iseed') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) iseed
         else if (arg(1:LongChaine(arg)).EQ.'-gamma') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) my_eta
            my_eta = 3.0d0-my_eta
         else if (arg(1:LongChaine(arg)).EQ.'-NbSE' .OR.
     $           arg(1:LongChaine(arg)).EQ.'-N') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) NbSE
         else if (arg(1:LongChaine(arg)).EQ.'-NomAmas' .OR.
     $           arg(1:LongChaine(arg)).EQ.'-Name') then
            iarg = iarg+1
            call getarg(iarg,NomAmas)
         else if (arg(1:LongChaine(arg)).EQ.'-Rmin') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) Rmin_cutoff
#if (_LIMITE_EXT_>0)
         else if (arg(1:LongChaine(arg)).EQ.'-TruncM') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) Mext_LimExt
            Mext_LimExt = 1 - Mext_LimExt ! Mass beyond the radial truncation radius
#endif
         else
            write(0,*) '!!! invalid option ',arg(1:LongChaine(arg)),
     $           ' !!!'
            call Meurt('Saisir_Parametres')
         end if
      end do
c
c---- saisie interactive des parametres manquants
c
      if (my_mu.EQ.-1.0d0) then
         write(0,*) ' >  fractional BH mass (mu): '
         read(*,*) my_mu
      end if
      if (my_mu.LT.0.0d0) then
         write(0,*) '!!! A positive (or null) BH mass should be used !!!'
         call Meurt('Saisir_Parametres')
      endif
      M_TN = my_mu
      M_ctr = my_mu

      if (my_eta.EQ.-1.0d0) then
         write(0,*) ' >  exponant (eta): '
         read(*,*) my_eta
      end if
      if (my_eta.GE.2.5d0 .OR. my_eta.LE.0.5d0) then
         write(0,*) '!!! Eta should be in ]0.5;2.5[ !!!'
         call Meurt('Saisir_Parametres')
      endif

      if (NbSE.EQ.-1) then
         write(0,*) ' >  Number of super-stars (NbSE): '
         read(*,*) NbSE
      end if
      if (NbSE.LE.0 .OR. NbSE.GT.iDimSE) then
         write(0,*) '!!! NbSE should be in [1;',iDimSE,'] !!!'
         call Meurt('Saisir_Param')
      end if

      if (NomAmas.EQ.'') then
         write(0,*) ' > Generic name for created files (NomAmas):'
         read(*,*) NomAmas
      end if

      write(0,90)
 90   format(/,
     $     '#########################################################',/,
     $     '## Parameters for eta-model cluster (Creer_Amas_EtaTN) ##',/,
     $     '#########################################################')
      write(0,100) NomAmas(1:LongChaine(NomAmas)),
     $     NbSE,my_eta,my_mu
 100  format(
     $     '> NomAmas   : ',A,/,
     $     '> NbSE      : ',I7,/,
     $     '> eta       : ',E12.5,/,
     $     '> M_TN (mu) : ',E12.5 )
      call SetParam_EtaModel(my_eta,my_mu)
c
c---- Compute various quantities of interest  
c
      Rbreak = 1.0d0/(2.0d0*my_eta-1.0d0) ! The stellar component of grav energy is -1/(2*(2*eta-1))
                                       ! See eq. 14 of Tremaine et al. (94)
      Sigma1D_break = Sigma1D_of_R(1.0d0)/sqrt(Rbreak)
      Rinfl = M_TN/Sigma1D_break**2
      write(0,120)
 120  format('>> Theoretical quantities',/,
     $       '   ----------------------')
      write(0,105) Rbreak
 105  format('> Break radius in N-body units (Rbreak): ',E12.5)
      write(0,108) Sigma1D_break
 108  format('> 1D dispersion at break radius (Sigma1D_break): ',E12.5)
      write(0,109) Rinfl
 109  format('> BH''s influence radius (Rinfl): ',E12.5)
      write(0,111) 1.0d0/2.0d0**my_eta
 111  format('> Stellar mass inside Rbreak (Mbreak): ',E12.5)
      write(0,112) (Rinfl/(1.0d0+Rinfl))**my_eta
 112  format('> Stellar mass inside Rinfl (Minfl): ',E12.5,/)
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_Amas
c----------------------------------------------------------------------
      include 'Creer_Amas_EtaTN_inc.f'
      include 'EvolAmas_LimExt_common.f'
      include 'EtaModel_common.f'
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision UnSurRac2
      integer iSE, itick, inext_tick
      double precision mi,Mr,R,V,F_VN,V2,f_V,VradSE,VtgSE,
     $     renR,renE,renJ,renV,Psi,Fmax,Vmax,DF_max
c
c functions :
c ^^^^^^^^^^^
      double precision random, R_of_Mr, Psi_of_R, DstFct_tbl
c......................................................................
c
c---- Global properties of particles beyong truncation radius
c
#if (_LIMITE_EXT_>0)
      MTrunc = 1.0d0 - Mext_LimExt
      if (Mext_LimExt) then
         R_LimExt = R_of_Mr(MTrunc) ! Truncation radius
         fNbEt_ext_LimExt = Mext_LimExt
         B_LimExt = Psi_of_R(R_LimExt)
         Egrav_ext_LimExt = XXXX ! Not including mass within truncation radius
         Ecin_ext_LimExt =  XXXX
      end if
#endif

      UnSurRac2 = 1.0d0/sqrt(2.0d0)
      write(0,90)
 90   format(/,
     $     '###########################',/,
     $     '## Creating cluster..... ##',/,
     $     '###########################',//,
     $     '0    10   20   30   40   50   60   70   80   90  100%',/,
     $     '|----|----|----|----|----|----|----|----|----|----|',/,
     $     $)
c
      call Init_random(iseed)
      M_amas = 1.0d0
      mi = 1.0d0/NbSE
      inext_tick=0
      itick=0
      do iSE=1,NbSE
         if (iSE.GE.inext_tick) then
            write(0,95)
 95         format('^',$)
            itick=itick+1
            inext_tick = 0.02d0*NbSE*itick
         end if
            
         M_SE(iSE)=mi
c
c---- Position
c
         R=-1.0d0
         do while (R.LT.Rmin_cutoff) 
                                ! pick random fractional mass
            Mr=random()
                                ! invert Mr(R) to find radius
            R = R_of_Mr(Mr)
         end do

         R_SE(iSE) = R
c
c---- Velocity
c
         Psi = Psi_of_R(R)
         DF_max = DstFct_tbl(Psi) 
         if (mu.GT.0.0d0 .AND. eta.GT.1.5d0 .AND. Psi.GT.0.75d0*e_max) then ! DF(eps) is not monotonically increasing
            DF_max = 1.2d0*f_max ! f_max is the max value of DF encountered during its tabulation
         end if
         Fmax = Psi * DF_max
c                \       \______________valeur max de la FD en R_i
c                 \
c                  \____________________valeur max de 0.5*V**2 (Eps=0)
c
         Vmax = sqrt(2.0*Psi)
 2       continue
         V  = Vmax*random()
         F_VN = Fmax*random()
         V2 = 0.5d0*V*V
         f_V = V2*DstFct_tbl(-V2+Psi)

         if (f_V.GT.Fmax) then
            write(0,*) '!!! f_V(=',f_V,') > Fmax(=',Fmax,') !!!'
            call Meurt('Calc_Amas')
         end if

         if (F_VN.GT.f_V) then
            goto 2
         end if

 3       continue
         VradSE = V*(1.0d0-2.0d0*random())
	 VtgSE = sqrt(2.0d0*V2-VradSE**2)
                                ! to avoid circular/radial orbits...
         if ( (VtgSE .LT.1.0d-6*V) .OR. (VradSE.LT.1.0d-6*V) ) then
            goto 3
         endif
         J_SE(iSE) = R*VtgSE
         T_SE(iSE) = V2 
      end do
      write(0,*) ''
c
c---- Adapt length scale to get total stellar gravitational energy = -1/2 (N-body units)
c
      call TriAmas
      call Calc_EnerGlob
      renR=-2.0d0*Epot_st
      renE=1.0d0/renR
      renJ=renR*sqrt(renE)
      write(0,110) Rbreak,renR
 110  format(
     $     '> factor on R to enforce stellar self-Egrav = -0.5: ',/,
     $     '                        theoretical : ',E12.5,/,
     $     '                        practical   : ',E12.5)
      
      do iSE=1,NbSE
         R_SE(iSE)=renR*R_SE(iSE)
         J_SE(iSE)=renJ*J_SE(iSE)
         T_SE(iSE)=T_SE(iSE)/renR
      end do
      Epot_st = Epot_st/renR
      Epot_bh = Epot_bh/renR
      Ecin = Ecin/renR

#if (_LIMITE_EXT_>0)
      if (Mext_LimExt.GT.0.0d0) then
         R_LimExt = R_LimExt/renR
         Egrav_ext_LimExt = Egrav_ext_LimExt/renR
         B_LimExt = B_LimExt/renR
         Ecin_ext_LimExt = Ecin_ext_LimExt/renE
      end if
#endif
c
c---- Adapt velocities to enforce virial equilibrium
c
      call Calc_EnerGlob

      renV = sqrt(-0.5d0*(Epot_st+Epot_bh)/Ecin)
      do iSE=1,NbSE
         J_SE(iSE)=renV*J_SE(iSE)
         T_SE(iSE)=renV**2*T_SE(iSE)
      end do
      Ecin = renV**2*Ecin

#if (_LIMITE_EXT_>0)
      if (Mext_LimExt.GT.0.0d0) then
         Ecin_ext_LimExt = Ecin_ext_LimExt*renV**2
      end if
#endif

      write(0,100) renV
 100  format(
     $     '> factor on V to enforce vir. equil.: ',E12.5)
      call Calc_EnerGlob
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_EnerGlob
c----------------------------------------------------------------------
      include 'Creer_Amas_EtaTN_inc.f' 
      include 'EvolAmas_LimExt_common.f'
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer i,iSE
      double precision M_i,R_i,Mint,Mprec
c......................................................................
c
      Ecin = 0.0d0
      Epot_st = 0.0d0
      Epot_bh = 0.0d0
      Mint = 0.0d0
      Mprec = 0.0d0
      do i=1,NbSE
         iSE = iListeRay(i)
         M_i = M_SE(iSE)
         R_i = R_SE(iSE)
         Ecin = Ecin + M_i*T_SE(iSE)
         Mint = Mint+Mprec
         Epot_st = Epot_st - M_i/R_i*(Mint+FracMasseInt*M_i)
         Epot_bh = Epot_bh - M_i*M_ctr/R_i
         Mprec = M_i
      end do

#if (_LIMITE_EXT_>0)
      if (Mext_LimExt.GT.0.0d0) then ! contribution of unmodelled stars (beyond the truncation)
         Ecin = Ecin + Ecin_ext_LimExt
         Epot_st = Epot_st + Egrav_ext_LimExt - Mint*B_LimExt
         Epot_bh = Epot_bh -M_ctr*B_LimExt
      end if
#endif
         
      write(0,100) Ecin,Epot_st,Epot_bh
 100  format(
     $     '> Energies : ',/,
     $     '           Ecin    = ',E12.5,/,
     $     '           Epot_st = ',E12.5,/,
     $     '           Epot_bh = ',E12.5
     $     )
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_Fichiers
c----------------------------------------------------------------------
      include 'Creer_Amas_EtaTN_inc.f' 
      include 'fxdr.inc'
c
c var locales :
c ^^^^^^^^^^^^^
      integer i,iXDR, ixdr_err, l
      double precision d
      character*256 NomFich
c
c fonctions :
c ^^^^^^^^^^^
      integer LongChaine
c......................................................................
c
      write(0,90)
 90   format(/,
     $     '###########################',/,
     $     '## Writing files...      ##',/,
     $     '###########################')

      l = LongChaine(NomAmas)
      if (l.GT.3 .AND. NomAmas(l-3:l).EQ.'.xdr') then
         l = l-4
      end if
      NomFich = NomAmas(1:l) // '.xdr'
      iXDR = initxdr(NomFich,'w',.TRUE.)
      call EcritAmas_XDR(iXDR)
      ixdr_err = ixdrclose(iXDR)
      if (M_TN.GT.0.0d0) then
         NomFich = NomAmas(1:l) // '_TN.xdr'
         iXDR = initxdr(NomFich,'w',.TRUE.)
         call EcritTN_XDR(iXDR)
         ixdr_err = ixdrclose(iXDR)
      end if
c 
      end
c
      

