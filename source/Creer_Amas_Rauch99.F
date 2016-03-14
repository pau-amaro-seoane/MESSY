
c======================================================================
c     Create initial cluster file for MC simulation according to 
c     cusp model from Rauch 99
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
      include 'Creer_Amas_Rauch99_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      character*20 arg
      integer iarg,narg,iseed
c
c fonctions :
c ^^^^^^^^^^^
      integer iargc, LongChaine
c......................................................................
c
      M_TN = -1.0d0
      NbSE = -1
      NomAmas = ''
      iseed = 67654
      expo_n = 7.0d0/4.0d0
      p_anis = 0.0d0
      a_min=1.0d-20
      a_max=1.0d0
      m_special=0.0d0
      a_special=-1.0d0
      e_special=-1.0d0
      iCorrec_SelfGrav = 2
c
c---- lecture des options sur la ligne de commande
c
      iarg = 0
      narg = iargc()
      do while(iarg.LT.narg)
         iarg = iarg+1
         call getarg(iarg,arg)
         if (arg(1:LongChaine(arg)).EQ.'-M_BH' .OR.
     $        arg(1:LongChaine(arg)).EQ.'-M_TN'.OR.
     $        arg(1:LongChaine(arg)).EQ.'-mu' ) then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) M_TN
         else if (arg(1:LongChaine(arg)).EQ.'-n' .OR.
     $           arg(1:LongChaine(arg)).EQ.'-gamma') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) expo_n
         else if (arg(1:LongChaine(arg)).EQ.'-p') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) p_anis
         else if (arg(1:LongChaine(arg)).EQ.'-a_min') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) a_min
         else if (arg(1:LongChaine(arg)).EQ.'-iseed') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) iseed
         else if (arg(1:LongChaine(arg)).EQ.'-NbSE' .OR.
     $           arg(1:LongChaine(arg)).EQ.'-N') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) NbSE
         else if (arg(1:LongChaine(arg)).EQ.'-NomAmas' .OR.
     $           arg(1:LongChaine(arg)).EQ.'-Name') then
            iarg = iarg+1
            call getarg(iarg,NomAmas)
         else if (arg(1:LongChaine(arg)).EQ.'-CorrectSelfGrav') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) iCorrec_SelfGrav
         else if (arg(1:LongChaine(arg)).EQ.'-SpecialStar') then 
                                ! Put a "special star" in orbit around MBH
                                ! with given mass, semi-major axis, eccentricity
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) m_special
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) a_special
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) e_special
            
         else
            write(0,*) '!!! invalid option ',arg(1:LongChaine(arg)),
     $           ' !!!'
            call Meurt('Saisir_Parametres')
         end if
      end do
c
c---- saisie interactive des parametres manquants
c
      if (M_TN.EQ.-1.0d0) then
         write(0,*) ' >  BH mass (>>1): '
         read(*,*) M_TN
      end if
      if (M_TN.LT.3.0d0) then ! quite arbirary...
         write(0,*)
     $        '!!! Warning: A large positive BH mass should be used !!!'
      end if
      M_ctr = M_TN
      if (p_anis.LT.0.0d0) then
         write(0,*) '!!! p_anis >= 0 should be used !!!'
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
     $     '###################################################',/,
     $     '## Parameters Rauch99 cluster (Creer_Amas_Rauch) ##',/,
     $     '###################################################')
      write(0,100) NomAmas(1:LongChaine(NomAmas)),
     $     NbSE,M_TN,p_anis,expo_n,a_min,a_max
 100  format(
     $     '> NomAmas   : ',A,/,
     $     '> NbSE      : ',I7,/,
     $     '> M_TN      : ',E12.5,/,
     $     '> p_anis    : ',E12.5,/,
     $     '> expo_n    : ',E12.5,/,
     $     '> a_min     : ',E12.5,/,
     $     '> a_max     : ',E12.5
     $     )
      if (a_special.GT.0.0d0) then
         if (m_special.LE.0.0d0) m_special=1.0d0/NbSE
         write(0,110) m_special,a_special,e_special
 110     format(
     $        '-- Special particle: ',/,
     $        ' > m_special :',E12.5,/,
     $        ' > a_special :',E12.5,/,
     $        ' > e_special :',E12.5
     $        )
      end if

      call init_random(iseed)
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_Amas
c----------------------------------------------------------------------
      include 'Creer_Amas_Rauch99_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      integer iSE,i,iseed,Nint,Npart,itry
      double precision amin_3n,amax_3n,One3n,MOneP,OneEp, e_mean,
     $     semi_a,Xrand,Mpart, ecc_ano, Rad,
     $     renR,renE,renJ,renV,e1,e2, M_int,
     $     Egrav, Ekin, PotShift, Energy, CoefPot, amax_over_amin
      logical GammaIs2, GammaIs3

      double precision mean_anom, ecc, newT
      common /param_f_ecc_anom/ mean_anom, ecc
c
c constants :
c ^^^^^^^^^^^
      double precision e_max, e_min
      parameter (e_min=1.0d-6, e_max=1.0d0-e_min)
      double precision prec_ecc_ano
      parameter (prec_ecc_ano=1.0d-6)
      integer Nmax_pos
      parameter (Nmax_pos=100)
      double precision tiny
      parameter (tiny=1.0d-10)
      double precision eps_bound
      parameter (eps_bound=1.d-12)
      integer ntry_max
      parameter (ntry_max=10 000)
c
c functions :
c ^^^^^^^^^^^
      double precision random, rtsafe
      external f_ecc_anom
c......................................................................
c
      Mpart=1.0d0/NbSE
      Npart=NbSE
      iSE=0

      if (a_special.GT.0.0d0) then ! Special particle first
         if (m_special.LE.0.0d0) m_special=Mpart
         iSE=1
         semi_a = a_special
         ecc = e_special
         mean_anom = r2Pi*random()
         e1=0.0d0
         e2=r2Pi
         ecc_ano=rtsafe(f_ecc_anom,e1,e2,prec_ecc_ano)
         Rad=semi_a*(1.0d0-ecc*cos(ecc_ano))
         R_SE(iSE) = Rad
         M_SE(iSE) = m_special 
         J_SE(iSE) = sqrt(M_TN*semi_a*(1.0d0-ecc**2))
         T_SE(iSE) = M_TN*(1.0d0/Rad-0.5d0/semi_a)
         
                                ! Masses of other particles
         Mpart=(1.0d0-m_special)/(NbSE-1)
         Npart=NbSE-1
      end if

      ! Values of the exponent requiring special care

      GammaIs3 = (abs(expo_n-3.d0).LT.1d-3)
      GammaIs2 = (abs(expo_n-2.d0).LT.1d-3)
      
      if (GammaIs3) then
         amax_over_amin = a_max/a_min
      else
         amin_3n=a_min**(3.0d0-expo_n)
         amax_3n=a_max**(3.0d0-expo_n)
         One3n = 1.0d0/(3.0d0-expo_n)
      end if
      if (p_anis.GT.tiny) then
         MOneP=-1.0d0/p_anis
         OneEp=1.0d0-exp(-p_anis)
      end if
      Nint=0
      if (GammaIs3) then
         CoefPot = 1.0d0/log(amax_over_amin)
         PotShift = CoefPot*(1.0d0/a_max-1.0d0/a_min)
      else
         if (GammaIs2) then
            PotShift = -1.0d0/a_max
         else
            CoefPot = 1.0/((2.0-expo_n)*a_max)
            PotShift = -(3.0d0-expo_n)*CoefPot
         end if
      end if

      e_mean = 0.0d0
      do while (iSE.LE.NbSE)
         Xrand = random()
         if (GammaIs3) then
            semi_a = a_min*amax_over_amin**Xrand
         else
            semi_a = (Xrand*amax_3n+(1.0d0-Xrand)*amin_3n)**One3n
         end if
         M_int = M_TN
         if (iCorrec_SelfGrav.GE.1) then
                                ! Naive "correction" for interior stellar mass
            M_int = M_TN+Xrand  ! Xrand turns out to be the interior stellar 
                                ! mass (in units where the totalk stellar mass is one)
         end if
         Xrand = random()
         if (p_anis.GT.tiny) then
            ecc = sqrt(MOneP*log(1.0d0-Xrand*OneEp))
         else
            ecc = sqrt(Xrand)
         end if
c         if (ecc.GT.e_min .AND. ecc.LT.e_max .AND.
c     $        semi_a*(1.0d0-ecc).LT.1.0d0) then
         if (ecc.GT.e_min .AND. ecc.LT.e_max) then

            do itry=1,ntry_max
c
c---- Pick position on orbit
c
               mean_anom = r2Pi*random()
               e1=0.0d0
               e2=r2Pi
               ecc_ano=rtsafe(f_ecc_anom,e1,e2,prec_ecc_ano)
               Rad=semi_a*(1.0d0-ecc*cos(ecc_ano))
c
c--- Check if bound (using approximate grav energy)
c
               Egrav = -M_TN/Rad
               if (Rad.LE.a_max) then
                  if (GammaIs3) then
                     if (Rad.LE.a_min) then
                        Egrav = Egrav+PotShift
                     else
                        Egrav = Egrav+CoefPot*(1.0d0/a_max
     $                       - 1.0d0/Rad*(log(Rad/a_min)+1.0d0))
                     end if
                  else
                     Egrav = Egrav + PotShift
                     if (.NOT.GammaIs2) then
                        Egrav = Egrav +
     $                       CoefPot*(Rad/a_max)**(2.0-expo_n)
                     else
                        Egrav = Egrav + log(Rad/a_max)/a_max
                     end if
                  end if
               else
                  Egrav = Egrav - 1.0d0/Rad
               end if
               Ekin = M_int*(1.0d0/Rad-0.5d0/semi_a)
               Energy = Ekin+Egrav
               if (Energy.LT.0.0d0 .OR.iCorrec_SelfGrav.LT.2) goto 99
            end do

            write(0,*) '!!! Couldn''t find a bound orbit with a,e = ',
     $           semi_a,ecc,' !!!'
            call exit(1)

 99         continue
            e_mean = e_mean + ecc

            ! mean eccentricity is supposed to be 
            ! [ sqrt(pi)/2*erf(sqrt(p))/sqrt(p) - exp(-p) ]/(1-exp(-p))

c            write (66,'(4(E12.5," "),I5)') semi_a,ecc,Rad,Egrav,itry !#############
            if (Rad.LT.a_max) Nint=Nint+1
            iSE = iSE+1
            R_SE(iSE) = Rad
            M_SE(iSE) = Mpart 
            J_SE(iSE) = sqrt(M_int*semi_a*(1.0d0-ecc**2))
            T_SE(iSE) = Ekin
         end if
      end do
      e_mean = e_mean/dble(NbSE)
      write(0,*) '> mean eccentricity : ',e_mean
c
c---- Transform lengths to N-body units to have stellar pot energy = -1/2
c
      call TriAmas
      call Calc_EnerGlob
      renR=-2.0d0*Epot_st
      write(0,110) renR,1.0d0*Nint/dble(NbSE)
 110  format(
     $     '> a_max in N-body units       :',E12.5,/,
     $     '> Mass frac interior to a_max :',E12.5)
      renE=1.0d0/renR
      renJ=sqrt(renR)
      do iSE=1,NbSE
         R_SE(iSE)=renR*R_SE(iSE)
         J_SE(iSE)=renJ*J_SE(iSE)
         T_SE(iSE)=renE*T_SE(iSE)
      end do
      Epot_st = Epot_st/renR
      Epot_bh = Epot_bh/renR
      Ecin = renE*Ecin
c
c---- DO NOT adapt velocities to enforce virial equilibrium
c     otherwise some particle close to the BH may end up with wrong semi-major axis
c     The procedure to create the cluster (assuming Keplerian orbits) assures
c     statistical dynamical equilibrium if the total stellar mass is much lower than
c     the mass of the central thingy
c
      call Calc_EnerGlob ! Just for information
      renV = sqrt(-0.5d0*(Epot_st+Epot_bh)/Ecin)
c$$$      do iSE=1,NbSE
c$$$         newT = renV**2*T_SE(iSE)
c$$$         J_SE(iSE)=renV*J_SE(iSE)
c$$$         T_SE(iSE)=newT
c$$$      end do
c$$$      Ecin = renV**2*Ecin
      write(0,100) renV
 100  format(
     $     '> factor on V to enforce vir. equil.: ',E12.5,
     $     ' (not applied)')
c$$$      call Calc_EnerGlob
         
c
      end
c
c----------------------------------------------------------------------
      subroutine f_ecc_anom(eta,f,df)
c----------------------------------------------------------------------
      implicit none

      double precision eta,f,df

      double precision mean_anom, ecc
      common /param_f_ecc_anom/ mean_anom, ecc
c......................................................................
c
      f=mean_anom-eta+ecc*sin(eta)
      df=ecc*cos(eta)-1.0d0
c
      end
c
c----------------------------------------------------------------------
      subroutine Calc_EnerGlob
c----------------------------------------------------------------------
      include 'Creer_Amas_Rauch99_inc.f'
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer i,iSE,iFirst
      double precision M_i,R_i,Mint,Mprec
c......................................................................
c
      Ecin = 0.0d0
      Epot_st = 0.0d0
      Epot_bh = 0.0d0
      Mint = 0.0d0
      Mprec = 0.0d0
      iFirst=1
                                ! do not take the "special particle"
                                ! into consideration to compute energies
      if (a_special.GT.0.0d0) then
         iFirst=2
      end if
      do i=iFirst,NbSE
         iSE = iListeRay(i)
         M_i = M_SE(iSE)
         R_i = R_SE(iSE)
         Ecin = Ecin + M_i*T_SE(iSE)
         Mint = Mint+Mprec
         Epot_st = Epot_st - M_i/R_i*(Mint+FracMasseInt*M_i)
         Epot_bh = Epot_bh - M_i*M_ctr/R_i
         Mprec = M_i
      end do
      write(0,100) Ecin,Epot_st,Epot_bh,
     $     (Epot_st+Epot_bh)/(Epot_st+Epot_bh+Ecin),
     $     (Epot_bh)/(Epot_bh+Ecin)
 100  format(
     $     '> Energies : ',/,
     $     '           Ecin    = ',E12.5,/,
     $     '           Epot_st = ',E12.5,/,
     $     '           Epot_bh = ',E12.5,/,
     $     '        Egrav/Etot = ',E12.5,/,
     $     '  Egrav_BH/Etot_BH = ',E12.5
     $     )
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_Fichiers
c----------------------------------------------------------------------
      include 'Creer_Amas_Rauch99_inc.f'
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
