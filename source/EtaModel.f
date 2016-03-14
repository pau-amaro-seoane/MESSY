c======================================================================
c     Properties of eta models, see Tremaine et al. 1994, AJ 107, 634
c
c     Use Eddington's formula to integrate the distribution function 
c     DstFct(eps) for eta-models with central BH
c======================================================================

c----------------------------------------------------------------------
      subroutine Init_DstFct_EtaModel
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
      include 'Table_DynFricCoefEtaModel_inc.f'
c......................................................................
c
      lTabDFCoef = .FALSE.
      if (mu.GT.0.0d0) then
         psi_max = 1.0d50
      else
         if (eta.GT.1.0d0) then
            psi_max = 1.0d0/(eta-1.0d0)
         else
            psi_max = 1.0d50
         end if
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine SetParam_EtaModel(my_eta,my_mu)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
      double precision my_eta,my_mu
c......................................................................
c
      mu = my_mu
      eta = my_eta
      call Init_DstFct_EtaModel
c
      end
c
c----------------------------------------------------------------------
      subroutine SetParam_EtaModelDrop(my_R,my_exp)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
      double precision my_R,my_exp
c......................................................................
c
      R_drop = my_R
      exp_drop = my_exp
c
      end
c
c----------------------------------------------------------------------
      double precision function d2RhodPsi2(u)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision u ! u = 1/r
c......................................................................
c
      if (u.LT.0.0d0) then
         write(0,*) '!!! d2RhodPsi2(u) called for negative u: ',u,' !!!'
         call Meurt('d2RhodPsi2')
      end if
      d2RhodPsi2 = eta *  u**2 * (1 + u)**(-3.0d0 + eta) *
     $     ( 2.0d0*(6.0d0 + u*(8.0d0 + 3.0d0*u - (2.0d0 + u)*eta)) + 
     $     (1.0d0 + u)**eta *
     $     (12.0d0 + u*(-8.0d0 + u*(-3.0d0 + eta))*(-2.0d0 + eta))*mu )
     $     / (4.0d0*Pi*(1.0d0 + (1.0d0 + u)**eta * mu)**3)
c
      end
c
c----------------------------------------------------------------------
      subroutine h_dh_u(u,h,dh)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision u, h, dh
c
c common block :
c ^^^^^^^^^^^^^^
      double precision Psi0
      common /Param_h/ Psi0
c......................................................................
c
      if (u.LT.0.0d0) then
         write(0,*) '!!! h_dh_u(u) called for negative u: ',u,' !!!'
         call Meurt('h_dh_u')
      end if
      if (eta.EQ.1.0d0) then
         h = log(1.0d0+u) + mu*u - Psi0
      else
         h = (1.0d0-(1.0d0+u)**(1.0d0-eta))/(eta-1.0d0) + mu*u - Psi0
      end if
      dh = (1.0d0+u)**(-eta) + mu
c
      end
c
c----------------------------------------------------------------------
      double precision function U_of_Psi(Psi)
c----------------------------------------------------------------------
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision Psi
c
c common block :
c ^^^^^^^^^^^^^^
      double precision Psi0
      common /Param_h/ Psi0
c
c functions :
c ^^^^^^^^^^^
      double precision rtsafe_rel
      external h_dh_u
c......................................................................
c
      if (Psi.LT.0.0d0) then
         write(0,*) '!!! U_of_Psi(Psi) called for negative Psi: ',
     $        Psi,' !!!'
         call Meurt('U_of_Psi')
      end if
      if (Psi.GT.psi_max)
     $     then
         write(0,*) '!!! U_of_Psi(Psi) called for Psi>Psi_max: ',
     $        Psi,psi_max,' !!!'
         call Meurt('U_of_Psi')
      end if
      if (Psi.LT.1.0d-7) then
         U_of_Psi = Psi/(1.0d0+mu)
         return
      end if

      Psi0 = Psi
      if (mu.EQ.0.0d0) then
         if (eta.EQ.1.0d0) then
            U_of_Psi = exp(Psi)-1.0d0
         else
            U_of_Psi = ( 1.0d0 - (eta-1.0d0)*Psi )**(1.0d0/(1.0d0-eta))
     $           - 1.0d0
         end if
      else
         U_of_Psi = rtsafe_rel(h_dh_u,Psi/(1.0d0+mu),Psi/mu,acc_u)
      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function DstFct(eps)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision eps
c
c constants :
c ^^^^^^^^^^^
      double precision coef_DstFct0,OneOverTwoPi52
      parameter (coef_DstFct0=0.03648844592d0,
     $     OneOverTwoPi52=0.01010532601d0)
c
c functions :
c ^^^^^^^^^^^
      double precision Integrand_DstFct, gammln
      external Integrand_DstFct, midsql, midpnt
c
c common block :
c ^^^^^^^^^^^^^^
      double precision eps0
      common /Param_Integrand_DstFct/ eps0
c......................................................................
c
      if (eps.LT.0.0d0) then
         write(0,*) '!!! DstFct(eps) called for negative eps: ',
     $        eps,' !!!'
         call Meurt('DstFct')
      end if
      if (eps.GT.psi_max) then
         write(0,*) '!!! DstFct(eps) called for eps>psi_max: ',
     $        eps,psi_max,' !!!'
         call Meurt('DstFct')
      end if
                                !---------------
                                ! central BH
                                !---------------
      if (mu.GT.0.0d0) then
         if (eps.GT.100000.0d0*mu) then ! deep in the BH's well
            DstFct = 0.005052663007d0*  ! Tremaine et al (94), Eq. 48
     $           exp(gammln(4.0d0-eta)-gammln(2.5d0-eta))*
     $           eta*eps**(1.5d0-eta)/mu**(3.0d0-eta)
         else
            eps0 = eps
            call qromo(Integrand_DstFct,0.0d0,sqrt(eps),acc_DstFct,
     $           DstFct,midpnt)
         end if

                                !---------------
      else                      ! no central BH
                                !---------------
         if (eps.LT.1.0d-7) then
            DstFct = Coef_DstFct0*eta*sqrt(eps)**5 ! Tremaine et al (94), Eq. 20 lacks factor eta

         else                          ! Tremaine et al (94), Eq. 19
            if (eta.LT.1.0d0 .AND. (eps.GT.1.0d9)) then
               DstFct = OneOverTwoPi52*eta*(3.0d0-eta)*
     $              (1.0d0-eta)**((1.0d0+eta)/(1.0d0-eta))*
     $              exp( gammln(2.0d0/(1.0d0-eta))-
     $              gammln(0.5d0+2.0d0/(1.0d0-eta)) )*
     $              eps**(0.5d0*(3.0d0+eta)/(1.0d0-eta))
            else if (eta.GE.1.0d0 .AND. (eps.GT.1.0d9)) then
               DstFct = 0.5d0*OneOverTwoPi52*exp(2.0d0*eps)
            else if (eta.GT.1.0d0 .AND.
     $              ((psi_max-eps).LT.1.0d-4*psi_max) ) then
               DstFct = OneOverTwoPi52*eta*(3.0d0-eta)*
     $              (eta-1.0d0)**((1.0d0+eta)/(1.0d0-eta))*
     $              exp( gammln((eta+1.0d0)/(eta-1.0d0)-0.5d0)-
     $              gammln((eta+1.0d0)/(eta-1.0d0)) )*
     $              (psi_max-eps)**(0.5d0*(3.0d0+eta)/(1.0d0-eta))
            else
               eps0 = eps
               call qromo(Integrand_DstFct,0.0d0,sqrt(eps),acc_DstFct,
     $              DstFct,midsql)
            end if
         end if
      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function DstFct_tbl(eps)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision eps
c
c constantes :
c ^^^^^^^^^^^^
      integer Ntbl
      parameter (Ntbl=10000)
      double precision alpha_x
      parameter (alpha_x=10.0d0)
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer i
      double precision X,Xmin,dX,Y_tbl(Ntbl),Y2_tbl(Ntbl),f,eps_
      logical lprem
      save lprem,Xmin,dX,Y_tbl,Y2_tbl
      data lprem /.TRUE./
      double precision eps_min, eps_max
      save eps_min, eps_max
      integer itick, inext_tick
c
c common block :
c ^^^^^^^^^^^^^^
      double precision eps0
      common /Param_Integrand_DstFct/ eps0
c
c functions :
c ^^^^^^^^^^^
      double precision DstFct
      double precision z, atanh
      atanh(z) = 0.5d0*log((1.0d0+z)/(1.0d0-z))
c......................................................................
c
      if (eps.LT.0.0d0) then
         write(0,*) '!!! DstFct_tbl(eps) called for negative eps: ',
     $        eps,' !!!'
         call Meurt('DstFct_tbl')
      end if

      if (eps.GT.psi_max) then
         DstFct_tbl = 0.0d0
         return
      end if
      
      if (lPrem) then           ! compute spline arrays
         f_max = -1.0d30

         write(0,90)
 90      format(/,
     $        '##################################################',/,
     $        '## Tabulating dist. function for eta model..... ##',/,
     $        '##################################################',//,
     $        '0    10   20   30   40   50   60   70   80   90  100%',/,
     $        '|----|----|----|----|----|----|----|----|----|----|',/,
     $        $)
         inext_tick=0
         itick=0

         if (psi_max.GT.1.0d30)
     $        then              ! eps in ]0,infinity[
            eps_max = 1.0d6
            eps_min = 1.0d-3
            Xmin=log(eps_min)
            dX=log(eps_max/eps_min)/(Ntbl-1)
            do i=1,Ntbl
               eps_ = eps_min*exp((i-1)*dX)
               f = DstFct(eps_)
               if (f.GT.f_max) then
                  f_max = f
                  e_max = eps_
               end if
               Y_tbl(i) = log(f)

               if (i.GE.inext_tick) then
                  write(0,95)
 95               format('^',$)
                  itick=itick+1
                  inext_tick = 0.02d0*Ntbl*itick
               end if

            end do
         else                   ! eps in ]0,psi_max[
            eps_min = 0.5d0*psi_max*(tanh(-alpha_x)+1.0d0)
            eps_max = 0.5d0*psi_max*(tanh( alpha_x)+1.0d0)
            dX = 2.0d0/(Ntbl-1)
            Xmin = -1.0d0
            do i=1,Ntbl
               eps_=0.5d0*psi_max*
     $              (tanh(alpha_x*(-1.0d0+(i-1)*dX))+1.0d0)
               f = DstFct(eps_)
               if (f.GT.f_max) then
                  f_max = f
                  e_max = eps_
               end if
               Y_tbl(i) = log(f)

               if (i.GE.inext_tick) then
                  write(0,95)
                  itick=itick+1
                  inext_tick = 0.02d0*Ntbl*itick
               end if

            end do
         end if

         write(0,*) ''
            
         call spline_reg(Y_tbl,Ntbl,1.0d30,1.0d30,Xmin,dX,Y2_tbl)
         lPrem = .FALSE.

      end if

      if ((eps.GT.eps_min).AND.(eps.LT.eps_max)) then
         if (psi_max.GT.1.0d30) then
            X = log(eps)
         else
            X = atanh(2.0d0*eps/psi_max-1.0d0)/alpha_x
         end if
         call splint_reg(Y_tbl,Y2_tbl,Ntbl,Xmin,dX,X,DstFct_tbl)
         DstFct_tbl=exp(DstFct_tbl)
      else
         DstFct_tbl = DstFct(eps)
      end if
c
      end
c     
c----------------------------------------------------------------------
      double precision function Integrand_DstFct(w) ! w = sqrt(eps-Psi)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision w
c
c constants :
c ^^^^^^^^^^^
      double precision norm
      parameter (norm = 2.0d0*0.03582244802d0) ! 2/(2**(3/2)*pi**2)
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision Psi
c
c functions :
c ^^^^^^^^^^^
      double precision d2RhodPsi2, U_of_Psi
c
c common block :
c ^^^^^^^^^^^^^^
      double precision eps0
      common /Param_Integrand_DstFct/ eps0
c......................................................................
c
      if (w.LT.0.0d0) then
         write(0,*) '!!! Integrand_DstFct(w) called for negative w: ',
     $        eps0,' !!!'
         call Meurt('Integrand_DstFct')
      end if
      Psi = eps0-w*w

      if (Psi.LT.0.0d0) then
         write(0,*) '!!! Integrand_DstFct(w) called with w>sqrt(eps0):',
     $        w,sqrt(eps0),' !!!'
         call Meurt('Integrand_DstFct')
      end if

      if (Psi.GT.psi_max) then
         Integrand_DstFct = 0.0d0
      else
         Integrand_DstFct = norm * d2RhodPsi2(U_of_Psi(Psi))
      end if
c
      end
c

c----------------------------------------------------------------------
      double precision function R_of_Mr(Mr)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision Mr
c......................................................................
c
      if (Mr.LT.0.0d0 .OR. Mr.GT.1.0d0) then
         write(0,*) '!!! R_of_Mr(Mr) called with Mr=',Mr,' !!!'
         call Meurt('R_of_Mr')
      end if

      R_of_Mr = 1.0d0/(Mr**(-1.0d0/eta)-1.0d0)
c
      end
c
c----------------------------------------------------------------------
      double precision function Mr_of_R(R)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R
c......................................................................
c
      if (R.LT.0.0d0) then
         write(0,*) '!!! Mr_of_R(R) called with R=',R,' !!!'
         call Meurt('Mr_of_R')
      end if

      Mr_of_R = (R/(1.0d0+R))**eta + mu
c
      end
c
c----------------------------------------------------------------------
      double precision function Psi_of_R(R)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      double precision R
c
c function :
c ^^^^^^^^^^
      double precision Psi_of_U
c......................................................................
c
      Psi_of_R = Psi_of_U(1.0d0/R)
c
      end
c
c----------------------------------------------------------------------
      double precision function Psi_of_U(u)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision u
c......................................................................
c
      if (u.LT.0.0d0) then
         write(0,*) '!!! Psi_of_U(u) called for negative u: ',u,' !!!'
         call Meurt('Psi_of_U')
      end if

      if (eta.EQ.1.0d0) then
         Psi_of_U = log(1.0d0+u) + mu*u
      else
         Psi_of_U = (1.0d0-(1.0d0+u)**(1.0d0-eta))/(eta-1.0d0) + mu*u
      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function Pot_NbodyUnits(R)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision u, Rb
c......................................................................
c
      Rb = 1.0d0/(2.0d0*eta-1.0d0)

      if (R.LE.0.0d0) then
         write(0,*)
     $        '!!! Pot_NBodyUnits called for negative R: ',R,' !!!'
         call Meurt('Pot_NBodyUnits')
      end if

      u = Rb/R

      if (eta.EQ.1.0d0) then
         Pot_NbodyUnits = log(1.0d0+u) + mu*u
      else
         Pot_NbodyUnits = (1.0d0-(1.0d0+u)**(1.0d0-eta))/(eta-1.0d0) +
     $        mu*u
      end if
      Pot_NbodyUnits = Pot_NbodyUnits/Rb
c
      end
c
c----------------------------------------------------------------------
      subroutine Get_RcircJmax(E,epsNR, Rcirc,Jmax)
c----------------------------------------------------------------------
c     compute Circular orbit radius and max angular momentum
c     for energy E
c
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision E, epsNR, Rcirc,Jmax ! epsNR is the required relative precison for Rcirc determination
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer i, Nmax
      parameter (Nmax=50)
      double precision Rmin,Rmax, g,dgdr, Phi, dPhidR, d2PhidR2
c
c functions :
c ^^^^^^^^^^^
      double precision rtsafe_rel
      external Get_gdg_Rcirc
c
c common bloc :
c ^^^^^^^^^^^^^
      double precision E_Rcirc
      common /Param_Rcirc/ E_Rcirc
c......................................................................
c
      E_Rcirc = E
c
c---- Find lower and upper bound on Rcirc
c
      Rmax = -(1.0d0+mu)/E_Rcirc ! upper bound
      Rmin = Rmax
      call Get_gdg_Rcirc(Rmin, g,dgdr)
      if (g.LT.0.0d0) then
         write(0,*) '!!! g(Rmax) = ',g,' !!!'
         call meurt('Get_RcircJmax')
      end if
      do i=1,Nmax
         Rmin = 0.25d0*Rmin
         call Get_gdg_Rcirc(Rmin, g,dgdr)
         if (g.LT.0.0d0) goto 9
      end do
      write(0,*) '!!! lower bound on Rcirc not found !!!'
      call meurt('Get_RcircJmax')

 9    continue
      Rcirc = rtsafe_rel(Get_gdg_Rcirc,Rmin,Rmax,epsNR)

      call Get_PhiAnd2Deriv_of_R(Rcirc, Phi, dPhidR, d2PhidR2)
      Jmax = 2.0d0*(E_Rcirc-Phi)
      if (Jmax.LT.0.0d0) then
         write(0,*) '!!! 2.0d0*(E_Rcirc-Phi) = ',
     $        2.0d0*(E_Rcirc-Phi),' !!!'
         call Meurt('Get_RcircJmax')
      end if
      Jmax = Rcirc*sqrt(Jmax)
c
      end
c
c----------------------------------------------------------------------
      subroutine Get_PhiAndDeriv_of_R(R, Phi, dPhidR)
c----------------------------------------------------------------------
c     compute potential (Psi) and first derivative of it
c
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R, Phi, dPhidR
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision InverseR
c......................................................................
c
      InverseR= 1.0d0/R
      if (eta.EQ.1.0d0) then
         Phi = -log(1.0d0 + InverseR) -mu*InverseR
      else
         Phi = ((R/(1.0d0+R))**(eta-1.0d0) -1.0d0)/(eta-1.0d0)
     $        -mu*InverseR
      end if

      dPhidR = ( (R/(1.0d0+R))**eta +mu )*InverseR**2
c
      end
c
c----------------------------------------------------------------------
      subroutine Get_PhiAnd2Deriv_of_R(R, Phi, dPhidR, d2PhidR2)
c----------------------------------------------------------------------
c     compute potential (Psi) and first and second derivative of it
c
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R, Phi, dPhidR, d2PhidR2
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision InverseR
c......................................................................
c
      InverseR= 1.0d0/R
      if (eta.EQ.1.0d0) then
         Phi = -log(1.0d0 + InverseR) -mu*InverseR
      else
         Phi = ((R/(1.0d0+R))**(eta-1.0d0) -1.0d0)/(eta-1.0d0)
     $        -mu*InverseR
      end if

      dPhidR = ( (R/(1.0d0+R))**eta +mu )*InverseR**2

      d2PhidR2 = ( (R/(1.0d0+R))**eta*(eta-2.0d0-eta*R*(1.0d0+R))
     $     - 2.0d0*mu )*InverseR**3
c
      end
c
c----------------------------------------------------------------------
      subroutine Get_gdg_Rcirc(R, g,dgdr)
c----------------------------------------------------------------------
c     g, dg/dR where g(Rc)=0 gives radius of circular orbit (for given 
c     energy.
c
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R, g,dgdr
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision Phi, dPhidR, d2PhidR2
c
c common bloc :
c ^^^^^^^^^^^^^
      double precision E_Rcirc
      common /Param_Rcirc/ E_Rcirc
c......................................................................
c
      call Get_PhiAnd2Deriv_of_R(R, Phi, dPhidR, d2PhidR2)

      g = 2.0d0*Phi -2.0d0* E_Rcirc + R*dPhidR

      dgdr = 3.0d0*dPhidR + R*d2PhidR2
c
      end
c
c----------------------------------------------------------------------
      double precision function Psi_stell_of_R(R)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R
c
c function :
c ^^^^^^^^^^
      double precision Psi_stell_of_U
c......................................................................
c
      if (R.EQ.0.0d0) then
         if (eta.LE.1.0d0) then
            write(0,*) '!!! Psi_stell_of_R(0)=infinity for eta<=1 !!!'
            call Meurt('Psi_stell_of_R')
         else
            Psi_stell_of_R = 1.0d0/(eta-1.0d0)
            return
         end if
      end if
      Psi_stell_of_R = Psi_stell_of_U(1.0d0/R)
c
      end
c
c----------------------------------------------------------------------
      double precision function Psi_stell_of_U(u)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision u
c......................................................................
c
      if (u.LT.0.0d0) then
         write(0,*) '!!! Psi_stell_of_U(u) called for negative u: ',
     $        u,' !!!'
         call Meurt('Psi_stell_of_U')
      end if

      if (eta.EQ.1.0d0) then
         Psi_stell_of_U = log(1.0d0+u)
      else
         Psi_stell_of_U = (1.0d0-(1.0d0+u)**(1.0d0-eta))/(eta-1.0d0)
      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function Sigma1D_of_R(R)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R
c
c functions :
c ^^^^^^^^^^^
      double precision Integrand_fv4,Integrand_fv2,Psi_of_R
      external Integrand_fv4,Integrand_fv2, midpnt, midsql
c
c local var:
c ^^^^^^^^^^
      double precision I1,I2
c
c common block :
c ^^^^^^^^^^^^^^
      double precision Psi0, Vmax
      common /Param_Integrand_Sigma1D/ Psi0, Vmax
c......................................................................
c
      if (R.LT.0.0d0) then
         write(0,*) '!!! Sigma1D_of_R(R) called for negative R: ',
     $        R,' !!!'
         call Meurt('Psi_of_R')
      end if

      Psi0 = Psi_of_R(R)
      Vmax = sqrt(2.0d0*Psi0)
      call qromo(Integrand_fv4,0.0d0,Vmax,acc_Sigma1D,I1,midpnt)
      call qromo(Integrand_fv2,0.0d0,Vmax,acc_Sigma1D,I2,midpnt)
      Sigma1D_of_R =  sqrt(0.333333333d0*I1/I2)
c
      end
c
c----------------------------------------------------------------------
      double precision function DynFricCoef_of_R(R)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R
c
c functions :
c ^^^^^^^^^^^
      double precision Integrand_fv2,Psi_of_R
      external Integrand_fv2, midpnt, midsql
c
c local var:
c ^^^^^^^^^^
      double precision I1,I2, Vcirc, Mencl
c
c common block :
c ^^^^^^^^^^^^^^
      double precision Psi0, Vmax
      common /Param_Integrand_Sigma1D/ Psi0, Vmax
c......................................................................
c
      if (R.LT.0.0d0) then
         write(0,*)
     $        '!!! DynFricIntegral_of_R(R) called for negative R: ',
     $        R,' !!!'
         call Meurt('DynFricIntegral_of_R')
      end if

      Psi0 = Psi_of_R(R)
      Vmax = sqrt(2.0d0*Psi0)
      Mencl = mu + (R/(1.0d0+R))**eta
      Vcirc = sqrt(Mencl/R)
      call qromo(Integrand_fv2,0.0d0,Vcirc,acc_Sigma1D,I1,midpnt)
      call qromo(Integrand_fv2,0.0d0,Vmax,acc_Sigma1D,I2,midpnt)
      !!write(0,*) '>>>>',R,I1,I2
      DynFricCoef_of_R = I1/I2
c
      end
c
c----------------------------------------------------------------------
      double precision function DynFricCoef_of_RandV(R,V)
c----------------------------------------------------------------------
c     R has to be in units of Rbreak and V in units of the local 1D 
c     velocity dispersion
c
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R, V
c
c functions :
c ^^^^^^^^^^^
      double precision Integrand_fv2,Psi_of_R,Sigma1D_of_R
      external Integrand_fv2, midpnt, midsql
c
c local var:
c ^^^^^^^^^^
      double precision I1,I2, sig1D
c
c common block :
c ^^^^^^^^^^^^^^
      double precision Psi0, Vmax
      common /Param_Integrand_Sigma1D/ Psi0, Vmax
c......................................................................
c
      if (R.LT.0.0d0) then
         write(0,*)
     $        '!!! DynFricIntegral_of_RandV(R) called for negative R: ',
     $        R,' !!!'
         call Meurt('DynFricIntegral_of_RandV')
      end if

      Psi0 = Psi_of_R(R)
      Vmax = sqrt(2.0d0*Psi0)
      sig1D = Sigma1D_of_R(R)
      call qromo(Integrand_fv2,0.0d0,min(V*sig1D,Vmax),
     $     acc_Sigma1D,I1,midpnt)
      call qromo(Integrand_fv2,0.0d0,Vmax,
     $     acc_Sigma1D,I2,midpnt)
      !!write(0,*) '>>>>',R,I1,I2
      DynFricCoef_of_RandV = I1/I2
c
      end
c
c----------------------------------------------------------------------
      double precision function Integrand_fv4(v)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision v
c
c functions :
c ^^^^^^^^^^^
      double precision DstFct_tbl
c
c common block :
c ^^^^^^^^^^^^^^
      double precision Psi0, Vmax
      common /Param_Integrand_Sigma1D/ Psi0, Vmax
c......................................................................
c
      if (v.LT.0.0d0) then
         write(0,*) '!!! Integrand_fv4(v) called for negative v: ',
     $        v,' !!!'
         call Meurt('Integrand_fv4')
      end if
      if (v.GT.Vmax) then
         write(0,*) '!!! Integrand_fv4(v) called with v > Vmax !!!'
         call Meurt('Integrand_fv4')
      end if

      Integrand_fv4 =  v**4 * DstFct_tbl(Psi0-0.5d0*v**2)
c
      end
c
c----------------------------------------------------------------------
      double precision function Integrand_fv2(v)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision v
c
c functions :
c ^^^^^^^^^^^
      double precision DstFct_tbl
c
c common block :
c ^^^^^^^^^^^^^^
      double precision Psi0, Vmax
      common /Param_Integrand_Sigma1D/ Psi0, Vmax
c......................................................................
c
      if (v.LT.0.0d0) then
         write(0,*) '!!! Integrand_fv2(v) called for negative v: ',
     $        v,' !!!'
         call Meurt('Integrand_fv2')
      end if
      if (v.GT.Vmax) then
         write(0,*) '!!! Integrand_fv2(v) called with v > Vmax !!!'
         call Meurt('Integrand_fv2')
      end if

      Integrand_fv2 =  v**2 * DstFct_tbl(Psi0-0.5d0*v**2)
c
      end
c
c----------------------------------------------------------------------
      double precision function Integrand_f(v)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision v
c
c functions :
c ^^^^^^^^^^^
      double precision DstFct_tbl
c
c common block :
c ^^^^^^^^^^^^^^
      double precision Psi0, Vmax
      common /Param_Integrand_Sigma1D/ Psi0, Vmax
c......................................................................
c
      if (v.LT.0.0d0) then
         write(0,*) '!!! Integrand_f(v) called for negative v: ',
     $        v,' !!!'
         call Meurt('Integrand_f')
      end if
      if (v.GT.Vmax) then
         write(0,*) '!!! Integrand_f(v) called with v > Vmax !!!'
         call Meurt('Integrand_f')
      end if

      Integrand_f =  DstFct_tbl(Psi0-0.5d0*v**2)
c
      end
c
c----------------------------------------------------------------------
      double precision function ProjDens_of_R(R)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R
c
c functions :
c ^^^^^^^^^^^
      double precision Integrand_ProjDens
      external Integrand_ProjDens, midpnt
c
c common block :
c ^^^^^^^^^^^^^^
      double precision R0
      common /Param_Integrand_PD/ R0
c......................................................................
c
      if (R.LT.0.0d0) then
         write(0,*) '!!! ProjDens_of_R(R) called for negative R: ',
     $        R,' !!!'
         call Meurt('ProjDens_of_R')
      end if
      R0 = R
      call qromo(Integrand_ProjDens,0.0d0,Pi2,acc_Proj,
     $     ProjDens_of_R,midpnt)
      ProjDens_of_R = OneOver2Pi* ProjDens_of_R * eta*R**(eta-2.0d0)
c
      end
c
c----------------------------------------------------------------------
      double precision function Integrand_ProjDens(u)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision u
c
c common block :
c ^^^^^^^^^^^^^^
      double precision R0
      common /Param_Integrand_PD/ R0
c......................................................................
c
c---- See Tremaine et al. (94) Eq. 23
c
      Integrand_ProjDens = sin(u)**2/(R0+sin(u))**(1.0d0+eta)
c
      end
c
c----------------------------------------------------------------------
      subroutine Compute_ProjDensSigma_of_R(R,PD,PS)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R
      double precision PD,PS ! LOS-projected density and velocity dispersion at R
c
c functions :
c ^^^^^^^^^^^
      double precision Integrand_ProjSigma, ProjDens_of_R
      external Integrand_ProjSigma, midpnt
c
c local var:
c ^^^^^^^^^^
      double precision Y, Z
c
c common block :
c ^^^^^^^^^^^^^^
      double precision R0, eta0
      common /Param_Integrand_PS/ R0, eta0
c......................................................................
c
      if (R.LT.0.0d0) then
         write(0,*) '!!! Compute_ProjDensSigma_of_R(R) called',
     $        ' for negative R: ',R,' !!!'
         call Meurt('Compute_ProjDensSigma_of_R')
      end if
c
c---- See Tremaine et al. (94) Eq. 49-51
c
      R0 = R
      eta0 = eta
      call qromo(Integrand_ProjSigma,0.0d0,Pi2,acc_Proj,Y,midpnt)
      Y = eta0/(2.0d0*Pi)*R**(2.0d0*eta0-3.0d0)*Y
      eta0 = 0.5d0*eta
      call qromo(Integrand_ProjSigma,0.0d0,Pi2,acc_Proj,Z,midpnt)
      Z = 2.0d0*eta0/(2.0d0*Pi)*R**(2.0d0*eta0-3.0d0)*Z
      PD = ProjDens_of_R(R)
      PS = sqrt((Y+mu*Z)/PD)
c
      end
c
c----------------------------------------------------------------------
      double precision function ProjSigma_of_R(R)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R
c
c local var:
c ^^^^^^^^^^
      double precision PD
c......................................................................
c
      call Compute_ProjDensSigma_of_R(R,PD,ProjSigma_of_R)
c
      end
c
c----------------------------------------------------------------------
      double precision function Integrand_ProjSigma(u)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      double precision u
c
c common block :
c ^^^^^^^^^^^^^^
      double precision R0, eta0
      common /Param_Integrand_PS/ R0, eta0
c......................................................................
c
c---- See Tremaine et al. (94) Eq. 31
c
      Integrand_ProjSigma = sin(u)**3*cos(u)**2 /
     $     (R0+sin(u))**(1.0d0+2.0d0*eta0)
c
      end
c
c----------------------------------------------------------------------
      subroutine Compute_PeriApo(E,J,R,eps, Rp,Ra)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $     E,J,R, ! Input: Energy, Angular momentum and a position on the orbit
     $     eps,   ! Tolerance for Newton-Raphson solver
     $     Rp,Ra  ! Output: Peri- and Apocenter distances
c
c common block :
c ^^^^^^^^^^^^^^
      double precision J_Vr2_dVr2, E_Vr2_dVr2
      common /param_Vr2_dVr2/ J_Vr2_dVr2, E_Vr2_dVr2
c
c constants :
c ^^^^^^^^^^^
      integer Nmax
      parameter (Nmax=100)
c
c local variables : 
c ^^^^^^^^^^^^^^^^^
      double precision Rmin, Rmax, Vr2,dVr2
      integer i
c
c routines :
c ^^^^^^^^^^
      double precision rtsafe_rel
      external Vr2_dVr2
c......................................................................
c
      J_Vr2_dVr2=J
      E_Vr2_dVr2=E

      call Vr2_dVr2(1.0d0/R,Vr2,dVr2)
      if (Vr2.LE.0.0d0) then
         write(0,*) '!!! Vr2(R)=',Vr2,' !!!'
         call Meurt('Compute_PeriApo')
      end if

      Rmin=R
      do i=1,Nmax
         Rmin=0.5d0*Rmin
         call Vr2_dVr2(1.0d0/Rmin,Vr2,dVr2)
         if (Vr2.LT.0.0d0) goto 9
      end do
      write(0,*) '!!! Rmin with Vr2(Rmin)<0 not found !!!'
      call Meurt('Compute_PeriApo')
 9    continue
      Rp = 1.0d0/rtsafe_rel(Vr2_dVr2,1.0d0/R,1.0d0/Rmin,eps)


      Rmax=R
      do i=1,Nmax
         Rmax=2.0d0*Rmax
         call Vr2_dVr2(1.0d0/Rmax,Vr2,dVr2)
         if (Vr2.LT.0.0d0) goto 19
      end do
      write(0,*) '!!! Rmax with Vr2(Rmin)<0 not found !!!'
      call Meurt('Compute_PeriApo')
 19   continue
      Ra = 1.0d0/rtsafe_rel(Vr2_dVr2,1.0d0/Rmax,1.0d0/R,eps)
c
      end
c
c----------------------------------------------------------------------
      subroutine Compute_Peri(E,J,R,eps, Rp)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $     E,J,R, ! Input: Energy, Angular momentum and a position on the orbit
     $     eps,   ! Tolerance for Newton-Raphson solver
     $     Rp     ! Output: Pericenter distance
c
c common block :
c ^^^^^^^^^^^^^^
      double precision J_Vr2_dVr2, E_Vr2_dVr2
      common /param_Vr2_dVr2/ J_Vr2_dVr2, E_Vr2_dVr2
c
c constants :
c ^^^^^^^^^^^
      integer Nmax
      parameter (Nmax=100)
c
c local variables : 
c ^^^^^^^^^^^^^^^^^
      double precision Rmin, Rmax, Vr2,dVr2
      integer i
c
c routines :
c ^^^^^^^^^^
      double precision rtsafe_rel
      external Vr2_dVr2
c......................................................................
c
      J_Vr2_dVr2=J
      E_Vr2_dVr2=E

      call Vr2_dVr2(1.0d0/R,Vr2,dVr2)
      if (Vr2.LE.0.0d0) then
         write(0,*) '!!! Vr2(R)=',Vr2,' !!!'
         call Meurt('Compute_Peri')
      end if

      Rmin=R
      do i=1,Nmax
         Rmin=0.5d0*Rmin
         call Vr2_dVr2(1.0d0/Rmin,Vr2,dVr2)
         if (Vr2.LT.0.0d0) goto 9
      end do
      write(0,*) '!!! Rmin with Vr2(Rmin)<0 not found !!!'
      call Meurt('Compute_Peri')
 9    continue
      Rp = 1.0d0/rtsafe_rel(Vr2_dVr2,1.0d0/R,1.0d0/Rmin,eps)
c
      end
c
c----------------------------------------------------------------------
      subroutine Compute_Apo(E,J,R,eps, Ra)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $     E,J,R, ! Input: Energy, Angular momentum and a position on the orbit
     $     eps,   ! Tolerance for Newton-Raphson solver
     $     Ra     ! Output: Apocenter distance
c
c common block :
c ^^^^^^^^^^^^^^
      double precision J_Vr2_dVr2, E_Vr2_dVr2
      common /param_Vr2_dVr2/ J_Vr2_dVr2, E_Vr2_dVr2
c
c constants :
c ^^^^^^^^^^^
      integer Nmax
      parameter (Nmax=100)
c
c local variables : 
c ^^^^^^^^^^^^^^^^^
      double precision Rmin, Rmax, Vr2,dVr2
      integer i
c
c routines :
c ^^^^^^^^^^
      double precision rtsafe_rel
      external Vr2_dVr2
c......................................................................
c
      J_Vr2_dVr2=J
      E_Vr2_dVr2=E

      call Vr2_dVr2(1.0d0/R,Vr2,dVr2)
      if (Vr2.LE.0.0d0) then
         write(0,*) '!!! Vr2(R)=',Vr2,' !!!'
         call Meurt('Compute_Apo')
      end if

      Rmax=R
      do i=1,Nmax
         Rmax=2.0d0*Rmax
         call Vr2_dVr2(1.0d0/Rmax,Vr2,dVr2)
         if (Vr2.LT.0.0d0) goto 19
      end do
      write(0,*) '!!! Rmax with Vr2(Rmin)<0 not found !!!'
      call Meurt('Compute_Apo')
 19   continue
      Ra = 1.0d0/rtsafe_rel(Vr2_dVr2,1.0d0/Rmax,1.0d0/R,eps)
c
      end
c
c----------------------------------------------------------------------
      subroutine Vr2_dVr2(U,Vr2,dVr2)
c----------------------------------------------------------------------
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $     U,        ! Input: 1/Distance from center
     $     Vr2,dVR2  ! Output: Vrad^2 and d/dU(Vrad^2)
c
c common block :
c ^^^^^^^^^^^^^^
      double precision J_Vr2_dVr2, E_Vr2_dVr2
      common /param_Vr2_dVr2/ J_Vr2_dVr2, E_Vr2_dVr2
c
c functions :
c ^^^^^^^^^^^
      double precision Psi_of_U
c......................................................................
c
      Vr2 = -(J_Vr2_dVr2*U)**2 + 2.0d0*(E_Vr2_dVr2+Psi_of_U(U))
      dVr2 = 2.0d0*(-J_Vr2_dVr2**2*U + (1.0d0+U)**(-eta) + mu)
c
      end
c
c----------------------------------------------------------------------
      subroutine Get_Rmax(E,epsNR, R_max)
c----------------------------------------------------------------------
c     compute Circular orbit radius and max angular momentum
c     for energy E
c
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision E, epsNR, R_max ! epsNR is the required relative precison for Rcirc determination
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer i, Nmax
      parameter (Nmax=50)
      double precision Rmin,Rmax, g,dgdr, Phi, dPhidR
c
c functions :
c ^^^^^^^^^^^
      double precision rtsafe_rel
      external Get_gdg_Rmax
c
c common bloc :
c ^^^^^^^^^^^^^
      double precision E_Rmax
      common /Param_Rmax/ E_Rmax
c......................................................................
c
      E_Rmax = E
c
c---- Find lower and upper bound on Rcirc
c
      Rmax = -(1.0d0+mu)/E_Rmax ! upper bound
      Rmin = Rmax
      call Get_gdg_Rmax(Rmin, g,dgdr)
      if (g.LT.0.0d0) then
         write(0,*) '!!! g(Rmax) = ',g,' !!!'
         call meurt('Get_Rmax')
      end if
      do i=1,Nmax
         Rmin = 0.25d0*Rmin
         call Get_gdg_Rmax(Rmin, g,dgdr)
         if (g.LT.0.0d0) goto 9
      end do
      write(0,*) '!!! lower bound on Rmax not found !!!'
      call meurt('Get_Rmax')

 9    continue
      R_max = rtsafe_rel(Get_gdg_Rmax,Rmin,Rmax,epsNR)
c
      end
c
c----------------------------------------------------------------------
      subroutine Get_gdg_Rmax(R, g,dgdr)
c----------------------------------------------------------------------
c     g, dg/dR where g(Rmax)=0 gives maximum radius atteinable by any
c     particle with energy E (radial orbit)
c
      implicit none
      include 'EtaModel_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision R, g,dgdr
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision Phi, dPhidR
c
c common bloc :
c ^^^^^^^^^^^^^
      double precision E_Rmax
      common /Param_Rmax/ E_Rmax
c......................................................................
c
      call Get_PhiAndDeriv_of_R(R, Phi, dPhidR)

      g = Phi -E_Rmax

      dgdr = dPhidR
c
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     Routine(s) to determine distribution of relative velocities
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c----------------------------------------------------------------------
      double precision function DistribVrel(R,Vrel) 
c----------------------------------------------------------------------
c     DistribVrel(R,Vrel)*dVrel = fractional number of encounters (at R) 
c     with relative velocities in [Vrel, Vrel+dVrel[
c     
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision R, Vrel
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision Dens
c
c common block to pass hidden parameters :
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      double precision LocalVmax, LocalPsi, Vrel_value, v1_value
      common /param_DistribVrel/ LocalVmax, LocalPsi,
     $     Vrel_value, v1_value
c
c functions :
c ^^^^^^^^^^^
      double precision Psi_of_R
      external Integrand_DstVrel, midpnt
c......................................................................
c
      LocalPsi = Psi_of_R(R)
      LocalVmax = sqrt(2.0d0*LocalPsi)

      if (Vrel.GE.2.0d0*LocalVmax) then
         DistribVrel =0.0d0
         return
      end if
      Vrel_value = Vrel

      Dens = eta/((4.0d0*Pi)*R**(3.0d0-eta)*(1.0d0+R)**(1.0d0+eta))
c
c---- Integrate product of distribution functions over velocity of first particle
c
c      call qromb(Integrand_DstVrel,0.0d0,LocalVmax,acc_DstVrel,
c     $     DistribVrel)
      call qromo(Integrand_DstVrel,0.0d0,LocalVmax,acc_DstVrel,
     $     DistribVrel,midpnt)
      DistribVrel = 8.0d0*(Pi*Vrel/Dens)**2 * DistribVrel
c
      end
c
c----------------------------------------------------------------------
      double precision function Integrand_DstVrel(v1)
c----------------------------------------------------------------------
c     
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision v1 ! v1 is the velocity of the first particle
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision u_max
c
c common block to pass hidden parameters :
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      double precision LocalVmax, LocalPsi, Vrel_value, v1_value
      common /param_DistribVrel/ LocalVmax, LocalPsi,
     $     Vrel_value, v1_value
c
c functions :
c ^^^^^^^^^^^
      double precision DstFct, DstFct_tbl
      external Integrand_Integr_DstVrel, midpnt
c......................................................................
c
      if (v1.GE.LocalVmax) then
         Integrand_DstVrel = 0.0d0
         return
      end if
      v1_value = v1
      u_max = 1.0d0
      if (v1*Vrel_value.GT.1.0d-10) then
         u_max = min(u_max,(LocalVmax**2-Vrel_value**2-v1**2)/
     $        (2.0d0*v1*Vrel_value))
         if (u_max.LE.-1.0d0) then
            Integrand_DstVrel = 0.0d0
            return
         end if
      end if
c      call qromb(Integrand_Integr_DstVrel,-1.0d0,u_max,acc_DstVrel,
c     $     Integrand_DstVrel)
      call qromo(Integrand_Integr_DstVrel,-1.0d0,u_max,acc_DstVrel,
     $     Integrand_DstVrel,midpnt)
      Integrand_DstVrel = v1**2*DstFct_tbl(LocalPsi-0.5d0*v1**2) *
     $     Integrand_DstVrel
c
      end
c
c----------------------------------------------------------------------
      double precision function Integrand_Integr_DstVrel(u)
c----------------------------------------------------------------------
c     
      implicit none
      include 'EtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision u ! u is the cosinus of the angle between the 
                         ! relative velocity and the velocity of the first particle
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision eps
c
c common block to pass hidden parameters :
c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      double precision LocalVmax, LocalPsi, Vrel_value, v1_value
      common /param_DistribVrel/ LocalVmax, LocalPsi,
     $     Vrel_value, v1_value
c
c functions :
c ^^^^^^^^^^^
      double precision DstFct, DstFct_tbl
c......................................................................
c
      eps = LocalPsi
     $     - 0.5d0*(v1_value**2+Vrel_value**2 +
     $              2.0d0*v1_value*Vrel_value*u)
      if (u.LE.-1.0d0 .OR. u.GE.1.0d0 .OR. eps.LE.0.0d0) then
         Integrand_Integr_DstVrel = 0.0d0
         return
      end if
      Integrand_Integr_DstVrel = DstFct_tbl(eps)
c
      end
c

c**********************************************************************
c     Routine to compute dynamical friction coefficient using a look-up
c     table (and to build, write and read it)
c**********************************************************************

c----------------------------------------------------------------------
      subroutine Determine_DFCoef
c----------------------------------------------------------------------
c     
      implicit none
      include 'EtaModel_inc.f'
      include 'Table_DynFricCoefEtaModel_inc.f'
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      character*512 file_name
      logical lFileDF
      integer iReadFileDFOK, l
c
c function :
c ^^^^^^^^^^
      integer LongChaine
c......................................................................
c
c---- Check for file containing Dynamical friction data
c
      write(file_name,100) eta,mu
 100  format('OmegaDynFricEtaModel_eta',F8.6,'_mu',F8.6,'.asc')
      l = LongChaine(PathFile_DFCoef)
      if (l.GT.0) then
         file_name = PathFile_DFCoef(1:l) // '/' // file_name
      end if
      l = LongChaine(file_name)
      write(0,*) '> Checking for file ',file_name(1:l)
      
      inquire(file=file_name,exist=lFileDF)
c
c---- Read data from file
c
      iReadFileDFOK = 99
      if (lFileDF) then
         write(0,*) '> Read tabulated dyn fric coef from ',
     $        file_name(1:l)
         open(11,file=file_name)
         call Read_Table_DFCoef(11,iReadFileDFOK)
         if (iReadFileDFOK.EQ.0) then
            write(0,*) '> Data read correctly'
         else
            write(0,*) '!!! Error encounter when reading data, code ',
     $           iReadFileDFOK,' !!!'
            call exit(100+iReadFileDFOK)
         end if
         close(11)
         return
      else
c
c---- Compute data a save it to disk for future reuse
c
         call Compute_Table_DFCoef
         open(11,file=file_name)
         call Write_Table_DFCoef(11)
         write(0,*) '> Data for tabulated dyn fric coef written to  ',
     $        file_name(1:l)
         close(11)
         return
      end if
c
      end 
c
c----------------------------------------------------------------------
      subroutine Compute_Table_DFCoef
c----------------------------------------------------------------------
c     
      implicit none
      include 'EtaModel_inc.f'
      include 'Table_DynFricCoefEtaModel_inc.f' ! Look there for explanations
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer i_R, i_nu
      integer i, itick, inext_tick
      double precision nu, R
c
c function :
c ^^^^^^^^^^
      double precision DynFricCoef_of_RandV, Sigma1D_of_R
c......................................................................
c
      lR_DFCoef_min = log10(R_DFCoef_min)
      lnu_DFCoef_min = log10(nu_DFCoef_min)
      dlR_DFCoef = log10(R_DFCoef_max/R_DFCoef_min)/
     $     dble(N_DFCoef_R-1)
      dlnu_DFCoef = log10(nu_DFCoef_max/nu_DFCoef_min)/
     $     dble(N_DFCoef_nu-1)

      write(0,90)
 90   format(/,
     $     '####################################################',/,
     $     '## Tabulating Dynamical Friction Coefficient..... ##',/,
     $     '####################################################',//,
     $     '0    10   20   30   40   50   60   70   80   90  100%',/,
     $     '|----|----|----|----|----|----|----|----|----|----|',/,
     $     $)
      inext_tick=0
      itick=0
      i=0
c
      do i_R=0,N_DFCoef_R-1
         R = R_DFCoef_min*10**(i_R*dlR_DFCoef)
         TabSig1D(i_R+1) = log10(Sigma1D_of_R(R))
         do i_nu=0,N_DFCoef_nu-1
            nu = nu_DFCoef_min*10**(i_nu*dlnu_DFCoef)
            TabDFCoef(i_nu+1,i_R+1) =
     $           log10(DynFricCoef_of_RandV(R,nu))

            i=i+1
            if (i.GE.inext_tick) then
               write(0,95)
 95            format('^',$)
               itick=itick+1
               inext_tick = 0.02d0*N_DFCoef_R*N_DFCoef_nu*itick
            end if

         end do
      end do
      write(0,*) ''
      lTabDFCoef = .TRUE.
c
      end

c----------------------------------------------------------------------
      subroutine Write_Table_DFCoef(iout)
c----------------------------------------------------------------------
c   
      implicit none
      include 'EtaModel_inc.f'
      include 'Table_DynFricCoefEtaModel_inc.f' ! Look there for explanations
c
c argument :
c ^^^^^^^^^^
      integer iout
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer i_R, i_nu
      character*64 format_str
c......................................................................
c
      write(iout,100) 
 100  format(
     $     '# Tabulated omega coefficient for dynamical friction ',
     $     'in eta model',/,
     $     '# Units are such that G = Mstars = Rbreak = 1',/,
     $     '#',/,
     $     '# TabDFCoef(j,i) contains the log10 of ',
     $     'omega(R_i,nu_j)',/,
     $     '# where R_i = Rmin*10**((i-1)*dlR)',
     $     ' is the radius in units of Rbreak (NOT N-body!)',/,
     $     '# nu_j = nu_min*10**((j-1)*dnu) is the velocity ',
     $     'in units of',/,
     $     '# the local 1D velocity dispersion',/,
     $     '#',/,
     $     '# omega (non dimensional) is defined through:',/,
     $     '# dV/dt = -4*pi*ln(Lambda)*G^2*Rho*M*(1+m/M)/V^2 *',/,
     $     '# omega(R/Rb,V/sigma_1D) * 1/V^2',/,
     $     '# TabSig1D(j) is log10 of the 1D velocity dispersion ',
     $     'at R_j, in units of sqrt(G*Mstars/Rbreak) (NOT N-body!)',/,
     $     '#',/,
     $     '# First data line: eta and mu',/,
     $     '# Second data line: Rmin,Rmax,N_R, nu_min,nu_max,N_nu',/,
     $     '# Rest of data: ',
     $     'TabSig1D(j), TabDFCoef(j,i); one line per j value',/,
     $     '#    -----> nu_i',/,
     $     '#    |',/,
     $     '#    |',/,
     $     '#    |',/,
     $     '#   \\/  R_j',/
     $     )

      write(iout,110) eta,mu
 110  format(2(E12.5,' '))
      write(iout,120)
     $     R_DFCoef_min, R_DFCoef_max, N_DFCoef_R,
     $     nu_DFCoef_min, nu_DFCoef_max, N_DFCoef_nu
 120  format(2(2(E12.5,' '),I3,' '))

      write(format_str,130) N_DFCoef_nu
 130  format('(E12.5,"   ",',I3,'(E12.5," "))')

      do i_R=1,N_DFCoef_R
         write(iout,format_str) TabSig1D(i_R),
     $        (TabDFCoef(i_nu,i_R),i_nu=1,N_DFCoef_nu)
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine Read_Table_DFCoef(ifile,iOK)
c----------------------------------------------------------------------
c   
      implicit none
      include 'EtaModel_inc.f'
      include 'Table_DynFricCoefEtaModel_inc.f' ! Look there for explanations
c
c argument :
c ^^^^^^^^^^
      integer ifile
      integer iOK ! set to 0 if the file was read without problems
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer NR_file, Nnu_file
      double precision Rmin_file, Rmax_file, numin_file, numax_file,
     $     eta_file, mu_file
      integer i_R, i_nu
      character*1 car
c......................................................................
c
      iOK = 0
c
c---- Skip over comments
c
      car = '#'
      do while (car.EQ.'#')
         read(ifile,*,ERR=991,END=991) car
      end do
c
c---- Read parameters and check that they are correct (or set horror flag)
c
      backspace(ifile) ! WHY ???
      read (ifile,*,ERR=992) eta_file, mu_file
      if (abs(eta_file-eta).GT.1e-4*eta) goto 994
      if (abs(mu_file-mu).GT.1d-4*mu) goto 994
      read (ifile,*,ERR=993) 
     $     Rmin_file, Rmax_file, NR_file,
     $     numin_file, numax_file, Nnu_file
      if (NR_file.NE.N_DFCoef_R) goto 995
      if (Nnu_file.NE.N_DFCoef_nu) goto 995
      if (abs(Rmin_file-R_DFCoef_min).GT.1d-4*Rmin_file) goto 996
      if (abs(Rmax_file-R_DFCoef_max).GT.1d-4*Rmax_file) goto 996
      if (abs(numin_file-nu_DFCoef_min).GT.1d-4*numin_file) goto 996
      if (abs(numax_file-nu_DFCoef_max).GT.1d-4*numax_file) goto 996
c
c---- Read array per se
c
      do i_R=1,N_DFCoef_R
         read(ifile,*,END=997,ERR=997) TabSig1D(i_R),
     $        (TabDFCoef(i_nu,i_R),i_nu=1,N_DFCoef_nu)
      end do
      lTabDFCoef = .TRUE.
      lR_DFCoef_min = log10(R_DFCoef_min)
      lnu_DFCoef_min = log10(nu_DFCoef_min)
      dlR_DFCoef = log10(R_DFCoef_max/R_DFCoef_min)/
     $     dble(N_DFCoef_R-1)
      dlnu_DFCoef = log10(nu_DFCoef_max/nu_DFCoef_min)/
     $     dble(N_DFCoef_nu-1)
      return
c
c---- set error flag
c
 991  continue
      iOK=1
      return
 992  continue
      iOK=2
      return
 993  continue
      iOK=3
      return
 994  continue
      iOK=4
      return
 995  continue
      iOK=5
      return
 996  continue
      iOK=6
      return
 997  continue
      iOK=7
      return
c
      end
c
c----------------------------------------------------------------------
      double precision function OmegaDFofRV(R,V) 
c----------------------------------------------------------------------
c     Returns the coefficient omega for dynamical friction in eta model.
c     We use the "real" local velocity distribution of the eta model.
c     See Table_DynFricCoefEtaModel_inc.f for explanantions
c
      implicit none
      include 'EtaModel_inc.f'
      include 'Table_DynFricCoefEtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision R, V     ! R is the radius in units of R_break (NOT N-body units!)
                                ! nu is velocity in units of sqrt(G*Mstars/R_break) (NOT N-body units!)
                                ! Mstars is the total mass in stars
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer iR_inf, iR_sup, inu_inf, inu_sup
      double precision lR, lnu, w_inf, w_Rinf, sigma_1D,
     $     value_Rinf, value_Rsup, value

c......................................................................
c
      if (.NOT.lTabDFCoef) call Compute_Table_DFCoef
c
c===  Interpolation in computed table of DF coefficients
c     using simple bilinear scheme (in log-log-log space)
c
c     Find R-index in table of DF coefficients
c
      lR = log10(R)
      iR_inf = (lR-lR_DFCoef_min)/dlR_DFCoef + 1
      iR_sup = iR_inf+1
c
c     Don't try to extrapolate in the R direction
c     just use the smallest or largest R value in the range
c
      if (iR_inf.LT.1) then
         iR_inf=1
         iR_sup=1
         w_Rinf = 1.0d0
      else if (iR_inf.GE.N_DFCoef_R) then
         iR_inf=N_DFCoef_R
         iR_sup=N_DFCoef_R
         w_Rinf = 1.0d0
      else
         w_Rinf = (lR_DFCoef_min+(iR_sup-1)*dlR_DFCoef-lR)/
     $        ((iR_sup-iR_inf)*dlR_DFCoef)
      end if
c
c     Get local velocity dispersion by interpolation
c
      sigma_1D = 10**(w_Rinf*TabSig1D(iR_inf) +
     $     (1.0d0-w_Rinf)*TabSig1D(iR_sup))

c
c     Find nu-index in table of DF coefficients
c     with nu = V/sigma_1D(R)
c
      lnu = log10(V/sigma_1D)
      inu_inf = (lnu-lnu_DFCoef_min)/dlnu_DFCoef + 1
      inu_sup = inu_inf+1
c
c     interpolation on nu for R_inf and R_sup
c
      value_Rsup = 0.0d0
      if (inu_inf.LT.1) then
         value_Rinf = TabDFCoef(1,iR_inf)-3.0d0*(lnu_DFCoef_min-lnu)
         if (iR_sup.GT.iR_inf) value_Rsup =
     $        TabDFCoef(1,iR_sup)-3.0d0*(lnu_DFCoef_min-lnu)
      else if (inu_sup.GT.N_DFCoef_nu) then
         value_Rinf = 0.0d0
      else
         w_inf = (lnu_DFCoef_min+(inu_sup-1)*dlnu_DFCoef-lnu)/
     $        ((inu_sup-inu_inf)*dlnu_DFCoef)
         value_Rinf =  w_inf*TabDFCoef(inu_inf,iR_inf) +
     $        (1.0d0-w_inf)*TabDFCoef(inu_sup,iR_inf)
         if (iR_sup.GT.iR_inf) value_Rsup =
     $        w_inf*TabDFCoef(inu_inf,iR_sup) +
     $        (1.0d0-w_inf)*TabDFCoef(inu_sup,iR_sup)
      end if  
c
c     interpolation on R
c
      value = w_Rinf*value_Rinf + (1.0d0-w_Rinf)*value_Rsup
c
c---- Return interpolated value
c
      OmegaDFofRV = 10**value
c
      end
c
c----------------------------------------------------------------------
      double precision function OmegaDFofRV_maxw(R,V) 
c----------------------------------------------------------------------
c     Returns the coefficient omega for dynamical friction in eta model.
c     Here we assume a locally Maxwellian velocity distribution so the
c     only things which actually matters is V/sigma1D(R)
c
      implicit none
      include 'EtaModel_inc.f'
      include 'Table_DynFricCoefEtaModel_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      double precision R, V     ! R is the radius in units of R_break (NOT N-body units!)
                                ! nu is velocity in units of sqrt(G*Mstars/R_break) (NOT N-body units!)
                                ! Mstars is the total mass in stars
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer iR_inf, iR_sup
      double precision lR, w_Rinf, sigma_1D
c
c function :
c ^^^^^^^^^^
      double precision OmegaDFCoef_maxw
c......................................................................
c
      if (.NOT.lTabDFCoef)
     $     call Compute_Table_DFCoef ! Quite an overkill for Maxwellian case!
                                     ! as only sigma1D(R) needs to be tabulated
c
c     Find R-index in table of DF coefficients
c
      lR = log10(R)
      iR_inf = (lR-lR_DFCoef_min)/dlR_DFCoef + 1
      iR_sup=iR_inf+1
c
c     Don't try to extrapolate in the R direction
c     just use the smallest or largest R value in the range
c
      if (iR_inf.LT.1) then
         iR_inf=1
         iR_sup=1
         w_Rinf = 1.0d0
      else if (iR_inf.GE.N_DFCoef_R) then
         iR_inf=N_DFCoef_R
         iR_sup=N_DFCoef_R
         w_Rinf = 1.0d0
      else
         w_Rinf = (lR_DFCoef_min+(iR_sup-1)*dlR_DFCoef-lR)/
     $        ((iR_sup-iR_inf)*dlR_DFCoef)
      end if
c
c     Get local velocity dispersion by interpolation
c
      sigma_1D = 10**(w_Rinf*TabSig1D(iR_inf) +
     $     (1.0d0-w_Rinf)*TabSig1D(iR_sup))
c
c     Compute omega value for Maxwellian distribution 
c
      OmegaDFofRV_maxw = OmegaDFCoef_maxw(V/sigma_1D)
c
      end
c
c----------------------------------------------------------------------
      double precision function OmegaDFCoef_maxw(nu) 
c----------------------------------------------------------------------
c     Returns the coefficient omega for dybamical friction for a 
c     background with a Maxwellian 
c
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      double precision nu       ! nu is velocity devided by (local) 1D velocity dispersion
      double precision sqrt2, sqrt_2overpi
c
c constants :
c ^^^^^^^^^^^
      parameter (sqrt2=0.141421356237310d1,
     $     sqrt_2overpi=0.797884560802865d0)
c
c function :
c ^^^^^^^^^^
      double precision erfcc
c......................................................................
c
      OmegaDFCoef_maxw  = 1.0d0-erfcc(nu/sqrt2)
     $     - sqrt_2overpi*nu*exp(-0.5d0*nu*nu)
c
      end
c

c----------------------------------------------------------------------
      block data Init_EtaModel
c----------------------------------------------------------------------
c
      implicit none
      include 'EtaModel_inc.f'
      include 'Table_DynFricCoefEtaModel_inc.f'

      data  acc_u,  acc_DstFct, acc_Sigma1D, acc_Proj, acc_DstVrel
     $     /1.0d-9, 1.0d-6,     1.0d-6,      1.0d-6,   1.0d-6/
      data lTabDFCoef /.FALSE./
      data PathFile_DFCoef /'.'/
      data R_drop /-1.0/ ! No central density depletion by default!
c
      end
