#define _DEFAULT_MR_ 1 /* 0: old relation; 1: new one (fall 2003) */ 
#if (_EVOL_STELL_<10) /* Simple stellar evol (MS-->remnant) or no stellar evolution */
#undef _NORMALISER_MR_ /* determine si il faut normaliser a Msol=1 */
                       /* Ce qui n'est pas forcement une bonne idee car le modele a 1 Msol semble bizarre */
c
c======================================================================
c     Routines implementant la relation M->R stellaire
c======================================================================
c----------------------------------------------------------------------
      subroutine Read_RelMR(iFile)
c----------------------------------------------------------------------
c     read M-R relation for MS stars from file
c
      include 'EvolAmas_RelMR_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFile
c
c var locales :
c ^^^^^^^^^^^^^
      character*80 line
c......................................................................
c
      NtabMR=0
      do while (.TRUE.)
         read(iFile,'(A)',END=99) line
         if (line(1:1).NE.'#') then
            if (NtabMR.GE.NtabMR_Max) then
               write(0,*) '!!! MS M--R relation can contain a maximum',
     $              ' of ',NtabMR_Max,' elements !!!'
               call Meurt('Read_RelMR')
            end if
            NtabMR = NtabMR+1
            read(line,*) TabM_MS(NtabMR), TabR_MS(NtabMR)
         end if
      end do
 99   continue
c
      end
c
c----------------------------------------------------------------------
      subroutine Initialiser_RelMR
c----------------------------------------------------------------------
      include 'EvolAmas_RelMR_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      double precision b,x_ext(Next),y_ext(Next)
      integer i
c......................................................................
c
c---- Relation M-R pour sequence principale
c
                                ! calcul des pentes pour l'extrapolation
      call MoindresCarres(TabM_MS,TabR_MS,Next,a1_MR,b)
      do i=1,Next
         x_ext(i)=TabM_MS(NTabMR-Next+i)
         y_ext(i)=TabR_MS(NTabMR-Next+i)
      end do
      call MoindresCarres(x_ext,y_ext,Next,a2_MR,b)
                                ! calcul des secondes derivees pour le spline
      call spline(TabM_MS,TabR_MS,NTabMR,a1_MR,a2_MR,y2_MR)
                                ! on determine la valeur du rayon d'une etoile de 1 Msol
      call splint(TabM_MS,TabR_MS,y2_MR,NTabMR,0.0d0,UnSurRsol)
      UnSurRsol = exp(-UnSurRsol)
c
      end
c
c----------------------------------------------------------------------
      double precision function RelMR(Masse,Age,iType)
c----------------------------------------------------------------------
      include 'EvolAmas_RelMR_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision
     $     Masse,               ! Masse en Msol
     $     Age                  ! Age en annees, compte a partir de la formation 
                                ! de l'etoile en tant que type iType
      integer*1 iType
c
c     Pour l'instant, on ne tient pas compte de l'age et il faut lui donner
c     une valeur negative pour montrer que l'on est conscient de cette limitation!
c
c constantes :
c ^^^^^^^^^^^^
      double precision UnTiers
      parameter (UnTiers=1.0d0/3.0d0)
c
c var locales :
c ^^^^^^^^^^^^^
      double precision lM,lR,m13,
     $     arg_sqrt
      logical lprem_MR
      data lprem_MR /.TRUE./
      save lprem_MR
c......................................................................
c
      RelMR = -1.0d0
                                ! initialiser les spline au premier appel
      if (lprem_MR) then
         call Initialiser_RelMR
         lprem_MR = .FALSE.
      end if

      if (Age.GE.0.0d0) then
         write(0,*) '!!! Impossible de tenir compte de',
     $        ' l''age des etoiles !!!'
         call Meurt('RelMR')
      end if

      if (Masse.LE.0.0d0) then
         RelMR = 0.0d0
      else

         if (iType.EQ.TypeEtoile_MS)
     $        then              ! SEQUENCE PRINCIPALE
                                ! ===================

            lM = log(Masse)
                                ! on traite a part les cas d'extrapolation :
                                ! dans ce cas on construit une relation afine
                                ! sur les 2 points extremes
            if (lM.LT.TabM_MS(1)) then
               lR = TabR_MS(1) + a1_MR*(lM-TabM_MS(1))
            else if (lM.GT.TabM_MS(NTabMR)) then
               ! lR = TabR_MS(NTabMR) + a2_MR*(lM-TabM_MS(NTabMR))
               ! Use M-R exponent for very massive stars from Goodman & Tan 03, Eq. 8
               lR = TabR_MS(NTabMR) + MR_exp_VMS*(lM-TabM_MS(NTabMR))
            else
               call splint(TabM_MS,TabR_MS,y2_MR,NTabMR,lM,lR)
            end if
#if (_NORMALISER_MR_>0)
            RelMR = UnSurRsol*exp(lR)
#else
            RelMR = exp(lR)
#endif
         else if (iType.EQ.TypeEtoile_WD)
     $           then           ! NAINE BLANCHE
                                ! =============
            if (Masse.GT.M_chandra) then
               write(0,*) '!!! Masse Naine blanche = ',Masse,
     $              ' > M_chandra = ',M_chandra,' !!!'
               !call Meurt('RelMR')
            end if
                                ! Relation M-R de Hurley et al.00, eq. 91
            m13 = (Masse/M_chandra)**UnTiers
#if (_SURVEILLE_SQRT_ > 0)
            arg_sqrt = 1.0d0-m13**4
            if (arg_sqrt.LT.0.0d0) then
               write(0,*) '!!! arg_sqrt = ',arg_sqrt,' !!!'
               call Meurt('RelMR')
            end if
            RelMR = max( R_NS, 0.0115d0*sqrt(arg_sqrt)/m13)
#else
            RelMR = max( R_NS, 0.0115d0*sqrt(1.0d0-m13**4)/m13)
#endif

         else if (iType.EQ.TypeEtoile_NS)
     $           then           ! ETOILE DE NEUTRONS
                                ! ==================c
            if (Masse.GT.Mmax_NS) then
               write(0,*) '!!! Masse Etoiles de neutrons = ',Masse,
     $              ' > Mmax_NS = ',Mmax_NS,' !!!'
               !call Meurt('RelMR')
            end if
            if (Masse.LT.Mmin_NS) then
               write(0,*) '!!! Masse Etoiles de neutrons = ',Masse,
     $              ' < Mmin_NS = ',Mmin_NS,' !!!'
               !call Meurt('RelMR')
            end if
            RelMR = R_NS        ! Relation M-R constante, cf Srinivasam
                                ! cours de Saas-Fee 25, p.153

         else if (iType.EQ.TypeEtoile_BH)
     $           then           ! TROU NOIR STELLAIRE
                                ! ===================
            RelMR = 4.24e-6*Masse ! Rayon de l'horizon

         else                   ! Autre type d'etoile
            write(0,*) '!!! Pas de relation M-R disponible pour Type ',
     $           iType,' !!!'
            call Meurt('RelMR')
         end if
            

      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function Rhalf_over_R_MS(M)
c----------------------------------------------------------------------
      include 'EvolAmas_RelMR_inc.f'
c
c argument :
c ^^^^^^^^^^
      double precision M ! ZAMS mass in M_sun
c
c var locales :
c ^^^^^^^^^^^^^
      double precision lM
      integer i
c......................................................................
c
      lM=log10(M)
      if (lM.LT.lM_RhR_MS(1)) then ! set to n=1.5 polytrope value at lower limit and n=3 at upper limit
         Rhalf_over_R_MS = RhR_MS(1)
      else if (lM.GT.lM_RhR_MS(NRhR)) then
         Rhalf_over_R_MS = RhR_MS(NRhR)
      else                      !  linear interpolation
         call locate(lM_RhR_MS,NRhR,lM,i)
         Rhalf_over_R_MS = RhR_MS(i) + (RhR_MS(i+1)-RhR_MS(i)) *
     $        (lM-lM_RhR_MS(i)) / (lM_RhR_MS(i+1)-lM_RhR_MS(i))
      end if
c
      end
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      block data init_MR
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      include 'EvolAmas_RelMR_inc.f'

#if (_DEFAULT_MR_==0)

      data NtabMR /45/
                                ! Relation ln(M)->ln(R) selon modeles
                                ! du groupe Maeder Meynet et al.
                                ! Faibles masses (<0.4 Msol) selon Chabrier+Baraffe 2000

      data TabM_MS /
     $     -6.882d+00,-6.120d+00,-5.524d+00,-4.932d+00,-4.428d+00,
     $     -3.880d+00,-3.293d+00,-2.966d+00,-2.809d+00,-2.680d+00,
     $     -2.519d+00,-2.261d+00,-2.109d+00,-1.925d+00,-1.798d+00,
     $     -1.635d+00,-1.462d+00,-1.262d+00,-1.149d+00,-1.038d+00,
     $     
     $     -0.916291d0,-0.693147d0,-0.510826d0,-0.356675d0,
     $     -0.223144d0,-0.105361d0,0d0,0.182322d0,0.405465d0,0.530628d0,
     $     0.693147d0,0.916291d0,1.09861d0,1.38629d0,1.60944d0,
     $     1.94591d0,2.19678d0,2.48157d0,2.6909d0,2.96119d0,3.17643d0,
     $     3.61674d0,3.94296d0,4.30695d0,4.61809d0,
     $     
     $     55*0.0d0/
      data TabR_MS /
     $     -2.235d+00,-2.244d+00,-2.244d+00,-2.273d+00,-2.303d+00,
     $     -2.375d+00,-2.465d+00,-2.513d+00,-2.538d+00,-2.501d+00,
     $     -2.303d+00,-2.048d+00,-1.897d+00,-1.778d+00,-1.671d+00,
     $     -1.556d+00,-1.427d+00,-1.262d+00,-1.171d+00,-1.097d+00,
     $
     $     -1.0106d0,-0.832409d0,-0.621757d0,-0.449417d0,
     $     -0.332679d0,-0.197232d0,0.0760347d0,0.354172d0,0.602128d0,
     $     0.71784d0,0.846726d0,0.965462d0,1.05918d0,1.19544d0,
     $     1.32601d0,1.49962d0,1.64577d0,1.82116d0,1.96389d0,
     $     2.17179d0,2.30308d0,2.72811d0,3.20153d0,3.05173d0,
     $     3.00924d0,
     $     
     $     55*0.0d0/

#elif (_DEFAULT_MR_==1)

c Mass-Radius relation for MS stars Z=0.001
c
c References:
c     * 1-120 M_sun: Schaller et al. 1992, Astronomy and Astrophysics Supl. 96, 269
c     * 0.4-1 M_sun: Charbonnel et al. 1999, Astronomy and Astrophysics Supl. 135, 405
c     * 0.001-0.4 M_sun: Chabrier & Baraffe 2000, Annual Review of Astronomy and Astrophysics 38, 337 (Fig 3.)
c
      data NtabMR /40/
      data  TabM_MS /
     &     -6.92774e+00, -5.86809e+00, -4.62850e+00, -3.92872e+00,
     &     -3.40889e+00, -3.08899e+00, -2.77911e+00, -2.68912e+00,
     &     -2.42920e+00, -2.21928e+00, -1.87939e+00, -1.67945e+00,
     &     -1.48951e+00, -1.29958e+00, -1.09964e+00, -9.16291e-01,
     &     -6.93147e-01, -5.10826e-01, -3.56675e-01, -2.23144e-01,
     &     -1.05361e-01, 0.00000e+00, 2.23144e-01, 4.05465e-01,
     &     5.30628e-01, 6.93147e-01, 9.16291e-01, 1.09861e+00,
     &     1.38629e+00, 1.60944e+00, 1.94590e+00, 2.19711e+00,
     &     2.48422e+00, 2.70657e+00, 2.99355e+00, 3.21672e+00,
     &     3.68722e+00, 4.09264e+00, 4.44019e+00, 4.78236e+00,
     $     
     $     60*0.0d0 /
      data TabR_MS /
     &     -2.21146e+00, -2.21359e+00, -2.28845e+00, -2.35197e+00,
     &     -2.43662e+00, -2.50964e+00, -2.52940e+00, -2.54893e+00,
     &     -2.23467e+00, -2.01849e+00, -1.74941e+00, -1.58829e+00,
     &     -1.44309e+00, -1.29388e+00, -1.14955e+00, -1.04410e+00,
     &     -8.23581e-01, -6.33927e-01, -4.86662e-01, -3.00765e-01,
     &     -1.23464e-01, -1.12594e-02, 6.61764e-02, 1.76821e-01,
     &     2.32817e-01, 3.15777e-01, 4.26743e-01, 5.41048e-01,
     &     7.25132e-01, 8.64423e-01, 1.08100e+00, 1.24329e+00,
     &     1.42723e+00, 1.57212e+00, 1.75049e+00, 1.88913e+00,
     &     2.16490e+00, 2.40621e+00, 2.59836e+00, 2.78126e+00,
     $     
     $     60*0.0d0 /
#endif

c
c     Relation M_R BC 2000 :
c     Ln(M) :
c     $     -6.882d+00,-6.120d+00,-5.524d+00,-4.932d+00,-4.428d+00,
c     $     -3.880d+00,-3.293d+00,-2.966d+00,-2.809d+00,-2.680d+00,
c     $     -2.519d+00,-2.261d+00,-2.109d+00,-1.925d+00,-1.798d+00,
c     $     -1.635d+00,-1.462d+00,-1.262d+00,-1.149d+00,-1.038d+00,
c     $     -9.233d-01,-8.543d-01,-7.691d-01,-6.954d-01,-5.872d-01,
c     $     -5.112d-01,-4.375d-01,-3.500d-01,-2.602d-01,-1.957d-01,
c     $     -1.174d-01, 1.151d-02
c     Ln(R) :
c     $     -2.235d+00,-2.244d+00,-2.244d+00,-2.273d+00,-2.303d+00,
c     $     -2.375d+00,-2.465d+00,-2.513d+00,-2.538d+00,-2.501d+00,
c     $     -2.303d+00,-2.048d+00,-1.897d+00,-1.778d+00,-1.671d+00,
c     $     -1.556d+00,-1.427d+00,-1.262d+00,-1.171d+00,-1.097d+00,
c     $     -9.943d-01,-9.339d-01,-8.723d-01,-7.787d-01,-6.463d-01,
c     $     -5.727d-01,-4.992d-01,-4.308d-01,-3.552d-01,-2.810d-01,
c     $     -1.863d-01, 0.000d+00
c

c
c---- This is some "fit" by hand of a M --> Rhalf/Rtot relation
c     obtained from Geneva models. BUT: Rhalf was obtained for 
c     the models used in SPH simulations (Z=0.04) while Rtot
c     is for Z=0.001 models at half-MS !!! See diag_Rhalf and
c     define_rap50_cursor in study_stell_tracks.smac
c
      data NRhR /64/

      data lM_RhR_MS / ! log10 is used here !
     &     -5.58434e-01, -5.08434e-01, -4.58434e-01, -4.08434e-01,
     &     -3.58434e-01, -3.08434e-01, -2.58434e-01, -2.08434e-01,
     &     -1.58434e-01, -1.08434e-01, -5.84340e-02, -8.43395e-03,
     &     4.15660e-02, 9.15660e-02, 1.41566e-01, 1.91566e-01,
     &     2.41566e-01, 2.91566e-01, 3.41566e-01, 3.91566e-01,
     &     4.41566e-01, 4.91566e-01, 5.41566e-01, 5.91566e-01,
     &     6.41566e-01, 6.91566e-01, 7.41566e-01, 7.91566e-01,
     &     8.41566e-01, 8.91566e-01, 9.41566e-01, 9.91566e-01,
     &     1.04157e+00, 1.09157e+00, 1.14157e+00, 1.19157e+00,
     &     1.24157e+00, 1.29157e+00, 1.34157e+00, 1.39157e+00,
     &     1.44157e+00, 1.49157e+00, 1.54157e+00, 1.59157e+00,
     &     1.64157e+00, 1.69157e+00, 1.74157e+00, 1.79157e+00,
     &     1.84157e+00, 1.89157e+00, 1.94157e+00, 1.99157e+00,
     &     2.04157e+00, 2.09157e+00, 2.14157e+00, 2.19157e+00,
     &     2.24157e+00, 2.29157e+00, 2.34157e+00, 2.39157e+00,
     &     2.44157e+00, 2.49157e+00, 2.54157e+00, 2.59157e+00 /

      data RhR_MS / ! set to n=1.5 polytrope value at lower limit and n=3 at upper limit
     &     5.21020e-01, 5.16224e-01, 5.10841e-01, 5.05122e-01,
     &     4.98558e-01, 4.86530e-01, 4.62237e-01, 4.25628e-01,
     &     3.85191e-01, 3.51476e-01, 2.89787e-01, 2.51520e-01,
     &     2.51674e-01, 2.57604e-01, 2.65348e-01, 2.74275e-01,
     &     2.83752e-01, 2.93146e-01, 3.01827e-01, 3.09343e-01,
     &     3.15721e-01, 3.21072e-01, 3.25506e-01, 3.29131e-01,
     &     3.32057e-01, 3.34393e-01, 3.36248e-01, 3.37733e-01,
     &     3.38955e-01, 3.39976e-01, 3.40775e-01, 3.41325e-01,
     &     3.41597e-01, 3.41562e-01, 3.41193e-01, 3.40460e-01,
     &     3.39336e-01, 3.37791e-01, 3.35802e-01, 3.33417e-01,
     &     3.30747e-01, 3.27906e-01, 3.25008e-01, 3.22168e-01,
     &     3.19470e-01, 3.16917e-01, 3.14496e-01, 3.12194e-01,
     &     3.09996e-01, 3.07891e-01, 3.05864e-01, 3.03903e-01,
     &     3.01994e-01, 3.00133e-01, 2.98313e-01, 2.96532e-01,
     &     2.94786e-01, 2.93069e-01, 2.91379e-01, 2.89711e-01,
     &     2.88062e-01, 2.86426e-01, 2.84801e-01, 2.83194e-01 /

      end
#endif /* (_EVOL_STELL_<10) */ /* Simple stellar evol (MS-->remnant) or no stellar evolution */
