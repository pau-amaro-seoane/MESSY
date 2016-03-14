c======================================================================
c     Routines to implement accretion of gas by stars, i.e. increase
c     of mass of the stars, in ME(SSY)**2.
c     Written by Marc Freitag for a part-3 project at the IoA
c     Started 23/10/07
c     Last update 6/12/07
c======================================================================
c
c---------------------------------------------------------------------------------------------------
c                              SOME EXPLANATIONS (6/12/07)
c---------------------------------------------------------------------------------------------------
c
c     Important subroutines to tinker with :
c   =========================================
c
c     * subroutine get_dMaccr_dt_star_MsunMyr(Mstar_Msun,time_Myr,dMdt)
c                  ^^^^^^^^^^^^^^^^^^^^^^^^^^
c     Computes the accretion rate for a star of given mass at given time
c     Unit of mass is Msun, unit of time Myr.
c     In the current implementation, it only depends on Mstar (and the total
c     accretion rate by the whole cluster, see below)
c     
c     * double precision function dMaccrdt_clust_MsunMyr(time_Myr)
c                                 ^^^^^^^^^^^^^^^^^^^^^^
c     Returns the total or maximum accretion rate for the whole cluster.
c     In the current implementation, it is just a constant.
c     If the value is negative, MESSY will interpret its absolute value
c     as a _maximum_ total rate. If it is positive, MESSY will try to 
c     enforce it "exactly". In order for the rate returned by 
c     `get_dMaccr_dt_star_MsunMyr' to comply with this constrain, the 
c     quantity `tAccrByStar_1Msun_Myr' is recomputed every time the routine
c     `Set_AccretionRateOnStars' is called, which is currently done every
c     time the timesteps are (re)computed by `Calc_Param_Tirage_Henon' 
c     (in EvolAmas_TirageHenon3.F).
c
c     Currently, the accretion-by-star business is parametrised by the 
c     following quantities, whose values can be set in `input_EvolAmas'.
c     Values in [] are default:
c
c     RateAccrPerStar_MsunMyr [1.0] : Average rate of accretion by stars in Msun/Myr; 
c                                     Used by `dMaccrdt_clust_MsunMyr'.
c                                     if negative, will be used as a maximum. So set it 
c                                     to, say, -1d30 if you don't want to put any constraint
c                                     on the total (whole-cluster) accretion rate.
c
c     tAccrByStar_1Msun_Myr [0.25] :  Parameter for accretion by stars. 
c                                     Used by `get_dMaccr_dt_star_MsunMyr'.
c                                     Formally, M/(dM/dt) for a 1 Msun star.
c                                     Recomputed and updated by MESSY if RateAccrPerStar_MsunMyr>0
c                                     or RateAccrPerStar_MsunMyr<0 and maximum total rate is reached.
c
c     exp_AccrByStar [0.0] :          Parameter for accretion by stars. 
c                                     Used by `get_dMaccr_dt_star_MsunMyr'.
c                                     Basically 
c                                                
c                                     dM/dt = 1Msun/tAccrByStar_1Msun_Myr * (M/1Msun)**exp_AccrByStar 
c
c                                     for M < Mmax_AccrByStar
c
c     Mmax_AccrByStar [80.0] :        Parameter for accretion by stars
c                                     Used by `get_dMaccr_dt_star_MsunMyr'.
c                                     Accretion rate becomes very small as this mass is
c                                     reached and exceeded.
c
c     WARNING: I haven't tested whether MESSY resumes correctly from a stop/kill 
c              (e.g., with gere_evolamas.pl -Continue)
c              when stellar accretion is enable. I guess it should. 
c---------------------------------------------------------------------------------------------------


#if (_ACCR_BY_STARS_>0)
#define _VERB_ACCR_BY_STARS_ 2 /* Set to >0 if you want a few messages from Set_AccretionRateOnStars */

c to compile test program (with g77):
c g77 -pipe  -ffixed-line-length-none -D_ACCR_BY_STARS_=1 -D_TEST_ACCR_BY_STARS_=1 EvolAmas_AccretionByStars.F RK2dbl.f

#   if (_TEST_ACCR_BY_STARS_>0) /* For test program */
c--------------------------------------------------
      implicit none

      double precision mass
      common /mass_common/ mass
      double precision m0,dt,t_end,t,m,dm,cdt, taccr, mEuler

      double precision Maccr_star_Msun

      logical lEuler
      common /TypeStep_Maccr/ lEuler

      write(0,*) 'Enter m(0),cdt,t_end'
      read(*,*) m,cdt,t_end
      mass = m
      mEuler = m
      
      write(6,110)
 110  format(
     $     '# 1: t_MYr 2: M_Msun 3: taccr_Myr 4: M_Euler 5: UsedEuler')
      t=0.0d0
      do while(t.LT.t_end)
         taccr = Maccr_star_Msun(1,t,-1.0d0)
         dt = cdt*taccr
         dm = Maccr_star_Msun(1,t,dt)
         m = m + dm
         mEuler = mEuler + mEuler/taccr*dt
         mass = m
         t=t+dt
         write(6,100) t,m, taccr, mEuler, lEuler
      end do
 100  format(4(e12.5,' '),L1)

      end

c--------------------------------------------------
      double precision function stellar_mass(i)
c--------------------------------------------------
      implicit none
      integer i
      double precision mass
      common /mass_common/ mass
      
      stellar_mass = mass
      end
c--------------------------------------------------
#   endif /* (_TEST_ACCR_BY_STARS_>0) */ /* END TEST SECTION */
c
c----------------------------------------------------------------------
      subroutine get_dMaccr_dt_star_MsunMyr(Mstar_Msun,time_Myr,dMdt)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_Param_common.f'

c
c arguments : 
c ^^^^^^^^^^^
      ! inputs:
      double precision Mstar_Msun ! Current mass in M_sun
      double precision time_Myr ! current time in Myr (in case we implement some explicit time-dependence)
      ! outputs:
      double precision dMdt     ! accretion rate in Msun/Myr
c
c constants :
c ^^^^^^^^^^^
      integer exp_quench
      parameter (exp_quench=10) ! random choice to get fast decrese of accretion
                                ! rate past the maximum mass
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision quench
c......................................................................
c
c---- Let's do something simple for now.
c     Only depends on initial mass. Not on position in cluster
c
c     Note that the value of tAccrByStar_1Msun_Myr might be manipulated by
c     Set_AccretionRateOnStars (in this file), in particular via
c     Calc_Tcarac (in EvolAmas_TirageHenon3.F), to enforce total accretion rate
c     for the whole cluster
c
      if (Mstar_Msun.LT.100.0d0*Mmax_AccrByStar) then
         dMdt = Mstar_Msun**exp_AccrByStar/tAccrByStar_1Msun_Myr
                                ! Stars more massive than a given limit do not accrete anymore
         quench = 1.0d0/((Mstar_Msun/Mmax_AccrByStar)**exp_quench+1.0d0)
         dMdt = quench*dMdt
      else
         dMdt = 1.0d-30
      end if
c
      end
c
c----------------------------------------------------------------------
      double precision function dMaccrdt_clust_MsunMyr(time_Myr)
c----------------------------------------------------------------------
c     Total accretion rate by whole cluster, in Msun/Myr.
c     If negative, the abs value will be used as a maximum.
c
      implicit none
c
c arguments : 
c ^^^^^^^^^^^
      double precision time_Myr ! current time in Myr (in case we implement some explicit time-dependence)
      include 'EvolAmas_Conserv_inc.f'
      include 'EvolAmas_Param_common.f'
c......................................................................
c
      dMaccrdt_clust_MsunMyr = RateAccrPerStar_MsunMyr*rNbEtoiles_ini
c
      end
c
c----------------------------------------------------------------------
      double precision function Maccr_star_Msun(iSE,tini_Myr,dt_Myr)
c----------------------------------------------------------------------
c     Returns the increase of stellar mass (in Msun) for particle iSE, 
c     due to accretion during the time interval
c     [tini_Myr,tini_Myr+dt_Myr]
c
c     If function is called with negative dt_Myr, return the current
c     accretion timescale in Myr
c
      implicit none
      include 'EvolAmas_Param_common.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iSE               ! the ID number of the particle 
                                ! (to get position, orbital & stellar properties, as required)
      double precision tini_Myr, dt_Myr ! initial time and time interval, in Myr
c
c constants :
c ^^^^^^^^^^^
      double precision frac_dtmax
      parameter (frac_dtmax=0.01d0)
      double precision coef_dt_Euler
      parameter (coef_dt_Euler=1.0d3)
      double precision accur
      parameter (accur=1.0d-3)
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision Mstar_Msun, dt_ini, taccr, Mend, dMdt
      integer Nsteps
      logical lEuler
      common /TypeStep_Maccr/ lEuler
c
c functions :
c ^^^^^^^^^^^
      double precision stellar_mass, dMaccr_dt_star_MsunMyr
      external get_dMaccr_dt_star_MsunMyr
c
c......................................................................
c
      Mstar_Msun = stellar_mass(iSE)
      lEuler = .FALSE.
ccc      R_SE_AccrByStar = R_SE(iSE) ! to pass "hidden" parameters (determining accretion rate) to 
ccc                                  ! "get_dMaccr_dt_star_MsunMyr" through common block (not used yet)
c
c---- Have we already reached a mass above
c     which no significant accretion should occur?
c
      if (Mstar_Msun.GT.10.0d0*Mmax_AccrByStar) then
         Maccr_star_Msun = 0.0d0
         return
      end if
c
c---- If time step very short compared to accretion time, 
c     just a simple explicit Euler step...
c
      call get_dMaccr_dt_star_MsunMyr(Mstar_Msun,tini_Myr,dMdt)
      taccr = Mstar_Msun/dMdt
      if (dt_Myr.LT.0.0d0) then ! Special case: return the current timescale
         Maccr_star_Msun = taccr
         return
      end if

      if (taccr.GT.coef_dt_Euler*dt_Myr) then
         Maccr_star_Msun = Mstar_Msun/taccr*dt_Myr
         lEuler = .TRUE.
      else
c
c---- Let's not be afraid of a bit of overkill and go for a full Runge-Kunta
c     integration for now
c
         dt_ini = min(dt_Myr,frac_dtmax*taccr)
         call IntegrateOverInterval_RKDOPRI54
     $        (Mstar_Msun,1,get_dMaccr_dt_star_MsunMyr,
     $        tini_Myr,tini_Myr+dt_Myr,dt_ini,accur,Mend,Nsteps)
         Maccr_star_Msun = Mend-Mstar_Msun
      end if
c
      end
c     
c----------------------------------------------------------------------
      subroutine Set_AccretionRateOnStars
c----------------------------------------------------------------------
c     Determine value of tAccrByStar_1Msun_Myr to enforce total 
c     accretion rate for the whole cluster
c
      implicit none
      include 'EvolAmas_Param_common.f'
      include 'Param_Amas.f'
      include 'VarSE.f'
      include 'EvolAmas_Dimensions_common.f'
      include 'EvolAmas_Pertes_inc.f'
      include 'EvolAmas_Temps_common.f'
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision 
     $     taccr, AccrByStars_TotRate, AccrByStars_Constr,
     $     Correction_AccrByStars
      double precision mstar_mean
      integer iSE, Nstars
c
c fonctions :
c ^^^^^^^^^^^
      double precision stellar_mass,
     $     Maccr_star_Msun, dMaccrdt_clust_MsunMyr
c......................................................................
c
c---- Compute "nominal" total accretion rate
c
      if (tAccrByStar_1Msun_Myr.LE.0.0d0) tAccrByStar_1Msun_Myr = 1.0d0

      Nstars = 0
      Mstar_mean = 0.0d0
      AccrByStars_TotRate = 0.0d0
      do iSE=1,NbSE
         if (R_SE(iSE).LT.R_perte) then
            Nstars = Nstars+1
            Mstar_mean = Mstar_mean + stellar_mass(iSE)
                                ! Maccr_star_Msun returns accretion time 
                                ! in Myr when 3rd argument is negatif!
            taccr = 1.0d6*UneAnnee*
     $           Maccr_star_Msun(iSE,1.0d-6*Tps_Amas/UneAnnee,-1.0d0)
            AccrByStars_TotRate = AccrByStars_TotRate +
     $           stellar_mass(iSE)/max(taccr,1.0d-10)
         end if
      end do 
      Mstar_mean = Mstar_mean/Nstars
c
c---- Get the value of the max/imposed value
c
      AccrByStars_Constr =
     $     dMaccrdt_clust_MsunMyr(1.0d-6*Tps_Amas/UneAnnee)

#if (_VERB_ACCR_BY_STARS_>1)
      write(0,*) ']] Accretion by stars : Mstar_mean = ',
     $     Mstar_mean
      write(0,*) ']] Accretion by stars : AccrByStars_TotRate = ',
     $     AccrByStars_TotRate
      write(0,*) ']] Accretion by stars : AccrByStars_Constr = ',
     $     AccrByStars_Constr
#endif
      if (AccrByStars_Constr.LT.0.0d0) then ! It is a maximum

         Correction_AccrByStars = min(1.0d0,
     $        -AccrByStars_Constr/AccrByStars_TotRate)

      else                                  ! It is a value to enforce

         Correction_AccrByStars = 
     $        AccrByStars_TotRate/AccrByStars_Constr

      end if
c
c---- Determine/correct tAccrByStar_1Msun_Myr to comply
c     with constrain on total rate
c
      if (abs(Correction_AccrByStars-1.0d0) .GT. 1d-6) then
         tAccrByStar_1Msun_Myr = Correction_AccrByStars *
     $        tAccrByStar_1Msun_Myr
      end if
#if (_VERB_ACCR_BY_STARS_>0)
      write(0,*) ']] Accretion by stars : tAccrByStar_1Msun_Myr = ',
     $     tAccrByStar_1Msun_Myr
#endif
c
      end
c

#endif /* (_ACCR_BY_STARS_>0) */
