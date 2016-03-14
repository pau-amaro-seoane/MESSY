c======================================================================
c     routines liees a la relaxation pour code EvolAmasXXX
c======================================================================
c
c----------------------------------------------------------------------
      subroutine CalcParam_Relax
c----------------------------------------------------------------------
      include 'EvolAmas_Relax_inc.f'
c......................................................................
c     
#if (_RELAX_NULLE_>0)
      write(0,110)
 110  format(
     $     ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,
     $     ' %% !!!!! RELAXATION NEUTRALISEE !!!!! %%',/,
     $     ' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/ )
#endif
      Correc_Log_Coul_kicks=1.0d0
#if (_KICKS_>0)
      if (fact_b0_kick.GT.1.0d0)
     $     Correc_Log_Coul_kicks=1.0d0/fact_b0_kick
#endif
      Log_Coulomb = LOG(Correc_Log_Coul_kicks*Gamma_relax*
     $     ((NbSE-NbSE_perte)*Net_par_SE))
                                ! Ceci n'est correct que pour un amas AUTOGRAVITANT
      Cste_Relax = 64.0d0/rPi*Log_Coulomb
      Met_moy = M_amas/((NbSE-NbSE_perte)*Net_par_SE)
c     
      end
c     
c----------------------------------------------------------------------
      subroutine CalcGlobRelax
c----------------------------------------------------------------------
      include 'EvolAmas_Relax_inc.f'
c
c constante :
c ^^^^^^^^^^^
      double precision r8Pi3
      parameter (r8Pi3=rPi*8.0d0/3.0d0)
      double precision Frac_Lim_Masse
      parameter (Frac_Lim_Masse=0.99d0)
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      integer iSE,lSE
      double precision M_r,Mi,Ri,Rc_prov,LC1,Demie_Masse,Lim_Masse
c......................................................................
c
      if (Log_Coulomb.EQ.0.0d0) call CalcParam_Relax
      Lim_Masse   = Frac_Lim_Masse*M_amas
      Demie_masse = 0.5d0*M_amas
      TrelMoy  = 0.0d0
      Trel_dt  = 0.0d0
      V2moy    = 0.0d0
      M_r      = 0.0d0
      Rh       = -1.0d0
      Rc       = -1.0d0
      do lSE=1,NbSE-NbSE_perte
         iSE = iListeRay(lSE)
         Mi = M_SE(iSE)
         Ri = R_SE(iSE)
         if (M_r.LT.Lim_Masse) then
            V2moy = V2moy + Mi*T_SE(iSE)
            M_r = M_r+Mi
            if (Rh.LE.0.0d0) then
               if (M_r.GE.0.5d0*M_amas) Rh = Ri
            end if
            if (Rc.LE.0.0d0) then
               Rc_prov = sqrt(V2moy*Ri**3/M_r**2)
               if (Rc_prov.LE.2.0d0*Ri) then
                  Rc = Rc_prov      ! AUDACIEUX !!!
                  V2c = V2moy/M_r
                  Rhoc = 3.0d0*M_r/(r4Pi*Ri**3)
               end if
            end if
         end if
      end do
      
      if (Rc.LE.0.0d0) then
         Rc = Rc_prov               ! AUDACIEUX !!!
         V2c = V2moy/M_r
         Rhoc = 3.0d0*M_r/(r4Pi*Ri**3)
      end if
      
      V2moy        = 2.0d0*V2moy/M_amas
      LC1          = 1.0d0/Log_Coulomb
      Trel_c       = a_spitzer * LC1 *
     $     V2c**1.5/(Rhoc*Met_moy)
      Trel_Rh      = a_spitzer * LC1 *
     $     V2moy**1.5*r8Pi3*Rh**3/(M_amas*Met_moy)
      Trel_Rh_spit = b_spitzer * LC1 *
     $     ((NbSE-NbSE_perte)*Net_par_SE)*sqrt(Rh**3/M_amas)
      if (Trel_Rh0.EQ.0.0d0) Trel_Rh0 = Trel_Rh
      Trel_dt  = Trel_dt/M_r
      q_Viriel = Rh*V2moy/M_r
      if (M_renc_tot.GT.0.0d0) then
         Trel_CM_moy  = Trel_CM_moy/(Frac_Trelax*M_renc_tot)
         Trel_tst_moy = Trel_tst_moy/(Frac_Trelax*M_renc_tot)
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_FichRelax_ASC(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Relax_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
      write(iFich,100)
 100  format('# Donnees globales de relaxation')
      write(iFich,110)
     $     '1:Rh','2:Rc','3:Tr_dt','4:Tr_CM','5:Tr_tst',
     $     '6:Tr_Rh','7:Tr_Rh_spit','8:Tr_c','9:V2c',
     $     '10:V2moy','11:Rhoc','12:q_Viriel' 
 110  format('#',12(A11," "))
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrire_Relax_ASC(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Relax_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c 
      write(iFich,100)
     $     Rh,Rc,Trel_dt,Trel_CM_moy,Trel_tst_moy,Trel_Rh,
     $     Trel_Rh_spit,Trel_c,V2c,V2moy,Rhoc,q_Viriel
 100  format(12(E11.4,' '))
c
      end
c
c----------------------------------------------------------------------
      double precision function TrelLoc(R)
c----------------------------------------------------------------------
      include 'EvolAmas_Relax_inc.f'
      include 'EvolAmas_Grille_param.f'
      include 'EvolAmas_Grille_common.f'
c
c argument :
c ^^^^^^^^^^
      double precision R
c
c constante :
c ^^^^^^^^^^^
      double precision Rac3
      parameter (Rac3=1.73205d0)
      double precision coefMax_Trel
      parameter (coefMax_Trel=100.0d0)
c
c var locales :
c ^^^^^^^^^^^^^
      integer iPG
      double precision Met_moy_loc
c
c fonctions :
c ^^^^^^^^^^^
      integer iCouche_PG
      double precision Rho_PG_lin, Sigma1D_PG_lin2
c......................................................................
c
      iPG = iCouche_PG(R)
      Met_moy_loc = M_PG(iPG)/(Net_par_SE*NSE_PG(iPG))
c
c     La valeur calculee ne tient pas compte d'une variation possible
c     du logarithme coulombien en fct du rayon. Ceci rend cette valeur
c     particulierement douteuse dans le cas d'un TN dominant
c
      TrelLoc = (Rac3*Sigma1D_PG_lin2(iPG,R))**3 /
     $     (Cste_Relax*Met_moy_loc*Rho_PG_lin(iPG,R))
c
      end
c
c----------------------------------------------------------------------
      double precision function Coef_Relax(iRang)
c----------------------------------------------------------------------
      include 'EvolAmas_Relax_inc.f'
      include 'EvolAmas_TN_common.f'
#if (_TYPE_LOG_COUL_==1) /* logarithme coulombien local */
      include 'EvolAmas_LambdaRelax_common.f'
#endif
c
c argument :
c ^^^^^^^^^^
      integer iRang
c
c constante :
c ^^^^^^^^^^^
      double precision coef_log
      parameter (coef_log=64.0d0/rPi)
c
c var locale :
c ^^^^^^^^^^^^
      integer i
      double precision arg_log
c......................................................................
c
#if (_TYPE_LOG_COUL_==1) /* logarithme coulombien local */
      i = iRang*real(N_Lambda_relax)/real(NbSE-NbSE_perte)
      i = min(max(1,i),N_Lambda_relax)
      arg_log = 1.0d0 + ( Correc_Log_Coul_kicks*Coef_Lambda*Gamma_relax*
     $     Lambda_relax(i) )**2 ! pour eviter des arguments du log < 1
      Coef_Relax = 0.5d0*coef_log*LOG(arg_log)
                                ! Lambda est le rapport Torb/Tmin ou Tmin
                                ! est une duree effective minimale pour une
                                ! rencontre 2 corps
#elif (_TYPE_LOG_COUL_==2) /* logarithme coulombien local comme Marchant et Shapiro 80 Eq. 8 */
      arg_log = Gamma_relax*iRang*Net_par_SE ! MS80 prennent le nb d'etoiles a l'interieur de l'apocentre !
     $     + M_ctr/Met_moy
      arg_log = Correc_Log_Coul_kicks*arg_log
      Coef_Relax = coef_log*LOG(max(1.0d0,arg_log))
#else /* logarithme coulombien global et constant */
      Coef_Relax = Cste_Relax
#endif
c
      end
c
c----------------------------------------------------------------------
      block data BD_Init_EvolAmas_Relax
      include 'EvolAmas_Relax_common.f'
      data
     $     Log_Coulomb, Cste_Relax, Met_moy, TSurTrel,
     $     Rh,Rc,TrelMoy,Trel_dt,Trel_Rh,Trel_Rh0,Trel_Rh_spit,Trel_c,
     $     V2c,V2moy,Rhoc,q_Viriel,M_renc_tot,Trel_CM_moy,Trel_tst_moy
     $     /19*0.0d0/
      end
c----------------------------------------------------------------------
