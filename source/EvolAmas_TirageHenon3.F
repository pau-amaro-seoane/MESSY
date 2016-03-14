#define _VERBOSITE_TH_ 10
#define _PONDERATION_
#define _TEST_TH_ 0 /* Test des procedures/ ecriture des donnees liees a ces procedures de tirage */
#define _TEST_PR_TH_ 0
#define _CALC_PR_ 1
#define _ECRIT_TADIAB_ 0
#define _SURV_T_EVAP_ 0
#define _TYPE_DT_TH_ 2
c======================================================================
c     routine pour faire les tirages des SE et de leurs positions
c     selon methode de Henon
c======================================================================
c
c----------------------------------------------------------------------
      subroutine Calc_Tcarac
     $     (Trel,Tcoll,Tkick,Torb,Tadiab,Taccr,Tevst,Tstaccr)
c----------------------------------------------------------------------
c     calcule les temps caracteristiques locaux pour chaque SE
c     en faisant des moyenne glissantes pour obtenir les proprietes
c     locales du milieu stellaire
c
c     In case of accretion by stars, also determines the coefficient
c     to comply with the constrain on total accretion rate by whole
c     cluster (5/12/07)
c
      include 'EvolAmas_TirageHenon3_inc.f'
      include 'EvolAmas_param.f'
      include 'EvolAmas_Relax_common.f'
      include 'EvolAmas_common.f'
      include 'EvolAmas_Param_common.f'
      include 'EvolAmas_TN_common.f'
      include 'EvolAmas_Conserv_inc.f'
#if (_TYPE_LOG_COUL_==1) /* logarithme coulombien local */
      include 'EvolAmas_LambdaRelax_common.f'
#endif
c
c arguments :
c ^^^^^^^^^^^
      double precision
     $     Trel(iDimSE),Tcoll(iDimSE),Tkick(iDimSE),Torb(iDimSE),
     $     Tadiab(iDimSE),Taccr(iDimSE),Tevst(iDimSE),Tstaccr(iDimSE)
c
c constantes :
c ^^^^^^^^^^^^
      double precision QuatreTiersPi,SeizeRacPi, DeuxPi, Rac3, UnSurRac3
      parameter ( QuatreTiersPi=4.0d0/3.0d0*rPi, SeizeRacPi=28.35926d0,
     $     DeuxPi=2.0d0*rPi, Rac3=1.7320508d0, UnSurRac3=1.0d0/Rac3 )
      double precision fact_Mint_min
      parameter (fact_Mint_min=3.0d-5)
c
c var locales :
c ^^^^^^^^^^^^^
      integer iRang,iSE,i,i_inf,i_sup,i_inf_,i_sup_,Nb_SE,NbSE_tranche,
     $     iPG
      double precision SigmaV,Dens,v3_sur_m2n,
     $     Mstel_moy,M2stel_moy,MRstel_moy,R2stel_moy,
     $     S_m,S_m2,S_mv2,
     $     S_Mstel_cpct,S_R2stel_cpct,S_MRstel_cpct,
     $     S_Mstel_non_cpct,S_R2stel_non_cpct,S_MRstel_non_cpct,
     $     mi,Mint, c_rel, c_kick,
     $     Ret, Met, A, B, Xi, Poti, Pprec, deltaT, deltaM_ctr, dM

      double precision Tprec, Mctr_prec
      save Tprec, Mctr_prec
      data Tprec, Mctr_prec /-1.0d0,-1.0d0/

      logical lPrem
      save lPrem
      data lPrem /.TRUE./

      integer i_Calc_Tcarac
      save i_Calc_Tcarac
      data i_Calc_Tcarac /0/

#if (_SURV_ADIAB_>0)
      double precision lRay(iDimSE), lPot(iDimSE), Y2spl(iDimSE), Xmaxspl
      integer Nspl
      parameter (Nspl=100)
      double precision Yspl(Nspl),dXspl,XminSpl
      save Yspl,dXspl,XminSpl
      character*32 FileName
#endif
c
c fonctions :
c ^^^^^^^^^^^
      double precision Coef_Relax, RelMR, T_MS,
     $     Rho_PG_lin2, Somme_MsurR_Arbre
#if (_TEST_TH_>0)
      double precision TrelLoc
#endif
#if (_FORCER_CROISSANCE_TN_>0)
      double precision M_TN_Croiss_Forc
#endif
#if (_ACCR_BY_STARS_>0)
      double precision Maccr_star_Msun
#endif
      double precision stellar_mass, stellar_radius,
     $     stellar_evoltimescale
      logical stellar_is_compact
c......................................................................
c
                                ! initialisations
#if (_TEST_TH_==1)
      write(0,*) '####### _TEST_TH_=1'
#elif (_TEST_TH_==2)
      write(0,*) '####### _TEST_TH_=2'
#elif (_TEST_TH_==3)
      write(0,*) '####### _TEST_TH_=3'
#endif
#if (_TEST_TH_==1)
      open(77,file='evolamas_TEST_TCARAC.asc',form='formatted')
      write(77,100)
#   if (_TYPE_LOG_COUL_==1) /* logarithme coulombien local */
 100  format('# 1: r 2: trel1 3: tcoll 4: torb 5: tevst 6: trel2',
     $     ' 7: sigmav 8: dens 9: Mstell_moy 10: Lambda',
     $     ' 11: coef_relax')
#   else
 100  format('# 1: r 2: trel1 3: tcoll 4: torb 5: tevst 6: trel2',
     $     ' 7: sigmav 8: dens 9: Mstell_moy 10:Tadiab 11:Taccr')
#   endif
#endif

#if (_ACCR_BY_STARS_>0)
      call Set_AccretionRateOnStars
#endif

#if (_FORCER_CROISSANCE_TN_>0)
      dM = -M_TN
      M_TN = M_TN_Croiss_Forc(Tps_amas)
      dM = dM + M_TN
      M_ctr = M_ctr + dM
#endif

      if (lPrem) then
         Tprec=Tps_amas
         Mctr_prec=M_ctr
      end if

      i_Calc_Tcarac = i_Calc_Tcarac + 1
      
      deltaT = Tps_amas - Tprec
      deltaM_ctr = M_ctr - Mctr_prec
      c_kick = rPi*fact_b0_kick**2
      NbSE_Tranche = Fac_SldAvg_Tscale*NbSECouchePot
      S_m = 0.0d0
      S_m2 = 0.0d0
      S_mv2 = 0.0d0
      S_R2stel_non_cpct = 0.0d0
      S_MRstel_non_cpct = 0.0d0
      S_R2stel_cpct = 0.0d0
      S_MRstel_cpct = 0.0d0
      S_Mstel_non_cpct = 0.0d0
      S_Mstel_cpct = 0.0d0
#if (_TYPE_LOG_COUL_==1) /* logarithme coulombien local */
      N_Lambda_relax = NbSE-NbSE_perte
#endif
                                ! initialisation des moyennes glissantes
      iRang=1
      iSE=iListeRay(iRang)
      i_inf = iRang-NbSE_Tranche/2
      i_sup = i_inf+NbSE_Tranche-1
      i_inf_ = max(1,i_inf)
      i_sup_ = min(NbSE-NbSE_perte,i_sup)
      Mint = M_ctr

      do i=i_inf_,i_sup_
         iSE = iListeRay(i)
         mi = M_SE(iSE)
         S_m   = S_m+mi
         S_m2  = S_m2+mi*mi
         S_mv2 = S_mv2+2.0d0*mi*T_SE(iSE)
         Ret = stellar_radius(iSE)
         Met = stellar_mass(iSE)
         if (.NOT.stellar_is_compact(iSE))
     $        then              ! etoile non compacte
            S_Mstel_non_cpct  = S_Mstel_non_cpct +Met
            S_MRstel_non_cpct = S_MRstel_non_cpct+Met*Ret
            S_R2stel_non_cpct = S_R2stel_non_cpct+Ret*Ret
         else                   ! etoile compacte
            S_Mstel_cpct  = S_Mstel_cpct +Met
            S_MRstel_cpct = S_MRstel_cpct+Met*Ret
            S_R2stel_cpct = S_R2stel_cpct+Ret*Ret
         end if
      end do
c
c==== Calcul des temps caracateristiques pour toutes les SE
c
      do iRang=1,NbSE-NbSE_perte+1
c
c---- temps carac pour la SE courante
c
                                ! Les moyennes glissantes

         Nb_SE = i_sup_-i_inf_+1 ! nb de SE dans la tranche employee pour la moyenne glissante
         SigmaV = UnSurRac3*sqrt(S_mv2/S_m)
         Dens = (i_sup_-i_inf_)*Net_par_SE / ( QuatreTiersPi*
     $        (R_SE(iListeRay(i_sup_))**3-R_SE(iListeRay(i_inf_))**3) )
         Mstel_moy  = S_m/(Nb_SE*Net_par_SE)
         M2stel_moy = S_m2/(Nb_SE*Net_par_SE**2)
         if (S_Mstel_non_cpct.GT.0.0d0) then
            MRstel_moy = Coef_ConvM*Coef_ConvR*S_MRstel_non_cpct/Nb_SE
            R2stel_moy = Coef_ConvR*Coef_ConvR*S_R2stel_non_cpct/Nb_SE
         else
            MRstel_moy = Coef_ConvM*Coef_ConvR*S_MRstel_cpct/Nb_SE
            R2stel_moy = Coef_ConvR*Coef_ConvR*S_R2stel_cpct/Nb_SE
         end if

         iSE = iListeRay(iRang)
         Mint = Mint + M_SE(iSE)
         Torb (iRang) = DeuxPi*
     $        sqrt(R_SE(iSE)**3/(Mint+FracMasseInt*M_SE(iSE)))
#if (_TYPE_LOG_COUL_==1) /* logarithme coulombien local */
         Lambda_relax(iRang) = SigmaV**3*Torb(iRang) / Mstel_moy
#endif
         c_rel = Coef_Relax(iRang)
         v3_sur_m2n = (Rac3*SigmaV)**3/(Dens*M2stel_moy)
c
c     RELAXATION
c     ----------
         if (c_rel.LT.1.0d-30) then
            Trel(iRang) = 1.0d50
         else
            Trel(iRang) = v3_sur_m2n/c_rel
         end if
c
c     LARGE ANGLE SCATTERING
c     ----------------------
         if (c_kick.LT.1.0d-30) then
            Tkick(iRang) = 1.0d50
         else
            Tkick(iRang) = v3_sur_m2n/c_kick
         end if
c
c     COLLISIONS
c     ----------
         Tcoll(iRang) = 1.0d0/( SeizeRacPi*Dens*SigmaV*R2stel_moy *
     $        (1.0d0+MRstel_moy/(2.0d0*SigmaV**2*R2stel_moy)) )
c
c     CENTRAL MASS INCREASE
c     ---------------------
         if (deltaM_ctr.GT.1.0d-20) then
            Taccr(iRang) = deltaT / deltaM_ctr *
     $           max(Mint,fact_Mint_min*Mamas_ini)
         else
#if (_FORCER_CROISSANCE_TN_>0)
            Taccr(iRang) = 1.0d-2*abs(T_TN_CF * Mint / M_TN_fin_CF)
#else
            Taccr(iRang) = 1.0d30 ! pas de limite effective
#endif
         end if
#if (_FORCER_CROISSANCE_TN_>0)
         Taccr(iRang) = min(Taccr(iRang),abs(T_TN_CF))
#endif
c
c     STELLAR EVOLUTION
c     -----------------

#if (_EVOL_STELL_>0)
         Tevst(iRang) = UneAnnee*stellar_evoltimescale(iSE)
#else
         Tevst(iRang) = 1.0d30
#endif
c
c     ACCRETION BY STARS
c     ------------------

#if (_ACCR_BY_STARS_>0)
                                ! Maccr_star_Msun returns accretion time 
                                ! in Myr when 3rd argument is negatif!
         Tstaccr(iRang) = 1.0d6*UneAnnee*
     $        Maccr_star_Msun(iSE,1.0d-6*Tps_Amas/UneAnnee,-1.0d0)
#else
         Tstaccr(iRang) = 1.0d30
#endif

#if (_TEST_TH_==1)
#   if (_TYPE_LOG_COUL_==1) /* logarithme coulombien local */
         write(77,110) R_SE(iListeRay(iRang)),Trel(iRang),
     $        Tcoll(iRang),Torb(iRang),Tevst(iRang),TrelLoc(R_SE(iSE)),
     $        SigmaV,Dens,
     $        Mstel_moy/Coef_ConvM, Lambda_relax(iRang), c_rel
 110     format(11(E12.5,' '))
#   else
         write(77,110) R_SE(iListeRay(iRang)),Trel(iRang),
     $        Tcoll(iRang),Torb(iRang),Tevst(iRang),TrelLoc(R_SE(iSE)),
     $        SigmaV,Dens,
     $        Mstel_moy/Coef_ConvM,Tadiab(iRang),Taccr(iRang)
 110     format(11(E12.5,' '))
#   endif
#endif

         if (iRang.EQ.(NbSE-NbSE_perte)) goto 9 ! c'est fini

         i_inf = i_inf+1
         i_sup = i_sup+1
         if (i_inf.GE.1) then
            i_inf_ = i_inf
            iSE = iListeRay(i_inf)
            mi = M_SE(iSE)
            S_m = S_m - mi
            S_m2 = S_m2 - mi*mi
            S_mv2 = S_mv2-2.0d0*mi*T_SE(iSE)
            Ret = stellar_radius(iSE)
            Met = stellar_mass(iSE)
            if (.NOT.stellar_is_compact(iSE))
     $           then           ! etoile non compacte
               S_Mstel_non_cpct  = S_Mstel_non_cpct -Met
               S_MRstel_non_cpct = S_MRstel_non_cpct-Met*Ret
               S_R2stel_non_cpct = S_R2stel_non_cpct-Ret*Ret
            else                ! etoile compacte
               S_Mstel_cpct  = S_Mstel_cpct -Met
               S_MRstel_cpct = S_MRstel_cpct-Met*Ret
               S_R2stel_cpct = S_R2stel_cpct-Ret*Ret
            end if
         end if
         if (i_sup.LE.(NbSE-NbSE_perte)) then
            i_sup_ = i_sup
            iSE = iListeRay(i_sup)
            mi = M_SE(iSE)
            S_m = S_m + mi
            S_m2 = S_m2 + mi*mi
            S_mv2 = S_mv2+2.0d0*mi*T_SE(iSE)
            Ret = stellar_radius(iSE)
            Met = stellar_mass(iSE)
            if (.NOT.stellar_is_compact(iSE))
     $           then           ! etoile non compacte
               S_Mstel_non_cpct  = S_Mstel_non_cpct +Met
               S_MRstel_non_cpct = S_MRstel_non_cpct+Met*Ret
               S_R2stel_non_cpct = S_R2stel_non_cpct+Ret*Ret
            else                ! etoile compacte
               S_Mstel_cpct  = S_Mstel_cpct +Met
               S_MRstel_cpct = S_MRstel_cpct+Met*Ret
               S_R2stel_cpct = S_R2stel_cpct+Ret*Ret
            end if
         end if
      end do

 9    continue
c
c---- Echelle de temps pour changement du potentiel
c     (pour assurer reaction adiabatique correcte) 
c
#if (_SURV_ADIAB_>0)
#   if (_ECRIT_TADIAB_>0)
      write(FileName,'("Tadiab",I5.5,".asc")') i_Calc_Tcarac
      open(iFich_Prov,file=FileName)
      write(iFich_Prov,200)iPas_Evol,Tps_Amas, deltaT
 200  format(
     $     '# tps de modification du pot pour iPas_Evol = ',I10,
     $     ' Tps_Amas = ',E12.5,/,
     $     '# deltaT = ',E12.5,/,
     $     '# 1:iSE 2:R 3:Pot_prec 4:dPot_rel 5:Tadiab 6:Trel') 
#   endif
                                ! spline d'interpolation pour le potentiel precedent
      call spline_reg(Yspl,Nspl,1.0d30,1.0d30,Xminspl,dXspl,Y2spl)
      Xmaxspl = Xminspl+(Nspl-1)*dXspl
      A = 0.0d0
      B = -Somme_MsurR_Arbre()
      do iRang=1,NbSE-NbSE_perte
         iSE = iListeRay(iRang)
         B = B + M_SE(iSE)/R_SE(iSE)
         Poti = (A-FracMasseInt*M_SE(iSE))/R_SE(iSE) + B
         Xi = log(R_SE(iSE))
         if (deltaT.LE.0.0d0) then
            Tadiab(iRang) = Tadiab_ini
            Pprec=Poti
         else
            if (Xi.LT.Xminspl .OR. Xi.GT.Xmaxspl) then
               Tadiab(iRang) = 1.0d30
            else
               call splint_reg(Yspl,Y2spl,Nspl,Xminspl,dXspl,Xi,Pprec)
               Pprec=-exp(Pprec)
               Tadiab(iRang) = -deltaT*Poti/
     $              max(-1.0d-10*Poti,abs(Poti-Pprec))
            end if
         end if
#if (_ECRIT_TADIAB_>0)
         write(iFich_Prov,'(I7,5(" ",E12.5))')
     $        iSE,R_SE(iSE),Pprec,(Poti-Pprec)/Pprec,
     $        Tadiab(iRang),Trel(iRang)
#endif
         lRay(iRang) = Xi
         lPot(iRang) = log(-Poti)
         A = A - M_SE(iSE)
      end do
                                ! Calcul de tableau de potentiel par interpolation (spline)
                                ! pour memoriser le pot en fct de R jusqu'au prochain appel
      Xminspl=lRay(1)
      dXspl=(lRay(NbSE-NbSE_perte)-Xminspl)/(Nspl-1)
      call spline(lRay,lPot,Nspl,1.0d30,1.0d30,Y2spl)
      do i=1,Nspl
         call splint(lRay,lPot,Y2spl,NbSE-NbSE_perte,
     $        Xminspl+(i-1)*dXspl,Yspl(i))
      end do    
#if (_ECRIT_TADIAB_>0)
      close(iFich_Prov)
#endif
#endif /* if (_SURV_ADIAB_>0) */

      lPrem = .FALSE.
      Tprec = Tps_amas
      Mctr_prec = M_ctr

#if (_TEST_TH_==1)
      close(77)
      stop
#endif
c
c---- Remarque :
c     
c     cette methode de calcul donne des valeurs de Trel un peu differentes de celles
c     plus exactes fournies par TrelLoc. En particulier les valeurs en R grand sont
c     farfelues car la densite est mal calculee. Cela est sans importance car c'est
c     l'integrale de 1/Trel qui nous importe ici et les grandes valeurs de Trel
c     n'y contribuent presque pas
c
      end         
c
c----------------------------------------------------------------------
      subroutine Calc_Param_Tirage_Henon
c----------------------------------------------------------------------
c     L'amas doit etre trie !
c
      include 'EvolAmas_TirageHenon3_inc.f'
      include 'EvolAmas_param.f'
      include 'EvolAmas_Compteurs_common.f'
      include 'EvolAmas_Param_common.f'
      include 'EvolAmas_TN_common.f'
      include 'EvolAmas_Conserv_inc.f'
c
c constante(s) :
c ^^^^^^^^^^^^^^
      double precision c_relx,c_coll,c_kick ! contribution des differents processus dans le calcul du pas de temps
      parameter (c_relx=1.0d0, c_coll=1.0d0, c_kick=1.0d0)
      double precision Fac_limite_dt_Torb
      parameter (Fac_limite_dt_Torb=2.0d0)
c
c var locales :
c ^^^^^^^^^^^^^
      double precision Pinf,Pr_tot,Pr_i,M_i,dt,c_tot,
     $     dt_max(iDimSE),
     $     dt_evp, dmrel_evp, T_evp
      integer iRang,i,i_TH,Nsuiv
      integer Nappels_TH
      double precision Tps_TH_prec, M_TH_Prec, Mctr_TH_Prec
      data Nappels_TH /0/
      save Nappels_TH, Tps_TH_prec, M_TH_Prec, Mctr_TH_Prec
      double precision
     $     Trel(iDimSE),Tcoll(iDimSE),Tkick(iDimSE),Torb(iDimSE),
     $     Tadiab(iDimSE),Taccr(iDimSE),Tevst(iDimSE),Tstaccr(iDimSE)
#if (_ECRIT_TSCALE_>0)
      double precision
     $     Trel_c,   Tcoll_c,  Tkick_c,
     $     Trel_01,  Trel_05,  Trel_09,
     $     Tcoll_01, Tcoll_05, Tcoll_09,
     $     Tkick_01, Tkick_05, Tkick_09,
     $     Tevst_01, Tevst_05, Tevst_09,
     $     Tadiab_01,Tadiab_05,Tadiab_09,
     $     Taccr_01,Taccr_05,Taccr_09
#endif
c
c fonction :
c ^^^^^^^^^^
      double precision P_TH, selip
#if (_TEST_TH_>0)
      double precision  PasTps_TH,h_TH
#endif

#if (_TEST_PR_TH_==1)
      integer Ntir_tot
      double precision UnSurNtot
      character*40 NomFichTestPr_TH
      integer Ntir_TH(Nmax_TH)
      common /common_test_Pr_TH/ Ntir_TH
#endif
c......................................................................
c
      Nappels_TH = Nappels_TH+1

#if (_TEST_PR_TH_==1) 
                                ! on ecrit les probabilites precedentes
                                ! et les frequences de tirage correspondantes
      if (Nappels_TH.GE.2) then
         Ntir_tot = 0 
         do i_TH=1,N_TH
            Ntir_tot=Ntir_tot+Ntir_TH(i_TH)
         end do
         if (Ntir_tot.GT.0) then
            UnSurNtot=1.0d0/Ntir_tot
         else
            UnSurNtot=1.0d0
         end if
         write(NomFichTestPr_TH,300) Nappels_TH
 300     format('evolamas_TEST_PR_TH_',I4.4,'.asc')
         open(11,file=NomFichTestPr_TH)
         write(11,310) Ntir_tot
 310     format(
     $        '# Test des probabilites de tirage',/,
     $        '# Nb de tirages : ',I6,/,
     $        '# 1: Probabilite 2: Frequence' )
         do i_TH=1,N_TH
            write(11,320) Pr_TH(i_TH),Ntir_TH(i_TH)*UnSurNtot
 320        format(2(E12.4,' '))
         end do
         close(11)
         do i_TH=1,Nmax_TH
            Ntir_TH(i_TH)=0
         end do
      end if
#endif
c
c---- Echelle de temps de la perte de masse
c
      if (Nappels_TH.LE.1) then
         Tps_TH_prec=Tps_Amas
         M_TH_prec=M_amas
      end if
      dt_evp=max(1d-20,Tps_Amas-Tps_TH_prec)
c$$$      dmrel_evp=2.0d0*abs(M_TH_Prec-M_amas)/(M_TH_Prec+M_amas)
      dmrel_evp=abs(M_TH_Prec-M_amas) / Mamas_ini ! Attempt Oct 12, 01
      if (dmrel_evp.GT.1.0d-10) then
         T_evp=dt_evp/dmrel_evp
      else
         T_evp=1.0d30
      end if
#if (_SURV_T_EVAP_>0)
      if (T_evp.LT.1.0d-6) then
         write(0,*) '!!! T_evp = ',T_evp,' !!!'
         write(0,*) '!!! ',Tps_TH_prec,Tps_Amas,
     $        M_TH_Prec,M_amas,Mamas_ini,' !!!'
         call Meurt('Calc_Param_Tirage_Henon')
      end if
#endif
      M_TH_prec=M_amas
      Tps_TH_prec=Tps_Amas
#if (_VERBOSITE_TH_>0)
      write(0,115) T_evp
 115  format(
     $     '}} Echelle tps de la perte de masse: T_evp = ',E12.4)
#endif
c
c---- calcul des temps caracteristiques
c
      call Calc_Tcarac(Trel,Tcoll,Tkick,Torb,Tadiab,Taccr,Tevst,Tstaccr)
c 
c---- on calcule les proba de tirage selon 1/Trelax + 1/Tcoll
c
      N_TH = min(Nmax_TH,NbSE-NbSE_perte)
      Pr_i = 0.0d0
      M_i = 0.0d0
      i_TH = 1
      Nsuiv=(NbSE-NbSE_perte-1)/N_TH+1
      do iRang = 1,NbSE-NbSE_perte
         dt=0.0d0
         c_tot=0.0d0
#if (_RELAXATION_>0)
         dt = dt + c_relx/(Frac_Trelax*Trel(iRang))
         c_tot = c_tot + c_relx
#endif
#if (_COLLISIONS_>0)
         dt = dt + c_coll/(Frac_Tcoll*Tcoll(iRang))
         c_tot = c_tot + c_coll
#endif
#if (_KICKS_>0)
         dt = dt + c_kick/(Frac_Tkick*Tkick(iRang))
         c_tot = c_tot + c_kick
#endif
         if (dt.GT.0.0d0) then
            dt = c_tot/dt
         else
            dt = 1.0d30
         end if

         dt = min(dt,Frac_Tevap*T_evp)        ! Pour resoudre la perte de masse (evaporation)

         dt = min(dt,Frac_Taccr*Taccr(iRang)) ! Pour resoudre la croissance du TN

#if (_EVOL_STELL_>0)
         dt = min(dt,Frac_Tevst*Tevst(iRang)) ! Pour resoudre l'evolution stellaire
#endif
#if (_ACCR_BY_STARS_>0)
         dt = min(dt,Frac_Tstaccr*Tstaccr(iRang)) ! To resolve mass accretion by stars (assumed adiabatic)
#endif
#if (_SURV_ADIAB_>0)
                                ! Pour s'assurer que l'on suit les modifications 
                                ! adiabatiques induites par l'evolution du pot
         dt = min(dt,Frac_Tadiab*Tadiab(iRang))
#endif
         dt_max(iRang) = dt     ! pas de tps maximal pour une particule de rang iRang
#if (_CALC_PR_==1) 
         Pr_i=max(Pr_i,1.0d0/dt)
#elif (_CALC_PR_==2)
         Pr_i=Pr_i + M_SE(iListeRay(iRang))/dt
         M_i = M_i + M_SE(iListeRay(iRang))
#else
         write(0,*) '!!! valeur de _CALC_PR_ incorrecte ',
     $        'dans EvolAmas_TirageHenon3.F !!!'
         call exit(1)
#endif
         if (iRang.GE.Nsuiv) then
#if (_CALC_PR_==1) 
            Pr_TH(i_TH)=Pr_i 
#elif (_CALC_PR_==2)
            Pr_TH(i_TH)=Pr_i/M_i
#endif
            Pr_i = 0.0d0
            i_TH=i_TH+1
            Nsuiv=(i_TH*(NbSE-NbSE_perte)-1)/N_TH+1
         end if
      end do
#if (_CALC_PR_==1) 
      if (i_TH.LE.N_TH) Pr_TH(N_TH) = Pr_i
#elif (_CALC_PR_==2)
      if (i_TH.LE.N_TH) Pr_TH(N_TH) = Pr_i/M_i
#endif
                                ! on fait en sorte que les probabilites soient
                                ! monotones decroissantes
      do i=N_TH-1,1,-1
         Pr_TH(i)=max(Pr_TH(i),Pr_TH(i+1))
      end do
                                ! on limite eventuellement la dynamique des proba
                                ! afin d'etre sur que toutes les SE vont evoluer
      if (Pr_TH(1).GT.FactTrelMax*Pr_TH(N_TH)) then
         Pinf=Pr_TH(1)/FactTrelMax
         do i=1,N_TH
            Pr_TH(i)=max(Pr_TH(i),Pinf)
         end do
      end if
                                ! normalisation des probabilites a 1
                                ! et calcul des probabilites integrees
      Pr_tot = 0.0d0
      do i=1,N_TH
         Pr_tot = Pr_tot+Pr_TH(i)
      end do
         
      Pinf = 0.0d0
      do i=1,N_TH
         Pr_TH(i)=Pr_TH(i)/Pr_tot
         Pr_int_TH(i)=Pinf+Pr_TH(i)
         Pinf=Pr_int_TH(i)
      end do
c
c---- on determine dt_TH de maniere a ce que le pas de tps soit tjs inferieur
c     au dt max determine ci-dessus
c
      dt_TH=rEnorme
      do iRang=1,NbSE-NbSE_perte
         dt_TH=min( dt_TH, dt_max(iRang)*P_TH(iRang) )
      end do

#if (_VERBOSITE_TH_>0)
      write(0,100) dt_TH,Pr_TH(1)/Pr_TH(N_TH)
 100  format(
     $     '}} Parametre pour tirage selon Henon : dt_TH = ',E12.4,/,
     $     '}}                 Rapport P_TH_max/P_TH_min = ',E12.4 )
#endif
      
#if (_ECRIT_TSCALE_>0)
      if (iPas_Evol.GT.0) then  ! sinon, le fichier n'a pas ete initialise 
#if (_RELAXATION_>0)
         Trel_01=selip( int(0.1d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Trel)
         Trel_05=selip( int(0.5d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Trel)
         Trel_09=selip( int(0.9d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Trel)
         Trel_c=Trel(1)
#else
         Trel_c=1.0d30
         Trel_01=1.0d30
         Trel_05=1.0d30
         Trel_09=1.0d30
#endif 
#if (_COLLISIONS_>0)
         Tcoll_01=selip( int(0.1d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Tcoll)
         Tcoll_05=selip( int(0.5d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Tcoll)
         Tcoll_09=selip( int(0.9d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Tcoll)
         Tcoll_c=Tcoll(1)
#else
         Tcoll_c=1.0d30
         Tcoll_01=1.0d30
         Tcoll_05=1.0d30
         Tcoll_09=1.0d30
#endif
#if (_KICKS_>0)
         Tkick_01=selip( int(0.1d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Tkick)
         Tkick_05=selip( int(0.5d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Tkick)
         Tkick_09=selip( int(0.9d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Tkick)
         Tkick_c=Tkick(1)
#else 
         Tkick_c=1.0d30
         Tkick_01=1.0d30
         Tkick_05=1.0d30
         Tkick_09=1.0d30
#endif
#if (_EVOL_STELL_>0)
         Tevst_01=selip( int(0.1d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Tevst)
         Tevst_05=selip( int(0.5d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Tevst)
         Tevst_09=selip( int(0.9d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Tevst)
#else
         Tevst_01=1.0d30
         Tevst_05=1.0d30
         Tevst_09=1.0d30
#endif
#if (_SURV_ADIAB_>0)
         Tadiab_01=selip( int(0.1d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Tadiab)
         Tadiab_05=selip( int(0.5d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Tadiab)
         Tadiab_09=selip( int(0.9d0*(NbSE-NbSE_perte)+0.5d0),
     $        NbSE-NbSE_perte, Tadiab)
#else
         Tadiab_01=1.0d30
         Tadiab_05=1.0d30
         Tadiab_09=1.0d30
#endif
         write(iFich_Tscale,510) iPas_Evol, Tps_Amas,
     $        Trel_c,  Trel_01,  Trel_05,  Trel_09,
     $        Tcoll_c, Tcoll_01, Tcoll_05, Tcoll_09,
     $        Tkick_c, Tkick_01, Tkick_05, Tkick_09,
     $        Tevst_01, Tevst_05, Tevst_09,
     $        Tadiab_01,Tadiab_05,Tadiab_09,
     $        Taccr_01, Taccr_05, Taccr_09,
     $        T_evp,
     $        dt_TH,Pr_TH(1)/Pr_TH(N_TH)
 510     format(I10,' ',E14.7,24(' ',E11.4))
         call flush_MD(iFich_Tscale)
      end if
#endif /* if (_ECRIT_TSCALE_>0) */

#if (_TEST_TH_==2)
      open(77,file='evolamas_TirageHenon3_TH.asc')
      write(77,220)
 220  format('# 1: Pr_TH 2: Pr_int_TH')
      do iRang=1,N_TH
         write(77,230) Pr_TH(iRang),Pr_int_TH(iRang)
 230     format(2(E12.4,' '))
      end do
      close(77)

      open(77,file='evolamas_TirageHenon3_dt.asc')
      write(77,200)
 200  format('# Pas de temps selon TirageHenon3')
      do iRang=1,NbSE-NbSE_perte
         write(77,210) PasTps_TH(iRang)
 210     format(E12.4)
      end do
      close(77)
      call Tester_TH
      stop
#endif
#if (_TEST_TH_==3)
      open(77,file='evolamas_TEST_Prob_TH.asc',form='formatted')
      write(77,120)
 120  format('# 1: Rang 2: iSE 3: R_SE 4: Trel 5: Tcoll',
     $     ' 6: Torb 7: P_TH 8: PasTps_TH 9: h_TH')
      do i=1,NbSE-NbSE_Perte
         write(77,125) i,iListeRay(i),R_SE(iListeRay(i)),Trel(i),
     $        Tcoll(i),Torb(i),P_TH(i),PasTps_TH(i),h_TH(i)
 125     format(2(I6,' '),7(E12.4,' '))
      end do
      stop
#endif
c
      if (dt_TH.LE.0.0d0) then
         write(0,*) '!!! dt_TH = ',dt_TH,' !!!'
         call Meurt('Calc_Param_Tirage_Henon')
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Ecrit_TH_bin(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_TirageHenon3_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer i
c......................................................................
c
      write(iFich) dt_TH,N_TH,(Pr_TH(i),i=1,N_TH)
c
      end
c
c----------------------------------------------------------------------
      subroutine Lit_TH_bin(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_TirageHenon3_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer i
c......................................................................
c
      read(iFich) dt_TH,N_TH,(Pr_TH(i),i=1,N_TH)
      Pr_int_TH(1)=Pr_TH(1)
      do i=2,N_TH
         Pr_int_TH(i)=Pr_int_TH(i-1)+Pr_TH(i)
      end do
c
      end
c
c----------------------------------------------------------------------
      double precision function h_TH(iRang)
c----------------------------------------------------------------------
c     renvoie la probabilite de tirer une SE de rang <= iRang
c      
      include 'EvolAmas_TirageHenon3_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iRang
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer i_TH
      double precision x_TH,Pinf
c......................................................................
c
      x_TH = dble(iRang*N_TH)/dble(NbSE-NbSE_perte)
      i_TH = ((iRang-1)*N_TH)/(NbSE-NbSE_perte)
      if (i_TH.EQ.0) then
         Pinf = 0.0d0
      else
         Pinf = Pr_int_TH(i_TH)
      end if
      h_TH =  Pinf + Pr_TH(i_TH+1)*(x_TH-i_TH)
c
      end
c
c----------------------------------------------------------------------
      integer function ihinv_TH(h)
c----------------------------------------------------------------------
c     Renvoie la fonction inverse de h_TH: rang de la SE _interieure_ 
c     de la paire
c      
      include 'EvolAmas_TirageHenon3_inc.f'
      include 'EvolAmas_Param_common.f'
c
c argument :
c ^^^^^^^^^^
      double precision h
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer i_TH
      double precision x

#if (_TEST_PR_TH_==1)
      integer Ntir_TH(Nmax_TH)
      common /common_test_Pr_TH/ Ntir_TH
#endif
c......................................................................
c
      if (h.LT.Pr_TH(1)) then
         x  = h/(N_TH*Pr_TH(1))
         i_TH = 0
      else
         call locate(Pr_int_TH,N_TH,h,i_TH)
         x  = ( i_TH +
     $        (h-Pr_int_TH(i_TH))/Pr_TH(i_TH+1)
     $        ) / N_TH ! x dans [0;1[
      end if

      if (x.GE.1.0d0) then
         write(0,*) '!!! x = ',x,' !!!'
         call Meurt('ihinv_TH')
      end if
      ihinv_TH = int((NbSE-NbSE_perte-DeltaRang_Paire)*x)+1 ! Modified on 11 june 2003 "-DeltaRang_Paire"
#if (_TEST_PR_TH_==1)
      Ntir_TH(i_TH+1) = Ntir_TH(i_TH+1)+1
#endif
c
      end
c
c----------------------------------------------------------------------
      double precision function P_TH(iRang)
c----------------------------------------------------------------------
c     renvoie la probabilite de tirer une SE de rang = iRang
c      
      include 'EvolAmas_TirageHenon3_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iRang
c
c var locales :
c ^^^^^^^^^^^^^
      integer i_TH
c......................................................................
c
      i_TH = min(N_TH,max(1,int((iRang-1)*N_TH/(NbSE-NbSE_perte))+1))
      P_TH = N_TH*Pr_TH(i_TH)/dble(NbSE-NbSE_perte)
c
      end
c

c----------------------------------------------------------------------
      double precision function PasTps_TH(iRang)
c----------------------------------------------------------------------
c     renvoie le pas de temps pour une SE de rang = iRang
c      
      include 'EvolAmas_TirageHenon3_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iRang
c
c fonction :
c ^^^^^^^^^^
      double precision P_TH
c......................................................................
c
      PasTps_TH = dt_TH/P_TH(iRang)
c
      end
c
c----------------------------------------------------------------------
      subroutine Tester_TH
c----------------------------------------------------------------------
c
      include 'EvolAmas_TirageHenon3_inc.f'
c
c constante(s) :
c ^^^^^^^^^^^^^^
      integer Ntest
      parameter (Ntest=100000)
c
c var locale(s) :
c ^^^^^^^^^^^^^^^
      integer iRang,i
      double precision h1,h2
c
c fonction(s) :
c ^^^^^^^^^^^^^
      double precision random, h_TH, P_TH
      integer ihinv_TH
c......................................................................
c
      do i=1,Ntest
         h1=random()
         iRang=ihinv_TH(h1)
         h2=h_Th(iRang)
         write(88,'(E14.6," ",I6," ",E14.6," ",E14.6))')
     $        h1,iRang,h2,P_TH(iRang)
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_FichTscale_ASC(iFich)
c----------------------------------------------------------------------
      implicit none
c
c argument:
c ^^^^^^^^^
      integer iFich
c......................................................................
c
      write(iFich,100)
 100  format(
     $     '# Echelles de temps pour evolution d''amas')
      call OutputFileHeader(iFich,
     $     'iPas_Evol Tps_Amas ' //
     $     'Trel_c Trel_01 Trel_05 Trel_09 ' //
     $     'Tcoll_c Tcoll_01 Tcoll_05 Tcoll_09 ' //
     $     'Tkick_c Tkick_01 Tkick_05 Tkick_09 ' //
     $     'Tevst_01 Tevst_05 Tevst_09 ' //
     $     'Tadiab_01 Tadiab_05 Tadiab_09 ' //
     $     'Taccr_01 Taccr_05 Taccr_09 ' //
     $     'T_evp dt_TH rap_PTH_maxmin')
c
      end
c
