c======================================================================
c    Blocs common des parametres de simulation pour EvolAmas.F
c    fichier genere automatiquement par EvolAmas_Param_rdb2common.sh
c======================================================================

      character*16
     &    Type_CtrObj_def       !  Default type of central object (TN or MS)

      common /Common_Param_char16/
     &    Type_CtrObj_def

      character*80
     &    NomSimul              !  Run name. Used to form the name of input/output files

      common /Common_Param_char80/
     &    NomSimul

      double precision
     &    FacNbPasSauv,         !  A complete save (snapshot and all) is done every FacNbPasSauv*NbSE steps (NbSE is the number of particles)
     &    FacNbPasSauvPart,     !  Same as FacNbPasSauv for partial save
     &    FacNbPasInfo,         !  General info about run will be given every FacNbPasInfo*NbSE steps
     &    FacNbPasRecArbre,     !  Binary tree rebuilt from scratch every FacNbPasRecArbre*NbSE steps
     &    FacNbPasRecRelax,     !  Relaxation parameters recomputed every FacNbPasRecRelax*NbSE steps
     &    FacNbPasRecTpsAmas,   !  Cluster (median) time recomputed every FacNbPasRecTpsAmas*NbSE steps
     &    FacNbPasDetTrelExtremes, !  Extreme values of relax times recomputed every FacNbPasDetTrelExtremes*NbSE steps
     &    FacNbPasDetEvap,      !  Check for evaporated stars performed every FacNbPasDetEvap*NbSE steps
     &    FacNbPasTestDuplic,   !  Test for particle duplication performed every FacNbPasTestDuplic*NbSE steps
     &    FacNbPasTestDemSauv,  !  Test for file asking for snapshot (_SAUV_DEMANDEE_) performed every FacNbPasTestDemSauv*NbSE steps
     &    T_fin,                !  Code termintes when time (in some units!?) reaches this value
     &    Frac_Trelax,          !  Requested maximum value for time_step/relaxation_time
     &    Frac_Tcoll,           !  Requested maximum value for time_step/collision_time
     &    Frac_Tkick,           !  Requested maximum value for time_step/large_angle_scattering_time
     &    Frac_Tadiab,          !  Requested maximum value for time_step/potential_change_time
     &    Frac_Taccr,           !  Requested maximum value for time_step/central_object_growth_time
     &    Frac_Tevap,           !  Requested maximum value for time_step/evaporation_time
     &    Frac_Tevst,           !  Requested maximum value for time_step/stellar_evolution_time
     &    Frac_Tstaccr,         !  Requested maximum value for time_step/accretion_by_star_time
     &    Tadiab_ini,           !  Initial value of time scale for potential evolutioin (just needs to be short)
     &    FactTrelMax,          !  Maximum ratio of time steps (misnomer!)
     &    FacMinPG,             !  Grid rebuilt if number of particle in any cell drops below FacMinPG*NbSECouchePot
     &    FacMaxPG,             !  Grid rebuilt if number of particle in any cell raises above FacMaxPG*NbSECouchePot
     &    MasseEtoileDef,       !  Default value for average stellar mass (overriden by info in "ETOILES.xdr" file if present)
     &    rNbEtoiles,           !  Total number of stars (independent of NbSE!)
     &    Frac_dt_max,          !  Obsolete, was used to determine value of sub-timestep
     &    TailleAmas_en_pc,     !  N-body length units in parsecs
     &    FacNbPasRecalcTH,     !  The parameters for the function used to select particle pairs are recomputed every FacNbPasRecalcTH*NbSE steps
     &    FacNbPasRecCB3c,      !  Not used. 3-body binary heating recomputed every FacNbPasRecCB3c*NbSE steps
     &    M_TN_ini_def,         !  Default initial mass of central object, in fraction of total stellar mass (overriden by info in "TN.xdr" file if present)
     &    R_Mar_ini_def,        !  Default initial value of tidal truncation radius in N-body units (overriden by info in "MAREE.xdr" file if present)
     &    Gamma_relax,          !  Proportionalite coefficient in Coulomb log. Lambda=Gamma_relax*N_star
     &    frac_accr_Coll,       !  Fraction of mass liberated in collisions accreted by MBH
     &    frac_accr_Dechir,     !  Fraction of mass liberated in tidal disruptions accreted by MBH
     &    frac_accr_EvolSt,     !  Fraction of mass liberated by stellar evolution accreted by MBH
     &    effic_conv_lum,       !  Efficiency for converting mass into light during accretion onto MBH
     &    Fac_dt_accr_M,        !  The amount of mass in the central reservoir (disk) that can be accreted onto the MBH is determined every Fac_dt_accr_M*NbSE steps
     &    Fac_dt_accr_Npart,    !  The amount of mass in the central reservoir (disk) that can be accreted onto the MBH is determined when the reservoir has accumulated more than Fac_dt_accr_Npart times the mass of an average particle
     &    Mu_mol_elec,          !  Average molecular mass per electron of stellar gas
     &    fact_b0_kick,         !  Encounters with impact parameter smaller than fact_b0_kick times the value for 90 degree deflection are treated as large-angle scatterings
     &    FracPas_AjustPassif,  !  Fraction of passive steps (change of orbital position only)
     &    FactR_CC,             !  Factor for distance to centre used in "central control". Probably obsolete
     &    StellMet_def,         !  Default stellar metallicity (0.02 is solar)
     &    Fac_SldAvg_Tscale,    !  NbSECouchePot*Fac_SldAvg_Tscale is the number of particles for sliding average to compute time scales in Calc_Tcarac
     &    CollMassLossFracDef,  !  Flat fractional mass loss in simplistic collisional treatment
     &    CollDminDestr_MS_vs_RMN, !  Distance below which the MS star is completely destroyed in encounters with compact stars (for _TRAITEMENT_COLL_MS_RMN_=1)
     &    CollAccrFrac_MS_on_WD, !  Fraction of mass of the MS star accreted onto WD if MS star tidally disrupted (for _TRAITEMENT_COLL_MS_RMN_=2)
     &    CollAccrFrac_MS_on_NS, !  Fraction of mass of the MS star accreted onto NS if MS star tidally disrupted (for _TRAITEMENT_COLL_MS_RMN_=2)
     &    CollAccrFrac_MS_on_BH, !  Fraction of mass of the MS star accreted onto BH if MS star tidally disrupted (for _TRAITEMENT_COLL_MS_RMN_=2)
     &    Coef_Trlx_GW_Capt,    !  Fudge factor on trelax when compared to GW-inspiral time for GW-capture
     &    Tmax_GW_Capt_Gyr,     !  Maximum inspiral time considered for GW captures, in Gyr
     &    Trlx_GW_def_Gyr,      !  Relaxation time when a constant value is used to look for captures (_TYPE_TREL_GW_=2; not recommended)
     &    R_tid_disr_def_NB,    !  Tidal disruption radius all stars in N-body units when _CONST_R_TID_DISRUP_>0
     &    RateAccrPerStar_MsunMyr, !  Average rate of accretion by stars in Msun/Myr (if _ACCR_BY_STARS_>0); if negative, will be used as a maximum
     &    tAccrByStar_1Msun_Myr, !  Parameter for accretion by stars (if _ACCR_BY_STARS_>0); overriden if RateAccrPerStar_MsunMyr>0
     &    exp_AccrByStar,       !  Parameter for accretion by stars
     &    Mmax_AccrByStar       !  Parameter for accretion by stars

      common /Common_Param_dp/
     &    FacNbPasSauv, FacNbPasSauvPart, FacNbPasInfo, FacNbPasRecArbre, 
     &    FacNbPasRecRelax, FacNbPasRecTpsAmas, FacNbPasDetTrelExtremes, 
     &    FacNbPasDetEvap, FacNbPasTestDuplic, FacNbPasTestDemSauv, T_fin, 
     &    Frac_Trelax, Frac_Tcoll, Frac_Tkick, Frac_Tadiab, Frac_Taccr, 
     &    Frac_Tevap, Frac_Tevst, Frac_Tstaccr, Tadiab_ini, FactTrelMax, 
     &    FacMinPG, FacMaxPG, MasseEtoileDef, rNbEtoiles, Frac_dt_max, 
     &    TailleAmas_en_pc, FacNbPasRecalcTH, FacNbPasRecCB3c, 
     &    M_TN_ini_def, R_Mar_ini_def, Gamma_relax, frac_accr_Coll, 
     &    frac_accr_Dechir, frac_accr_EvolSt, effic_conv_lum, 
     &    Fac_dt_accr_M, Fac_dt_accr_Npart, Mu_mol_elec, fact_b0_kick, 
     &    FracPas_AjustPassif, FactR_CC, StellMet_def, Fac_SldAvg_Tscale, 
     &    CollMassLossFracDef, CollDminDestr_MS_vs_RMN, 
     &    CollAccrFrac_MS_on_WD, CollAccrFrac_MS_on_NS, 
     &    CollAccrFrac_MS_on_BH, Coef_Trlx_GW_Capt, Tmax_GW_Capt_Gyr, 
     &    Trlx_GW_def_Gyr, R_tid_disr_def_NB, RateAccrPerStar_MsunMyr, 
     &    tAccrByStar_1Msun_Myr, exp_AccrByStar, Mmax_AccrByStar

      integer
     &    iVerbeux,             !  Determines how verbose the code will be (on the standard error)
     &    NbSECouchePot,        !  Requested number of particles per cell in density radial grid (misnomer!)
     &    TypeEtoileDef,        !  Default stellar Type (overriden by info in "ETOILES.xdr" file if present; 1 for MS stars)
     &    iRand_Seed,           !  Seed for random number generator
     &    NbSECoucheCB3c,       !  Not used. Number of grid elements used to determine 3-body binary heating. Never worked well.
     &    NbSauv_Conserve,      !  Only 1 snapshot every NbSauv_Conserve is kept. Others are erased as the simulation proceeds to save disk space
     &    Rang_CC,              !  Number of particles subject to special "central control". Probably obsolete
     &    DeltaRang_Paire       !  Rank difference for two particles in an interacting pair

      common /Common_Param_int/
     &    iVerbeux, NbSECouchePot, TypeEtoileDef, iRand_Seed, 
     &    NbSECoucheCB3c, NbSauv_Conserve, Rang_CC, DeltaRang_Paire

      integer*8
     &    i_ini,                !  Initial step value (Code will attempt to read corresponding snapshot)
     &    i_fin                 !  Code terminates when step counter reaches this value

      common /Common_Param_int8/
     &    i_ini, i_fin
