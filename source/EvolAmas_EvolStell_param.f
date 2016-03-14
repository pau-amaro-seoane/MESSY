      double precision Mzams_max_WD_def, Mzams_max_NS_def,
     $     Mzams_max_fallbackBH_def,
     $     M_max_WD_def, M_max_NS_def, Mremn_min_BH_def,
     $     M_WD_def, M_NS_def, M_BH_def,
     $     NatKick_Sig1D_1_def, NatKick_Sig1D_2_def, NatKick_Prob_1_def
      parameter (
     $     Mzams_max_WD_def = 8.0d0,  ! Masse limite (ZAMS) conduisant a la formation d'une WD
     $     Mzams_max_NS_def = 30.0d0, ! Masse limite (ZAMS) conduisant a la formation d'une NS
     $     Mzams_max_fallbackBH_def = 42.0d0, ! Masse limite (ZAMS) conduisant a la formation d'un BH par fall-back
     $     M_max_WD_def = 1.44d0,     ! Masse maximal d'une WD (Chandrasekhar mass)
     $     M_max_NS_def = 3.0d0,      ! Masse maximal d'une NS
     $     Mremn_min_BH_def = 3.0d0,  ! Masse minimale d'un TN (Fryer & Kalogera 01) 
     $     M_WD_def = 0.6d0,      ! Masse naine blanche (approx constante)
     $     M_NS_def = 1.4d0,      ! Masse etoile de neutrons (approx constante)
     $     M_BH_def = 7.0d0)      ! Masse trou noir (approx constante)

c
c---- Default parameters for natal kicks (single Maxwellian from Hobbs et al. 2005, MNRAS 360, 973)
c
      parameter (
     $     NatKick_Sig1D_1_def = 265.0d0, ! 1D dispersion in km/s
     $     NatKick_Sig1D_2_def = 0.0d0,
     $     NatKick_Prob_1_def = 1.0d0
     $     )
