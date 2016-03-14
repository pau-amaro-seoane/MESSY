c======================================================================
c     Variables for simple-minded stellar evolution
c     See EvolAmas_EvolStell_param.f for default values
c======================================================================

      double precision
     $     Mzams_max_WD, Mzams_max_NS, Mzams_max_fallbackBH,
     $     M_max_WD, M_max_NS, Mremn_min_BH,
     $     M_WD, M_NS, M_BH,
     $     NatKick_Sig1D_1, NatKick_Sig1D_2,
     $     NatKick_Prob_1

      common /EvolStellCommon/
     $     Mzams_max_WD, Mzams_max_NS, Mzams_max_fallbackBH,
     $     M_max_WD, M_max_NS, Mremn_min_BH,
     $     M_WD, M_NS, M_BH,
     $     NatKick_Sig1D_1, NatKick_Sig1D_2,
     $     NatKick_Prob_1

      logical lInitEvolStell
      common /EvolStellControl/ lInitEvolStell

      double precision NatalKickLastEvol 
      common /EvolStell_NatalKick_Common/
     $     NatalKickLastEvol

      character*64 StringTypeEvolStell
      common /EvolStell_Type_Common/
     $     StringTypeEvolStell
