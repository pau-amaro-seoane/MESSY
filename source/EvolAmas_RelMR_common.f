c======================================================================
c     commons pour EvolAmas_RelMR.F
c======================================================================

      double precision tabM_MS(NtabMR_Max),tabR_MS(NtabMR_Max),
     $     y2_MR(NtabMR_Max),a1_MR,a2_MR,UnSurRsol
      integer NTabMR
      common /Relation_MR_common/ tabM_MS,tabR_MS,y2_MR,a1_MR,a2_MR,
     $     UnSurRsol, NTabMR
      double precision lM_RhR_MS(N_RhR_max), RhR_MS(N_RhR_max)
      integer NRhR
      common /Rel_MRhalf_common/ lM_RhR_MS, RhR_MS, NRhR
