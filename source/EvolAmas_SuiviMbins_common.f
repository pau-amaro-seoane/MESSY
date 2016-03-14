c======================================================================
c     commons pour EvolAmas_SuivMbins.F
c======================================================================
c
c constante :
c ^^^^^^^^^^^
      integer N_MbinsMax
      parameter (N_MbinsMax=50)
      integer N_MFracR_MbinMax
      parameter (N_MFracR_MbinMax=20)
c
c common :
c ^^^^^^^^
      integer N_Mbins, N_MFracR_Mbin
      integer iMode_Mbins ! Indicate how Inf_Mbin and Sup_Mbin have to be interpreted
                          ! iMode_Mbins=1 : inf,sup limits of each bin in M_sun
                          !             2 : Minf-Msup parametrization 0 for (current) Minf of the MF, 1 for Msup of the MF
                          !             3 : limits given in fraction of the total number of stars less massive than the mass
                          !             4 : limits given in fraction of the total mass   of stars less massive than the mass
      double precision Inf_Mbin(N_MbinsMax), Sup_Mbin(N_MbinsMax)
      double precision MFracR_Mbin(N_MFracR_MbinMax)
      common /Mbins_Common/ Inf_Mbin, Sup_Mbin, MFracR_Mbin,
     $     N_Mbins, N_MFracR_Mbin, iMode_Mbins
