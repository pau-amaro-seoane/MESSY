      implicit none

      integer N_MostMassivePart
      parameter (N_MostMassivePart=10) ! the number of most massive particles to be tracked

      integer N_MassSpec_Max
      parameter (N_MassSpec_Max=50)
      integer N_MassSpec
      double precision M_MassSpec(N_MassSpec_Max) ! array to contain the fractional masses for the cumulative mass spectrum

      common /common_MassSpec/ M_MassSpec, N_MassSpec
