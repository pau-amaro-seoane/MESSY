c======================================================================
c     Parameters for Peri/Apo routines to exchange diagnostic
c     information with outside world
c======================================================================

      integer
     $     iOK_err_PeriApo,
     $     iRsupGTPeri_err_PeriApo,
     $     iRinfLTPeri_err_PeriApo,
     $     iRinfGTApo_err_PeriApo,
     $     iRsupLTApo_err_PeriApo,
     $     iRinfNotFound_err_PeriApo,
     $     iRsupNotFound_err_PeriApo,
     $     iRsupMaxBissec_err_PeriApo,
     $     iNegDelta_err_PeriApo
      parameter (
     $     iOK_err_PeriApo            = 0,
     $     iRsupGTPeri_err_PeriApo    = 1,
     $     iRinfLTPeri_err_PeriApo    = 2,
     $     iRinfGTApo_err_PeriApo     = 3,
     $     iRsupLTApo_err_PeriApo     = 4,
     $     iRinfNotFound_err_PeriApo  = 5,
     $     iRsupNotFound_err_PeriApo  = 6,
     $     iRsupMaxBissec_err_PeriApo = 7,
     $     iNegDelta_err_PeriApo      = 8
     $     )
