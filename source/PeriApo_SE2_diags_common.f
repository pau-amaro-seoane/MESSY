c======================================================================
c     Common block for Peri/Apo routines to exchange diagnostic
c     information with outside world
c======================================================================
      double precision
     $     Mse_err_PeriApo, Ese_err_PeriApo, Jse_err_PeriApo,
     $     Rloc_err_PeriApo, Qloc_err_PeriApo, dQdRloc_err_PeriApo,
     $     A_err_PeriApo, B_err_PeriApo
      integer iErr_PeriApo
      logical lStop_err_PeriApo

      common /PeriApo_SE2_diags_common/ 
     $     Mse_err_PeriApo, Ese_err_PeriApo, Jse_err_PeriApo,
     $     Rloc_err_PeriApo, Qloc_err_PeriApo, dQdRloc_err_PeriApo,
     $     A_err_PeriApo, B_err_PeriApo,
     $     iErr_PeriApo, lStop_err_PeriApo
