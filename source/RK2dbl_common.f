
      double precision Prec__RK, dt_done_RK ! dt_done_RK indicates the value of last time step 
      integer iRej__RK          ! iRej__RK compte le nombre de pas 
                                ! d integration sans rejet    
      common /common__RK/ Prec__RK, dt_done_RK, iRej__RK
