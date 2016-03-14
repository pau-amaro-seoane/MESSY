c======================================================================
c     grandeurs liees a la limitation par champ de maree pour evolution
c     d'amas
c======================================================================

      double precision
     $     R_Mar,               ! rayon de maree
     $     C_Mar,               ! defint le rayon de maree par R_mar=C_Mar*M**1/3
                                ! avec M = M_amas+M_ctr-Mstel_TN - M_evap
                                !          \-------\/---------/
                                !              constant au cours de l'evolution
     $     E_mar
      
      common /common_Maree/ R_Mar,C_Mar,E_mar


      character*(*) EnTeteFichMaree_XDR
      parameter (EnTeteFichMaree_XDR='%%% XDR Cluster Tide File %%%')

