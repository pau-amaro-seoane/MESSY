c======================================================================
c     parametres decrivant la fct des masses 
c======================================================================
c
c constantes :
c ^^^^^^^^^^^^
      integer Nmax_morc
      parameter (Nmax_morc=20)
      integer N_morc
      double precision
     $     Mlim(0:Nmax_morc),   ! Les limites en M_sol des "morceaux"
     $     alpha(Nmax_morc)     ! Les exposants pour chaque morceau

      double precision Age_amas, Tsf
      common /Param_FM_common/
     $     Age_amas,Tsf,
     $     Mlim,alpha,N_morc


      logical lRenorm
      common /Control_FM_common/ lRenorm
