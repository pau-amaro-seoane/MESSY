c======================================================================
c     Physical and astronomical constants in cgs units
c     Source: Binney & Tremaine 87, Appendix 1.A
c======================================================================

      double precision Grav_G, Light_c, Planck_h, Boltzmann_k,
     $     Electron_q, Electron_q_inCoulomb, Proton_m, Electron_m,
     $     StefanBoltz_Sigma, Thomson_CrossSec
      parameter (
     $     Grav_G               = 6.672d-8,
     $     Light_c	        = 2.99792458d10,
     $     Planck_h             = 6.62618d-27,
     $     Boltzmann_k          = 1.38066d-16,
     $     Electron_q           = 4.80324d-10,
     $     Electron_q_inCoulomb = 1.602189d-19,
     $     Proton_m             = 1.672649d-24,
     $     Electron_m           = 9.10953d-28,
     $     StefanBoltz_Sigma    = 5.6703d-5,
     $     Thomson_CrossSec     = 6.65245d-25
     $     )
      double precision AstroUnit, Parsec, SideralYear, Sun_M, Sun_GM,
     $     Sun_R, Sun_L, Sun_Vesc, Sun_MagV, Sun_MagB, Earth_M,
     $     Sun_R_inPc
      parameter (
     $     AstroUnit   = 1.49597892d13,
     $     Parsec      = 3.08567802d18,
     $     SideralYear = 3.1558149984d7,
     $     Sun_M       = 1.989d33,
     $     Sun_GM      = 1.32712497d26,
     $     Sun_R       = 6.9599d10,
     $     Sun_R_inPc  = Sun_R/Parsec,
     $     Sun_L       = 3.826d33,
     $     Sun_Vesc    = 6.175d7,
     $     Sun_MagV    = 4.83d0,
     $     Sun_MagB    = 5.48d0,
     $     Earth_M     = 5.976d27
     $     )
