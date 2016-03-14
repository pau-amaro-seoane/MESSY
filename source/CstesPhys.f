c======================================================================
c     Constantes physiques, conversion d'unites
c======================================================================
c     conversions d'unites
c
      double precision pc_en_cm, Rsol_en_cm, Rsol_en_pc, pc_en_Rsol,
     $     AU_en_cm, AU_en_pc, Msol_en_g, yr_en_sec
      parameter (
     $     pc_en_cm = 3.085678d18,
     $     Rsol_en_cm = 6.9599d10,
     $     Rsol_en_pc = Rsol_en_cm/pc_en_cm,
     $     pc_en_Rsol = pc_en_cm/Rsol_en_cm,
     $     AU_en_cm   = 1.49597892d13,
     $     AU_en_pc   = AU_en_cm/pc_en_cm,
     $     Msol_en_g  = 1.989d33,
     $     yr_en_sec  = 3.1558149984d7
     $     )
c
c     constantes physiques
c
      double precision c_en_cgs, c_en_Vesc_sol, mp_en_cgs,
     $     sigma_th_en_cgs, G_en_cgs
      parameter (
     $     c_en_cgs = 2.99792458d10,      ! c : vitesse de la lumiere
     $     c_en_Vesc_sol = 486.d0,        ! Vesc_sol=sqrt(2*G*Msol/Rsol)=617.4 km/s
     $     mp_en_cgs = 1.672645d-24,      ! masse du proton en grammes
     $     sigma_th_en_cgs = 6.65245d-25, ! section efficace de Thompson, en cm^2
     $     G_en_cgs =  6.672d-8           ! Constante de la gravitation
     $     )
