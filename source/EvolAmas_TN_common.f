c=========================================================================
c     Quantities linked to central object (MBH or VMS: very massive star)
c=========================================================================
      
      integer id_CtrObj
      parameter (id_CtrObj=0)

      double precision
     $     M_TN,                ! Mass of central object (in code units)
     $     M_ctr,               ! Total central mass: central VMS or MBH plus gas resevoir, i.e. accretion disk (M_rsrv=M_ctr-M_TN)
     $     BirthDate_CtrObj     ! Effective date of birth of central object, in code units (important if it is not a MBH but a VMS)
      integer*1 iType_CtrObj
      
      common /common_TN/
     $     M_TN, M_ctr, BirthDate_CtrObj, iType_CtrObj
c
c---- parametres pour la croissance forcee (adiabatique) du TN
c
      integer Nmax_expl_TN_CF
      parameter (Nmax_expl_TN_CF=1000) ! Nb max de points dans la relation explicite T->M_TN
      integer EXPL, PARAM
      parameter (PARAM=1, EXPL=2)
      integer N_expl_TN_CF, iType_TN_CF
      double precision
     $     T_TN_CF,M_TN_fin_CF,
     $     T_expl_TN_CF(Nmax_expl_TN_CF), M_expl_TN_CF(Nmax_expl_TN_CF),
     $     Y2_expl_TN_CF(Nmax_expl_TN_CF),
     $     T_TN_CF_ext, M_TN_fin_CF_ext
      common /common_TN_CF/
     $     T_TN_CF,M_TN_fin_CF,
     $     T_expl_TN_CF, M_expl_TN_CF, Y2_expl_TN_CF,
     $     T_TN_CF_ext, M_TN_fin_CF_ext,
     $     N_expl_TN_CF, iType_TN_CF

      character*(*) EnTeteFichTN_XDR
      parameter (EnTeteFichTN_XDR=
     $     '%%% XDR Cluster Central BH File %%%')

