c
c constante(s) :
c ^^^^^^^^^^^^^^
      double precision C4_1D,C4_3D
      parameter (C4_1D=90.0d0, C4_3D=4.21d3) ! coefficient de chauffage par les binaires
c
c blocs common :
c ^^^^^^^^^^^^^^
      double precision Coef_CB3c
      common /common_CB3c/ Coef_CB3c

      character*(*) EnTeteFichCB3c_XDR
      parameter (EnTeteFichCB3c_XDR=
     $     '%%% XDR Cluster binary heating File %%%')
