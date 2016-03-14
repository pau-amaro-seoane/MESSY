c======================================================================
c     commons pour le suivi des differents types stellaires
c======================================================================
c
c constante :
c ^^^^^^^^^^^
      integer N_RayLag_STS_Max
      parameter (N_RayLag_STS_Max=100)
c
c common :
c ^^^^^^^^
      integer N_RayLag_STS
      double precision
     $     f_RayLag_STS(N_RayLag_STS_Max),
     $                          ! Masse dans l'amas stellaire
     $     M_STS(Nb_TypesEtoile),
     $                          ! Masse accretee au centre
     $     dMaccr_evst_STS(Nb_TypesEtoile),
     $     dMaccr_dechir_STS(Nb_TypesEtoile),
     $     dMaccr_disp_STS(Nb_TypesEtoile),
     $     dMaccr_captGW_STS(Nb_TypesEtoile),
     $     dMaccr_coll_STS(Nb_TypesEtoile),
     $                          ! Masse ejectee de l'amas
     $     dMejec_evst_STS(Nb_TypesEtoile),
     $     dMejec_dechir_STS(Nb_TypesEtoile),
     $     dMejec_coll_STS(Nb_TypesEtoile),
     $     dMejec_evap_STS(Nb_TypesEtoile)
      integer 
     $     N_STS(Nb_TypesEtoile)
      common /STS_Common/
     $     f_RayLag_STS,M_STS,
     $     dMaccr_evst_STS,dMaccr_dechir_STS,
     $     dMaccr_disp_STS,dMaccr_captGW_STS,
     $     dMaccr_coll_STS,
     $     dMejec_evst_STS,dMejec_dechir_STS,
     $     dMejec_coll_STS,dMejec_evap_STS,
     $     N_RayLag_STS,N_STS
c
c tableaux "constants" (initialises dans block data)
c
      character*32 TypeFich_STS(Nb_TypesEtoile)
      integer iFich_STS(Nb_TypesEtoile)
      
      common /STS_common_csts/ iFich_STS, TypeFich_STS
