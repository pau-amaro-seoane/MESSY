c======================================================================
c     constantes et blocs common utilises pour "suivre" les SE
c     perdues par evaporation ou destruction
c======================================================================

      double precision Tps_perte,R_perte,
     $     R_evap,R_dest_coll,R_dechir,R_disp,R_captGW,
     $     R_perte_inconnu,
     $     R_vide
      parameter (
     $     Tps_perte=1.0d30,
     $     R_perte=1.0d30,
     $     R_evap=1*R_perte,R_dest_coll=2*R_perte,R_dechir=3*R_perte,
     $     R_disp=4*R_perte, R_captGW=5*R_perte,
     $     R_perte_inconnu=9.9d0*R_perte,
     $     R_vide=1.0d40        ! For particles that should not be counted at all
     $                          ! even to compute initial mass and such global quantities
     $     )

      integer
     $     NbSE_perte,NbSE_evap,NbSE_dest_coll,NbSE_dest_dechir,
     $     NbSE_dest_disp, NbSE_dest_captGW
      common /SE_Perdues/
     $     NbSE_perte,NbSE_evap,NbSE_dest_coll,NbSE_dest_dechir,
     $     NbSE_dest_disp, NbSE_dest_captGW
