c======================================================================
c     Types d'etoiles
c======================================================================

      integer*1 TypeEtoile_Unspec,
     $     TypeEtoile_MS,TypeEtoile_RG,TypeEtoile_WD,
     $     TypeEtoile_NS,TypeEtoile_BH
      parameter (
     $     TypeEtoile_Unspec=0,
     $     TypeEtoile_MS=1,
     $     TypeEtoile_RG=2,
     $     TypeEtoile_WD=3,
     $     TypeEtoile_NS=4,
     $     TypeEtoile_BH=5 )
      integer Nb_TypesEtoile
      parameter (Nb_TypesEtoile=5)

      double precision solar_metal, gc_metal
      parameter (solar_metal=0.02d0, gc_metal=1.0d0-4)
