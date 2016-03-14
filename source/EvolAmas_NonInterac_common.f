c======================================================================
c     common blocks containing information about which stellar type
c     interact with which one (for relaxation)
c======================================================================

      logical lRelaxesWith(Nb_TypesEtoile,Nb_TypesEtoile),
     $     lNonInterac, lNonInteracAsymmetrical
      common /NonInterac_common/
     $     lRelaxesWith, lNonInterac, lNonInteracAsymmetrical
c
c     lRelaxesWith(Type1,Type2) = .TRUE. means that Type1 is affected by relaxation with Type2
c     By default, all elements are .TRUE.
c
c     lNonInterac = .TRUE. means that there are non-interacting combinations (.FALSE. by default)
c
c     lNonInteracAsymmetrical = .TRUE. means that there exist Type1,Type2 with
c     lRelaxesWith(Type1,Type2) = .TRUE. but lRelaxesWith(Type2,Type1) = .FALSE.
c     Therefore, energy cannot be strictly conserved. .FALSE. by default
c
