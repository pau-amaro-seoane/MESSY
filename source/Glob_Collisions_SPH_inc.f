c
c---- parametres globaux pour les collisions selon simulations SPH
c
      integer Nvar,Nfct
      parameter (
     $     Nvar=4,              ! Nb de variables independantes
     $     Nfct=4               ! Nb de "fcts": grandeurs decrivant le resultat d'une collision
     $     )
c
c---- Les noms des grandeurs
c
      character*64 NomsVar_ResColl(Nvar), NomsFct_ResColl(Nfct)
      common /ResColl_Noms/ NomsVar_ResColl, NomsFct_ResColl
