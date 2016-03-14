c======================================================================
c     includes pour Creer_Amas_FdeE.F
c======================================================================

      implicit none

      include 'CstesMath.f'
      include 'Param_Amas.f'
      include 'VarSE.f'
      include 'PotSE.f'
c
c---- fichiers :
c
      integer iFich,iErr
      parameter (iFich=11, iErr=0)
      
      character*(80) NomOut
      common /NomsFich/ NomOut
c
c---- integration de RK
c
      double precision PrecRK
      parameter (PrecRK = 3.0e-10) ! precision de l integration
c
c---- tableaux de profils de l'amas
c
      integer DimTab
      parameter (DimTab=100000)
      
      double precision
     $     tabR(0:DimTab),   tabM(0:DimTab),
     $     tabRho(0:DimTab), tabPsi(0:DimTab)
      double precision Rho0
      integer
     $     NbElemTab
      logical
     $     lRayonFini
      common /TabProf/
     $     tabR, tabM, tabRho, tabPsi, Rho0, NbElemTab,
     $     lRayonFini
c
c---- parametres complementaires de l'amas
c
      double precision Psi0, Dist_typ, Dist_min, Rayon_Amas
      common /ParamComp/ Dist_typ, Dist_min, Rayon_Amas
c
c---- fct de distribution
c
      integer NbParamMax
      parameter (NbParamMax=4)
      integer NbParamFD
      double precision ParamFD(NbParamMax), ParamAux1
      common /FctDist/ ParamFD, ParamAux1, NbParamFD
      character*40 NomFD
      common /Nom_FctDist/ NomFD 
c
c--- Energie(s)
c
      double precision tol_Vir
      parameter (tol_Vir=0.03)  ! tolerance sur le viriel (pour NbSE=10000)
      double precision
     $     Omega,Viriel         ! energie pot totale, viriel
      common /Ener_tot/
     $     Omega,Viriel
c
c---- Facteur de normalisation pour passer aux unites N-corps
c
      double precision renR, renT, renJ
      common /normalisation/ renR, renT, renJ

c
c---- Graine pour le generateur aleatoire
c
      integer iRand_Seed
      common /seed_common/ iRand_Seed

