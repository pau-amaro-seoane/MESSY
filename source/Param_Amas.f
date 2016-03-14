c
c---- Nombre Max de Super Etoiles dans la simulation d'amas
c
      integer iDimSE 
      parameter (iDimSE =  2 000 000)
c
c---- Nb max de noeuds dans l'arbre
c
      integer iDimArbre
      parameter (iDimArbre = iDimSE) 
c
c---- Fraction de masse 'auto-gravitante' d'une SE
c
      double precision FracMasseInt, FracMasseInt1 
      parameter (FracMasseInt=0.0d0, FracMasseInt1=FracMasseInt-1.0d0)
c
c---- grandeurs de "garage" pour les SE evaporees
c
      integer iPG_evap
      parameter (iPG_evap=1000000)
      double precision Revap,Tps_evap,dT_evap
      parameter (Revap = 1.0d30, Tps_evap=1.0d30, dT_evap=1.0d30)
c
c---- grandeurs de "garage" pour les SE tombees dans le TN
c
      integer iPG_TN
      parameter (iPG_TN=iPG_evap)
      double precision R_TN,Tps_TN,dT_TN
      parameter (R_TN=2.0d0*Revap, Tps_TN=Tps_evap, dT_TN=dT_evap)
