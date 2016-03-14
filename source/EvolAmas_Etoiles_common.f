c======================================================================
c     grandeurs "stellaires" des SE pour amas a la Henon
c======================================================================
c
      double precision Met_SE(iDimSE), DNet_SE(iDimSE)
      integer*1 iTet_SE(iDimSE)

      common /common_Etoiles/
     $     Met_SE,              ! Masse des etoiles dans chaque SE (en unites de masse solaire)
     $     DNet_SE,             ! "Date" de naissance des etoiles de la SE (unites de tps du code)
     $     iTet_SE              ! Type des etoiles de la SE
c  
      character*(*) EnTeteFichEtoiles_XDR
      parameter (EnTeteFichEtoiles_XDR='%%% XDR Cluster Stars File %%%')
