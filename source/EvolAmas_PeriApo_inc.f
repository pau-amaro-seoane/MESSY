c======================================================================
c     blocs common des informations complementaires utiles concernant
c     les position peri/apocentriques d'une SE
c======================================================================

      double precision A_apo, A_peri, B_apo, B_peri
      integer iRang_peri,iRang_apo
      common /Potperiapo/ A_apo, A_peri, B_apo, B_peri,
     $     iRang_peri,iRang_apo
      
