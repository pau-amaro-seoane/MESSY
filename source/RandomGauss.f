
c----------------------------------------------------------------------
      double precision function RandomGauss()
c----------------------------------------------------------------------
c     
      implicit none
c     
      integer iset
      double precision fac,gset,rsq,v1,v2
      double precision random
c     
      save iset,gset
      data iset /0/
c     
      if (iset.eq.0) then 
 1       v1 = 2.0d0*RANDOM()-1.0d0
         v2 = 2.0d0*RANDOM()-1.0d0
         rsq = v1*v1 + v2*v2
c     
c---- tester si (v1,v2) est dans cercle unite
c     
         if (rsq.GE.1.0d0 .OR. rsq.EQ.0.0d0) then
            goto 1
         endif
c     
         fac = SQRT(-2.0d0*LOG(rsq)/rsq)
         gset = v1*fac
         RandomGauss = v2*fac
c     
c---- iset est mis a 1 pour signaler que gset contient deja un tirage
c     pour le prochain appel de "TirageGauss"
c     
         iset = 1
      else
         RandomGauss = gset
         iset = 0
      end if
c
      end
c
