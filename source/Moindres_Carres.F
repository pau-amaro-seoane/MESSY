
c----------------------------------------------------------------------
      subroutine MoindresCarres(X,Y,N,a,b)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      integer N
      double precision X(N),Y(N),a,b
c
c var locales :
c ^^^^^^^^^^^^^
      double precision sX,sX2,sXY,sY,delta
      integer i
c......................................................................
c
      sX  = 0.0d0
      sX2 = 0.0d0
      sXY = 0.0d0
      sY  = 0.0d0
      do i=1,N
         sX  = sX  + X(i)
         sX2 = sX2 + X(i)*X(i)
         sXY = sXY + X(i)*Y(i)
         sY  = sY  + Y(i)
      end do
      delta = N*sX2-sX**2
      a = (N*sXY-sX*sY)  /Delta
      b = (sX2*sY-sXY*sX)/Delta
c
      end
c
