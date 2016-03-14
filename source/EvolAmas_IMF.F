c======================================================================
c     Tirage de masses stellaires selon IMF en loi de puissances par 
c     morceaux
c======================================================================

c----------------------------------------------------------------------
      subroutine Tirage_IMF_LPPM(N_morc, Mlim, alpha, Net, Met)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      integer N_morc
      double precision
     $     Mlim(0:N_morc),      ! Les limites en M_sol des "morceaux"
     $     alpha(N_morc)        ! Les exposants pour chaque morceau
      integer Net               ! Le nb de masses a tirer (Salpeter : 2.35)
      double precision Met(Net) ! Les masses a tirer
c
c constantes :
c ^^^^^^^^^^^^
      integer Nmax_morc
      parameter (Nmax_morc=20)
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      double precision Ncum(0:Nmax_morc), C(Nmax_morc), a1, Y
      integer i,j
c
c fonctions :
c ^^^^^^^^^^^
      double precision random
c......................................................................
c
      if (N_morc.GT.Nmax_morc) then
         write(0,*) '!!! L''IMF doit etre composee de moins de ',
     $        Nmax_morc+1,' morceaux !!!'
         call Meurt('Tirage_IMF_LPPM')
      end if
c
c---- Coefficients de continuite de l'IMF
c
      C(1)=1.0d0
      do i=2,N_morc
         C(i) = C(i-1)*Mlim(i-1)**(alpha(i)-alpha(i-1))
      end do
c
c---- Nb cumulatifs aux limites des morceaux
c
      Ncum(0)=0.0d0
      do i=1,N_morc
         a1=1.0d0-alpha(i)
         Ncum(i) = Ncum(i-1)+C(i)/a1*(Mlim(i)**a1-Mlim(i-1)**a1)
      end do
                                ! normalisation a 1
      do i=1,N_morc
         C(i) = 1.00000001d0*C(i)/Ncum(N_morc)
         Ncum(i) = 1.00000001d0*Ncum(i)/Ncum(N_morc)
      end do
c
c---- Tirages par inversion de la distribution cumulative
c
      do j=1,Net
         Y=random()
         call locate0(Ncum,N_morc+1,Y,i)
         i=i-1
         a1=1.0d0-alpha(i+1)
         if (abs(a1).LT.1.0d-4) then
            Met(j)=Mlim(i)*exp((Y-Ncum(i))/C(i+1))
         else
            Met(j)=(Mlim(i)**a1+a1*(Y-Ncum(i))/C(i+1))**(1.0d0/a1)
         end if
      end do
c
      end
c
c----------------------------------------------------------------------
      SUBROUTINE locate0(xx,n,x,j)
      INTEGER j,n
      DOUBLE PRECISION x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END
