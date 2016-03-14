c======================================================================
c     Solve Roche problem.
c     A star of radius R and mass m is encountering an object of mass
c     M. We assume pure (Newtonian) gravitational encounter and tha
c     the star keeps spherical shape. 
c     What is the critical value of d_min/R ?
c======================================================================
#define _TEST_ROUTINE_ 0
#if (_TEST_ROUTINE_>0)
      implicit none
c
c local variables :
c ^^^^^^^^^^^^^^^^^ 
      double precision q,nu

      do while(.TRUE.)
         read(5,*,END=99) q
         call Solve_Roche(q,nu)
         write(6,100) q,nu
 100     format(2(E12.5,' ' ))
      end do
 99   continue
c
      end
#endif
c----------------------------------------------------------------------
      subroutine Solve_Roche(q,nu)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      double precision q        ! input: q = m/M
      double precision nu       ! output: nu = d_min/R (critical)
c
c functions :
c ^^^^^^^^^^^
      double precision rtsafe_rel
c
c local variables :
c ^^^^^^^^^^^^^^^^^ 
      double precision x1,x2
c
c various :
c ^^^^^^^^^
      double precision prec
      parameter (prec=1.0d-6)
      external fdf_roche
      double precision q_fdf
      common /param_fdf_Roche/ q_fdf
c......................................................................
c
      q_fdf=q
      x1=1.0d0
      x2=1.0d0+1/sqrt(q)
      nu = rtsafe_rel(fdf_roche,x1,x2,prec)
c
      end
c
c----------------------------------------------------------------------
      subroutine fdf_roche(nu,f,df)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      double precision nu, f,df
c
c various :
c ^^^^^^^^^
      double precision q_fdf
      common /param_fdf_Roche/ q_fdf
c----------------------------------------------------------------------
c
      f = q_fdf*(nu*(nu-1.0d0))**2 -2.0d0*nu +1.0d0
      df = 2.0d0*(q_fdf*nu*(nu-1.0d0)*(2.0d0*nu-1.0d0)-1.0d0)
c
      end
c
