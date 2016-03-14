c
c---- abscisses pour methode de quadrature de Gauss-Chebyshev
c
c constantes :
c ^^^^^^^^^^^^
      integer n_GC_max,k_GC_max ! ordre max pour Gauss-Chebyshev
      parameter (n_GC_max=100,k_GC_max=((n_GC_max)*(n_GC_max+1))/2)
      
      double precision x_GC(k_GC_max)

      common /common_GC/ x_GC
