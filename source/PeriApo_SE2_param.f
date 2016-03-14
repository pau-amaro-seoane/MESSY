      double precision cR_inf,cR_sup
      parameter (cR_inf=0.5d0,cR_sup=2.0d0)

      integer Max_cR
      parameter (Max_cR=100)

      double precision epsR,epsR1
      parameter (epsR=1.0d-6,epsR1=1.0d0-epsR)

      double precision Tol_R
      parameter (Tol_R=1.0d-9) ! determine la tolerance  pour R<Rperi ou R>Rapo

      double precision Tol_circul
      parameter (Tol_circul=1.0d-9) ! determine la tolerance sur la
                                    ! detection des orbites ~circulaires
      double precision eps_dQdR
      parameter (eps_dQdR=1.0d-10)  ! tolerance pour la resolution de dQdR(Rc)=0

      integer Max_Bissec            ! nb max de pas de bissection
      parameter (Max_Bissec=1000)   ! pour la resolution de dQdR(Rc)=0
