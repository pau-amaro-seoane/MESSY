c
c======================================================================
c     Calcul numerique d'integrale 1D selon Runge Kutta
c     Dormand-Prince 5(4)
c======================================================================
c     
c----------------------------------------------------------------------
      double precision function Int_RKDOPRI54(Fct,tmin,tmax,dt0,Prec)
c----------------------------------------------------------------------
c     Fct          fonction reelle a integrer (Nom de la routine)
c     tmin,tmax    bornes d'integration
c     dt0          pas de temps initial ("first guess")
c     Prec         precision requise
c----------------------------------------------------------------------
c     
      implicit none
c     
      double precision tmin,tmax,dt0,Prec
      double precision Fct
      external Fct
c     
      double precision c2,c3,c4,c5,c6,c7,
     &     a21, a31,a32, a41,a42,a43, a51,a52,a53,a54, a61,a62,a63,a64,a65,
     &     a71,a72,a73,a74,a75,a76,
     &     bc1,bc2,bc3,bc4,bc5,bc6,bc7
c     
      parameter ( c2=1.d0/5.d0, c3=3.d0/10.d0, c4=4.d0/5.d0, c5=8.d0/9.d0, c6=1.d0, c7=1.d0 )
      parameter (
     &     a21=1.d0/5.d0,
     &     a31=3.d0/40.d0, a32=9.d0/40.d0,
     &     a41=44.d0/45.d0, a42=-56.d0/15.d0, a43=32.d0/9.d0,
     &     a51=19372.d0/6561.d0, a52=-25360.d0/2187.d0, a53=64448.d0/6561.d0, a54=-212.d0/729.d0,
     &     a61=9017.d0/3168.d0, a62=-355.d0/33.d0, a63=46732.d0/5247.d0, a64=49.d0/176.d0, a65=-5103.d0/18656.d0,
     &     a71=35.d0/384.d0, a72=0.d0, a73=500.d0/1113.d0, a74=125.d0/192.d0, a75=-2187.d0/6784.d0, a76=11.d0/84.d0 )
      parameter ( bc1= 5179.d0/57600.d0, bc2=0.d0, bc3=7571.d0/16695.d0, bc4=393.d0/640.d0, 
     &            bc5=-92097.d0/339200.d0, bc6=187.d0/2100.d0, bc7=1.d0/40.d0 )
c     
      double precision p
      parameter (p=5.d0)
c---- ordre de la methode
c     
      double precision Fac, FacMin, FacMax
      parameter (Fac = 0.85d0, FacMin = 0.05d0, FacMax = 1.5d0)
c---- facteur de correction du pas dt
c
      double precision dt,dtnew, tRK, FacMaxEff,
     &     k1,k2,k3,k4,k5,k6,k7,
     &     XRK, Xordrep,Xordrep1, err, DivErr, X, T
      integer iRej

      integer Npas,Nacc
      common /Nb_Pas/ Npas,Nacc
c......................................................................
c

      if (tmax.LT.tmin) then
         write(0,*) "!!! La borne inf doit etre inf. a la borne sup !!!"
         stop
      end if 
c
      Npas = 0
      Nacc = 0
      iRej = 0
      dt  = dt0
      T = tmin                  ! variable d integration
      X   = 0.0                 ! valeur de l integrale
c     
c---- Calcul du pas selon RK
c     ======================
c
 1    continue
c!!!!!
c!!!!!      write(0,*)'Npas=',Npas,' T=',T,' dt=',dt
c!!!!!
      Npas=Npas+1
c     
c---- Calcul k1
c     
      tRK = T
      XRK = X
      k1=Fct(tRK)
c     
c---- Calcul k2
c     
      tRK = T + c2*dt
      XRK = X + dt*a21*k1
      k2=Fct(tRK)
c     
c---- Calcul k3
c     
      tRK = T + c3*dt
      XRK = X + dt* (a31*k1 + a32*k2)
      k3=Fct(tRK)
c     
c---- Calcul k4
c     
      tRK = T + c4*dt
      XRK = X + dt* (a41*k1 + a42*k2 + a43*k3)
      k4=Fct(tRK)
c     
c---- Calcul k5
c     
      tRK = T + c5*dt
      XRK = X + dt* (a51*k1 + a52*k2 + a53*k3 + a54*k4)
      k5=Fct(tRK)
c     
c---- Calcul k6
c     
      tRK = T + c6*dt
      XRK = X + dt* (a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5) 
      k6=Fct(tRK)
c     
c---- Calcul k7 et Calcul de X a l'ordre p-1
c     
      tRK = T + c7*dt
      XOrdrep1 = X + dt* (a71*k1 + a73*k3 + a74*k4 + a75*k5 + a76*k6) 
c
c !!!! Terme a72 supprime (car nul)
c
      k7=Fct(tRK)
c     
c---- Calcul de X a l'ordre p
c     
      XOrdrep = X + dt * (
     &      bc1*k1 + bc3*k3 + bc4*k4 + bc5*k5 + bc6*k6 + bc7*k7 ) 
c
c !!!! Terme bc2 supprime (car nul)
c     
c---- Calcul de l'erreur (relative) par comparaison des estimateurs d'ordre p-1 et p
c     
      DivErr = MAX( ABS(XOrdrep) , 1.0e-8 ) 
      err = 1.e-20
      err = MAX( ABS(XOrdrep1-XOrdrep)/DivErr , err )
      err = err / (2.0**(p-1.)-1.0)      
c     
c---- Calcul du dt ideal en fct de la precision requise
c     
      FacMaxEff = FacMax
      if (iRej.LT.3) FacMaxEff = 1.0
      dtnew = dt * MIN( FacMaxEff, MAX( FacMin , Fac*(Prec/err)**(1.0/(p+1.0)) ) )
c     
c---- Acceptation ou rejet du pas ?
c     
      if (err.LT.Prec) then
c     
c---- EN CAS D'ACCEPTATION
c
         Nacc = Nacc+1
         iRej = iRej+1
         T = T + dt
         X = XOrdrep 
         dtnew = MIN(dtnew,tmax-T) ! Faire en sorte de ne pas depasser la borne sup.
         if (T.GE.tmax) goto 999 ! integration terminee
      else
         iRej = 0
      end if
c
c---- On repart au debut
c
      dt = dtnew
      goto 1
c
c---- Fin de l'integration
c
 999  continue
      Int_RKDOPRI54 = X
c
      end
c





