c%%%%%%% CONTIENT DES ERREURS : TROP INEFFICACE !!!
c
c======================================================================
c     integration selon Runge Kutta
c     Version 1 : Selon Dormand-Prince 5(4)
c======================================================================
c     
c----------------------------------------------------------------------
      subroutine RKDOPRI54
     &           (X,DimX,NX,T,NC,RoutDeriv,tmin,tmax,dt0,Prec)
c----------------------------------------------------------------------
c     X(DimX,0:NX) Suite de vecteurs de dim DimX : Solution de l'ED
c                  de dimension DimX (ordre 1)
c                  X(:,0) doit contenir la condition initiale
c     DimX         Dimensionnalite de l'ED d'ordre 1
c     NX           Taille de X (nb de pas max.)
c     T(0:NX)      Suite des temps (var indep) correspondants aux X
c                  T(0) est l'instant initial
c     NC           En sortie : Nb de pas calcules
c                  Si Nc = 0 -> Erreur
c                  Si Nc < 0 -> NX trop faible pour integrer jusqu'en tmax
c     RoutDeriv    Routine de calcul de la derivee en X au temps t
c                  est appellee par "call RoutDer(X,t,Xp)" ou X et
c                  Xp sont des vect de dim DimX. Xp est la derivee 
c                  calculee
c     tmin,tmax    bornes d'integration
c     dt0          pas de temps initial ("first guess")
c     Prec         precision requise
c----------------------------------------------------------------------
c     
      implicit none
      include 'RK2dbl_param.f'
c     
      integer DimX,NX,NC
      double precision X(DimX,0:NX),T(0:NX),tmin,tmax,dt0,Prec
      external RoutDeriv
c     
      double precision minuscule
      parameter (minuscule=1.0d-30)
      double precision c2,c3,c4,c5,c6,c7,
     &     a21, a31,a32, a41,a42,a43, a51,a52,a53,a54,
     $     a61,a62,a63,a64,a65,
     &     a71,a72,a73,a74,a75,a76,
     &     bc1,bc2,bc3,bc4,bc5,bc6,bc7
c     
      parameter ( c2=1.d0/5.d0, c3=3.d0/10.d0, c4=4.d0/5.d0,
     $     c5=8.d0/9.d0, c6=1.d0, c7=1.d0 )
      parameter (
     &     a21=1.d0/5.d0,
     &     a31=3.d0/40.d0, a32=9.d0/40.d0,
     &     a41=44.d0/45.d0, a42=-56.d0/15.d0, a43=32.d0/9.d0,
     &     a51=19372.d0/6561.d0, a52=-25360.d0/2187.d0,
     $     a53=64448.d0/6561.d0, a54=-212.d0/729.d0,
     &     a61=9017.d0/3168.d0, a62=-355.d0/33.d0, a63=46732.d0/5247.d0,
     $     a64=49.d0/176.d0, a65=-5103.d0/18656,
     &     a71=35.d0/384.d0, a72=0.d0, a73=500.d0/1113.d0,
     $     a74=125.d0/192.d0, a75=-2187.d0/6784.d0, a76=11.d0/84.d0 )
      parameter (
     $     bc1= 5179.d0/57600.d0, bc2=0.d0, bc3=7571.d0/16695.d0,
     $     bc4=393.d0/640.d0,bc5=-92097.d0/339200.d0,
     $     bc6=187.d0/2100.d0,bc7=1.d0/40.d0 )
c     
      double precision p,ExpoDim,ExpoAug
      parameter (p=4.d0, ExpoDim=-0.25d0,ExpoAug=-0.20d0)
      
      double precision CoefSecur, FacMin, FacMax
      parameter (CoefSecur = 0.9d0, FacMin = 0.05d0, FacMax = 5.0d0)    
      integer MemRej            ! nb pas pdt lesquel la "memoire" d'un rejet persiste
      parameter (MemRej=3)
c
      double precision dt,dtnew, Tps, SigneTps,
     &     k1(MaxDim),k2(MaxDim),k3(MaxDim),k4(MaxDim),
     $     k5(MaxDim),k6(MaxDim),k7(MaxDim),
     &     XRK(MaxDim), Xordrep(MaxDim),Xordrep1(MaxDim), err,
     $     DivErr(MaxDim)
      integer iDim, iRej, nRej
c......................................................................
c
c TEST COEF
c$$$c\\\\\
c$$$      write(0,*) a31+a32-c3
c$$$      write(0,*) a41+a42+a43-c4
c$$$      write(0,*) a51+a52+a53+a54-c5
c$$$      write(0,*) a61+a62+a63+a64+a65-c6
c$$$      write(0,*) a71+a72+a73+a74+a75+a76-c7
c$$$      write(0,*) bc1+bc2+bc3+bc4+bc5+bc6+bc7
c$$$      stop
c$$$c/////
      if (tmax.GT.tmin) then
         SigneTps = 1.0
      else
         SigneTps = -1.0
      end if 
c
      nRej = 0
      iRej = 0
      NC = 0
      dt = dt0
      Tps = tmin
      T(0) = tmin
c
c---- calcul de la valeur initiale de la derivee
c
      do iDim=1,DimX
         XRK(iDim) = X(iDim,NC)
      end do 
      call RoutDeriv(XRK,Tps,k1)
c     
c---- Calcul du pas selon RK
c     ======================
c
 1    continue
c     
c---- Tester si il faut et si on peut continuer
c     
      if ((SigneTps*Tps).GE.(SigneTps*tmax)) then
         goto 999
      else if (NC.GE.NX) then
         NC = -1
         goto 999
      end if

      if (Tps.EQ.Tps+dt) then
         write(0,*) '!!! La procedure RKDOPRI54 stagne !!!'
         stop
      end if
c     
c---- Calcul k1
c    
      if (iRej.GT.0) then       ! si il y a eu rejet, ce n'est pas la peine
                                ! de recalculer k1
         do iDim=1,DimX
            XRK(iDim) = X(iDim,NC)
         end do 
         call RoutDeriv(XRK,Tps,k1)
      end if
c     
c---- Calcul k2
c     
      do iDim=1,DimX
         XRK(iDim) = X(iDim,NC) + dt*a21*k1(iDim)
      end do 
      call RoutDeriv(XRK,Tps+c2*dt,k2)
c     
c---- Calcul k3
c     
      do iDim=1,DimX
         XRK(iDim) = X(iDim,NC) + dt* (a31*k1(iDim) + a32*k2(iDim))
      end do 
      call RoutDeriv(XRK,Tps+c3*dt,k3)      
c     
c---- Calcul k4
c     
      do iDim=1,DimX
         XRK(iDim) = X(iDim,NC) + dt* (a41*k1(iDim) + a42*k2(iDim) +
     &        a43*k3(iDim))
      end do 
      call RoutDeriv(XRK,Tps+c4*dt,k4)
c     
c---- Calcul k5
c     
      do iDim=1,DimX
         XRK(iDim) = X(iDim,NC) + dt* (a51*k1(iDim) + a52*k2(iDim) +
     &        a53*k3(iDim) + a54*k4(iDim))
      end do 
      call RoutDeriv(XRK,Tps+c5*dt,k5)
c     
c---- Calcul k6
c     
      do iDim=1,DimX
         XRK(iDim) = X(iDim,NC) + dt* (a61*k1(iDim) + a62*k2(iDim) +
     &        a63*k3(iDim) + a64*k4(iDim) + a65*k5(iDim)) 
      end do 
      call RoutDeriv(XRK,Tps+c6*dt,k6)
c     
c---- Calcul k7 et Calcul de X a l'ordre p-1
c     
      do iDim=1,DimX
         XOrdrep1(iDim) = X(iDim,NC) + dt* (a71*k1(iDim) + 
     &        a73*k3(iDim) + a74*k4(iDim) + a75*k5(iDim) + a76*k6(iDim))
c
c !!!! Terme a72 supprime (car nul)
c
      end do 
      call RoutDeriv(XOrdrep1,Tps+c7*dt,k7)
c     
c---- Calcul de X a l'ordre p
c     
      do iDim=1,DimX
         XOrdrep(iDim) = X(iDim,NC) + dt * (
     &        bc1*k1(iDim) + bc3*k3(iDim) + bc4*k4(iDim) +  bc5*k5(iDim) +  
     &        bc6*k6(iDim) + bc7*k7(iDim) )
      end do 
c
c !!!! Terme bc2 supprime (car nul)
c     
c---- Calcul de l'erreur (relative) par comparaison des estimateurs d'ordre p-1 et p
c     
      do iDim=1,DimX
         DivErr(iDim) = ABS(XOrdrep(iDim))+ABS(dt*k1(iDim))+minuscule
      end do 
      err = 1.d-20
      do iDim=1,DimX
         err = MAX( ABS(XOrdrep1(iDim)-XOrdrep(iDim))/DivErr(iDim) , err )
      end do 
c     
c---- Calcul du dt ideal en fct de la precision requise
c     Acceptation ou rejet du pas ?
c
      err = err/Prec
      if (err.GT.1.0d0) then    ! il faut diminuer le pas de temps et
                                ! rejeter le pas d'integration
         dtnew=max(FacMin,CoefSecur*err**ExpoDim) * dt
         nRej = nRej+1
         iRej = 0
      else                      ! on peut augmenter le pas de temps
                                ! et accepter le pas d'integration
         if (iRej.LT.MemRej) then
            dtnew=dt
         else
            dtnew=min(FacMax,CoefSecur*err**ExpoAug) * dt
         end if
         iRej = iRej+1
         NC = NC + 1
         Tps = Tps + dt
         T(NC) = Tps
         do iDim=1,DimX
            X(iDim,NC) = XOrdrep(iDim)
         end do  
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
      write(0,*) 'NbOK,NbRej : ',Nc,nRej 
c
      end
c
c----------------------------------------------------------------------
      subroutine init_RKDOPRI54_1pas (Prec)
c----------------------------------------------------------------------
c     Fait les initialisation necessaires a l'emploi de la routine 
c     ci-dessous
c     Fixe la valeur de la precision relative exigee
c
      implicit none
      include 'RK2dbl_common.f'
c
      double precision Prec
c
c......................................................................
c
      Prec__RK = Prec
      iRej__RK = 0
c
      end
c   
c----------------------------------------------------------------------
      subroutine RKDOPRI54_1pas (X,DimX,T,dt,RoutDeriv)
c----------------------------------------------------------------------
c     Integre une ED d'ordre 1 sur 1 pas de temps
c
c     X(DimX) Valeur du vecteur-variable au debut du pas
c             est modifie pour devenir la valeur a la fin du pas
c
c     DimX    Dimensionnalite de l'ED d'ordre 1
c
c     T       Valeur du temps (var indep.) au debut du pas
c             est modifie pour devenir la valeur a la fin du pas
c             
c     RoutDeriv    Routine de calcul de la derivee en X au temps t
c                  est appellee par "call RoutDer(X,t,Xp)" ou X et
c                  Xp sont des vect de dim DimX. Xp est la derivee 
c                  calculee
c
c     dt           pas de temps initial ("first guess")
c                  est modifie pour devenir l'evaluation initiale du 
c                  pas de temps suivant
c----------------------------------------------------------------------
c     
      implicit none
      include 'RK2dbl_common.f'
      include 'RK2dbl_param.f'
c     
      integer DimX
      double precision X(DimX),T,dt
      external RoutDeriv
c     
      double precision c2,c3,c4,c5,c6,c7,
     &     a21, a31,a32, a41,a42,a43, a51,a52,a53,a54,
     $     a61,a62,a63,a64,a65,
     &     a71,a72,a73,a74,a75,a76,
     &     bc1,bc2,bc3,bc4,bc5,bc6,bc7
c     
      parameter ( c2=1.d0/5.d0, c3=3.d0/10.d0, c4=4.d0/5.d0,
     $     c5=8.d0/9.d0, c6=1.d0, c7=1.d0 )
      parameter (
     &     a21=1.d0/5.d0,
     &     a31=3.d0/40.d0, a32=9.d0/40.d0,
     &     a41=44.d0/45.d0, a42=-56.d0/15.d0, a43=32.d0/9.d0,
     &     a51=19372.d0/6561.d0, a52=-25360.d0/2187.d0,
     $     a53=64448.d0/6561.d0, a54=-212.d0/729.d0,
     &     a61=9017.d0/3168.d0, a62=-355.d0/33.d0, a63=46732.d0/5247.d0,
     $     a64=49.d0/176.d0, a65=-5103.d0/18656,
     &     a71=35.d0/384.d0, a72=0.d0, a73=500.d0/1113.d0,
     $     a74=125.d0/192.d0, a75=-2187.d0/6784.d0, a76=11.d0/84.d0 )
      parameter (
     $     bc1= 5179.d0/57600.d0, bc2=0.d0, bc3=7571.d0/16695.d0,
     $     bc4=393.d0/640.d0,bc5=-92097.d0/339200.d0,
     $     bc6=187.d0/2100.d0,bc7=1.d0/40.d0 )
c     c     
      double precision p
      parameter (p=5.)
c---- ordre de la methode
c     
      double precision Fac, FacMin, FacMax
      parameter (Fac = 0.85, FacMin = 0.05, FacMax = 1.5)
c---- facteur de correction du pas dt
c
      double precision tRK, FacMaxEff, SigneTps, dtold, 
     &     k1(MaxDim),k2(MaxDim),k3(MaxDim),k4(MaxDim),k5(MaxDim),k6(MaxDim),k7(MaxDim),
     &     XRK(MaxDim), Xordrep(MaxDim),Xordrep1(MaxDim), err, DivErr(MaxDim)
      integer iDim
c
c--- nombre de rejets du pas
c
      integer NbRej_RK
      common /RK_cmpl/ NbRej_RK
c......................................................................
c     
      NbRej_RK = 0
      dt_done_RK = 0.0d0

      if (dt.GT.0) then
         SigneTps = 1.0
      else
         SigneTps = -1.0
      end if 
c
c---- Calcul du pas selon RK
c     ======================
c
 1    continue
c     
c---- Calcul k1
c     
      tRK = T
      do iDim=1,DimX
         XRK(iDim) = X(iDim)
      end do 
      call RoutDeriv(XRK,tRK,k1)
c     
c---- Calcul k2
c     
      tRK = T+ c2*dt
      do iDim=1,DimX
         XRK(iDim) = X(iDim) + dt*a21*k1(iDim)
      end do 
      call RoutDeriv(XRK,tRK,k2)
c     
c---- Calcul k3
c     
      tRK = T+ c3*dt
      do iDim=1,DimX
         XRK(iDim) = X(iDim) + dt* (a31*k1(iDim) + a32*k2(iDim))
      end do 
      call RoutDeriv(XRK,tRK,k3)      
c     
c---- Calcul k4
c     
      tRK = T + c4*dt
      do iDim=1,DimX
         XRK(iDim) = X(iDim) + dt* (a41*k1(iDim) + a42*k2(iDim) +
     &        a43*k3(iDim))
      end do 
      call RoutDeriv(XRK,tRK,k4)
c     
c---- Calcul k5
c     
      tRK = T + c5*dt
      do iDim=1,DimX
         XRK(iDim) = X(iDim) + dt* (a51*k1(iDim) + a52*k2(iDim) +
     &        a53*k3(iDim) + a54*k4(iDim))
      end do 
      call RoutDeriv(XRK,tRK,k5)
c     
c---- Calcul k6
c     
      tRK = T + c6*dt
      do iDim=1,DimX
         XRK(iDim) = X(iDim) + dt* (a61*k1(iDim) + a62*k2(iDim) +
     &        a63*k3(iDim) + a64*k4(iDim) + a65*k5(iDim)) 
      end do 
      call RoutDeriv(XRK,tRK,k6)
c     
c---- Calcul k7 et Calcul de X a l'ordre p-1
c     
      tRK = T + c7*dt
      do iDim=1,DimX
         XOrdrep1(iDim) = X(iDim) + dt* (a71*k1(iDim) + 
     &        a73*k3(iDim) + a74*k4(iDim) + a75*k5(iDim) + a76*k6(iDim)) 
c
c !!!! Terme a72 supprime (car nul)
c
      end do 
      call RoutDeriv(XOrdrep1,tRK,k7)
c     
c---- Calcul de X a l'ordre p
c     
      do iDim=1,DimX
         XOrdrep(iDim) = X(iDim) + dt * (
     &        bc1*k1(iDim) + bc3*k3(iDim) + bc4*k4(iDim) +  bc5*k5(iDim) +  
     &        bc6*k6(iDim) + bc7*k7(iDim) )
      end do 
c
c !!!! Terme bc2 supprime (car nul)
c     
c---- Calcul de l'erreur (relative) par comparaison des estimateurs d'ordre p-1 et p
c     
      do iDim=1,DimX
         DivErr(iDim) = ABS(XOrdrep(iDim))+1.0e-8
      end do 
      err = 1.e-20
      do iDim=1,DimX
         err = MAX( ABS(XOrdrep1(iDim)-XOrdrep(iDim))/DivErr(iDim) , err )
      end do 
      !!err = err !!/ (2.0**(p-1.)-1.0)     
c     
c---- Calcul du dt ideal en fct de la precision requise
c     
      FacMaxEff = FacMax
      if (iRej__RK.LT.3) FacMaxEff = 1.0
      dtold = dt
      dt = dt * MIN( FacMaxEff, MAX( FacMin , Fac*(Prec__RK/err)**0.20d0 ) ) 
c
c---- Acceptation ou rejet du pas ?
c     
      if (err.LT.Prec__RK) then
c     
c---- EN CAS D'ACCEPTATION
c     
         iRej__RK = iRej__RK+1
         T = T + dtold
         dt_done_RK = dtold
         do iDim=1,DimX
            X(iDim) = XOrdrep(iDim)
         end do  
      else
c     
c---- EN CAS DE REJET
c 
         NbRej_RK = NbRej_RK+1
         iRej__RK = 0
         goto 1
      end if
c
c---- Fin de l'integration
c
      end
c

c----------------------------------------------------------------------
      subroutine IntegrateOverInterval_RKDOPRI54
     $     (X_ini,DimX,RoutDeriv,t_ini,t_end,dt_guess,accur,
     $     X_end,Nsteps)
c----------------------------------------------------------------------
      implicit none
      include 'RK2dbl_param.f'
c
c arguments :
c ^^^^^^^^^^^
      integer DimX
      double precision X_ini(DimX) ! initial conditions (at t_ini)
      external RoutDeriv
      double precision t_ini,t_end
      double precision dt_guess ! guess of time step
      double precision accur    ! (relative) accuracy 
      double precision X_end(DimX) ! Final value (result of integration)
      integer Nsteps            ! Number of steps that were required
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision t,dt,t_max
      integer i
c......................................................................
c
      if (t_ini.GE.t_end) then
         Nsteps = -1
         write(0,*) '!!! IntegrateOverInterval_RKDOPRI54:',
     $        ' needs t_ini < t_end !!!'
         return
      end if

      Nsteps = 0
      call init_RKDOPRI54_1pas(accur)

      dt = dt_guess
      t = t_ini
      do i=1,DimX
         X_end(i)=X_ini(i)
      end do
      !t_max = t_ini+(1.0d0-eps)*(t_end+t_ini) 
     
      do while (t.LT.t_end) 
         dt = max(min(dt,t_end-t),1.0d-30)
         call RKDOPRI54_1pas(X_end,DimX,t,dt,RoutDeriv)
         Nsteps = Nsteps+1
      end do
c
      end
c

      

