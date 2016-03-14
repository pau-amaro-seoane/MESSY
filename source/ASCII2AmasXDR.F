c======================================================================
c     Read an ascii file and convert it into a *AMAS.xdr file
c     for the Monte Carlo code. Read from stdin. Write file specified
c     as last command-line option
c
c     With obtion -atakan the format must be R Vr Vt m
c     With obtion -atakan-nomass the format must be R Vr Vt
c     With obtion -nbody the format must be  m x y z vx vy vz
c     With obtion -nbodyp the format must be  m x y z vx vy vz pot
c     With obtion -messy the format must be  m R T J 
c     With obtion -messy2 the format must be  R m T J 
c
c     By default '-nbody' format is assumed
c     
c     Use '-NoMassInfo' if there is no column for the masses (single-mass cluster)
c
c     Option '-Sort' produces and output file sorted in radius
c            '-N Nsub' select (at random) a sub-sample of Nsub particles
c            '-mu M' includes a central MBH with mass mu (M_BH/M_stars)
c                    Warning: A file *TN.xdr will NOT be created!
c            '-MBH' or '-FirstIsMBH' indicates that the first particle
c                   of the input file has to be treated as a central MBH
c                   (works for format -nbody only)
c     In the case of option -nbodyp, the grav potential is also read 
c     and used to determine if particles are bound. If the particle is 
c     not bound, a very large value is given to its radius (R_perte in 
c     EvolAmas_Pertes_inc.f, should be 1e30)
c======================================================================

      implicit none
      include 'Param_Amas.f'
      include 'EvolAmas_Pertes_inc.f'
      include 'VarSE.f'
      include 'fxdr.inc'
c
c local variables:
c ^^^^^^^^^^^^^^^^
      double precision Ecin, Egrav_stel, Egrav_bh, Mclust,
     $     mi, Ri, Vri, Vti,
     $     xi,yi,zi, vxi,vyi,vzi, Jxi,Jyi,Jzi, Ti, Ji, Pi,
     $     convM, convR, convT, convJ, mu  
      double precision M_SE_copy(iDimSE), R_SE_copy(iDimSE),
     $     T_SE_copy(iDimSE), J_SE_copy(iDimSE)
      integer i,l,idum,iXDR, ixdr_err, N_to_remove
      integer iarg, narg
      character*64 arg
      character*256 FileName
      integer iFormat, iSort
      logical lFirstPartIsMBH, lConvertUnits, lSpecifiedMasses,
     $     lSpecifiedPot
      double precision mMBH, xMBH,yMBH,zMBH, vxMBH,vyMBH,vzMBH
      integer Nsub
      integer iRandomSeed
c
c constants :
c ^^^^^^^^^^^
      integer iFormAto, iFormAtoNoMass, iFormNbody, iFormMESSY
      parameter (iFormAto=10, iFormAtoNoMass=11, iFormNbody=20,
     $     iFormMESSY=30,  iFormMESSY2=31)
c
c functions :
c ^^^^^^^^^^^
      integer LongChaine, iargc
      double precision random
c......................................................................
c
c---- Read command line arguments
c
      iRandomSeed=978463
      iFormat=iFormNbody
      iSort=0
      lFirstPartIsMBH=.FALSE.
      lConvertUnits=.TRUE.
      lSpecifiedMasses=.TRUE.
      lSpecifiedPot=.FALSE.
      Nsub=0

      mu = 0.0d0
      mMBH=0.0d0
      xMBH=0.0d0
      yMBH=0.0d0
      zMBH=0.0d0
      vxMBH=0.0d0
      vyMBH=0.0d0
      vzMBH=0.0d0

      iarg = 1
      narg = iargc()
      call getarg(iarg,arg)

      do while(arg(1:1).EQ.'-')
         if (arg .EQ. '-mu' .OR. arg .EQ. '-AddMBH') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) mu
         elseif (arg .EQ. '-nbody') then
            iFormat=iFormNbody
         elseif (arg .EQ. '-nbodyp') then
            iFormat=iFormNbody
            lSpecifiedPot=.TRUE.
         elseif (arg .EQ. '-atakan') then
            iFormat=iFormAto
         elseif (arg .EQ. '-atakan-nomass') then
            iFormat=iFormAto
            lSpecifiedMasses=.FALSE.
         elseif (arg .EQ. '-messy') then
            iFormat=iFormMessy
         elseif (arg .EQ. '-messy2') then
            iFormat=iFormMessy2
         elseif (arg .EQ. '-MBH' .OR.
     $           arg .EQ. '-FirstIsMBH') then
            lFirstPartIsMBH=.TRUE.
         elseif (arg .EQ. '-DontConvertUnits') then
            lConvertUnits=.FALSE.
         elseif (arg .EQ. '-NoMassInfo') then
            lSpecifiedMasses=.FALSE.
         elseif (arg .EQ. '-Sort') then
            iSort=1
         elseif (arg .EQ. '-N') then ! Only output N particles (at random)
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) Nsub
         elseif (arg .EQ. '-RandomSeed') then
            iarg = iarg+1
            call getarg(iarg,arg)
            read(arg,*) iRandomSeed
         else
            write(0,*) '!!! option ',arg(1:LongChaine(arg)),
     $           ' non comprise !!!'
            call exit(1)
         end if
         iarg = iarg+1
         if (iarg.LE.narg) call getarg(iarg,arg)
      end do

      call init_random(iRandomSeed)
c
c---- Read particle data
c
      mi = -1.0d0
      if (iFormat.EQ.iFormAto) then
         NbSE=0
         do while (.TRUE.)
            if (lSpecifiedMasses) then
               read(5,*,END=9) Ri, Vri, Vti, mi
            else
               read(5,*,END=9) Ri, Vri, Vti
            end if
            NbSE=NbSE+1
            if (NbSE.GT.iDimSE) then
               write(0,*) '!!! NbSE=',NbSE,' exceeds maximum number of',
     $              ' particles, iDimSE=',iDimSE,' !!!'
               call exit(1)
            end if
            M_SE(NbSE)=mi
            R_SE(NbSE)=Ri
            T_SE(NbSE)=0.5d0*(Vri**2+Vti**2)
            J_SE(NbSE)=Ri*Vti
         end do

      elseif (iFormat.EQ.iFormNbody) then
         NbSE=0
         if (lFirstPartIsMBH) then
            read(5,*,END=9) mMBH, xMBH,yMBH,zMBH, vxMBH,vyMBH,vzMBH
            write(0,*) '> central MBH of mass ',mMBH
         end if
         do while (.TRUE.)
            if (lSpecifiedMasses) then
               if (lSpecifiedPot) then
                  read(5,*,END=9) mi, xi,yi,zi, vxi,vyi,vzi, pi
               else
                  read(5,*,END=9) mi, xi,yi,zi, vxi,vyi,vzi
               end if
            else
               if (lSpecifiedPot) then
                  read(5,*,END=9) xi,yi,zi, vxi,vyi,vzi, pi
               else
                  read(5,*,END=9) xi,yi,zi, vxi,vyi,vzi
               end if
            end if
            if (lFirstPartIsMBH) then
               xi=xi-xMBH
               yi=yi-yMBH
               zi=zi-zMBH
               vxi=vxi-vxMBH
               vyi=vyi-vyMBH
               vzi=vzi-vzMBH
            end if
            NbSE=NbSE+1
            if (NbSE.GT.iDimSE) then
               write(0,*) '!!! NbSE=',NbSE,' exceeds maximum number of',
     $              ' particles, iDimSE=',iDimSE,' !!!'
               call exit(1)
            end if
            M_SE(NbSE)=mi
            R_SE(NbSE)=sqrt(xi**2+yi**2+zi**2)
            T_SE(NbSE)=0.5d0*(vxi**2+vyi**2+vzi**2)
            Jxi=yi*vzi-zi*vyi
            Jyi=zi*vxi-xi*vzi
            Jzi=xi*vyi-yi*vxi
            J_SE(NbSE)=sqrt(Jxi**2+Jyi**2+Jzi**2)
            if (lSpecifiedPot.AND.(pi+T_SE(NbSE).GT.0.0d0)) then
               R_SE(NbSE)=R_perte
            end if
         end do

      elseif (iFormat.EQ.iFormMessy .OR. iFormat.EQ.iFormMessy2) then
         NbSE=0
         do while (.TRUE.)
            if (lSpecifiedMasses) then
               if (iFormat.EQ.iFormMessy2) then
                  read(5,*,END=9) ri, mi, Ti, Ji
               else
                  read(5,*,END=9) mi, ri, Ti, Ji
               end if
            else
               read(5,*,END=9) ri, Ti, Ji
            end if
            NbSE=NbSE+1
            if (NbSE.GT.iDimSE) then
               write(0,*) '!!! NbSE=',NbSE,' exceeds maximum number of',
     $              ' particles, iDimSE=',iDimSE,' !!!'
               call exit(1)
            end if
            M_SE(NbSE)=mi
            R_SE(NbSE)=ri
            T_SE(NbSE)=Ti
            J_SE(NbSE)=Ji
         end do
      end if

 9    continue


      if (.NOT.lSpecifiedMasses) then
         mi=1.0d0/NbSE
         do i=1,NbSE
            M_SE(i)=mi
         end do
      end if

      write(0,110) NbSE
 110  format(' > Ascii file contains ',I7,' particles')
c
c---- Select a smaller number of particles if required
c
      if (Nsub.GT.0 .AND. Nsub.LT.NbSE) then
         N_to_remove=NbSE-Nsub
         do while (N_to_remove.GT.0)
            i=1+int(random()*NbSE)
            if (M_SE(i).GT.0.0d0) then
               M_SE(i)=-1.0d0
               N_to_remove=N_to_remove-1
            end if
         end do
         i=1
         NbSE=0
         do while (NbSE.LT.Nsub)
            if (M_SE(i).GT.0.d0) then
               NbSE=NbSE+1
               M_SE(NbSE) = M_SE(i)
               R_SE(NbSE) = R_SE(i)
               J_SE(NbSE) = J_SE(i)
               T_SE(NbSE) = T_SE(i)
            end if
            i=i+1
         end do
      end if
c
c---- Convert to N-body units in which Egrav=-0.5 & Ecin=0.25
c
      call TriAmas

      Mclust=0.0d0
      do l=1,NbSE
         i=iListeRay(l)
         Mclust=Mclust+M_SE(i)
      end do

      if (lFirstPartIsMBH) then
         mu=mMBH/Mclust
      else
         mMBH=mu*Mclust
      end if

      Mclust=0.0d0
      Egrav_bh=0.0d0
      Egrav_stel=0.0d0
      Ecin=0.0d0
      do l=1,NbSE
         i=iListeRay(l)
         Egrav_stel=Egrav_stel -
     $        (Mclust+FracMasseInt*M_SE(i))*M_SE(i)/R_SE(i)
         Egrav_bh=Egrav_bh - mMBH*M_SE(i)/R_SE(i)
         Ecin=Ecin+M_SE(i)*T_SE(i)
         Mclust=Mclust+M_SE(i)
      end do
      write(0,120) Mclust,Ecin,Egrav_stel,Egrav_bh,
     $     (Egrav_stel+Egrav_bh)/Ecin
         
c
c     N-body units: set Mcluster to one by scaling masses
c
      convM = 1.0d0/Mclust
      Egrav_bh = convM**2*Egrav_bh
      Egrav_stel = convM**2*Egrav_stel
      Ecin = convM*Ecin
      mMBH = convM*mMBH
c
c     N-body units: set Egrav_stel to -1/2 by scaling distances
c
      convR = -2.0d0*Egrav_stel
      Egrav_bh = Egrav_bh/convR
      Egrav_stel = Egrav_stel/convR
c
c     Global virial : set Ecin_tot to -1/2 * Egrav_tot by scaling velocities
c
      convT= -0.5d0*(Egrav_bh+Egrav_stel)/Ecin
      Ecin = convT * Ecin
      convJ= convR*sqrt(convT)

      write(0,100) convM, convR, convT, convJ
 100  format(
     $     '> conversion factors to impose N-body units and ',
     $     'virial equilibrium:',/,
     $     '     convM, convR, convT, convJ : ',4(E12.5,' '))

      if (lConvertUnits) then
         do i=1,NbSE
            M_SE(i)=convM*M_SE(i)
            R_SE(i)=convR*R_SE(i)
            T_SE(i)=convT*T_SE(i)
            J_SE(i)=convJ*J_SE(i)
         end do
         
         Egrav_bh=0.0d0
         Egrav_stel=0.0d0
         Mclust=0.0d0
         Ecin=0.0d0
         do l=1,NbSE
            i=iListeRay(l)
            Egrav_stel=Egrav_stel -
     $           (Mclust+FracMasseInt*M_SE(i))*M_SE(i)/R_SE(i)
            Egrav_bh=Egrav_bh - mMBH*M_SE(i)/R_SE(i)
            Mclust=Mclust+M_SE(i)
            Ecin=Ecin+M_SE(i)*T_SE(i)
         end do
         write(0,120) Mclust,Ecin,Egrav_stel,Egrav_bh,
     $        (Egrav_stel+Egrav_bh)/Ecin
 120     format('> Mclust,Ecin,Egrav_stel,Egrav_bh,Egrav_tot/Ecin : ',
     $        5(E12.5,' '))
      else
         write (0,*) 'WARNING: No unit conversion/adjustment done'
      endif
c
c---- Sort if required
c
      if (iSort.EQ.1) then 
         do i=1,NbSE
            M_SE_copy(i)=M_SE(i)
            R_SE_copy(i)=R_SE(i)
            T_SE_copy(i)=T_SE(i)
            J_SE_copy(i)=J_SE(i)
         end do
         do l=1,NbSE
            i=iListeRay(l)
            M_SE(l) = M_SE_copy(l)
            R_SE(l) = R_SE_copy(l)
            T_SE(l) = T_SE_copy(l)
            J_SE(l) = J_SE_copy(l)
         end do
      end if
c
c---- Write *AMAS.xdr file
c
      call getarg(narg,FileName)
      iXDR = initxdr(FileName,'w',.TRUE.)
      call EcritAmas_XDR(iXDR)
      ixdr_err = ixdrclose(iXDR)
c
      end
c
