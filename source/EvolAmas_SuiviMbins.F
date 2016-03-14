c======================================================================
c     Routines to follow distribution of stars binned my mass
c======================================================================

#if (_MONITOR_MBINS_>0)
c----------------------------------------------------------------------
      subroutine Read_Mbins(iIn)
c----------------------------------------------------------------------
c     read the specifications of mass-bins and the fractional masses
c     for radii to be monitored
c
      include 'EvolAmas_SuiviMbins_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iIn
c
c var locales :
c ^^^^^^^^^^^^^
      integer iColNb
      character*64 line
      double precision value
c......................................................................
c
c---- Read specifications of mass-bins
c
      read(iIn,*) iMode_Mbins, iColNb

      if (iColNb.LE.1) then
         N_Mbins=-1
      else
         N_Mbins=0
      end if

      do while(.TRUE.) 
         read(iIn,*,END=9) line
         if (line(1:1).EQ.'#') goto 9
         N_Mbins = N_Mbins + 1
         if (N_Mbins.GT.N_MbinsMax) then
            write(0,*) '!!! Too many mass bins specified (max is ',
     $           N_MbinsMax,' !!!'
            call Meurt('Read_Mbins')
         end if
         if (iColNb.LE.1) then
            if (N_Mbins.GE.1) then
               Inf_Mbin(N_Mbins) = value
            end if
            read(line,*) value
            if (N_Mbins.GE.1) then
               Sup_Mbin(N_Mbins) = value
            end if
         else
            read(line,*) Inf_Mbin(N_Mbins), Sup_Mbin(N_Mbins)
         end if
      end do
 9    continue
c
c---- Read fractional masses
c
      if (line(1:1).EQ.'#') then
         
         N_MFracR_Mbin = 0
         do while(.TRUE.) 
            read(iIn,*,END=19) value
            N_MFracR_Mbin = N_MFracR_Mbin+1
            if (N_MFracR_Mbin.GT.N_MFracR_MbinMax) then
               write(0,*) '!!! Too many frac masses specified (max is ',
     $              N_MFracR_MbinMax,' !!!'
               call Meurt('Read_Mbins')
            end if
            MFracR_Mbin(N_MFracR_Mbin)=value
         end do
      end if
 19   continue
c
      end
c
c----------------------------------------------------------------------
      subroutine Init_MbinsMonitor_ASC(iOut)
c----------------------------------------------------------------------
c     read the specifications of mass-bins and the fractional masses
c     for radii to be monitored
c
      include 'EvolAmas_SuiviMbins_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iOut
c
c var locales :
c ^^^^^^^^^^^^^
      integer i,j, expo_coef10, k
      double precision coef10
      character*256 frm
c......................................................................
c
      expo_coef10 = -( int(log10(MFracR_Mbin(1))+100.0d0)-100 )
      expo_coef10 = max(expo_coef10,3)
      coef10 = 10.0d0**expo_coef10
      write(frm,200) expo_coef10+1,expo_coef10+1
 200  format('(I3,": R",I',I1,'.',I1,',"_",I2.2," ",$)')
      write(iOut,100) N_Mbins, N_MFracR_Mbin
 100  format(
     $     '# Monitoring of distribution of particles binned by',
     $     ' stellar mass',/,
     $     '# Number of bins : ',I3,/,
     $     '# Number of Lagr. radii followed : ',I3,/,
     $     '# 1: iPas_Evol 2: Time_FPunits 3: NbBoundPart 4: Mbound ',$)
      k=4
      do i=1,N_Mbins
         k=k+1
         write(iOut,110)k,i
 110     format( I3,': Mmin_',I2.2,' ',$)
         k=k+1
         write(iOut,115)k,i
 115     format( I3,': Mmax_',I2.2,' ',$)
         k=k+1
         write(iOut,120)k,i
 120     format( I3,': Mav_',I2.2,' ',$)
         k=k+1
         write(iOut,125)k,i
 125     format( I3,': fracM_',I2.2,' ',$)
         k=k+1
         write(iOut,130)k,i
 130     format( I3,': fracN_',I2.2,' ',$)

         do j=1,N_MFracR_Mbin
            k=k+1
            write(iOut,frm) k,NINT(coef10*MFracR_Mbin(j)),i
         end do
      end do

      write(iOut,*)
c
      end
c
c----------------------------------------------------------------------
      subroutine Write_MbinsMonitor_ASC(iOut)
c----------------------------------------------------------------------
      include 'EvolAmas_SuiviMbins_inc.f'
      include 'EvolAmas_common.f'
      include 'EvolAmas_Temps_common.f'
c
c argument :
c ^^^^^^^^^^
      integer iOut
c
c var locales :
c ^^^^^^^^^^^^^
      double precision Ri(iDimSE), Mi(iDimSE),
     $     Msort(iDimSE), Mcum(iDimSE), Mstar(iDimSE)
      double precision Mmin,Mmax,Mtot, Mmin_rq_bin,Mmax_rq_bin,
     $     Mmin_bin,Mmax_bin,Mtot_bin, Msum, Mlim, Mav_bin, Rlagr,
     $     Mtot_NB, R_km1
      integer ibin, iray, NbBoundPart, indx(iDimSE), Npart_bin,i,k,
     $     imin, imax, indx_bin(iDimSE)
c
c functions :
c ^^^^^^^^^^^
      double precision stellar_mass
c......................................................................
c
c---- Sort particles by stellar mass; Find min/max stellar masses
c
      NbBoundPart = 0
      Mtot=0.0d0
      Mtot_NB=0.0d0
      do k=1,NbSE
         Mstar(k) = stellar_mass(k)
      end do
      call indexx(NbSE,Mstar,indx)
      do k=1,NbSE
         if (R_SE(indx(k)).LT.R_perte) then
            NbBoundPart = NbBoundPart + 1
            indx(NbBoundPart) = indx(k)
            Mtot=Mtot+Mstar(indx(k))
            Mtot_NB=Mtot_NB+M_SE(indx(k))
            Mcum(NbBoundPart)=Mtot
            Msort(NbBoundPart)=Mstar(indx(k))
         end if
      end do
      Mmin=Msort(1)
      Mmax=Msort(NbBoundPart)

      write(iOut,100) iPas_evol, Tps_amas/UT_FP, NbBoundPart,
     $     Mtot_NB
 100  format(I10,' ',E12.5,' ',I7,' ',E12.5,' ',$)

      do ibin=1,N_Mbins
c
c---- Find required min/max masses of bin
c
         if (iMode_Mbins.EQ.1) then
            Mmin_rq_bin = Inf_Mbin(ibin)
            Mmax_rq_bin = Sup_Mbin(ibin)
         else if (iMode_Mbins.EQ.2) then
            Mmin_rq_bin = Mmin+Inf_Mbin(ibin)*(Mmax-Mmin)
            Mmax_rq_bin = Mmin+Sup_Mbin(ibin)*(Mmax-Mmin)
         else if (iMode_Mbins.EQ.3) then
            Mmin_rq_bin = Mstar(indx(max(
     $           1,int(Inf_Mbin(ibin)*NbBoundPart) )))
            Mmax_rq_bin = Mstar(indx(min(
     $           NbBoundPart,int(Sup_Mbin(ibin)*NbBoundPart)+1 )))
         else if (iMode_Mbins.EQ.4) then
            call locate(Mcum,NbBoundPart,Inf_Mbin(ibin)*Mtot,i)
            if (i.LE.0) then
               Mmin_rq_bin = Msort(1)
            else
               Mmin_rq_bin = 0.5d0*( Msort(min(i,NbBoundPart)) +
     $              Msort(min(i+1,NbBoundPart)) )
            end if
            call locate(Mcum,NbBoundPart,Sup_Mbin(ibin)*Mtot,i)
            if (i.LE.0) then
               Mmax_rq_bin = Msort(1)
            else
               Mmax_rq_bin = 0.5d0*( Msort(min(i,NbBoundPart)) +
     $              Msort(min(i+1,NbBoundPart)) )
               if (i.EQ.NbBoundPart) then
                  Mmax_rq_bin = Mmax+1d-6*(Mmax-Mmin)
               end if

            end if
         else
            write(0,*) '!!! Invalid value of iMode_Mbins: ',
     $           iMode_Mbins,' !!!'
            call Meurt('Write_MbinsMonitor_ASC')
         end if
c
c---- Extract stars in required bin...
c
         Mtot_bin=0.0d0
         Npart_bin=0
         call locate(Msort,NbBoundPart,Mmin_rq_bin,imin)
         call locate(Msort,NbBoundPart,Mmax_rq_bin,imax)
         imin=imin+1
         Mmin_bin=Msort(imin)
         Mmax_bin=Msort(imax)
         !write(0,*) '>>>>**',Mmin,Mmax
         !write(0,*) '>>>>>>>',ibin,Inf_Mbin(ibin),Sup_Mbin(ibin),Mmin_rq_bin,Mmax_rq_bin,Mmin_bin,Mmax_bin !##########
         if (imax.GT.imin) then
            do i=imin,imax
               Npart_bin=Npart_bin+1
               Ri(Npart_bin)=R_SE(indx(i))
               Mi(Npart_bin)=Msort(i)
               Mtot_bin=Mtot_bin+Msort(i)
            end do
            Mav_bin=Mtot_bin/dble(Npart_bin)
            write(iOut,110) Mmin_bin, Mmax_bin, Mav_bin,
     $           Mtot_bin/Mtot, dble(Npart_bin)/dble(NbBoundPart)
 110        format(5(E12.5,' '),$)
c
c---- ...Sort them according to distance from center...
c
            call indexx(Npart_bin,Ri,indx_bin)
c
c---- ...And determine the Lagrange radii
c
            Msum=0.0d0
            k=0
            do iray=1,N_MFracR_Mbin
               Mlim = min(MFracR_Mbin(iray),0.9999d0)*Mtot_bin
               do while (Msum.LT.Mlim .AND. k.LT.Npart_bin)
                  k=k+1
                  i=indx_bin(k)
                  Msum=Msum+Mi(i)
               end do
               if (k.EQ.Npart_bin .AND. Msum.LT.Mlim) then
                  Rlagr = Ri(indx_bin(Npart_bin))
               else
                  if (k.EQ.1) then
                     R_km1=0.0d0
                  else
                     R_km1=Ri(indx_bin(k-1))
                  end if
                  Rlagr = Ri(indx_bin(k)) - 
     $                 (Ri(indx_bin(k))-R_km1)*
     $                 (Msum-Mlim)/Mi(indx_bin(k))
               end if
               write(iOut,120) Rlagr
 120           format(E12.5,' ',$)
            end do

         else                   ! empty bin
            write(iOut,110) Mmin_rq_bin,Mmax_rq_bin,
     $           sqrt(Mmin_rq_bin*Mmax_rq_bin),0.0d0,0.0d0
            do iray=1,N_MFracR_Mbin
               write(iOut,120) -1.0d0
            end do
         end if

      end do

      write(iOut,*)            

c
      end
c
c----------------------------------------------------------------------
      block data BD_Init_SuivMbins
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_SuiviMbins_common.f'
      data iMode_Mbins /4/
      data N_Mbins, N_MFracR_Mbin /12,8/
      data Inf_Mbin /
     $     0.0d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0,
     $     0.8d0, 0.9d0,  0.95d0, 0.99d0,
     $     38*-1.0d0/
      data Sup_Mbin /
     $     0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, 0.7d0, 0.8d0,
     $     0.9d0, 0.95d0, 0.99d0, 1.0001d0,
     $     38*-1.0d0/
      data MFracR_Mbin /
     $     1d-3, 1d-2, 0.1d0, 0.25d0, 0.5d0, 0.75d0, 0.9d0, 0.99d0,
     $     12*-1.0d0/
      end
c----------------------------------------------------------------------
#endif
