c======================================================================
c     convertit un fichier d'amas au format XDR en fichier ascii
c     au format /rdb (sortie standard)
c======================================================================
c
      implicit none
c......................................................................
c
      include 'fxdr.inc'
c
c var locales :
c ^^^^^^^^^^^^^
      character*80 NomFich
      integer iXDR, ixdr_err
      integer iSeed
      data iSeed /26167/
c
c fonctions :
c ^^^^^^^^^^^
      integer iargc
c......................................................................
c
      if (iargc().LT.1) then
         write(0,*) '> Nom du fichier d''entree (format XDR) :'
         read(5,'(A)') NomFich
      else 
         call getarg(1,NomFich)
      end if

      iXDR = initxdr(NomFich,'r',.TRUE.)
      call LitAmas_XDR(iXDR)
      ixdr_err = ixdrclose(iXDR)

      call WriteCluster_Nbody(iSeed,6)
c
      end
c
c----------------------------------------------------------------------
      subroutine WriteCluster_Nbody(iRandSeed,iFich)
c----------------------------------------------------------------------
      include 'Param_Amas.f'
      include 'VarSE.f'
      include 'EvolAmas_Pertes_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iRandSeed,iFich
c
c constants :
c ^^^^^^^^^^^
      double precision TwoPi
      parameter (TwoPi=2.0d0*3.141592653589793116d0)
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      double precision
     $      x(iDimSE), y(iDimSE), z(iDimSE),
     $     vx(iDimSE),vy(iDimSE),vz(iDimSE),
     $      x_cm, y_cm, z_cm,
     $     vx_cm,vy_cm,vz_cm,
     $     M_tot,
     $     Xrand, cphi,sphi,ctheta,stheat,stau,ctau,
     $     Vtg,Vr,Vphi,Vtheta
      double precision iSE
c
c functions :
c ^^^^^^^^^^^
      double precision random
c......................................................................
c
      call init_random(iRandSeed)

      x_cm  = 0.0d0
      y_cm  = 0.0d0
      z_cm  = 0.0d0
      vx_cm = 0.0d0
      vy_cm = 0.0d0
      vz_cm = 0.0d0
      M_tot = 0.0d0

      do iSE=1,NbSE
         if (R_SE(iSE).LT.R_perte) then
            M_tot = M_tot + M_SE(iSE)
            Xrand=random()
            ctheta = 1.0d0-2.0d0*Xrand
            stheta = sqrt(1-ctheta**2)
            Xrand=random()
            cphi   = cos(TwoPi*Xrand)
            sphi   = sin(TwoPi*Xrand)
            x(iSE) = stheta*cphi*R_SE(iSE)
            y(iSE) = stheta*sphi*R_SE(iSE)
            z(iSE) = ctheta     *R_SE(iSE)

            Vtg = J_SE(iSE)/R_SE(iSE)
            Vr  = sqrt(2.0d0*T_SE(iSE)-Vtg**2)
            if (random().LT.0.5d0) Vr=-Vr
            Xrand=random()
            ctau = cos(TwoPi*Xrand)
            stau = sin(TwoPi*Xrand)
            Vphi   = ctau*Vtg
            Vtheta = stau*Vtg

            vx(iSE) = stheta*cphi*Vr - sphi*Vphi - ctheta*cphi*Vtheta
            vy(iSE) = stheta*sphi*Vr + cphi*Vphi - ctheta*sphi*Vtheta
            vz(iSE) = ctheta     *Vr             + stheta     *Vtheta

            x_cm  = x_cm  + M_SE(iSE)* x(iSE)
            y_cm  = y_cm  + M_SE(iSE)* y(iSE)
            z_cm  = z_cm  + M_SE(iSE)* z(iSE)
            vx_cm = vx_cm + M_SE(iSE)*vx(iSE)
            vy_cm = vy_cm + M_SE(iSE)*vy(iSE)
            vz_cm = vz_cm + M_SE(iSE)*vz(iSE)
         else
           x(iDimSE) = 0.0d0
           y(iDimSE) = 0.0d0
           z(iDimSE) = R_SE(iSE)
           vx(iDimSE) = 0.0d0
           vy(iDimSE) = 0.0d0
           vz(iDimSE) = 0.0d0
        end if
      end do
  
      ! adjust to CM reference frame

      x_cm  =  x_cm/M_tot
      y_cm  =  y_cm/M_tot
      z_cm  =  z_cm/M_tot
      vx_cm = vx_cm/M_tot
      vy_cm = vy_cm/M_tot
      vz_cm = vz_cm/M_tot
      write(0,105) 
     $      x_cm, y_cm, z_cm,
     $     vx_cm,vy_cm,vz_cm
 105  format(
     $     '> CM position  x, y, z = ',3(E12.5,' '),/,
     $     '> CM velocity vx,vy,vz = ',3(E12.5,' '))

      Nbound = 0
      do iSE=1,NbSE
         if (R_SE(iSE).LT.R_perte) then
            Nbound = Nbound+1
            x(iSE) = x(iSE) - x_cm
            y(iSE) = y(iSE) - y_cm
            z(iSE) = z(iSE) - z_cm
            vx(iSE) = vx(iSE) - vx_cm
            vy(iSE) = vy(iSE) - vy_cm
            vz(iSE) = vz(iSE) - vz_cm
            write(iFich,100) M_SE(iSE),x(iSE),y(iSE),z(iSE),
     $           vx(iSE),vy(iSE),vz(iSE)
 100        format(7(E20.13,' '))
         end if
      end do

      write(0,110) Nbound
 110  format('> ASCII N-body file written; contains ',I8,
     $     ' particles')
c
      end
c
