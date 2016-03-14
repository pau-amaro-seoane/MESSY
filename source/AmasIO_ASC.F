c======================================================================
c     routines d'ecriture de fichier d'amas en ascii
c======================================================================
#define _HAUTE_PRECISION_ 1
c----------------------------------------------------------------------
      subroutine EcritAmas_ASC(iFich,lArbreComplet)
c----------------------------------------------------------------------
      include 'AmasIO_ASC_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
      logical lArbreComplet     ! indique si il faut ecrire tout l'arbre (.TRUE.)
                                ! ou seulement les SE a proprement parler (.FALSE.)
c
c var locale :
c ^^^^^^^^^^^^
      integer i
c......................................................................
c
      NbNoeuds=MAX(NbNoeuds,NbSE) ! au cas ou NbNoeuds serait encore a 0
c
c--- Header
c
      write(iFich,100)NbSE

 100  format(
     +     '# Data for spherical star cluster (N-body units) ',/,
     +     '#================================================',/,
     +     '#',/,
     +     '# NbSE           = ',I8,' # Number of particles')
      if (lArbreComplet) then
         write(iFich,102)NbNoeuds
 102     format(
     +     '# NbNoeuds       = ',I8,' # Number of nodes in binary tree')
      end if
      if (R_mar.GT.0.0d0 .AND. R_mar.LT.1d30) then
         write(iFich,104) R_mar
 104     format(
     +     '# R_mar          = ',E22.15,' # Tidal truncation radius')
      end if
      if (M_TN.GT.0.0d0 .OR. M_ctr.GT.0.0d0) then
         write(iFich,106) M_TN, M_ctr, iType_CtrObj, BirthDate_CtrObj
 106     format(
     $     '# M_TN           = ',E22.15,' # Mass of central object',/,
     $     '# M_Tctr         = ',E22.15,' # Total central mass ',
     $        '(with unresolved, unaccreted gas)',/,
     $     '# iType_CtrObj   = ',I2,' # Type of central object',/,
     $     '# BirthDate_CtrObj = ',E22.15,
     $        ' # Birth date of central object '
     $        )
      end if

      write(iFich,108)
 108  format(
     +     '# 1: R_SE 2: M_SE 3: T_SE 4: J_SE' )
c
c--- Particle data
c
      do i=1,NbSE
         write(iFich,110) R_SE(i),M_SE(i),T_SE(i),J_SE(i)
#if (_HAUTE_PRECISION_>0)
 110     format(4(E22.15,' '))
#else
 110     format(4(E15.8,' '))
#endif
      end do
c
c--- ecriture des donnees des noeuds vides
c
      if (lArbreComplet) then
         do i=NbSE+1,NbNoeuds
            write(iFich,120) R_SE(i),M_SE(i)
#if (_HAUTE_PRECISION_>0)
 120     format(2(E20.13,' '),'0.0 0.0')
#else
 120     format(2(E15.8,' '),'0.0 0.0')
#endif
         end do
      end if
c
      end
c
c#if defined(sun)
c----------------------------------------------------------------------
      subroutine EcritAmas_RDB(iFich,lArbreComplet)
c----------------------------------------------------------------------
      include 'AmasIO_ASC_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFich
      logical lArbreComplet     ! indique si il faut ecrire tout l'arbre (.TRUE.)
                                ! ou seulement les SE a proprement parler (.FALSE.)
c
c var locale :
c ^^^^^^^^^^^^
      integer i
      character*1 t_            ! tabulation
c......................................................................
c
      NbNoeuds=MAX(NbNoeuds,NbSE) ! au cas ou NbNoeuds serait encore a 0
      t_ = CHAR(9)
c
c--- ecriture en-tete
c
      write(iFich,100) 'R_SE',t_,'M_SE',t_,'T_SE',t_,'J_SE'
      write(iFich,100) '----',t_,'----',t_,'----',t_,'----'
 100  format(3(A4,A1),A4)
c
c--- ecriture des donnees des SE
c
      do i=1,NbSE
         write(iFich,110) R_SE(i),t_,M_SE(i),t_,T_SE(i),t_,J_SE(i)
#if (_HAUTE_PRECISION_>0)
 110     format(3(E22.15,A1),E22.15)
#else
 110     format(3(E15.8,A1),E15.8)
#endif
      end do
c
c--- ecriture des donnees des noeuds vides
c
      if (lArbreComplet) then
         do i=NbSE+1,NbNoeuds
            write(iFich,110) R_SE(i),t_,M_SE(i),t_,0.0d0,t_,0.0d0
         end do
      end if
c
      end
c
c#endif
c
c----------------------------------------------------------------------
      subroutine EcritAmas_Nbody(iFich,iRandSeed)
c----------------------------------------------------------------------
c
c     Write each SE as a (m;x,y,z;vx,vy,vz) particle.
c
c     Do not do anything clever or fancy, just select position at random 
c     on sphere and select direction of tangential velocity and sign of
c     radial velocity at random too.
c
c     For refined MC --> Nbody conversion, see Amas_XDR2NCadavre.F
c
      include 'AmasIO_ASC_inc.f'
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
     $     Mass(iDimSE+1),
     $     M_tot, mi,
     $     Xrand, cphi,sphi,ctheta,stheta,stau,ctau,
     $     Vtg,Vr,Vphi,Vtheta
      integer iSE, iCadavre, nCadavre
c
c functions :
c ^^^^^^^^^^^
      double precision random
c......................................................................
c
      call init_random(iRandSeed)

      M_tot = 0.0d0
      nCadavre = 0
      do iSE=1,NbSE
         if (R_SE(iSE).LT.R_perte) then
            nCadavre = nCadavre+1
            Mi = M_SE(iSE)
            Mass(nCadavre) = Mi
            M_tot = M_tot + Mi
            Xrand=random()
            ctheta = 1.0d0-2.0d0*Xrand
            stheta = sqrt(1-ctheta**2)
            Xrand=random()
            cphi   = cos(TwoPi*Xrand)
            sphi   = sin(TwoPi*Xrand)
            x(nCadavre) = stheta*cphi*R_SE(iSE)
            y(nCadavre) = stheta*sphi*R_SE(iSE)
            z(nCadavre) = ctheta     *R_SE(iSE)

            Vtg = J_SE(iSE)/R_SE(iSE)
            Vr  = sqrt(2.0d0*T_SE(iSE)-Vtg**2)
            if (random().LT.0.5d0) Vr=-Vr
            Xrand=random()
            ctau = cos(TwoPi*Xrand)
            stau = sin(TwoPi*Xrand)
            Vphi   = ctau*Vtg
            Vtheta = stau*Vtg

            vx(nCadavre) = stheta*cphi*Vr - sphi*Vphi
     $           - ctheta*cphi*Vtheta
            vy(nCadavre) = stheta*sphi*Vr + cphi*Vphi
     $           - ctheta*sphi*Vtheta
            vz(nCadavre) = ctheta     *Vr
     $           + stheta     *Vtheta
        end if
      end do
c
c---- Write particle data in ascii format
c
      write(iFich,200) nCadavre
 200  format(
     $     '# N-body cluster data ',/,
     $     '# Npart = ',I7,/,
     $     '# 1:m 2:x 3:y 4:z 5:vx 6:vy 7:vz')

      do iCadavre=1,nCadavre
         write(iFich,100) Mass(iCadavre),
     $         x(iCadavre), y(iCadavre), z(iCadavre),
     $        vx(iCadavre),vy(iCadavre),vz(iCadavre)
 100     format(7(E22.15,' '))
      end do
c
      end
c
