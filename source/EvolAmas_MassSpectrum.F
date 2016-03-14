c======================================================================
c     subroutines to monitor the evolution of the stellar mass spetrum
c======================================================================

c----------------------------------------------------------------------
      subroutine Init_MassSpec_ASC(iOut)
c----------------------------------------------------------------------
      include 'EvolAmas_MassSpectrum_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iOut
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer i,j,ibase
c......................................................................
c
      write(iOut,100)
 100  format(
     $     '# Stellar cumulative mass spectrum',/,
     $     '# Masses are in  M_sun, Ages in yrs',/,
     $     '# 1: iPas_Evol 2: N_part 3: N_star 4: M_cluster ',$)
      j=4
      do i=1,N_MassSpec
         write(iOut,110) i+j, int(100.0d0*M_MassSpec(i)+0.5d0)
 110     format(I2,': M',I3.3,' ',$)
      end do
      j=j+N_MassSpec
      do i=1,N_MostMassivePart
         ibase = 4*(i-1)+j+1
         write(iOut,120) ibase, i, ibase+1, i, ibase+2, i, ibase+3, i
 120     format(I2,': iSE_max',I2.2,' ',I2,': iTet_max',I2.2,' ',
     $        I2,': Met_max',I2.2,' ',I2,': Age_max',I2.2,' ',$)
      end do
      write(iOut,*) ''
c
      end
c

c----------------------------------------------------------------------
      subroutine Write_MassSpec_ASC(iOut)
c----------------------------------------------------------------------
      include 'EvolAmas_MassSpectrum_inc.f'
      include 'Param_Amas.f'
      include 'EvolAmas_common.f'
      include 'VarSE.f'
      include 'EvolAmas_Pertes_inc.f'
      include 'EvolAmas_Temps_common.f'
c
c argument :
c ^^^^^^^^^^
      integer iOut
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer i,iSE,i_sp,N_part, iType
      double precision N_star,M_cluster,M_cum,M_sp, Mstar(iDimSE)
      integer indx(iDimSE)
c
c functions :
c ^^^^^^^^^^^
      double precision stellar_mass, stellar_dateofbirth
      integer*1 stellar_generaltype
c......................................................................
c
c---- total mass, particle number
c
      M_cluster=M_amas/Coef_ConvM
      N_part=NbSE-NbSE_perte
      N_star=(NbSE-NbSE_perte)*Net_par_SE
      write(iOut,100) iPas_Evol, N_part, N_star, M_cluster
 100  format(I11,' ',I9,' ',E12.5,' ',E12.5,' ',$)
c
c---- sort particles according to their stellar mass
c
      do iSE=1,NbSE
         Mstar(iSE) = stellar_mass(iSE)
      end do
         
      call indexx(NbSE,Mstar,indx)
c
c---- write cumulative mass spectrum
c
      M_cum=0.0d0
      i_sp=1
      M_sp=M_MassSpec(i_sp)*M_amas
      do i=1,NbSE
         iSE=indx(i)
         if (R_SE(iSE).LT.R_perte) then
            M_cum=M_cum+M_SE(iSE)
            do while (M_cum.GE.M_sp)
               write(iOut,110) Mstar(iSE)
 110           format(E12.5,' ',$)
               i_sp = i_sp+1
               M_sp=M_MassSpec(i_sp)*M_amas
               if (i_sp.GT.N_MassSpec) goto 9
            end do
         end if
      end do
 9    continue
c
c---- write properies of most massive particles
c
      i_sp = 1
      do i=NbSE,1,-1
         iSE=indx(i)
         if (R_SE(iSE).LT.R_perte) then
            iType = stellar_generaltype(iSE)
            write(iOut,120)iSE,iType,Mstar(iSE),
     $           Tps_amas/UneAnnee-stellar_dateofbirth(iSE)
 120        format(I9,' ',I2,' ',E12.5,' ',E12.5,' ',$)
            i_sp = i_sp+1
            if (i_sp.GT.N_MostMassivePart) goto 19
         end if
      end do
 19   continue
      write(iOut,*) ''

      call flush_MD(iOut)
c
      end
c
c----------------------------------------------------------------------
      block data BD_Init_MassSpec
c----------------------------------------------------------------------
      include 'EvolAmas_MassSpectrum_inc.f'

      data N_MassSpec /24/
      data M_MassSpec /
     $     0.00d0,0.05d0,0.10d0,0.15d0,0.20d0,0.25d0,0.30d0,0.35d0,
     $     0.40d0,0.45d0,0.50d0,0.55d0,0.60d0,0.65d0,0.70d0,0.75d0,
     $     0.80d0,0.85d0,0.90d0,0.95d0,
     $     0.96d0,0.97d0,0.98d0,0.99d0,26*-1.0d0 /
      end
c----------------------------------------------------------------------
