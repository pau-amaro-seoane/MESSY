c======================================================================
c     routines accessoires liees au temps propres et pas de temps
c     pour EvolAmas.F
c======================================================================

c----------------------------------------------------------------------
      subroutine EcrireTempsBIN(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Temps_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c
c var locales :
c ^^^^^^^^^^^^^
      integer i
c......................................................................
c
      write(iFich)NbSE,Tps_amas,
     $     (Tps_SE(i),i=1,NbSE)
c
      end
c
c----------------------------------------------------------------------
      subroutine EcrireTempsRDB(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Temps_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c
c var locales :
c ^^^^^^^^^^^^^
      integer i
      character*1 t_            ! tabulation
c......................................................................
c
      t_ = CHAR(9)
c
c--- ecriture en-tete
c
      write(iFich,100) 'Tps_SE'
      write(iFich,100) '------'
 100  format(A)
c
c--- ecriture des donnees
c
      do i=1,NbSE
         write(iFich,110) Tps_SE(i)
 110     format(E15.8)
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine EcrireTempsASC(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Temps_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c
c var locales :
c ^^^^^^^^^^^^^
      integer i
c......................................................................
c
      write(iFich,100)NbSE,Tps_amas
 100  format(
     $     '# Fichier des temps pour Amas',/,
     $     '# NbSE = ',I7,/,
     $     '# Tps_amas = ',E12.4,/,
     $     '# 1: Tps_SE' )
      do i=1,NbSE
         write(iFich,110) Tps_SE(i)
 110     format(E14.6)
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine LireTempsBIN(iFich,lCalc_Tps_Amas)
c----------------------------------------------------------------------
      include 'EvolAmas_Temps_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
      logical lCalc_Tps_Amas    ! indique si il faut calculer le temps
                                ! global (median) de l'amas
c
c var locales :
c ^^^^^^^^^^^^^
      integer i
c......................................................................
c
      read(iFich,end=999,err=999)NbSE,Tps_amas,
     $     (Tps_SE(i),i=1,NbSE)
      if (lCalc_Tps_Amas) call Det_TpsAmas

      return

 999  continue
      write(0,*)'!!! Fichier ''Temps'' defectueux !!!'
      call Meurt('LireTempsBIN')
c
      end
c
c----------------------------------------------------------------------
      subroutine EcritTemps_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_Temps_inc.f'
      include 'fxdr.inc'
c
c argument :
c ^^^^^^^^^^
      integer iFich_XDR
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err
c......................................................................
c
      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichTemps_XDR)
c
      call XDR_Ecrit_d (iFich_XDR,'Tps_Amas',Tps_Amas)
      call XDR_Ecrit_td(iFich_XDR,'Tps_SE'  ,NbSE,Tps_SE)
c
      end
c
c----------------------------------------------------------------------
      subroutine LitTemps_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_Temps_inc.f'
      include 'fxdr.inc'
c
c argument :
c ^^^^^^^^^^
      integer iFich_XDR
c
c variables locales :
c ^^^^^^^^^^^^^^^^^^^
      character*20 NomType,NomVar
      character*80 EnTete
      integer idum, ixdr_err
c
c fonction :
c ^^^^^^^^^^
      integer ixdr_read_string
c......................................................................
c
      ixdr_err = ixdr_read_string(iFich_XDR,EnTete)
      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichTemps_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     $        ' de temps XDR valable !!!'
         call Meurt('LitTemps_XDR')
      end if
    
      do while (.TRUE.)

         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         if (ixdr_err.LT.0) goto 99

         if (NomVar.EQ.'Tps_SE') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,NbSE)
            ixdr_err = ixdrdmat(iFich_XDR,NbSE,Tps_SE)
         else if (NomVar.EQ.'Tps_Amas') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Tps_Amas)
         else
            write(0,*) '!!! Variable inconnue : ',NomVar,' !!!'
            call Meurt('LitTemps_XDR')
         end if

      end do

 99   continue
      if (NbSE.GT.iDimSE) then
         write(0,*) '!!! NbSE(',NbSE,') > iDimSE(',iDimSE,') !!!'
         call Meurt('LitTemps_XDR')
      end if
c
      end
c
c----------------------------------------------------------------------
      subroutine Det_TpsAmas
c----------------------------------------------------------------------
      include 'EvolAmas_Temps_inc.f'
c
c var locales :
c ^^^^^^^^^^^^^
      integer i,N
c
c fonction(s) :
c ^^^^^^^^^^^^^
      double precision selip
c......................................................................
c
      Tps_amas = selip((NbSE-NbSE_perte)/2,NbSE,Tps_SE)
      dTps_amas_inf = Tps_amas - selip((NbSE-NbSE_perte)/6,NbSE,Tps_SE)
      dTps_amas_sup = selip(5*(NbSE-NbSE_perte)/6,NbSE,Tps_SE)
     $     - Tps_amas
c
c--- temps moyen
c
      Tps_amas_moy = 0.0d0
      Sigma_Tps_amas = 0.0d0
      N = 0
      do i=1,NbSE
         if (R_SE(i).LT.R_perte) then
            Tps_amas_moy = Tps_amas_moy+Tps_SE(i)
            N=N+1
         end if
      end do
      Tps_amas_moy = Tps_amas_moy/N
c
c--- dispertion
c
      do i=1,NbSE
         if (R_SE(i).LT.R_perte) then
            Sigma_Tps_amas = Sigma_Tps_amas +
     $           (Tps_SE(i)-Tps_amas_moy)**2
         end if
      end do
      Sigma_Tps_amas = sqrt(Sigma_Tps_amas/N)
c
      end
c
c----------------------------------------------------------------------
      block data BD_Init_Temps
c----------------------------------------------------------------------
      implicit none
      logical lCalc_Tps_Amas
      common /Flag_TempsIO/ lCalc_Tps_Amas
      
      data lCalc_Tps_Amas /.TRUE./ ! indique qu'il faut calculer le temps de
                                   ! l'amas lors de la lecture d'un fichier de tps
      end
c----------------------------------------------------------------------
