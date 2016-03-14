c======================================================================
c     Subroutines to keep track of cumulative probabilities
c     of collision and large-angle scatterings
c======================================================================

#if (_WRITE_CUMUL_COLL_PROB_>0)
c----------------------------------------------------------------------
      subroutine Init_CumulPcoll
c----------------------------------------------------------------------
      implicit none
      include 'Param_Amas.f'
      include 'VarSE.f'
      include 'EvolAmas_CumulProb_common.f'
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer iSE
c......................................................................
c
      do iSE=1,NbSE
         CumulPcoll_SE(iSE) = 0.0d0
         CumulPkick_SE(iSE) = 0.0d0
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine Write_CumulPcoll_ASC(iFile)
c----------------------------------------------------------------------
      implicit none
      include 'Param_Amas.f'
      include 'VarSE.f'
      include 'EvolAmas_CumulProb_common.f'
      include 'EvolAmas_Pertes_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFile
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer iSE
c......................................................................
c
      write(iFile,100)
 100  format(
     $     '# Cumulative probabilities for:',/,
     $     '#          * collisions',/,
     $     '#          * large-angle scatterings (kicks)',/,
     $     '# Probabilities are accumulated from t=0',/,
     $     '# Negative values are used for particles not',
     $     ' in the cluster anymore',/,
     $     '#',/,
     $     '# 1: CumPcoll 2: CumPkick')
      do iSE=1,NbSE
         if (R_SE(iSE).LT.R_perte) then
            write(iFile,110) CumulPcoll_SE(iSE), CumulPkick_SE(iSE)
         else
            write(iFile,110) -CumulPcoll_SE(iSE), -CumulPkick_SE(iSE)
         end if
 110     format(2(' ',E12.5))
      end do
c
      end
c
c----------------------------------------------------------------------
      subroutine Read_CumulPcoll_ASC(iFile,iOK) ! iOK=0 if OK
c----------------------------------------------------------------------
      implicit none
      include 'Param_Amas.f'
      include 'VarSE.f'
      include 'EvolAmas_CumulProb_common.f'
      include 'EvolAmas_Pertes_inc.f'
c
c arguments :
c ^^^^^^^^^^^
      integer iFile,iOK
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer nSE
      character*1024 OneLine
      double precision p1,p2
c......................................................................
c
      iOK = 0

c---- Skip header

      OneLine='###'
      do while (OneLine(1:1).EQ.'#')
         read(iFile,'(A)') OneLine
      end do

c---- First particle
      
      nSE=1
      read(OneLine,*,ERR=999) p1,p2
      CumulPcoll_SE(nSE) = abs(p1) 
      CumulPkick_SE(nSE) = abs(p2)

c---- All other particles

      do while(.TRUE.)
         read(iFile,*,END=9,ERR=999) p1,p2
         nSE=nSE+1
         if (nSE.GT.NbSE) then
            write(0,*) '!!! Too many lines in cumul prob file !!!'
            iOK = 1
            return
         end if
         CumulPcoll_SE(nSE) = abs(p1) 
         CumulPkick_SE(nSE) = abs(p2)
      end do

 9    continue
      if (nSE.NE.NbSE) then
         write(0,*) '!!! Wrong number of lines in cumul prob file !!!'
         iOK = 2
         return
      end if

      iOK = 0
      return

 999  continue
      write(0,*) '!!! Error reading data for iSE = ',nSE,' !!!'
      iOK = 3
c
      end
c      
         
#endif
