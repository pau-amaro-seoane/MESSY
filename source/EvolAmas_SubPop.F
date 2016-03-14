#if (_FOLLOW_SUBPOP_>0)

c----------------------------------------------------------------------
      subroutine Read_ListSubPop(iFich)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_SubPop_common.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer iread, j
c......................................................................
c
      N_SubPop=0
      do while (.TRUE.)
         read(iFich,*,END=99) iread
         if (N_SubPop.EQ.Nmax_SubPop) then
            write(0,*) '!!! The max number of member of a',
     $           ' subpopulation to be tracked is ',Nmax_SubPop,' !!!'
            call Meurt('Read_SpecificationSubPop')
         end if
         j=N_SubPop
         N_SubPop=N_SubPop+1
                                ! Stupid bubble sort
         do while (iread.LT.iSE_SubPop(j) .AND. j.GE.1)
            iSE_SubPop(j+1) = iSE_SubPop(j)
            j=j-1
         end do
         iSE_SubPop(j+1)=iread
      end do
 99   continue
c
      end
c
c----------------------------------------------------------------------
      subroutine write_SubPop(iFich)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_SubPop_common.f'
      include 'Param_Amas.f'
      include 'VarSE.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      integer i,iSE
c......................................................................
c
      write(iFich,100)
 100  format(
     $     '# Subpopulation tracking',/,
     $     '# 1: iSE 2: M_SE 3: R_SE 4: T_SE 5: J_SE')
      do i=1,N_SubPop
         iSE=iSE_SubPop(i)
         write(iFich,110) iSE, M_SE(iSE),R_SE(iSE),T_SE(iSE),J_SE(iSE)
 110     format(I8,4(' ',E12.5))
      end do
c   
      end
c
#endif
