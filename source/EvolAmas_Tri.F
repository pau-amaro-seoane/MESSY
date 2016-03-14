c
c----------------------------------------------------------------------
      subroutine TriAmas
c----------------------------------------------------------------------
c     Procedure de tri Quicksort pour amas a la Henon
c
      include 'EvolAmas_Tri_inc.f'
c
      integer m, nstack
      parameter (m=7,nstack=200)
c
      integer istack(nstack),i,j,k,l,ir,jstack,itemp,iListeRayt
      double precision a
c
      do j=1,NbSE
         iListeRay(j)=j
      enddo
c
c---- placer les "sentinelles"
c
      iListeRay(0) = 0
      iListeRay(NbSE+1) = NbSE+1
c
      jstack=0
      l=1
      ir=NbSE
c
    1 continue
      if(ir-l.lt.m)then
         do j=l+1,ir
            iListeRayt=iListeRay(j)
            a=R_SE(iListeRayt)
            do i=j-1,1,-1
               if(R_SE(iListeRay(i)).le.a)go to 2
               iListeRay(i+1)=iListeRay(i)
            enddo
            i=0
    2       iListeRay(i+1)=iListeRayt
         enddo
         if(jstack.eq.0) goto 99
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack - 2
      else
         k=(l+ir)/2
         itemp=iListeRay(k)
         iListeRay(k)=iListeRay(l+1)
         iListeRay(l+1)=itemp
         if(R_SE(iListeRay(l+1)).gt.R_SE(iListeRay(ir)))then
            itemp=iListeRay(l+1)
            iListeRay(l+1)=iListeRay(ir)
            iListeRay(ir)=itemp
         end if
         if(R_SE(iListeRay(l)).gt.R_SE(iListeRay(ir)))then
            itemp=iListeRay(l)
            iListeRay(l)=iListeRay(ir)
            iListeRay(ir)=itemp
         end if
         if(R_SE(iListeRay(l+1)).gt.R_SE(iListeRay(l)))then
            itemp=iListeRay(l+1)
            iListeRay(l+1)=iListeRay(l)
            iListeRay(l)=itemp
         end if
c
         i=l+1
         j=ir
         iListeRayt=iListeRay(l)
         a=R_SE(iListeRayt)
c
    3    continue
         i=i+1
         if(R_SE(iListeRay(i)).lt.a)go to 3
    4    continue
         j=j-1
         if(R_SE(iListeRay(j)).gt.a)go to 4
         if(j.lt.i)go to 5
c
         itemp=iListeRay(i)
         iListeRay(i)=iListeRay(j)
         iListeRay(j)=itemp
         go to 3
c
    5    continue
         iListeRay(l)=iListeRay(j)
         iListeRay(j)=iListeRayt
         jstack=jstack + 2
         if(jstack.gt.nstack)then
            write(*,*)'nstack too small in indexx!'
            stop
         end if
c
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         end if
      end if
      go to 1

 99   continue
      i_TriAmas = i_TriAmas+1
      iSE_Rmin_CC = iListeRay(1)
      Rmin_CC = R_SE(iSE_Rmin_CC)
c
      end
c
c----------------------------------------------------------------------
      subroutine VerifTriAmas
c----------------------------------------------------------------------
c     Verification du tri des SE a des fins de deboguage
c     
      include 'EvolAmas_Tri_inc.f'

      integer i
c
      do i=1,NbSE-1
         if (R_SE(iListeRay(i+1)).LT.R_SE(iListeRay(i))) then
            write(0,*) '!!! Arbre Pot mal trie !!!'
            call Meurt('TrierAmas')
            stop
         end if
      end do
c      write(*,'(/,A)')" TRI CORRECT."
c
      end
c
c----------------------------------------------------------------------
      block data bd_TriAmas
c----------------------------------------------------------------------
      implicit none

      integer i_TriAmas
      common /common_triAmas/ i_TriAmas
      data i_TriAmas /0/
c
      end
