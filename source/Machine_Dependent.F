c _COMPIL_

c 1: old sun f77 (199X)
c 2: g77 (all platforms?)
c 3: pgf77 (linux?)
c 4: intel fortran 95 (linux?)
c 5: sun f95 (2007)
c 6: g95 & gfortran (2007)

c----------------------------------------------------------------------
      integer function getpid_MD()
c----------------------------------------------------------------------
      implicit none
      integer getpid
#if (_CONDOR_==0)
#if (_COMPIL_<=6 )
      getpid_MD = getpid() ! works with g77 on linux PCs and with f77 on SUNs
                           ! and with 'ifc -Vaxlib' (intel compiler)
      return
#endif
#else
      getpid_MD = -999
      return
#endif
      write(0,*) '!!! getpid_MD',
     $     ' not implemented for this compiler/machine !!!'
      call exit(1)
c
      end
c
c----------------------------------------------------------------------
      integer function time_MD()
c----------------------------------------------------------------------
      implicit none
      integer time
#if (_CONDOR_==0)
#if (_COMPIL_==1 || _COMPIL_==2 || _COMPIL_==3) 
      time_MD = time() ! works with g77 on linux PCs and with f77 on SUNs
      return
#endif 
#if ( _COMPIL_==4 || _COMPIL_==5 || _COMPIL_==6) 
      call SYSTEM_CLOCK(COUNT=time_MD)
      return
#endif
#else
      time_MD = 0.0d0
      return
#endif
      write(0,*) '!!! time_MD',
     $     ' not implemented for this compiler/machine !!!'
      call exit(1)
c
      end
c
c----------------------------------------------------------------------
      subroutine flush_MD(iUnit)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      integer iUnit
c......................................................................
c
#if (_NO_FLUSH_==0)
#   if (_COMPIL_<=6) /* version sparc, g77 et pgf77 */
c      write(0,*) '>>> calling flush for unit ', iUnit !#########
      call flush(iUnit)         ! works with g77 on linux PCs and with f77 on SUNs
                                ! and with 'ifc -Vaxlib' (intel compiler)
c      write(0,*) '>>> done calling flush ' !#########
      return
#   endif
      write(0,*) '!!! flush_MD',
     $     ' not implemented for this compiler/machine !!!'
      call exit(1)
      
#endif
      continue
c
      end
c
c----------------------------------------------------------------------
      subroutine QuelleHeure_MD(ih,im,is)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      integer ih,im,is          ! heure,minute,secondes
c
c var locale :
c ^^^^^^^^^^^^
      integer iarray(3), ivalues(8)
c......................................................................
c
#if (_CONDOR_==0)
#if (_COMPIL_==1 || _COMPIL_==2 || _COMPIL_==3)
      call itime(iarray)
      ih=iarray(1)
      im=iarray(2)
      is=iarray(3)
      return
#endif
#if ( _COMPIL_==4 || _COMPIL_==5  || _COMPIL_==6)
      call DATE_AND_TIME(VALUES=ivalues)
      ih=ivalues(5)
      im=ivalues(6)
      is=ivalues(7)
      return
#endif
#else
      ih=99
      im=99
      is=99
      return
#endif
      write(0,*) '!!! QuelleHeure_MD',
     $     ' not implemented for this compiler/machine !!!'
      call exit(1)
c
      end
c
c----------------------------------------------------------------------
      subroutine QuelleHeure_str(str)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      character*(*) str
c
c var locales :
c ^^^^^^^^^^^^^
      integer ih,im,is 
c......................................................................
c
      call QuelleHeure_MD(ih,im,is)
      write(str,100) ih,im,is
 100  format(I2.2,'h',I2.2,'m',I2.2,'s')
c
      end
c
c----------------------------------------------------------------------
      subroutine QuelleDate_MD(ij,im,ia)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      integer ij,im,ia          ! jour,mois,annee
c
c var locale :
c ^^^^^^^^^^^^
      integer iarray(3), ivalues(8)
c......................................................................
c
#if (_CONDOR_==0)
#if (_COMPIL_==1 || _COMPIL_==2) /* version sparc, g77 */
      call idate (iarray) ! version sparc et g77
      ij = iarray(1)
      im = iarray(2)
      ia = iarray(3)
      return
#endif
#if ( _COMPIL_==3) /* vesion pgf77 */ 
      call idate (im,ij,ia)
      if (ia.LT.50) then 
         ia=2000+ia
      else
         ia=1900+ia
      end if
      return
#endif
#if ( _COMPIL_==4 || _COMPIL_==5 || _COMPIL_==6 ) /* fortran 95 */ 
      call DATE_AND_TIME(VALUES=ivalues)
      ia=ivalues(1)
      im=ivalues(2)
      ij=ivalues(3)
      return
#endif
#else
      ia=9999
      im=99
      ij=99
      return
#endif
      write(0,*) '!!! QuelleDate_MD',
     $     ' not implemented for this compiler/machine !!!'
      call exit(1)
c
      end
c
c----------------------------------------------------------------------
      subroutine QuelleDate_str(str)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      character*(*) str
c
c var locales :
c ^^^^^^^^^^^^^
      integer ij,im,ia
c......................................................................
c
      call QuelleDate_MD(ij,im,ia)
      write(str,100) ij,im,ia
 100  format(I2.2,'.',I2.2,'.',I4.4)
c
      end
c

c----------------------------------------------------------------------
      subroutine Hostname_MD(petit_nom)
c----------------------------------------------------------------------
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      character*(*) petit_nom
c
c var locale :
c ^^^^^^^^^^^^
      integer istat
c......................................................................
c
#if (_CONDOR_==0)
#if (_COMPIL_==1 || _COMPIL_==2 || _COMPIL_==3) /* version sparc, g77 et pgf77 */
      call HostNm(petit_nom,istat)
      if (istat.NE.0) petit_nom='_unknown_'
      return
#endif
#if ( _COMPIL_==4 || _COMPIL_==5) /* version ifc (fortran 95) */ 
      call HostNm(petit_nom)
      return
#endif
#if ( _COMPIL_==6) /* g95 */ 
      istat=HostNm(petit_nom)
      if (istat.NE.0) petit_nom='_unknown_'
      return
#endif
#else
      petit_nom='XXX_condor_xxx'
      return
#endif
      write(0,*) '!!! Hostname_MD',
     $     ' not implemented for this compiler/machine !!!'
      call exit(1)
c
      end
c
#if (_COMPIL_<=6) /* version sparc, g77, pgf77 et ifc */
c----------------------------------------------------------------------
      real function TpsUtil_MD() 
c----------------------------------------------------------------------
      implicit none
      real ttot,tarray(2)
      real etime
c
#if (_CONDOR_==0)
      ttot = etime(tarray)
      TpsUtil_MD = tarray(1)
#else
      TpsUtil_MD = 0.0d0
#endif
c
      end
c----------------------------------------------------------------------
      real function TpsSyst_MD() 
c----------------------------------------------------------------------
      implicit none
      real ttot,tarray(2)
      real etime
c
#if (_CONDOR_==0)
      ttot = etime(tarray)
      TpsSyst_MD = tarray(2)
#else
      TpsSyst_MD = 0.0d0
#endif
c
      end
c
#endif
#if ( _COMPIL_==666) /* version ifc (fortran 95) */ 
c----------------------------------------------------------------------
      subroutine Set_TrapFPE_MD
c----------------------------------------------------------------------
      implicit none
      integer result(4), EnableMask(4)

      INTERFACE 
          INTEGER(4) FUNCTION FOR_SET_FPE(EnableMask) 
              INTEGER(4) EnableMask 
          END FUNCTION 
      END INTERFACE

      EnableMask
      result = FOR_SET_FPE(EnableMask)
c
      end
c
#endif

#if (_CONDOR_==0)
c----------------------------------------------------------------------
      integer function df_MD(path,Filesystem,OneK_Blocks,Used,Available,
     $     UsePercent,Mounted_on)
c----------------------------------------------------------------------
c     gives disk usage information for directory path, see 'man df'
c     returns 0 if OK
c
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      character*(*) path
      character*(*) Filesystem,Mounted_on
      integer*8 OneK_Blocks,Used,Available
      integer UsePercent
c
c local variables :
c ^^^^^^^^^^^^^^^^^
      character*1 dum
      character*256 direc, line
      character*64 tmp_file
      integer i1,i2, i
      double precision x
c
c functions :
c ^^^^^^^^^^^
      integer getpid_MD
c......................................................................
c
      write(tmp_file,'(I8.8)') getpid_MD()
      direc=path
      tmp_file = '/tmp/_df_MD_' // tmp_file

      call system_MD('df -k ' // direc(1:index(direc,' ')-1) // ' > ' //
     $     tmp_file,df_MD)

      if (df_MD.EQ.0) then

         open(11,file=tmp_file)
         read(11,'(A1)') dum
         read(11,'(A)') line
         close(11,STATUS='DELETE')

         if (df_MD.EQ.0) then
            i1=1
            do while(line(i1:i1).NE.' ')
               i1=i1+1
            end do
            Filesystem=line(1:i1-1) 
            do while(line(i1:i1).EQ.' ')
               i1=i1+1
            end do
            i2=i1
            do while(line(i2:i2).NE.' ')
               i2=i2+1
            end do
            x=456.78e4
            do i=1,100000
               x=sin(x)
            end do
            !! Code chokes here if compiled with g77 in 64 bits abd run on the Sun grid!!
            read(line(i1:i2-1),*) OneK_Blocks

            i1=i2     
            do while(line(i1:i1).EQ.' ')
               i1=i1+1
            end do
            i2=i1
            do while(line(i2:i2).NE.' ')
               i2=i2+1
            end do
            read(line(i1:i2-1),*) Used
            
            i1=i2
            do while(line(i1:i1).EQ.' ')
               i1=i1+1
            end do
            i2=i1
            do while(line(i2:i2).NE.' ')
               i2=i2+1
            end do
            read(line(i1:i2-1),*) Available
            
            i1=i2     
            do while(line(i1:i1).EQ.' ')
               i1=i1+1
            end do
            i2=i1
            do while(line(i2:i2).NE.' ')
               i2=i2+1
            end do
            read(line(i1:i2-2),*) UsePercent
            
            i1=i2     
            do while(line(i1:i1).EQ.' ')
               i1=i1+1
            end do
            i2=len(line)    
            do while(line(i2:i2).EQ.' ')
               i2=i2-1
            end do
            Mounted_on=line(i1:i2)
         end if
      end if
c
      end
c
#endif /* if (_CONDOR_==0) */

#if (_CONDOR_==0)
c----------------------------------------------------------------------
      subroutine system_MD(command,exit_status)
c----------------------------------------------------------------------
c
      implicit none
c
c arguments :
c ^^^^^^^^^^^
      character*(*) command
      integer exit_status
c
c functions :
c ^^^^^^^^^^^
#if ( _COMPIL_==3 || _COMPIL_==4 || _COMPIL_==5) /* Portland Group & Intel compilers */ 
      integer system
#endif
c......................................................................
c
#if ( _COMPIL_==3 || _COMPIL_==4 || _COMPIL_==5) /* Portland Group & Intel compilers */
      exit_status = system(command)
#else
c      write(0,*) '>>> calling system with ',command !########
      call system(command,exit_status)
c      write(0,*) '>>> done calling system with ',command !########
#endif
c
      end
#endif /* if (_CONDOR_==0) */

      
