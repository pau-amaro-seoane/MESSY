#define _TYPE_RAND_ 2
#define _PLTFM_ 1
c----------------------------------------------------------------------
      subroutine init_random(iseed)
c----------------------------------------------------------------------
      implicit none
c
c argument :
c ^^^^^^^^^^
      integer iseed
c
c bloc common :
c ^^^^^^^^^^^^^
      integer idum_random
      logical linit_random
      common /BC_random/ idum_random, linit_random
c
c variable locale :
c ^^^^^^^^^^^^^^^^^
      double precision dum
c
c fonctions :
c ^^^^^^^^^^^
#if (_PLTFM_==1) /* Sparc */
      double precision DRAND
#elif (_PLTFM_==2) /* Silicon Graphics */
      double precision SRAND
#endif
      double precision ran2, ran3
c......................................................................
c
#if   (_TYPE_RAND_==0) /* Generateur systeme : mauvais mais rapide */
#if (_PLTFM_==1) /* Sparc */
      dum = drand(1)
#elif (_PLTFM_==2) /* Silicon Graphics */
      dum = srand(iseed)
#endif
#elif (_TYPE_RAND_==2) /* ran2 des numerical recipes */
      idum_random = -abs(iseed)
      dum = ran2(idum_random)
#elif (_TYPE_RAND_==3) /* ran3 des numerical recipes */
      idum_random = -abs(iseed)
      dum = ran3(idum_random)
#endif
      end
c
c----------------------------------------------------------------------
      double precision function random()
c----------------------------------------------------------------------
      implicit none
c
c bloc common :
c ^^^^^^^^^^^^^
      integer idum_random
      logical linit_random
      common /BC_random/ idum_random, linit_random
c
c fonctions :
c ^^^^^^^^^^^
#if (_PLTFM_==1) /* Sparc */
      double precision DRAND
#elif (_PLTFM_==2) /* Silicon Graphics */
      double precision RAND
#endif
      double precision ran2, ran3
c......................................................................
c
#if (_TYPE_RAND_==0) /* Generateur systeme : mauvais mais rapide */
#if (_PLTFM_==1) /* Sparc */
      random = drand(0)
#elif (_PLTFM_==2) /* Silicon Graphics */
      random = rand()
#else
      write(0,*) '!!! _PLTFM_ not defined in random.F !!!'
      stop 
#endif
#elif (_TYPE_RAND_==2) /* ran2 des numerical recipes */
      random = ran2(idum_random)
#elif (_TYPE_RAND_==3) /* ran3 des numerical recipes */
      random = ran3(idum_random)
#else
      write(0,*) '!!! _TYPE_RAND_ not defined in random.F !!!'
      stop
#endif
c
      end
c
c----------------------------------------------------------------------
      subroutine save_random(iOut)
c----------------------------------------------------------------------
c     Sauvegrade les variables internes du generateur de nb aleatoires
c     dans le fichier binaire de numero iOut
c     Permet ainsi de continuer une sequence donnee de nb aleatoires
c
      implicit none
c
c argument :
c ^^^^^^^^^^
      integer iOut
c
c blocs common :
c ^^^^^^^^^^^^^^
      integer idum_random
      logical linit_random
      common /BC_random/ idum_random, linit_random

#if (_TYPE_RAND_==2)    
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     variables d'etat du generateur
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer NTAB
      parameter (NTAB=32)
      integer iv(NTAB),iy,idum2
      common /common_ran2/ iv,iy,idum2
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c......................................................................
c
      write(iOut,'(35(I10," "))') idum_random,iv,iy,idum2
#elif (_TYPE_RAND_==3)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     variables d'etat du generateur
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer iff,inext,inextp
      integer ma(55)
      common /common_ran3/  iff,inext,inextp,ma
c......................................................................
c
      write(iOut,'(59(I10," "))') idum_random,iff,inext,inextp,ma
#else
      write(0,*) '!!! la routine "save_random" n''est pas disponible',
     $     ' pour ce generateur aleatoire !!!'
      stop
#endif
      end
c
c----------------------------------------------------------------------
      subroutine read_random(iIn)
c----------------------------------------------------------------------
c     Lit les variables internes du generateur de nb aleatoires
c     dans le fichier binaire de numero iIn
c     Permet ainsi de continuer une sequence donnee de nb aleatoires
c
      implicit none
c
c argument :
c ^^^^^^^^^^
      integer iIn
c
c blocs common :
c ^^^^^^^^^^^^^^
      integer idum_random
      logical linit_random
      common /BC_random/ idum_random, linit_random

#if (_TYPE_RAND_==2)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     variables d'etat du generateur
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer NTAB
      parameter (NTAB=32)
      integer iv(NTAB),iy,idum2
      common /common_ran2/ iv,iy,idum2
c......................................................................
c
      read(iIn,*) idum_random,iv,iy,idum2
#elif (_TYPE_RAND_==3)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     variables d'etat du generateur
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer iff,inext,inextp
      integer ma(55)
      common /common_ran3/  iff,inext,inextp,ma
c......................................................................
c
      read(iIn,*)  idum_random,iff,inext,inextp,ma
#else
      write(0,*) '!!! la routine "read_random" n''est pas disponible',
     $     ' pour ce generateur aleatoire !!!'
      stop
#endif
      end
c
#if (_XDR_>0)
c----------------------------------------------------------------------
      subroutine EcritRandom_XDR(iFich_XDR)
c----------------------------------------------------------------------
c     Sauvegrade les variables internes du generateur de nb aleatoires
c     dans fichier binaire au format XDR
c     Permet ainsi de continuer une sequence donnee de nb aleatoires
c
      implicit none
      include 'fxdr.inc'
c
c argument :
c ^^^^^^^^^^
      integer iFich_XDR
c
c var locale :
c ^^^^^^^^^^^^
      integer ixdr_err
c
c blocs common :
c ^^^^^^^^^^^^^^
      integer idum_random
      logical linit_random
      common /BC_random/ idum_random, linit_random

#if (_TYPE_RAND_==2)    
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     variables d'etat du generateur
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer NTAB
      parameter (NTAB=32)
      integer iv(NTAB),iy,idum2
      common /common_ran2/ iv,iy,idum2

      character*(*) EnTeteFichRand_XDR
      parameter (EnTeteFichRand_XDR='%%% XDR ran 2 File %%%')
c......................................................................
c
      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichRand_XDR)

      call XDR_Ecrit_i (iFich_XDR,'idum_random',idum_random)
      call XDR_Ecrit_ti(iFich_XDR,'iv',NTAB,iv)
      call XDR_Ecrit_i (iFich_XDR,'iy',iy)
      call XDR_Ecrit_i (iFich_XDR,'idum2',idum2)

#elif (_TYPE_RAND_==3)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     variables d'etat du generateur
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer iff,inext,inextp
      integer ma(55)
      common /common_ran3/  iff,inext,inextp,ma

      character*(*) EnTeteFichRand_XDR
      parameter (EnTeteFichRand_XDR='%%% XDR ran 3 File %%%')
c......................................................................
c
      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichRand_XDR)

      call XDR_Ecrit_i (iFich_XDR,'idum_random',idum_random)
      call XDR_Ecrit_i (iFich_XDR,'iff',iff)
      call XDR_Ecrit_i (iFich_XDR,'inext',inext)
      call XDR_Ecrit_i (iFich_XDR,'inextp',inextp)
      call XDR_Ecrit_ti(iFich_XDR,'ma',55,ma)
#else
      write(0,*) '!!! la routine "EcritRandom_XDR" n''est pas disponible',
     $     ' pour ce generateur aleatoire !!!'
      stop
#endif
      end
c
c----------------------------------------------------------------------
      subroutine LitRandom_XDR(iFich_XDR)
c----------------------------------------------------------------------
c     Relit les variables internes du generateur de nb aleatoires
c     dans fichier binaire au format XDR
c     Permet ainsi de continuer une sequence donnee de nb aleatoires
c
      implicit none
      include 'fxdr.inc'
c
c argument :
c ^^^^^^^^^^
      integer iFich_XDR
c
c var locales :
c ^^^^^^^^^^^^^
      character*20 NomType,NomVar
      character*80 EnTete
      integer idum, ixdr_err
c
c fonction :
c ^^^^^^^^^^
      integer ixdr_read_string
c
c blocs common :
c ^^^^^^^^^^^^^^
      integer idum_random
      logical linit_random
      common /BC_random/ idum_random, linit_random

#if (_TYPE_RAND_==2)    
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     variables d'etat du generateur
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer NTAB
      parameter (NTAB=32)
      integer iv(NTAB),iy,idum2
      common /common_ran2/ iv,iy,idum2

      character*(*) EnTeteFichRand_XDR
      parameter (EnTeteFichRand_XDR='%%% XDR ran 2 File %%%')
c......................................................................
c
      ixdr_err = ixdr_read_string(iFich_XDR,EnTete)
      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichRand_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     $        ' ran2 XDR valable !!!'
         call Meurt('LitRandom_XDR')
      end if
    
      do while (.TRUE.)

         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         if (ixdr_err.LT.0) goto 99

         if (NomVar.EQ.'idum_random') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,idum_random)
         else if (NomVar.EQ.'iv') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrimat(iFich_XDR,NTAB,iv)
         else if (NomVar.EQ.'iy') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,iy)
         else if (NomVar.EQ.'idum2') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,idum2)
         else
            write(0,*) '!!! Variable inconnue : ',NomVar,' !!!'
            call Meurt('LitRandom_XDR')
         end if

      end do

 99   continue
c
#elif (_TYPE_RAND_==3)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     variables d'etat du generateur
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer iff,inext,inextp
      integer ma(55)
      common /common_ran3/  iff,inext,inextp,ma

      character*(*) EnTeteFichRand_XDR
      parameter (EnTeteFichRand_XDR='%%% XDR ran 3 File %%%')
c......................................................................
c
      ixdr_err = ixdr_read_string(iFich_XDR,EnTete)
      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichRand_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     $        ' ran3 XDR valable !!!'
         call Meurt('LitRandom_XDR')
      end if
    
      do while (.TRUE.)

         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         if (ixdr_err.LT.0) goto 99

         if (NomVar.EQ.'idum_random') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,idum_random)
         else if (NomVar.EQ.'iff') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,iff)
         else if (NomVar.EQ.'inext') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,inext)
         else if (NomVar.EQ.'inextp') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrint(iFich_XDR,inextp)
         else if (NomVar.EQ.'ma') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrmat(iFich_XDR,55,ma)
         else
            write(0,*) '!!! Variable inconnue : ',NomVar,' !!!'
            call Meurt('LitRandom_XDR')
         end if

      end do

 99   continue
#else
      write(0,*) '!!! la routine "LitRandom_XDR" n''est pas disponible',
     $     ' pour ce generateur aleatoire !!!'
      stop
#endif
      end
c
#endif /* if (_XDR_>0) */

#if (_TYPE_RAND_==2) /* ran2 des numerical recipes */
c----------------------------------------------------------------------
      FUNCTION ran2(idum)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     variables d'etat du generateur
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer NTAB
      parameter (NTAB=32)
      integer idum2,iv(NTAB),iy
      common /common_ran2/ iv,iy,idum2
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NDIV
      DOUBLE PRECISION ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NDIV=1+IMM1/NTAB,EPS=1.2d-16,RNMX=1.d0-EPS)
      INTEGER j,k

      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
#endif

#if (_TYPE_RAND_==3) /* ran3 des numerical recipes */
      FUNCTION ran3(idum)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     variables d'etat du generateur
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer iff,inext,inextp
      integer ma(55)
      common /common_ran3/  iff,inext,inextp,ma
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     DOUBLE PRECISION MBIG,MSEED,MZ
      DOUBLE PRECISION ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,ii,k
      INTEGER mj,mk
C     DOUBLE PRECISION mj,mk,ma(55)

      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END
#endif
c----------------------------------------------------------------------
      block data random_BD
      integer idum_random
      logical linit_random
      common /BC_random/ idum_random, linit_random
      data idum_random  /12345677/
      data linit_random /.FALSE./
      end

#if (_TYPE_RAND_==2)
      block data ran2_BD
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     variables d'etat du generateur ran2
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer NTAB
      parameter (NTAB=32)
      integer idum2,iv(NTAB),iy
      common /common_ran2/ iv,iy,idum2
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      end
#elif (_TYPE_RAND_==3)
      block data ran3_BD
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     variables d'etat du generateur ran3
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer iff,inext,inextp
      integer ma(55)
      common /common_ran3/  iff,inext,inextp,ma
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      end
#endif
c----------------------------------------------------------------------
