c======================================================================
c     routines liees a la limite radiale 'reflechissante'
c     pour code d'evolution d'amas
c======================================================================

      subroutine moins_que_rien
      end
c
#if (_LIMITE_EXT_>0)
c----------------------------------------------------------------------
      subroutine EcritLimExt_XDR(iFich_XDR)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_LimExt_common.f'
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
      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichLimExt_XDR)
c
      call XDR_Ecrit_d (iFich_XDR,'R_LimExt',R_LimExt)
      call XDR_Ecrit_d (iFich_XDR,'B_LimExt',B_LimExt)
      call XDR_Ecrit_d (iFich_XDR,'Mext_LimExt',Mext_LimExt)
      call XDR_Ecrit_d (iFich_XDR,'Egrav_ext_LimExt',Egrav_ext_LimExt)
      call XDR_Ecrit_d (iFich_XDR,'Ecin_ext_LimExt' ,Ecin_ext_LimExt )
      call XDR_Ecrit_d (iFich_XDR,'fNbEt_ext_LimExt',fNbEt_ext_LimExt)
c
      end
c
c----------------------------------------------------------------------
      subroutine LitLimExt_XDR(iFich_XDR)
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_LimExt_common.f'
      include 'EvolAmas_TN_common.f'
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
      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichLimExt_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     $        ' de reflection radiale valable !!!'
         call Meurt('LitLimExt_XDR')
      end if
      M_ctr=0.0d0
      M_TN=0.0d0
      write(0,*) '???? Why do I put M_TN=0 in  LitLimExt_XDR ???'
      call exit(1)
      do while (.TRUE.)

         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         if (ixdr_err.LT.0) goto 99

         if (NomVar.EQ.'R_LimExt') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,R_LimExt)
         else if (NomVar.EQ.'B_LimExt') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,B_LimExt)
         else if (NomVar.EQ.'Mext_LimExt') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Mext_LimExt)
         else if (NomVar.EQ.'Egrav_ext_LimExt') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Egrav_ext_LimExt)
         else if (NomVar.EQ.'Ecin_ext_LimExt') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,Ecin_ext_LimExt)
         else if (NomVar.EQ.'fNbEt_ext_LimExt') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,fNbEt_ext_LimExt)
         else
            write(0,*) '!!! Variable inconnue : ',NomVar,' !!!'
            call Meurt('LitLimExt_XDR')
         end if

      end do

 99   continue
c
      end
c
#endif  /* (_LIMITE_EXT_>0) */

c----------------------------------------------------------------------
      block data bd_init_LimExt
c----------------------------------------------------------------------
      implicit none
      include 'EvolAmas_LimExt_common.f'
      
      data R_LimExt /1.0d50/
      data B_LimExt /0.0d0/
c
      end
c----------------------------------------------------------------------
