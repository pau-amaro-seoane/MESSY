c======================================================================
c     routines liees a la limitation par champ de maree pour evolution
c     d'amas
c======================================================================

c----------------------------------------------------------------------
      subroutine read_Maree(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Maree_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
      read(iFich) R_Mar,C_Mar
      call regulariser_CR_Maree
c
      end
c
c----------------------------------------------------------------------
      subroutine write_Maree(iFich)
c----------------------------------------------------------------------
      include 'EvolAmas_Maree_inc.f'
c
c argument :
c ^^^^^^^^^^
      integer iFich
c......................................................................
c
      call regulariser_CR_Maree
      write(iFich) R_Mar,C_Mar   
c
      end
c
#ifndef DONT_USE_XDR
c----------------------------------------------------------------------
      subroutine EcritMaree_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_Maree_inc.f'
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
      call regulariser_CR_Maree
c
      ixdr_err = ixdrstring(iFich_XDR,EnTeteFichMaree_XDR)
c
      call XDR_Ecrit_d (iFich_XDR,'C_Mar',C_Mar)
      call XDR_Ecrit_d (iFich_XDR,'R_Mar',R_Mar)
c
      end
cc----------------------------------------------------------------------
      subroutine LitMaree_XDR(iFich_XDR)
c----------------------------------------------------------------------
      include 'EvolAmas_Maree_inc.f'
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
      if (ixdr_err.LT.0 .OR. EnTete.NE.EnTeteFichMaree_XDR) then
         write(0,*) '!!! Le fichier n''est pas un fichier',
     $        ' de maree XDR valable !!!'
         call Meurt('LitMaree_XDR')
      end if
    
      do while (.TRUE.)

         ixdr_err = ixdr_read_string(iFich_XDR,NomVar)
         if (ixdr_err.LT.0) goto 99

         if (NomVar.EQ.'C_Mar') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,C_Mar)
         else if (NomVar.EQ.'R_Mar') then
            ixdr_err = ixdr_read_string(iFich_XDR,NomType)
            ixdr_err = ixdrint(iFich_XDR,idum)
            ixdr_err = ixdrdouble(iFich_XDR,R_Mar)
         else
            write(0,*) '!!! Variable inconnue : ',NomVar,' !!!'
            call Meurt('LitMaree_XDR')
         end if

      end do

 99   continue
      call regulariser_CR_Maree
c
      end
c
#endif
c----------------------------------------------------------------------
      double precision function R_tronc_mar()
c----------------------------------------------------------------------
      include 'EvolAmas_Maree_inc.f'
      include 'EvolAmas_Conserv_inc.f'
c......................................................................
c
      if (C_mar.GE.1.0d20 .OR. C_mar.LE.0.0d0) then
         R_mar = 1.0d20
      else
#if (_FIGER_RMAR_>0)
         R_mar = min(1.0d20,C_Mar*(Mamas_ini+Mtn_ini)**0.333333333333d0)
#else
         R_mar = min(1.0d20,C_Mar*(M_amas+M_ctr)**0.333333333333d0)
#endif
      end if
      R_tronc_mar = R_mar
c
      end
c
c----------------------------------------------------------------------
      subroutine Regulariser_CR_Maree
c----------------------------------------------------------------------
      include 'EvolAmas_Maree_inc.f'
c
c var locale :
c ^^^^^^^^^^^^
      double precision M
c......................................................................
c
      if (M_amas.LE.0.0d0) then
         M = 1.0d0
      else
         M = M_amas+M_ctr
      end if
c$$$
c$$$      if (R_mar.LE.0.0d0) then
c$$$         if (C_Mar.GT.0.0d0) then
c$$$            R_Mar = C_Mar*M**0.3333333333d0
c$$$         else
c$$$            R_mar = 1.0d20
c$$$         end if
c$$$      end if
c$$$      R_mar=min(R_mar,1.0d20)

      if (C_mar.LE.0.0d0) then
         if (R_Mar.GT.0.0d0) then
            C_Mar = R_Mar/M**0.3333333333d0
         else
            C_mar = 1.0d20
         end if
      end if
      C_mar=min(C_mar,1.0d20)
      R_Mar = C_Mar*M**0.3333333333d0
      R_mar=min(R_mar,1.0d20)
c
      end
c
c----------------------------------------------------------------------
      block data bd_init_Maree
c----------------------------------------------------------------------
      include 'EvolAmas_Maree_common.f'
      
      data R_Mar,C_Mar /2*1.0d20/ ! valeurs par defaut : pas de troncature de maree
c
      end
c----------------------------------------------------------------------
